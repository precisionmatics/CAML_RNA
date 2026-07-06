"""
Step 10 — Multiple Kernel Learning (MKL)

Combines:
  K_ph  = RBF kernel on ES_CS+WCh PH features (normalized)
  K_lig = Tanimoto kernel on Morgan ECFP4 fingerprints
  K_fm  = RBF kernel on RNA-FM embeddings

K_combined = alpha*K_ph + beta*K_lig + gamma*K_fm
where alpha+beta+gamma=1, tuned in inner CV.

SVR with precomputed kernel on K_combined.
Then applies hybrid (aptamer + ribosomal_asite subtype override).
"""

import numpy as np
import pandas as pd
import os, warnings
from scipy.stats import pearsonr, spearmanr
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, ParameterGrid, LeaveOneOut, GridSearchCV
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.metrics import mean_squared_error, r2_score

warnings.filterwarnings('ignore')

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

def load_feat(fname):
    d = np.load(os.path.join(FEAT_DIR, fname), allow_pickle=True)
    return d['X'].astype(np.float64), d['y'].astype(np.float64)

def metrics(y_true, y_pred):
    r,  _ = pearsonr(y_true, y_pred)
    rho,_ = spearmanr(y_true, y_pred)
    rmse  = np.sqrt(mean_squared_error(y_true, y_pred))
    r2    = r2_score(y_true, y_pred)
    return r, rho, rmse, r2

def bootstrap_ci(y_true, y_pred, n=10000):
    rng = np.random.default_rng(SEED)
    ns  = len(y_true)
    rs  = [pearsonr(y_true[idx:=rng.integers(0,ns,ns)],
                    y_pred[idx])[0] for _ in range(n)]
    return np.percentile(rs, 2.5), np.percentile(rs, 97.5)

def rbf_kernel(X, Y=None, gamma=None):
    if Y is None: Y = X
    if gamma is None: gamma = 1.0 / X.shape[1]
    XX = np.sum(X**2, axis=1, keepdims=True)
    YY = np.sum(Y**2, axis=1, keepdims=True)
    D2 = XX + YY.T - 2*X@Y.T
    return np.exp(-gamma * np.maximum(D2, 0))

def tanimoto_kernel(X, Y=None):
    """Tanimoto (Jaccard) kernel for binary fingerprints."""
    if Y is None: Y = X
    XY = X @ Y.T
    XX = X.sum(axis=1, keepdims=True)
    YY = Y.sum(axis=1, keepdims=True)
    denom = XX + YY.T - XY
    denom = np.maximum(denom, 1e-10)
    return XY / denom

def normalize_kernel(K):
    d = np.sqrt(np.diag(K))
    d[d == 0] = 1.0
    return K / np.outer(d, d)

def build_combined_kernel(K_ph, K_lig, K_fm, alpha, beta, gamma_fm):
    return alpha*K_ph + beta*K_lig + gamma_fm*K_fm

def svr_precomputed(K_train, y_train, K_test, C=10.0):
    svr = SVR(kernel='precomputed', C=C)
    svr.fit(K_train, y_train)
    return svr.predict(K_test)

def nested_mkl_cv(K_ph, K_lig, K_fm, y):
    """Nested 5-fold CV with MKL weight + SVR-C tuning in inner."""
    n = len(y)
    outer = KFold(5, shuffle=True, random_state=SEED)
    inner = KFold(5, shuffle=True, random_state=SEED)
    oof   = np.zeros(n)

    alpha_grid = [0.6, 0.7, 0.8, 0.9]
    beta_grid  = [0.05, 0.10, 0.20]
    C_grid     = [1, 10, 100]

    for tr, te in outer.split(np.arange(n)):
        best_r_in, best_params = -99, None

        for alpha in alpha_grid:
            for beta in beta_grid:
                gf = max(1 - alpha - beta, 0)
                for C in C_grid:
                    # inner CV
                    inner_oof = np.zeros(len(tr))
                    for tr2, te2 in inner.split(tr):
                        tri, tei = tr[tr2], tr[te2]
                        K_comb_tr = build_combined_kernel(
                            K_ph[np.ix_(tri,tri)],
                            K_lig[np.ix_(tri,tri)],
                            K_fm[np.ix_(tri,tri)],
                            alpha, beta, gf)
                        K_comb_te = build_combined_kernel(
                            K_ph[np.ix_(tei,tri)],
                            K_lig[np.ix_(tei,tri)],
                            K_fm[np.ix_(tei,tri)],
                            alpha, beta, gf)
                        inner_oof[te2] = svr_precomputed(
                            K_comb_tr, y[tri], K_comb_te, C)
                    try:
                        r_in = pearsonr(y[tr], inner_oof)[0]
                    except Exception:
                        r_in = -99
                    if r_in > best_r_in:
                        best_r_in  = r_in
                        best_params = (alpha, beta, gf, C)

        a, b, g, C = best_params
        K_comb_tr = build_combined_kernel(
            K_ph[np.ix_(tr,tr)], K_lig[np.ix_(tr,tr)],
            K_fm[np.ix_(tr,tr)], a, b, g)
        K_comb_te = build_combined_kernel(
            K_ph[np.ix_(te,tr)], K_lig[np.ix_(te,tr)],
            K_fm[np.ix_(te,tr)], a, b, g)
        oof[te] = svr_precomputed(K_comb_tr, y[tr], K_comb_te, C)
        print(f"    fold best params: α={a:.2f} β={b:.2f} γ={g:.2f} C={C}")

    return oof

def svr_subtype(X_ph, y, use_loo=False):
    """Reuse step07 subtype SVR logic."""
    from sklearn.pipeline import Pipeline
    from sklearn.feature_selection import VarianceThreshold
    pipe = Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])
    grid = {'model__C':[0.1,1,10,100,500], 'model__gamma':['scale','auto']}
    cv_out = LeaveOneOut() if use_loo else KFold(5, shuffle=True, random_state=SEED)
    cv_in  = KFold(min(5, len(y)-1), shuffle=True, random_state=SEED)
    oof = np.zeros(len(y))
    for tr, te in cv_out.split(X_ph):
        if len(tr) < 5:
            oof[te] = y.mean()
            continue
        gs = GridSearchCV(pipe, grid, cv=cv_in, scoring='r2', n_jobs=-1)
        gs.fit(X_ph[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X_ph[te])
    return oof

def main():
    df  = pd.read_csv(DATA_CSV)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values

    # Load features
    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X_ph = np.hstack([X_es, X_wch])  # 10260
    X_lig, _ = load_feat('step08_ligand_rna.npz')
    X_lig = X_lig[:, :2048]          # Morgan only
    X_fm,  _ = load_feat('step09_rnafm.npz')

    # Normalize PH features before building kernel
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.pipeline import Pipeline
    from sklearn.feature_selection import VarianceThreshold
    pipe_ph = Pipeline([
        ('vt', VarianceThreshold(1e-4)),
        ('sc', StandardScaler()),
        ('pc', PCA(n_components=100, random_state=SEED)),
    ])
    X_ph_red = pipe_ph.fit_transform(X_ph)

    pipe_fm = Pipeline([
        ('sc', StandardScaler()),
        ('pc', PCA(n_components=50, random_state=SEED)),
    ])
    X_fm_red = pipe_fm.fit_transform(X_fm)

    # Build normalized kernels
    print("Building kernels...")
    gamma_ph = 1.0 / X_ph_red.shape[1]
    gamma_fm = 1.0 / X_fm_red.shape[1]

    K_ph  = normalize_kernel(rbf_kernel(X_ph_red, gamma=gamma_ph))
    K_lig = normalize_kernel(tanimoto_kernel(X_lig))
    K_fm  = normalize_kernel(rbf_kernel(X_fm_red, gamma=gamma_fm))

    print(f"K_ph range:  [{K_ph.min():.3f}, {K_ph.max():.3f}]")
    print(f"K_lig range: [{K_lig.min():.3f}, {K_lig.max():.3f}]")
    print(f"K_fm range:  [{K_fm.min():.3f}, {K_fm.max():.3f}]")

    print("\nRunning nested MKL CV (α*K_ph + β*K_lig + γ*K_fm)...")
    oof_mkl = nested_mkl_cv(K_ph, K_lig, K_fm, y)
    r, rho, rmse, r2 = metrics(y, oof_mkl)
    lo, hi = bootstrap_ci(y, oof_mkl)
    print(f"\nGlobal MKL:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}  "
          f"95%CI [{lo:.3f},{hi:.3f}]")
    np.save(os.path.join(MODEL_DIR, 'oof_pred_mkl.npy'), oof_mkl)

    # Hybrid: override aptamer + ribosomal_asite with subtype SVR
    oof_hybrid = oof_mkl.copy()

    for st, use_loo in [('aptamer', False), ('ribosomal_asite', True)]:
        idx = np.where(sub == st)[0]
        oof_st = svr_subtype(X_ph[idx], y[idx], use_loo=use_loo)
        r_st = pearsonr(y[idx], oof_st)[0]
        r_gl = pearsonr(y[idx], oof_mkl[idx])[0]
        print(f"\n{st}: subtype r={r_st:.4f}  mkl r={r_gl:.4f}")
        if r_st > r_gl:
            oof_hybrid[idx] = oof_st
            print(f"  → using subtype model")
        else:
            print(f"  → keeping MKL")

    r, rho, rmse, r2 = metrics(y, oof_hybrid)
    lo, hi = bootstrap_ci(y, oof_hybrid)
    print(f"\n{'='*55}")
    print(f"MKL HYBRID:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*55}")

    # Per-subtype
    print("\nPer-subtype breakdown:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        try:
            rs, *_ = metrics(y[idx], oof_hybrid[idx])
            print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")
        except Exception:
            print(f"  {st:20s}  n={len(idx):2d}  r=N/A")

    np.save(os.path.join(MODEL_DIR, 'oof_pred_step10_mkl_hybrid.npy'), oof_hybrid)
    pd.DataFrame({'pdb': df['pdb'].values, 'y_true': y,
                  'y_pred': oof_hybrid, 'subtype': sub}).to_csv(
        os.path.join(RES_DIR, 'step10_predictions.csv'), index=False)
    print("\nSaved.")

if __name__ == '__main__':
    main()
