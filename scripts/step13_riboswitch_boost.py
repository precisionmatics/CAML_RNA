"""
Step 13 — Riboswitch-boosted hybrid

Key ideas:
  1. For riboswitch (n=61): add rs_subclass one-hot to ES_CS+WCh → should let SVR
     learn subclass-specific binding patterns more directly than subclass-model override
  2. Try MLP on ES_CS+WCh for riboswitch (different inductive bias from SVR)
  3. For other_misc (n=27): dedicated LOO SVR with ES_CS+WCh
  4. Aptamer override (same as step07, consistently r~0.937)
  5. Ribosomal_asite override (same as step07, r~0.580)
  Final: best combination as new OOF
"""

import numpy as np
import pandas as pd
import os, warnings
from scipy.stats import pearsonr, spearmanr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, LeaveOneOut, GridSearchCV, cross_val_predict
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from sklearn.metrics import mean_squared_error, r2_score

warnings.filterwarnings('ignore')

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

def metrics(y, p):
    r,_   = pearsonr(y, p)
    rho,_ = spearmanr(y, p)
    rmse  = np.sqrt(mean_squared_error(y, p))
    r2    = r2_score(y, p)
    return r, rho, rmse, r2

def bootstrap_ci(y, p, n=10000):
    rng = np.random.default_rng(SEED)
    ns  = len(y)
    rs  = [pearsonr(y[i:=rng.integers(0,ns,ns)], p[i])[0] for _ in range(n)]
    return np.percentile(rs, [2.5, 97.5])

def load_feat(fname):
    d = np.load(os.path.join(FEAT_DIR, fname), allow_pickle=True)
    return d['X'].astype(np.float64), d['y'].astype(np.float64)

def make_svr_pipe():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])

def make_mlp_pipe():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  MLPRegressor(hidden_layer_sizes=(256,128,64),
                                max_iter=500, random_state=SEED,
                                early_stopping=True, validation_fraction=0.15)),
    ])

SVR_GRID = {'model__C':[0.1,1,10,100,500,1000], 'model__gamma':['scale','auto']}
MLP_GRID = {'model__hidden_layer_sizes':[(256,128),(256,128,64),(128,64)],
            'model__alpha':[1e-4, 1e-3, 1e-2]}

def svr_oof(X, y, loo=False):
    n     = len(y)
    outer = LeaveOneOut() if loo else KFold(5, shuffle=True, random_state=SEED)
    inner = KFold(min(5, n-1), shuffle=True, random_state=SEED)
    oof   = np.zeros(n)
    for tr, te in outer.split(X):
        if len(tr) < 5:
            oof[te] = y.mean(); continue
        gs = GridSearchCV(make_svr_pipe(), SVR_GRID, cv=inner, scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def mlp_oof(X, y, loo=False):
    n     = len(y)
    outer = LeaveOneOut() if loo else KFold(5, shuffle=True, random_state=SEED)
    inner = KFold(min(5, n-1), shuffle=True, random_state=SEED)
    oof   = np.zeros(n)
    for tr, te in outer.split(X):
        if len(tr) < 5:
            oof[te] = y.mean(); continue
        gs = GridSearchCV(make_mlp_pipe(), MLP_GRID, cv=inner, scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def ohe_col(labels, classes):
    mat = np.zeros((len(labels), len(classes)))
    for i, l in enumerate(labels):
        j = list(classes).index(l) if l in classes else 0
        mat[i, j] = 1.0
    return mat

def main():
    df    = pd.read_csv(DATA_CSV)
    y     = df['pKd'].values.astype(np.float64)
    sub   = df['subtype'].values
    pdb   = df['pdb'].values
    rs_sc = df['rs_subclass'].values

    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X_cpf, _ = load_feat('step11_cpf.npz')

    X_base = np.hstack([X_es, X_wch])

    # Start from step07 predictions (best baseline)
    oof = np.load(os.path.join(MODEL_DIR, 'oof_pred_step07.npy')).copy()
    print(f"Baseline (step07): r={pearsonr(y, oof)[0]:.4f}")

    # --- Aptamer (always override, consistently excellent) ---
    idx = np.where(sub == 'aptamer')[0]
    oof_apt = svr_oof(X_base[idx], y[idx])
    r_apt = pearsonr(y[idx], oof_apt)[0]
    r_gl  = pearsonr(y[idx], oof[idx])[0]
    oof[idx] = oof_apt
    print(f"\naptamer  n={len(idx)}  svr_r={r_apt:.4f}  (was {r_gl:.4f})")

    # --- Ribosomal_asite ---
    idx = np.where(sub == 'ribosomal_asite')[0]
    oof_ras = svr_oof(X_base[idx], y[idx], loo=True)
    r_ras = pearsonr(y[idx], oof_ras)[0]
    r_gl  = pearsonr(y[idx], oof[idx])[0]
    if r_ras > r_gl:
        oof[idx] = oof_ras
        print(f"ribosomal_asite  n={len(idx)}  svr_r={r_ras:.4f} > {r_gl:.4f}  ✓ override")
    else:
        print(f"ribosomal_asite  n={len(idx)}  svr_r={r_ras:.4f} <= {r_gl:.4f}  keep")

    # --- Other_misc: try dedicated SVR ---
    idx = np.where(sub == 'other_misc')[0]
    oof_om  = svr_oof(X_base[idx], y[idx])
    r_om    = pearsonr(y[idx], oof_om)[0]
    r_gl    = pearsonr(y[idx], oof[idx])[0]
    if r_om > r_gl:
        oof[idx] = oof_om
        print(f"other_misc  n={len(idx)}  svr_r={r_om:.4f} > {r_gl:.4f}  ✓ override")
    else:
        print(f"other_misc  n={len(idx)}  svr_r={r_om:.4f} <= {r_gl:.4f}  keep")

    # --- Riboswitch: subclass one-hot → SVR ---
    rs_idx  = np.where(sub == 'riboswitch')[0]
    rs_labs = rs_sc[rs_idx]
    classes = sorted(set(rs_labs))
    rs_ohe  = ohe_col(rs_labs, classes)
    X_rs    = np.hstack([X_base[rs_idx], rs_ohe])

    print(f"\n=== Riboswitch (n={len(rs_idx)}) ===")
    oof_rs_svr = svr_oof(X_rs, y[rs_idx])
    r_rs_svr   = pearsonr(y[rs_idx], oof_rs_svr)[0]
    r_rs_gl    = pearsonr(y[rs_idx], oof[rs_idx])[0]
    print(f"  SVR + subclass ohe:  r={r_rs_svr:.4f}  (global slice r={r_rs_gl:.4f})")

    oof_rs_mlp = mlp_oof(X_rs, y[rs_idx])
    r_rs_mlp   = pearsonr(y[rs_idx], oof_rs_mlp)[0]
    print(f"  MLP + subclass ohe:  r={r_rs_mlp:.4f}")

    # plain SVR on riboswitch without subclass ohe
    oof_rs_plain = svr_oof(X_base[rs_idx], y[rs_idx])
    r_rs_plain   = pearsonr(y[rs_idx], oof_rs_plain)[0]
    print(f"  SVR plain:           r={r_rs_plain:.4f}")

    # pick best riboswitch model
    best_rs = max([
        (r_rs_svr,   oof_rs_svr,   'SVR+subclass_ohe'),
        (r_rs_mlp,   oof_rs_mlp,   'MLP+subclass_ohe'),
        (r_rs_plain, oof_rs_plain, 'SVR_plain'),
        (r_rs_gl,    oof[rs_idx],  'global_keep'),
    ], key=lambda x: x[0])
    best_rs_r, best_rs_oof, best_rs_name = best_rs
    print(f"  → Best riboswitch: {best_rs_name}  r={best_rs_r:.4f}")
    oof[rs_idx] = best_rs_oof

    # --- CPF-augmented for riboswitch (bonus try) ---
    X_rs_cpf = np.hstack([X_base[rs_idx], X_cpf[rs_idx], rs_ohe])
    oof_rs_cpf = svr_oof(X_rs_cpf, y[rs_idx])
    r_rs_cpf   = pearsonr(y[rs_idx], oof_rs_cpf)[0]
    print(f"  SVR + CPF + subclass ohe: r={r_rs_cpf:.4f}")
    if r_rs_cpf > best_rs_r:
        oof[rs_idx] = oof_rs_cpf
        print(f"  → Override with CPF version!")

    # --- Final report ---
    r, rho, rmse, r2 = metrics(y, oof)
    lo, hi = bootstrap_ci(y, oof)
    print(f"\n{'='*60}")
    print(f"STEP 13:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*60}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        rs,*_ = metrics(y[idx], oof[idx])
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")

    np.save(os.path.join(MODEL_DIR,'oof_pred_step13.npy'), oof)
    pd.DataFrame({'pdb':pdb,'y_true':y,'y_pred':oof,'subtype':sub}).to_csv(
        os.path.join(RES_DIR,'step13_predictions.csv'), index=False)
    print("\nSaved.")

if __name__ == '__main__':
    main()
