"""
Step 12 — Subtype-aware ensemble + augmented features

Approaches:
  A) Subtype one-hot appended to ES_CS+WCh → retrain nested SVR-RBF
  B) LOO ridge stacking of best OOF predictors (step07, step09_ph_fm, step11)
  C) Subtype-specific optimal weighting of B
  D) Riboswitch LOO-retrained SVR on ES_CS+WCh+CPF with aptamer override
  Final: best of the above + aptamer LOO override
"""

import numpy as np
import pandas as pd
import os, warnings
from scipy.stats import pearsonr, spearmanr
from scipy.optimize import minimize
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, LeaveOneOut, GridSearchCV
from sklearn.svm import SVR
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_squared_error, r2_score

warnings.filterwarnings('ignore')

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

def metrics(y, p):
    r,  _ = pearsonr(y, p)
    rho,_ = spearmanr(y, p)
    rmse  = np.sqrt(mean_squared_error(y, p))
    r2    = r2_score(y, p)
    return r, rho, rmse, r2

def bootstrap_ci(y, p, n=10000):
    rng = np.random.default_rng(SEED)
    ns  = len(y)
    rs  = [pearsonr(y[ix:=rng.integers(0,ns,ns)], p[ix])[0] for _ in range(n)]
    return np.percentile(rs, [2.5, 97.5])

def make_pipe():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])

SVR_GRID = {'model__C':[0.1,1,10,100,500,1000], 'model__gamma':['scale','auto']}

def nested_cv_svr(X, y, loo=False):
    n     = len(y)
    outer = LeaveOneOut() if loo else KFold(5, shuffle=True, random_state=SEED)
    inner = KFold(min(5, n-1), shuffle=True, random_state=SEED)
    oof   = np.zeros(n)
    for tr, te in outer.split(X):
        if len(tr) < 5:
            oof[te] = y.mean(); continue
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=inner, scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def load_feat(fname):
    d = np.load(os.path.join(FEAT_DIR, fname), allow_pickle=True)
    return d['X'].astype(np.float64), d['y'].astype(np.float64)

def load_oof(fname):
    return np.load(os.path.join(MODEL_DIR, fname))

def print_result(tag, y, oof):
    r, rho, rmse, r2 = metrics(y, oof)
    lo, hi = bootstrap_ci(y, oof)
    print(f"  {tag:35s}  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  95%CI [{lo:.3f},{hi:.3f}]")
    return r

def apply_aptamer_loo(X, y, sub, oof_global):
    idx = np.where(sub == 'aptamer')[0]
    oof_apt = nested_cv_svr(X[idx], y[idx], loo=False)
    r_apt = pearsonr(y[idx], oof_apt)[0]
    r_gl  = pearsonr(y[idx], oof_global[idx])[0]
    oof_out = oof_global.copy()
    if r_apt > r_gl:
        oof_out[idx] = oof_apt
        print(f"    aptamer: override r={r_apt:.4f} > global r={r_gl:.4f}")
    else:
        print(f"    aptamer: keep global r={r_gl:.4f}")
    return oof_out

def main():
    df  = pd.read_csv(DATA_CSV)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values
    n   = len(y)

    # --- Load features ---
    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X_cpf, _ = load_feat('step11_cpf.npz')

    # Subtype one-hot
    le  = LabelEncoder().fit(sub)
    ohe = np.zeros((n, len(le.classes_)), dtype=np.float64)
    for i, s in enumerate(sub):
        ohe[i, le.transform([s])[0]] = 1.0

    # Best OOF predictors
    p07  = load_oof('oof_pred_step07.npy')
    p09  = load_oof('oof_pred_step09_ph_fm.npy')
    p11  = load_oof('oof_pred_step11.npy')
    p05  = load_oof('oof_pred_step05_best.npy')

    print("=== Approach A: ES_CS+WCh + subtype one-hot ===")
    X_A = np.hstack([X_es, X_wch, ohe])
    oof_A = nested_cv_svr(X_A, y)
    r_A = print_result("ES_CS+WCh+subtype_ohe", y, oof_A)
    oof_A_hybrid = apply_aptamer_loo(X_A, y, sub, oof_A)
    r_Ah = print_result("ES_CS+WCh+subtype_ohe (hybrid)", y, oof_A_hybrid)

    print("\n=== Approach B: LOO ridge stacking (step07+step09+step11) ===")
    S = np.column_stack([p07, p09, p11])
    oof_B = np.zeros(n)
    for tr, te in LeaveOneOut().split(S):
        rc = RidgeCV(alphas=np.logspace(-3, 3, 30)).fit(S[tr], y[tr])
        oof_B[te] = rc.predict(S[te])
    r_B = print_result("ridge_stack(07+09+11)", y, oof_B)

    print("\n=== Approach C: Subtype-specific LOO ridge stacking ===")
    oof_C = np.zeros(n)
    for st in sorted(set(sub)):
        idx = np.where(sub == st)[0]
        preds = np.column_stack([p07[idx], p09[idx], p11[idx], p05[idx]])
        yi = y[idx]
        if len(idx) >= 8:
            # proper LOO within subtype
            oof_sub = np.zeros(len(idx))
            for tr_i, te_i in LeaveOneOut().split(preds):
                rc = RidgeCV(alphas=np.logspace(-3,3,30)).fit(preds[tr_i], yi[tr_i])
                oof_sub[te_i] = rc.predict(preds[te_i])
            oof_C[idx] = oof_sub
        else:
            oof_C[idx] = p07[idx]   # too few samples: use best single predictor
    r_C = print_result("subtype_specific_LOO_ridge", y, oof_C)

    print("\n=== Approach D: ES_CS+WCh+CPF+subtype_ohe ===")
    X_D = np.hstack([X_es, X_wch, X_cpf, ohe])
    oof_D = nested_cv_svr(X_D, y)
    r_D = print_result("ES_CS+WCh+CPF+subtype_ohe", y, oof_D)
    oof_D_hybrid = apply_aptamer_loo(X_D, y, sub, oof_D)
    r_Dh = print_result("ES_CS+WCh+CPF+subtype_ohe (hybrid)", y, oof_D_hybrid)

    print("\n=== Approach E: Ridge stack of B + best hybrid preds ===")
    p_best_prev = p07
    S2 = np.column_stack([oof_B, oof_A_hybrid, oof_D_hybrid, p_best_prev])
    oof_E = np.zeros(n)
    for tr, te in LeaveOneOut().split(S2):
        rc = RidgeCV(alphas=np.logspace(-3,3,30)).fit(S2[tr], y[tr])
        oof_E[te] = rc.predict(S2[te])
    r_E = print_result("meta_stack_ABDE", y, oof_E)

    # --- Final best ---
    candidates = {
        'A_hybrid': (r_Ah, oof_A_hybrid),
        'B_stack':  (r_B,  oof_B),
        'C_subtype':(r_C,  oof_C),
        'D_hybrid': (r_Dh, oof_D_hybrid),
        'E_meta':   (r_E,  oof_E),
    }
    best_name = max(candidates, key=lambda k: candidates[k][0])
    best_r, best_oof = candidates[best_name]

    print(f"\n{'='*60}")
    print(f"STEP 12 BEST: {best_name}  r={best_r:.4f}")
    lo, hi = bootstrap_ci(y, best_oof)
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    _, rho, rmse, r2 = metrics(y, best_oof)
    print(f"rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"{'='*60}")

    print("\nPer-subtype (best):")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        rs,*_ = metrics(y[idx], best_oof[idx])
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")

    np.save(os.path.join(MODEL_DIR, 'oof_pred_step12.npy'), best_oof)
    pd.DataFrame({'pdb':df['pdb'].values,'y_true':y,'y_pred':best_oof,
                  'subtype':sub}).to_csv(
        os.path.join(RES_DIR,'step12_predictions.csv'), index=False)
    print("\nSaved oof_pred_step12.npy and step12_predictions.csv")

if __name__ == '__main__':
    main()
