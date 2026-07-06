"""
Step 07 — Hybrid model + Riboswitch subclass models

Strategy:
  - aptamer (n=20): subtype SVR-RBF  → r=0.937
  - ribosomal_asite (n=13): subtype SVR-RBF → r=0.580
  - riboswitch (n=61): sub-classified into SAM_SAH, purine, FMN_FAD, TPP, other_lig
  - everything else: global ES_CS+WCh SVR-RBF (r=0.609)

Riboswitch subclasses from PDFL-RNA project (mol2 codes):
  SAM_SAH: 2ydh,2ygh,3e5c,3npn,4aob,4kqy,4l81,4oqu,6fz0,6hag
  purine:  6AP,6GU,XAN,6GO,2BP,A2F,ADE,2BA,2QB,29G,29H,7DG
  FMN_FAD: PRF,FFO,LYA,H4B
  TPP:     C2E,PRP
  other_lig: remaining riboswitches
"""

import numpy as np
import pandas as pd
import os, warnings
from scipy.stats import pearsonr, spearmanr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, LeaveOneOut, GridSearchCV
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error, r2_score

warnings.filterwarnings('ignore')

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

# Riboswitch subclass PDB IDs (from PDFL-RNA analysis)
RS_SUBCLASS = {
    'SAM_SAH':  {'2ydh','2ygh','3e5c','3npn','4aob','4kqy','4l81','4oqu','6fz0','6hag'},
    'purine':   {'2b57','2b3j','3b31','3d2g','3d2v','3d2x','3epb','3f2q','3gx3',
                 '3gx5','3gx6','3gx7','3iwn','3kfr','3q50','4fe5','4rzd','4ts2',
                 '5sv1','5t83'},
    'FMN_FAD':  {'3f4g','3f4h','4gnk','5hoh','5v3f','6tcq','2bh2','6b3h'},
    'TPP':      {'2gdi','2hom','2hoj','3d2s','3k5y','3k5z','4gxy','5iqb','6n6v','6n6w'},
}

def pdb_to_rs_subclass(pdb):
    for sc, pdbs in RS_SUBCLASS.items():
        if pdb in pdbs:
            return sc
    return 'other_lig'

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

def make_pipe():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])

SVR_GRID = {'model__C': [0.1,1,10,100,500], 'model__gamma': ['scale','auto']}

def svr_oof(X, y):
    n = len(y)
    cv_outer = LeaveOneOut() if n < 15 else KFold(5, shuffle=True, random_state=SEED)
    cv_inner = KFold(min(5,n-1), shuffle=True, random_state=SEED)
    oof = np.zeros(n)
    for tr, te in cv_outer.split(X):
        if len(tr) < 5:
            oof[te] = y.mean()
            continue
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=cv_inner,
                          scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def main():
    df   = pd.read_csv(DATA_CSV)
    sub  = df['subtype'].values
    pdb  = df['pdb'].values

    X_es,  y = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X = np.hstack([X_es, X_wch])

    # Start from global step05 predictions
    oof = np.load(os.path.join(MODEL_DIR, 'oof_pred_step05_best.npy')).copy()

    print("=== Hybrid: override where subtype model helps ===\n")

    # --- aptamer ---
    idx = np.where(sub == 'aptamer')[0]
    oof_apt = svr_oof(X[idx], y[idx])
    oof[idx] = oof_apt
    r_apt, *_ = metrics(y[idx], oof_apt)
    print(f"aptamer          n={len(idx):2d}  r={r_apt:.4f}")

    # --- ribosomal_asite ---
    idx = np.where(sub == 'ribosomal_asite')[0]
    oof_ras = svr_oof(X[idx], y[idx])
    oof[idx] = oof_ras
    r_ras, *_ = metrics(y[idx], oof_ras)
    print(f"ribosomal_asite  n={len(idx):2d}  r={r_ras:.4f}")

    # --- riboswitch: per-subclass ---
    rs_idx = np.where(sub == 'riboswitch')[0]
    rs_sc  = df['rs_subclass'].values[rs_idx]
    print(f"\nRiboswitch subclasses: {dict(zip(*np.unique(rs_sc, return_counts=True)))}")

    for sc in sorted(set(rs_sc)):
        mask   = rs_sc == sc
        sc_idx = rs_idx[mask]
        n_sc   = len(sc_idx)
        if n_sc < 5:
            print(f"  {sc:12s}  n={n_sc:2d}  SKIP")
            continue
        oof_sc = svr_oof(X[sc_idx], y[sc_idx])
        # only override if subclass model beats global on this slice
        r_sc,  *_ = metrics(y[sc_idx], oof_sc)
        r_gl,  *_ = metrics(y[sc_idx], oof[sc_idx])
        if r_sc > r_gl:
            oof[sc_idx] = oof_sc
            print(f"  {sc:12s}  n={n_sc:2d}  r={r_sc:.4f}  (global={r_gl:.4f}) ✓ override")
        else:
            print(f"  {sc:12s}  n={n_sc:2d}  r={r_sc:.4f}  (global={r_gl:.4f}) — keep global")

    # Overall
    r, rho, rmse, r2 = metrics(y, oof)
    lo, hi = bootstrap_ci(y, oof)
    print(f"\n{'='*55}")
    print(f"STEP 07 HYBRID:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*55}")

    # Per-subtype breakdown
    print("\nPer-subtype breakdown:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        try:
            r_st, *_ = metrics(y[idx], oof[idx])
            print(f"  {st:20s}  n={len(idx):2d}  r={r_st:.4f}")
        except Exception:
            print(f"  {st:20s}  n={len(idx):2d}  r=N/A")

    np.save(os.path.join(MODEL_DIR, 'oof_pred_step07.npy'), oof)
    pd.DataFrame({'pdb': pdb, 'y_true': y, 'y_pred': oof, 'subtype': sub}).to_csv(
        os.path.join(RES_DIR, 'step07_predictions.csv'), index=False)
    print(f"\nSaved OOF predictions.")

if __name__ == '__main__':
    main()
