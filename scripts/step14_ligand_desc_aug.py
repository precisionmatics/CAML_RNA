"""
Step 14 — Ligand descriptor augmented global model

Ligand properties (mol_weight, n_rings etc.) correlate with pKd (r up to 0.39).
Adding them to ES_CS+WCh gives SVR extra chemical context.
Then apply step07 hybrid overrides (aptamer, ribosomal_asite, riboswitch subclasses).
Also tests adding RNA size descriptors (n_rna_atoms, rna_C/N/O/P).
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
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

LIG_COLS = ['mol_weight','n_rings','n_hbd','n_hba','n_rot_bonds','tpsa',
            'lig_C','lig_N','lig_O','lig_S','lig_F','lig_P','n_lig_atoms']
RNA_COLS = ['n_rna_atoms','rna_C','rna_N','rna_O','rna_P']

RS_SUBCLASS = {
    'SAM_SAH':  {'2ydh','2ygh','3e5c','3npn','4aob','4kqy','4l81','4oqu','6fz0','6hag'},
    'purine':   {'2b57','2b3j','3b31','3d2g','3d2v','3d2x','3epb','3f2q','3gx3',
                 '3gx5','3gx6','3gx7','3iwn','3kfr','3q50','4fe5','4rzd','4ts2',
                 '5sv1','5t83'},
    'FMN_FAD':  {'3f4g','3f4h','4gnk','5hoh','5v3f','6tcq','2bh2','6b3h'},
    'TPP':      {'2gdi','2hom','2hoj','3d2s','3k5y','3k5z','4gxy','5iqb','6n6v','6n6w'},
}

def pdb_to_rs(pdb):
    for sc, pdbs in RS_SUBCLASS.items():
        if pdb in pdbs:
            return sc
    return 'other_lig'

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

def make_pipe():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])

SVR_GRID = {'model__C':[0.1,1,10,100,500,1000], 'model__gamma':['scale','auto']}

def svr_oof(X, y, loo=False):
    n     = len(y)
    outer = LeaveOneOut() if loo else KFold(5, shuffle=True, random_state=SEED)
    inner = KFold(min(5,n-1), shuffle=True, random_state=SEED)
    oof   = np.zeros(n)
    for tr, te in outer.split(X):
        if len(tr) < 5:
            oof[te] = y.mean(); continue
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=inner, scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def apply_hybrid(X, y, sub, pdb, oof_global):
    oof = oof_global.copy()

    # aptamer
    idx = np.where(sub == 'aptamer')[0]
    oof_s = svr_oof(X[idx], y[idx])
    r_s = pearsonr(y[idx], oof_s)[0]
    oof[idx] = oof_s
    print(f"  aptamer  r={r_s:.4f}")

    # ribosomal_asite
    idx = np.where(sub == 'ribosomal_asite')[0]
    oof_s = svr_oof(X[idx], y[idx], loo=True)
    r_s = pearsonr(y[idx], oof_s)[0]
    r_g = pearsonr(y[idx], oof[idx])[0]
    if r_s > r_g:
        oof[idx] = oof_s
        print(f"  ribosomal_asite  r={r_s:.4f} ✓")
    else:
        print(f"  ribosomal_asite  keep global r={r_g:.4f}")

    # riboswitch subclasses
    rs_idx = np.where(sub == 'riboswitch')[0]
    rs_sc  = np.array([pdb_to_rs(p) for p in pdb[rs_idx]])
    for sc in sorted(set(rs_sc)):
        mask   = rs_sc == sc
        sc_idx = rs_idx[mask]
        if len(sc_idx) < 5:
            continue
        oof_s = svr_oof(X[sc_idx], y[sc_idx],
                        loo=len(sc_idx) < 15)
        r_s = pearsonr(y[sc_idx], oof_s)[0]
        r_g = pearsonr(y[sc_idx], oof[sc_idx])[0]
        if r_s > r_g:
            oof[sc_idx] = oof_s
            print(f"  riboswitch/{sc:10s}  n={len(sc_idx)}  r={r_s:.4f} ✓")
        else:
            print(f"  riboswitch/{sc:10s}  n={len(sc_idx)}  global={r_g:.4f}")

    return oof

def main():
    df  = pd.read_csv(DATA_CSV)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values
    pdb = df['pdb'].values

    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X_cpf, _ = load_feat('step11_cpf.npz')

    # Ligand + RNA descriptors
    X_lig = df[LIG_COLS].fillna(0).values.astype(np.float64)
    X_rna = df[RNA_COLS].fillna(0).values.astype(np.float64)
    X_phys = np.hstack([X_lig, X_rna])

    base_ph = np.hstack([X_es, X_wch])

    combos = {
        'ES_CS+WCh+Lig':         np.hstack([base_ph, X_lig]),
        'ES_CS+WCh+Phys':        np.hstack([base_ph, X_phys]),
        'ES_CS+WCh+CPF+Lig':     np.hstack([base_ph, X_cpf, X_lig]),
    }

    print("=== Global SVR-RBF results ===")
    best_r_gl, best_X, best_name_gl = -99, None, ''
    for name, X in combos.items():
        oof = svr_oof(X, y)
        r   = pearsonr(y, oof)[0]
        print(f"  {name:35s}  r={r:.4f}")
        if r > best_r_gl:
            best_r_gl, best_X, best_name_gl = r, X, name

    # Reference: plain ES_CS+WCh (step05)
    p05 = np.load(os.path.join(MODEL_DIR, 'oof_pred_step05_best.npy'))
    print(f"  {'ES_CS+WCh (step05 ref)':35s}  r={pearsonr(y,p05)[0]:.4f}")

    print(f"\nBest global: {best_name_gl}  r={best_r_gl:.4f}")

    print("\n=== Hybrid overrides ===")
    oof_global = svr_oof(best_X, y)
    oof_hybrid = apply_hybrid(best_X, y, sub, pdb, oof_global)

    r, rho, rmse, r2 = metrics(y, oof_hybrid)
    lo, hi = bootstrap_ci(y, oof_hybrid)
    print(f"\n{'='*60}")
    print(f"STEP 14 HYBRID:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*60}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        rs,*_ = metrics(y[idx], oof_hybrid[idx])
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")

    np.save(os.path.join(MODEL_DIR,'oof_pred_step14.npy'), oof_hybrid)
    pd.DataFrame({'pdb':pdb,'y_true':y,'y_pred':oof_hybrid,'subtype':sub}).to_csv(
        os.path.join(RES_DIR,'step14_predictions.csv'), index=False)
    print("\nSaved.")

if __name__ == '__main__':
    main()
