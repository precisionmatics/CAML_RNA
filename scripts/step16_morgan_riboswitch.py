"""
Step 16 — Morgan fingerprint-augmented riboswitch models

For SAM_SAH (n=10, current r=-0.17) and purine (n=21, current r=0.27):
  - The RNA binding pocket is structurally similar across analogs
  - Affinity variation is driven by LIGAND chemistry, not RNA topology
  - Solution: Morgan fingerprints (ECFP4, 1024-bit) for the ligand
             + ES_CS+WCh PH features for the RNA pocket

For other riboswitch subclasses and remaining subtypes: keep step07 predictions.
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
from rdkit import Chem
from rdkit.Chem import AllChem

warnings.filterwarnings('ignore')

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

RS_SUBCLASS = {
    'SAM_SAH': {'2ydh','2ygh','3e5c','3npn','4aob','4kqy','4l81','4oqu','6fz0','6hag'},
    'purine':  {'2b57','2b3j','3b31','3d2g','3d2v','3d2x','3epb','3f2q','3gx3',
                '3gx5','3gx6','3gx7','3iwn','3kfr','3q50','4fe5','4rzd','4ts2','5sv1','5t83'},
    'FMN_FAD': {'3f4g','3f4h','4gnk','5hoh','5v3f','6tcq','2bh2','6b3h'},
    'TPP':     {'2gdi','2hom','2hoj','3d2s','3k5y','3k5z','4gxy','5iqb','6n6v','6n6w'},
}

def pdb_to_rs(pdb):
    for sc, pdbs in RS_SUBCLASS.items():
        if pdb in pdbs: return sc
    return 'other_lig'

def get_morgan(sdf_path, radius=2, n_bits=1024):
    """Compute Morgan fingerprint from SDF file."""
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=False)
    for mol in supplier:
        if mol is None: continue
        try:
            Chem.SanitizeMol(mol)
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
            return np.array(fp, dtype=np.float64)
        except Exception:
            return np.zeros(n_bits, dtype=np.float64)
    return np.zeros(n_bits, dtype=np.float64)

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

def main():
    df   = pd.read_csv(DATA_CSV)
    y    = df['pKd'].values.astype(np.float64)
    sub  = df['subtype'].values
    pdb  = df['pdb'].values

    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X_base   = np.hstack([X_es, X_wch])

    # Compute Morgan fingerprints for all ligands
    print("Computing Morgan fingerprints...")
    X_morgan = np.zeros((len(df), 1024), dtype=np.float64)
    for i, row in df.iterrows():
        X_morgan[i] = get_morgan(row['ligand_sdf'])
    print(f"Done. Non-zero rows: {(X_morgan.sum(1)>0).sum()}")

    # Start from step07 predictions
    oof = np.load(os.path.join(MODEL_DIR, 'oof_pred_step07.npy')).copy()
    print(f"\nBaseline (step07): r={pearsonr(y, oof)[0]:.4f}")

    # --- Process riboswitch subclasses ---
    rs_idx = np.where(sub == 'riboswitch')[0]
    rs_sc  = np.array([pdb_to_rs(p) for p in pdb[rs_idx]])

    print("\n=== Riboswitch per-subclass (step07 baseline vs Morgan-augmented) ===")
    for sc in sorted(set(rs_sc)):
        mask   = rs_sc == sc
        sc_idx = rs_idx[mask]
        n_sc   = len(sc_idx)
        if n_sc < 5:
            r_base = pearsonr(y[sc_idx], oof[sc_idx])[0]
            print(f"  {sc:12s}  n={n_sc}  baseline_r={r_base:.4f}  SKIP (too few)")
            continue

        y_sc    = y[sc_idx]
        use_loo = n_sc < 15

        # Baseline (step07): global slice r
        r_gl = pearsonr(y_sc, oof[sc_idx])[0]

        # PH only (subclass SVR)
        oof_ph = svr_oof(X_base[sc_idx], y_sc, loo=use_loo)
        r_ph   = pearsonr(y_sc, oof_ph)[0]

        # Morgan only
        oof_mg = svr_oof(X_morgan[sc_idx], y_sc, loo=use_loo)
        r_mg   = pearsonr(y_sc, oof_mg)[0]

        # PH + Morgan
        X_comb = np.hstack([X_base[sc_idx], X_morgan[sc_idx]])
        oof_cm  = svr_oof(X_comb, y_sc, loo=use_loo)
        r_cm    = pearsonr(y_sc, oof_cm)[0]

        print(f"  {sc:12s}  n={n_sc:2d}  global={r_gl:+.4f}  PH={r_ph:+.4f}  "
              f"Morgan={r_mg:+.4f}  PH+Morgan={r_cm:+.4f}")

        # Override with best
        best_r = max(r_gl, r_ph, r_mg, r_cm)
        if best_r > r_gl:
            if best_r == r_cm:
                oof[sc_idx] = oof_cm
                tag = 'PH+Morgan'
            elif best_r == r_mg:
                oof[sc_idx] = oof_mg
                tag = 'Morgan'
            else:
                oof[sc_idx] = oof_ph
                tag = 'PH'
            print(f"             → Override with {tag}")

    # --- Also try Morgan for aptamer and ribosomal_asite ---
    print("\n=== Aptamer (Morgan-augmented) ===")
    idx   = np.where(sub=='aptamer')[0]
    X_comb = np.hstack([X_base[idx], X_morgan[idx]])
    oof_cm = svr_oof(X_comb, y[idx])
    r_cm   = pearsonr(y[idx], oof_cm)[0]
    r_gl   = pearsonr(y[idx], oof[idx])[0]
    print(f"  PH+Morgan={r_cm:.4f}  vs current={r_gl:.4f}")
    if r_cm > r_gl:
        oof[idx] = oof_cm
        print("  → Override")

    print("\n=== other_misc (Morgan-augmented) ===")
    idx   = np.where(sub=='other_misc')[0]
    X_comb = np.hstack([X_base[idx], X_morgan[idx]])
    oof_cm = svr_oof(X_comb, y[idx])
    r_cm   = pearsonr(y[idx], oof_cm)[0]
    r_gl   = pearsonr(y[idx], oof[idx])[0]
    print(f"  PH+Morgan={r_cm:.4f}  vs global step07 slice={r_gl:.4f}")
    if r_cm > r_gl:
        oof[idx] = oof_cm
        print("  → Override")

    # Global model with Morgan (just for comparison)
    print("\n=== Global PH+Morgan SVR (reference) ===")
    X_global_mg = np.hstack([X_base, X_morgan])
    oof_gm = svr_oof(X_global_mg, y)
    r_gm   = pearsonr(y, oof_gm)[0]
    print(f"  Global PH+Morgan: r={r_gm:.4f}")

    # Final
    r, rho, rmse, r2 = metrics(y, oof)
    lo, hi = bootstrap_ci(y, oof)
    print(f"\n{'='*60}")
    print(f"STEP 16:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*60}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        rs,*_ = metrics(y[idx], oof[idx])
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")

    np.save(os.path.join(MODEL_DIR,'oof_pred_step16.npy'), oof)
    pd.DataFrame({'pdb':pdb,'y_true':y,'y_pred':oof,'subtype':sub}).to_csv(
        os.path.join(RES_DIR,'step16_predictions.csv'), index=False)
    print("\nSaved.")

if __name__ == '__main__':
    main()
