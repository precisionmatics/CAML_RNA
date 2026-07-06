"""
Step 08 — Ligand Morgan FP + RNA pocket composition features

Ligand features (from SDF):
  - Morgan ECFP4 (2048-bit)
  - MACCS keys (167-bit)
  - Physicochemical: MW, logP, TPSA, HBD, HBA, rings, rotbonds (7-dim)

RNA pocket features (from dataset_clean.csv):
  - Element counts: C,N,O,P (4-dim)
  - Structural composition: pocket size, purine/pyrimidine fraction (3-dim)

Saves:
  - features/step08_ligand_features.npz  (2222-dim)
  - features/step08_combined.npz         (ES_CS + WCh + ligand = 12482-dim)

Then runs nested CV SVR-RBF on combined and reports result.
"""

import numpy as np
import pandas as pd
import os, warnings
from scipy.stats import pearsonr, spearmanr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error, r2_score
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, Descriptors

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

def compute_ligand_features(sdf_path):
    """Returns 2222-dim vector: Morgan(2048) + MACCS(167) + physchem(7)."""
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=False)
    for mol in supplier:
        if mol is None:
            continue
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            pass
        try:
            morgan = np.array(
                AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048),
                dtype=np.float32)
        except Exception:
            morgan = np.zeros(2048, dtype=np.float32)
        try:
            maccs = np.array(MACCSkeys.GenMACCSKeys(mol), dtype=np.float32)
        except Exception:
            maccs = np.zeros(167, dtype=np.float32)
        try:
            phys = np.array([
                Descriptors.MolWt(mol),
                Descriptors.MolLogP(mol),
                Descriptors.TPSA(mol),
                Descriptors.NumHDonors(mol),
                Descriptors.NumHAcceptors(mol),
                Descriptors.RingCount(mol),
                Descriptors.NumRotatableBonds(mol),
            ], dtype=np.float32)
        except Exception:
            phys = np.zeros(7, dtype=np.float32)
        return np.concatenate([morgan, maccs, phys])
    return np.zeros(2222, dtype=np.float32)

def compute_rna_features(row):
    """7-dim RNA pocket composition from dataset columns."""
    total = max(row['n_rna_atoms'], 1)
    return np.array([
        row['n_rna_atoms'],
        row['n_lig_atoms'],
        row['rna_C'] / total,
        row['rna_N'] / total,
        row['rna_O'] / total,
        row['rna_P'] / total,
        (row['rna_C'] + row['rna_N']) / total,  # base fraction proxy
    ], dtype=np.float32)

def make_pipe():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])

SVR_GRID = {'model__C': [0.1,1,10,100,500], 'model__gamma': ['scale','auto']}

def nested_cv_svr(X, y):
    outer = KFold(5, shuffle=True, random_state=SEED)
    inner = KFold(5, shuffle=True, random_state=SEED)
    oof   = np.zeros(len(y))
    for tr, te in outer.split(X):
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=inner,
                          scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def main():
    df = pd.read_csv(DATA_CSV)
    n  = len(df)

    print(f"Computing ligand + RNA features for {n} complexes...")
    lig_feats = np.zeros((n, 2222), dtype=np.float32)
    rna_feats = np.zeros((n, 7),    dtype=np.float32)

    for i, row in df.iterrows():
        lig_feats[i] = compute_ligand_features(row['ligand_sdf'])
        rna_feats[i] = compute_rna_features(row)
        if (i+1) % 20 == 0:
            print(f"  [{i+1}/{n}]")

    extra = np.hstack([lig_feats, rna_feats])  # 2229-dim
    y = df['pKd'].values.astype(np.float64)

    np.savez(os.path.join(FEAT_DIR, 'step08_ligand_rna.npz'),
             X=extra, y=y, pdbs=df['pdb'].values)
    print(f"Saved ligand+RNA features: {extra.shape}")

    # Load PH features and combine
    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X_ph = np.hstack([X_es, X_wch])           # 10260
    X_all = np.hstack([X_ph, extra.astype(np.float64)])  # 12489

    print(f"\nRunning nested CV on combined features: {X_all.shape}")

    # Ligand-only baseline
    print("\nLigand-only SVR-RBF:")
    oof_lig = nested_cv_svr(extra.astype(np.float64), y)
    r, rho, rmse, r2 = metrics(y, oof_lig)
    lo, hi = bootstrap_ci(y, oof_lig)
    print(f"  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  95%CI [{lo:.3f},{hi:.3f}]")

    # PH + ligand combined
    print("\nPH (ES_CS+WCh) + Ligand+RNA SVR-RBF:")
    oof_comb = nested_cv_svr(X_all, y)
    r, rho, rmse, r2 = metrics(y, oof_comb)
    lo, hi = bootstrap_ci(y, oof_comb)
    print(f"  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}  "
          f"95%CI [{lo:.3f},{hi:.3f}]")

    np.save(os.path.join(MODEL_DIR, 'oof_pred_step08_combined.npy'), oof_comb)
    np.save(os.path.join(MODEL_DIR, 'oof_pred_step08_ligand.npy'),   oof_lig)

    print(f"\n{'='*55}")
    print(f"Step 08 best (PH+Ligand+RNA): r={r:.4f}")
    print(f"{'='*55}")

if __name__ == '__main__':
    main()
