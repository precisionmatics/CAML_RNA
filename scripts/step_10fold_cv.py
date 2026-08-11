"""
10-Fold Cross-Validation — fair comparison with literature methods

Reviewer concern: CAML-RNA uses LOO-CV while competing methods use 10-fold CV.
This script evaluates every feature combination under nested 10-fold CV
(outer for OOF predictions, inner 5-fold for hyperparameter tuning)
so results are directly comparable with RLASIF, RLaffinity, AffiGrapher etc.

Also reports LOO-CV for each combination so both protocols are documented.
"""

import numpy as np
import pandas as pd
import os, warnings, time
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

SVR_GRID = {'model__C': [0.1, 1, 10, 100, 500],
            'model__gamma': ['scale', 'auto']}

def make_pipe():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])

def metrics(y, p):
    r,   _ = pearsonr(y, p)
    rho, _ = spearmanr(y, p)
    rmse   = np.sqrt(mean_squared_error(y, p))
    r2     = r2_score(y, p)
    return dict(r=r, rho=rho, rmse=rmse, r2=r2)

def bootstrap_ci(y, p, n=10000):
    rng = np.random.default_rng(SEED)
    ns  = len(y)
    rs  = [pearsonr(y[i := rng.integers(0, ns, ns)], p[i])[0] for _ in range(n)]
    return np.percentile(rs, 2.5), np.percentile(rs, 97.5)

def nested_cv(X, y, n_outer=10):
    outer = KFold(n_splits=n_outer, shuffle=True, random_state=SEED)
    inner = KFold(n_splits=5, shuffle=True, random_state=SEED)
    oof   = np.zeros(len(y))
    for tr, te in outer.split(X):
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=inner,
                          scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def nested_loo(X, y):
    loo   = LeaveOneOut()
    inner = KFold(n_splits=5, shuffle=True, random_state=SEED)
    oof   = np.zeros(len(y))
    for tr, te in loo.split(X):
        if len(tr) < 5:
            oof[te] = y.mean(); continue
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=inner,
                          scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def load_npz(fname, key='X'):
    d = np.load(os.path.join(FEAT_DIR, fname), allow_pickle=True)
    return d[key].astype(np.float64)

def main():
    df  = pd.read_csv(DATA_CSV)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values

    # Load features (use new step02 outputs where available)
    d02 = np.load(os.path.join(FEAT_DIR, 'step02_es_cs_features.npz'), allow_pickle=True)
    X_ph   = d02['X'].astype(np.float64)          # PH only (bipartite corrected)
    X_psrt = d02['X_psrt'].astype(np.float64)     # PSRT f/h-vectors
    X_full = d02['X_full'].astype(np.float64)     # PH + PSRT combined

    X_cpf = load_npz('step11_cpf.npz')
    X_fm  = load_npz('step09_rnafm.npz')
    X_wch = load_npz('step03_wch.npz')

    # Morgan fingerprints + physicochemical (from step33)
    from rdkit import Chem
    from rdkit.Chem import AllChem
    def safe_morgan(row):
        for col in ['ligand_sdf', 'lig_sdf']:
            p = str(row.get(col, ''))
            if p and os.path.exists(p):
                suppl = Chem.SDMolSupplier(p, removeHs=True, sanitize=False)
                for m in suppl:
                    if m:
                        try:
                            m.UpdatePropertyCache(strict=False)
                            return np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024))
                        except:
                            pass
        return np.zeros(1024)

    X_morgan = np.array([safe_morgan(row) for _, row in df.iterrows()])

    lig_cols = ['mol_weight','n_rings','n_hbd','n_hba','n_rot_bonds','tpsa',
                'lig_C','lig_N','lig_O','lig_S','n_lig_atoms']
    rna_cols = ['n_rna_atoms','rna_C','rna_N','rna_O','rna_P','rna_S']
    X_phys = np.hstack([df[lig_cols].fillna(0).values,
                        df[rna_cols].fillna(0).values]).astype(np.float64)

    combos = {
        # Individual modalities
        'PH_bipartite':        X_ph,
        'PSRT_fvec':           X_psrt,
        'PH+PSRT':             X_full,
        'CPF':                 X_cpf,
        'RNA-FM':              X_fm,
        'Morgan':              X_morgan,
        'Physicochemical':     X_phys,
        'WCh':                 X_wch,
        # Key combinations
        'PH+CPF':              np.hstack([X_ph, X_cpf]),
        'PH+PSRT+CPF':         np.hstack([X_full, X_cpf]),
        'PH+PSRT+CPF+FM':      np.hstack([X_full, X_cpf, X_fm]),
        'PH+PSRT+CPF+FM+Phys': np.hstack([X_full, X_cpf, X_fm, X_phys]),
        'Full_ensemble':       np.hstack([X_full, X_cpf, X_fm, X_morgan, X_phys]),
    }

    rows = []
    print(f"{'Combination':<30} {'10fold-R':>9} {'10fold-RMSE':>12} {'LOO-R':>8} {'95%CI':>16}")
    print("=" * 80)

    for name, X in combos.items():
        t0 = time.time()
        print(f"  {name:<28} ... ", end='', flush=True)

        oof_10 = nested_cv(X, y, n_outer=10)
        m10    = metrics(y, oof_10)
        lo, hi = bootstrap_ci(y, oof_10)

        oof_loo = nested_loo(X, y)
        mloo    = metrics(y, oof_loo)

        elapsed = time.time() - t0
        print(f"R={m10['r']:.4f}  RMSE={m10['rmse']:.4f}  "
              f"LOO-R={mloo['r']:.4f}  [{lo:.3f},{hi:.3f}]  {elapsed:.0f}s")

        row = {'combination': name,
               'r_10fold':   m10['r'],   'rho_10fold': m10['rho'],
               'rmse_10fold': m10['rmse'], 'r2_10fold': m10['r2'],
               'ci_lo_10fold': lo, 'ci_hi_10fold': hi,
               'r_loo':  mloo['r'],  'rho_loo':  mloo['rho'],
               'rmse_loo': mloo['rmse'], 'r2_loo': mloo['r2'],
               'n': len(y)}

        # Per-subtype breakdown
        for st in sorted(set(sub)):
            idx = np.where(sub == st)[0]
            if len(idx) >= 3:
                row[f'r_10fold_{st}'] = pearsonr(y[idx], oof_10[idx])[0]
                row[f'r_loo_{st}']    = pearsonr(y[idx], oof_loo[idx])[0]

        rows.append(row)

        np.save(os.path.join(MODEL_DIR, f"oof_10fold_{name.replace(' ','_')}.npy"),
                oof_10)
        np.save(os.path.join(MODEL_DIR, f"oof_loo_{name.replace(' ','_')}.npy"),
                oof_loo)

    out_csv = os.path.join(RES_DIR, 'step_10fold_cv_results.csv')
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")

    # Print summary table
    print("\n" + "=" * 80)
    print("SUMMARY — 10-fold CV vs LOO-CV")
    print("=" * 80)
    for row in sorted(rows, key=lambda x: x['r_10fold'], reverse=True):
        print(f"  {row['combination']:<30}  "
              f"10fold={row['r_10fold']:.4f}  LOO={row['r_loo']:.4f}")

if __name__ == '__main__':
    main()
