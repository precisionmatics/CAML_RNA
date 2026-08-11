"""
Full Ablation Study — contribution of each feature modality

Addresses Reviewer 1 (Comment 3) and Reviewer 3 (Comment 2):
  - Ablation of PH vs PSRT vs each modality (per reviewer 1)
  - Ablation of PH alone vs PSRT alone vs PH+PSRT (per reviewer 3)
  - All combinations evaluated under both 10-fold CV and LOO-CV
  - Per-subtype performance for each combination

Modalities:
  PH    = bipartite persistent homology (corrected step02 β₀+β₁)
  PSRT  = persistent f/h-vector (step02 PSRT features)
  CPF   = contact pair fingerprint (step11, 660-dim)
  FM    = RNA-FM embeddings (step09, 640-dim)
  MFP   = Morgan fingerprints (ECFP4, 1024-dim)
  PHY   = physicochemical (17-dim)
  WCH   = weighted chromatic PH (step03, 5940-dim)
"""

import numpy as np
import pandas as pd
import os, warnings, time
from itertools import chain, combinations
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
    return dict(r=round(r,4), rho=round(rho,4),
                rmse=round(rmse,4), r2=round(r2,4))

def bootstrap_ci(y, p, n=10000):
    rng = np.random.default_rng(SEED)
    ns  = len(y)
    rs  = [pearsonr(y[i := rng.integers(0, ns, ns)], p[i])[0] for _ in range(n)]
    return round(np.percentile(rs, 2.5), 3), round(np.percentile(rs, 97.5), 3)

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

def safe_morgan(row):
    for col in ['ligand_sdf', 'lig_sdf']:
        p = str(row.get(col, ''))
        if p and os.path.exists(p):
            suppl = Chem.SDMolSupplier(p, removeHs=True, sanitize=False)
            for m in suppl:
                if m:
                    try:
                        m.UpdatePropertyCache(strict=False)
                        return np.array(
                            AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024))
                    except Exception:
                        pass
    return np.zeros(1024)

def powerset_nonempty(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

def main():
    df  = pd.read_csv(DATA_CSV)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values
    n   = len(y)

    print("Loading features...")
    d02    = np.load(os.path.join(FEAT_DIR, 'step02_es_cs_features.npz'), allow_pickle=True)
    X_ph   = d02['X'].astype(np.float64)
    X_psrt = d02['X_psrt'].astype(np.float64)

    X_cpf  = np.load(os.path.join(FEAT_DIR, 'step11_cpf.npz'),   allow_pickle=True)['X'].astype(np.float64)
    X_fm   = np.load(os.path.join(FEAT_DIR, 'step09_rnafm.npz'), allow_pickle=True)['X'].astype(np.float64)
    X_wch  = np.load(os.path.join(FEAT_DIR, 'step03_wch.npz'),   allow_pickle=True)['X'].astype(np.float64)

    X_morgan = np.array([safe_morgan(row) for _, row in df.iterrows()])

    lig_cols = ['mol_weight','n_rings','n_hbd','n_hba','n_rot_bonds','tpsa',
                'lig_C','lig_N','lig_O','lig_S','n_lig_atoms']
    rna_cols = ['n_rna_atoms','rna_C','rna_N','rna_O','rna_P','rna_S']
    X_phys = np.hstack([df[lig_cols].fillna(0).values,
                        df[rna_cols].fillna(0).values]).astype(np.float64)

    modalities = {
        'PH':   X_ph,
        'PSRT': X_psrt,
        'CPF':  X_cpf,
        'FM':   X_fm,
        'MFP':  X_morgan,
        'PHY':  X_phys,
    }

    subtypes = sorted(set(sub))
    rows = []

    # Priority ablations first (these are the ones for the paper)
    priority_combos = [
        ('PH',               ['PH']),
        ('PSRT',             ['PSRT']),
        ('PH+PSRT',          ['PH', 'PSRT']),
        ('CPF',              ['CPF']),
        ('FM',               ['FM']),
        ('MFP',              ['MFP']),
        ('PHY',              ['PHY']),
        ('PH+CPF',           ['PH', 'CPF']),
        ('PSRT+CPF',         ['PSRT', 'CPF']),
        ('PH+PSRT+CPF',      ['PH', 'PSRT', 'CPF']),
        ('PH+PSRT+FM',       ['PH', 'PSRT', 'FM']),
        ('PH+PSRT+CPF+FM',   ['PH', 'PSRT', 'CPF', 'FM']),
        ('PH+PSRT+CPF+FM+MFP+PHY', ['PH', 'PSRT', 'CPF', 'FM', 'MFP', 'PHY']),
    ]

    print(f"\n{'Combination':<35} {'10fold-R':>9} {'RMSE':>7} {'LOO-R':>7}  subtype Rs")
    print("=" * 100)

    for name, mod_keys in priority_combos:
        Xs = [modalities[k] for k in mod_keys]
        X  = np.hstack(Xs) if len(Xs) > 1 else Xs[0]

        t0 = time.time()
        print(f"  {name:<33} ", end='', flush=True)

        oof_10  = nested_cv(X, y, n_outer=10)
        oof_loo = nested_loo(X, y)

        m10  = metrics(y, oof_10)
        mloo = metrics(y, oof_loo)
        lo, hi = bootstrap_ci(y, oof_10)

        st_rs = {}
        for st in subtypes:
            idx = np.where(sub == st)[0]
            if len(idx) >= 3:
                st_rs[st] = round(pearsonr(y[idx], oof_10[idx])[0], 3)

        elapsed = time.time() - t0
        print(f"R10={m10['r']:.4f}  RMSE={m10['rmse']:.4f}  "
              f"LOO={mloo['r']:.4f}  [{lo:.3f},{hi:.3f}]  {elapsed:.0f}s")

        row = {'combination': name,
               'modalities': '+'.join(mod_keys),
               'n_modalities': len(mod_keys),
               'dim': X.shape[1],
               'r_10fold':    m10['r'],
               'rho_10fold':  m10['rho'],
               'rmse_10fold': m10['rmse'],
               'r2_10fold':   m10['r2'],
               'ci_lo':  lo, 'ci_hi': hi,
               'r_loo':  mloo['r'],
               'rmse_loo': mloo['rmse'],
               }
        for st, rv in st_rs.items():
            row[f'r_{st}'] = rv
        rows.append(row)

        np.save(os.path.join(MODEL_DIR, f"oof_ablation_{name}.npy"), oof_10)

    out_csv = os.path.join(RES_DIR, 'step_ablation_results.csv')
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    print(f"\nSaved: {out_csv}")

    # Print formatted table for paper
    print("\n" + "=" * 80)
    print("ABLATION TABLE — 10-fold CV (R ± notes)")
    print("=" * 80)
    header = f"{'Combination':<35}  {'R':>6}  {'RMSE':>6}  {'LOO-R':>6}  {'95%CI':>14}"
    print(header)
    print("-" * 80)
    for row in sorted(rows, key=lambda x: x['r_10fold'], reverse=True):
        ci = f"[{row['ci_lo']:.3f},{row['ci_hi']:.3f}]"
        print(f"  {row['combination']:<33}  "
              f"{row['r_10fold']:>6.4f}  {row['rmse_10fold']:>6.4f}  "
              f"{row['r_loo']:>6.4f}  {ci:>14}")

if __name__ == '__main__':
    main()
