"""
Resume both 10-fold CV and ablation jobs — only runs missing combos.
Skips any combo whose OOF file already exists.
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
    return dict(r=round(r,4), rho=round(rho,4), rmse=round(rmse,4), r2=round(r2,4))

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
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=inner, scoring='r2', n_jobs=-1)
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
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=inner, scoring='r2', n_jobs=-1)
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
                        return np.array(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024))
                    except Exception:
                        pass
    return np.zeros(1024)

def load_oof(prefix, name):
    safe = name.replace(' ', '_').replace('+', '+')
    p = os.path.join(MODEL_DIR, f"{prefix}{safe}.npy")
    if os.path.exists(p):
        return np.load(p)
    return None

def save_oof(prefix, name, oof):
    safe = name.replace(' ', '_')
    np.save(os.path.join(MODEL_DIR, f"{prefix}{safe}.npy"), oof)

def main():
    df  = pd.read_csv(DATA_CSV)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values

    print("Loading features...")
    d02    = np.load(os.path.join(FEAT_DIR, 'step02_es_cs_features.npz'), allow_pickle=True)
    X_ph   = d02['X'].astype(np.float64)
    X_psrt = d02['X_psrt'].astype(np.float64)
    X_full = d02['X_full'].astype(np.float64)
    X_cpf  = np.load(os.path.join(FEAT_DIR, 'step11_cpf.npz'),   allow_pickle=True)['X'].astype(np.float64)
    X_fm   = np.load(os.path.join(FEAT_DIR, 'step09_rnafm.npz'), allow_pickle=True)['X'].astype(np.float64)
    X_wch  = np.load(os.path.join(FEAT_DIR, 'step03_wch.npz'),   allow_pickle=True)['X'].astype(np.float64)
    X_morgan = np.array([safe_morgan(row) for _, row in df.iterrows()])
    lig_cols = ['mol_weight','n_rings','n_hbd','n_hba','n_rot_bonds','tpsa',
                'lig_C','lig_N','lig_O','lig_S','n_lig_atoms']
    rna_cols = ['n_rna_atoms','rna_C','rna_N','rna_O','rna_P','rna_S']
    X_phys = np.hstack([df[lig_cols].fillna(0).values,
                        df[rna_cols].fillna(0).values]).astype(np.float64)

    # ── 10-fold CV — full combo set ───────────────────────────────────────────
    combos_10 = {
        'PH_bipartite':           X_ph,
        'PSRT_fvec':              X_psrt,
        'PH+PSRT':                X_full,
        'CPF':                    X_cpf,
        'RNA-FM':                 X_fm,
        'Morgan':                 X_morgan,
        'Physicochemical':        X_phys,
        'WCh':                    X_wch,
        'PH+CPF':                 np.hstack([X_ph, X_cpf]),
        'PH+PSRT+CPF':            np.hstack([X_full, X_cpf]),
        'PH+PSRT+CPF+FM':         np.hstack([X_full, X_cpf, X_fm]),
        'PH+PSRT+CPF+FM+Phys':    np.hstack([X_full, X_cpf, X_fm, X_phys]),
        'Full_ensemble':          np.hstack([X_full, X_cpf, X_fm, X_morgan, X_phys]),
    }

    print("\n=== 10-fold CV (resuming) ===")
    rows_10 = []
    for name, X in combos_10.items():
        oof_10 = load_oof('oof_10fold_', name)
        oof_loo = load_oof('oof_loo_', name)
        t0 = time.time()

        if oof_10 is None:
            print(f"  Running 10fold: {name} ...", flush=True)
            oof_10 = nested_cv(X, y, n_outer=10)
            save_oof('oof_10fold_', name, oof_10)
        else:
            print(f"  [cached] {name}", flush=True)

        if oof_loo is None:
            print(f"  Running LOO:    {name} ...", flush=True)
            oof_loo = nested_loo(X, y)
            save_oof('oof_loo_', name, oof_loo)

        m10  = metrics(y, oof_10)
        mloo = metrics(y, oof_loo)
        lo, hi = bootstrap_ci(y, oof_10)

        st_rs = {}
        for st in sorted(set(sub)):
            idx = np.where(sub == st)[0]
            if len(idx) >= 3:
                st_rs[f'r_10fold_{st}'] = round(pearsonr(y[idx], oof_10[idx])[0], 4)
                st_rs[f'r_loo_{st}']    = round(pearsonr(y[idx], oof_loo[idx])[0], 4)

        row = {'combination': name,
               'r_10fold': m10['r'], 'rho_10fold': m10['rho'],
               'rmse_10fold': m10['rmse'], 'r2_10fold': m10['r2'],
               'ci_lo': lo, 'ci_hi': hi,
               'r_loo': mloo['r'], 'rmse_loo': mloo['rmse'],
               'n': len(y), **st_rs}
        rows_10.append(row)
        print(f"    R10={m10['r']:.4f}  LOO={mloo['r']:.4f}  [{lo:.3f},{hi:.3f}]  {time.time()-t0:.0f}s")

    pd.DataFrame(rows_10).to_csv(os.path.join(RES_DIR, 'step_10fold_cv_results.csv'), index=False)
    print(f"\nSaved: step_10fold_cv_results.csv")

    # ── Ablation — priority combos ────────────────────────────────────────────
    modalities = {
        'PH':   X_ph,
        'PSRT': X_psrt,
        'CPF':  X_cpf,
        'FM':   X_fm,
        'MFP':  X_morgan,
        'PHY':  X_phys,
    }
    priority_combos = [
        ('PH',                     ['PH']),
        ('PSRT',                   ['PSRT']),
        ('PH+PSRT',                ['PH', 'PSRT']),
        ('CPF',                    ['CPF']),
        ('FM',                     ['FM']),
        ('MFP',                    ['MFP']),
        ('PHY',                    ['PHY']),
        ('PH+CPF',                 ['PH', 'CPF']),
        ('PSRT+CPF',               ['PSRT', 'CPF']),
        ('PH+PSRT+CPF',            ['PH', 'PSRT', 'CPF']),
        ('PH+PSRT+FM',             ['PH', 'PSRT', 'FM']),
        ('PH+PSRT+CPF+FM',         ['PH', 'PSRT', 'CPF', 'FM']),
        ('PH+PSRT+CPF+FM+MFP+PHY', ['PH', 'PSRT', 'CPF', 'FM', 'MFP', 'PHY']),
    ]

    print("\n=== Ablation (resuming) ===")
    rows_ab = []
    for name, mod_keys in priority_combos:
        Xs = [modalities[k] for k in mod_keys]
        X  = np.hstack(Xs) if len(Xs) > 1 else Xs[0]

        oof_10 = load_oof('oof_ablation_', name)
        t0 = time.time()

        if oof_10 is None:
            print(f"  Running: {name} ...", flush=True)
            oof_10 = nested_cv(X, y, n_outer=10)
            save_oof('oof_ablation_', name, oof_10)
        else:
            print(f"  [cached] {name}", flush=True)

        oof_loo = load_oof('oof_loo_ablation_', name)
        if oof_loo is None:
            oof_loo = nested_loo(X, y)
            save_oof('oof_loo_ablation_', name, oof_loo)

        m10  = metrics(y, oof_10)
        mloo = metrics(y, oof_loo)
        lo, hi = bootstrap_ci(y, oof_10)

        st_rs = {}
        for st in sorted(set(sub)):
            idx = np.where(sub == st)[0]
            if len(idx) >= 3:
                st_rs[f'r_{st}'] = round(pearsonr(y[idx], oof_10[idx])[0], 4)

        row = {'combination': name,
               'modalities': '+'.join(mod_keys),
               'n_modalities': len(mod_keys),
               'dim': X.shape[1],
               'r_10fold': m10['r'], 'rho_10fold': m10['rho'],
               'rmse_10fold': m10['rmse'], 'r2_10fold': m10['r2'],
               'ci_lo': lo, 'ci_hi': hi,
               'r_loo': mloo['r'], 'rmse_loo': mloo['rmse'],
               **st_rs}
        rows_ab.append(row)
        print(f"    R10={m10['r']:.4f}  LOO={mloo['r']:.4f}  [{lo:.3f},{hi:.3f}]  {time.time()-t0:.0f}s")

    pd.DataFrame(rows_ab).to_csv(os.path.join(RES_DIR, 'step_ablation_results.csv'), index=False)
    print(f"\nSaved: step_ablation_results.csv")

    # ── Final summary ─────────────────────────────────────────────────────────
    print("\n" + "=" * 80)
    print("ABLATION TABLE (sorted by 10-fold R)")
    print(f"{'Combination':<35}  {'R10':>6}  {'RMSE':>6}  {'LOO-R':>6}  {'95%CI':>14}")
    print("-" * 80)
    for r in sorted(rows_ab, key=lambda x: x['r_10fold'], reverse=True):
        ci = f"[{r['ci_lo']:.3f},{r['ci_hi']:.3f}]"
        print(f"  {r['combination']:<33}  {r['r_10fold']:>6.4f}  "
              f"{r['rmse_10fold']:>6.4f}  {r['r_loo']:>6.4f}  {ci:>14}")

    print("\n" + "=" * 80)
    print("10-FOLD CV TABLE (sorted by R10)")
    print(f"{'Combination':<35}  {'R10':>6}  {'RMSE':>6}  {'LOO-R':>6}")
    print("-" * 80)
    for r in sorted(rows_10, key=lambda x: x['r_10fold'], reverse=True):
        print(f"  {r['combination']:<33}  {r['r_10fold']:>6.4f}  "
              f"{r['rmse_10fold']:>6.4f}  {r['r_loo']:>6.4f}")

if __name__ == '__main__':
    main()
