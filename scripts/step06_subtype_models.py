"""
Step 06 — Per-subtype SVR-RBF models on ES_CS+WCh features

Uses the subtype column from dataset.csv to train separate SVR-RBF
models per RNA subtype. Global SVR-RBF used for subtypes with n < 8.
Final OOF predictions stitched together.
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
from sklearn.linear_model import Ridge
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
    rs  = []
    for _ in range(n):
        idx = rng.integers(0, ns, ns)
        try:
            rs.append(pearsonr(y_true[idx], y_pred[idx])[0])
        except Exception:
            pass
    return np.percentile(rs, 2.5), np.percentile(rs, 97.5)

def make_pipe(model):
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  model),
    ])

def svr_cv(X, y, use_loo=False):
    """Returns OOF predictions using nested CV or LOO for small n."""
    param_grid = {'model__C': [0.1, 1, 10, 100, 500],
                  'model__gamma': ['scale', 'auto']}
    if use_loo or len(y) < 15:
        # LOO outer, 5-fold inner
        loo = LeaveOneOut()
        oof = np.zeros(len(y))
        for tr, te in loo.split(X):
            if len(tr) < 6:
                oof[te] = y.mean()
                continue
            inner = min(5, len(tr))
            pipe  = make_pipe(SVR(kernel='rbf'))
            gs    = GridSearchCV(pipe, param_grid,
                                 cv=KFold(inner, shuffle=True, random_state=SEED),
                                 scoring='r2', n_jobs=-1)
            gs.fit(X[tr], y[tr])
            oof[te] = gs.best_estimator_.predict(X[te])
        return oof
    else:
        outer = KFold(n_splits=5, shuffle=True, random_state=SEED)
        inner = KFold(n_splits=5, shuffle=True, random_state=SEED)
        oof   = np.zeros(len(y))
        for tr, te in outer.split(X):
            pipe = make_pipe(SVR(kernel='rbf'))
            gs   = GridSearchCV(pipe, param_grid, cv=inner,
                                scoring='r2', n_jobs=-1)
            gs.fit(X[tr], y[tr])
            oof[te] = gs.best_estimator_.predict(X[te])
        return oof

def main():
    df = pd.read_csv(DATA_CSV)

    # Load best combined features from step05
    X_es,  y = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X = np.hstack([X_es, X_wch])   # 10260-dim, same as step05 best

    # Subtype column
    subtypes = df['subtype'].values
    unique   = sorted(df['subtype'].unique())
    print(f"Subtypes: {unique}")
    print(f"Counts:   {df['subtype'].value_counts().to_dict()}\n")

    oof_global = np.load(os.path.join(MODEL_DIR, 'oof_pred_step05_best.npy'))
    oof_sub    = oof_global.copy()   # start from global, override per subtype

    rows = []
    for st in unique:
        idx   = np.where(subtypes == st)[0]
        n_sub = len(idx)
        X_sub, y_sub = X[idx], y[idx]

        if n_sub < 5:
            print(f"  {st:20s}  n={n_sub:3d}  SKIP (too small, use global)")
            continue

        oof_s = svr_cv(X_sub, y_sub, use_loo=(n_sub < 15))
        oof_sub[idx] = oof_s

        try:
            r, rho, rmse, r2 = metrics(y_sub, oof_s)
            print(f"  {st:20s}  n={n_sub:3d}  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}")
            rows.append({'subtype': st, 'n': n_sub, 'r': r, 'rho': rho,
                         'rmse': rmse, 'r2': r2})
        except Exception as e:
            print(f"  {st:20s}  n={n_sub:3d}  metrics failed: {e}")

    # Overall after subtype models
    r, rho, rmse, r2 = metrics(y, oof_sub)
    lo, hi = bootstrap_ci(y, oof_sub)
    print(f"\n{'='*55}")
    print(f"COMBINED (subtype models):  r={r:.4f}  rho={rho:.4f}  "
          f"rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*55}")
    print(f"vs step05 global:  r={pearsonr(y, oof_global)[0]:.4f}")

    pd.DataFrame(rows).to_csv(
        os.path.join(RES_DIR, 'step06_subtype_results.csv'), index=False)
    np.save(os.path.join(MODEL_DIR, 'oof_pred_step06.npy'), oof_sub)
    print(f"\nSaved results and OOF predictions.")

if __name__ == '__main__':
    main()
