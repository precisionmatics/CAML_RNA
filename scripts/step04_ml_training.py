"""
Step 04 — ML Training and Evaluation (CAML-RNA)

Trains Ridge, Lasso, ElasticNet, SVR-RBF, SVR-Linear, RandomForest,
GradientBoosting, XGBoost on all 4 feature sets:
  - ES+CS (step02): 4320-dim
  - Chromatic (step03): 5940-dim
  - Weighted (step03): 4320-dim
  - WCh (step03): 5940-dim

Evaluation: nested 5-fold CV (outer for OOF predictions,
inner 5-fold grid search for hyperparameters).
Metrics: Pearson r, RMSE, R², Spearman rho.
Bootstrap 95% CI (10000 resamples) on best model per feature set.
Saves: results/step04_results.csv, models/step04_best_models.pkl
"""

import numpy as np
import pandas as pd
import os, joblib, warnings
from scipy.stats import pearsonr, spearmanr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.linear_model import Ridge, Lasso, ElasticNet
from sklearn.svm import SVR
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
import xgboost as xgb

warnings.filterwarnings('ignore')

FEAT_DIR   = "/home/stalin/Desktop/CAML/features"
RES_DIR    = "/home/stalin/Desktop/CAML/results"
MODEL_DIR  = "/home/stalin/Desktop/CAML/models"
os.makedirs(RES_DIR, exist_ok=True)
os.makedirs(MODEL_DIR, exist_ok=True)

SEED = 42
N_OUTER = 5
N_INNER = 5

# ── Feature sets ──────────────────────────────────────────────────────────────

FEATURE_SETS = {
    'ES_CS':      'step02_es_cs_features.npz',
    'Chromatic':  'step03_chromatic.npz',
    'Weighted':   'step03_weighted.npz',
    'WCh':        'step03_wch.npz',
}

# ── Model grid ────────────────────────────────────────────────────────────────

def get_models():
    return {
        'Ridge': (
            Ridge(),
            {'model__alpha': [0.1, 1, 10, 100, 1000, 5000]}
        ),
        'Lasso': (
            Lasso(max_iter=5000),
            {'model__alpha': [0.001, 0.01, 0.1, 1.0, 10.0]}
        ),
        'ElasticNet': (
            ElasticNet(max_iter=5000),
            {'model__alpha': [0.01, 0.1, 1.0], 'model__l1_ratio': [0.2, 0.5, 0.8]}
        ),
        'SVR_RBF': (
            SVR(kernel='rbf'),
            {'model__C': [0.1, 1, 10, 100], 'model__gamma': ['scale', 'auto']}
        ),
        'SVR_Lin': (
            SVR(kernel='linear'),
            {'model__C': [0.01, 0.1, 1, 10]}
        ),
        'RF': (
            RandomForestRegressor(n_estimators=200, random_state=SEED, n_jobs=-1),
            {'model__max_depth': [3, 5, None], 'model__min_samples_leaf': [2, 4]}
        ),
        'GBR': (
            GradientBoostingRegressor(n_estimators=300, random_state=SEED),
            {'model__learning_rate': [0.05, 0.1], 'model__max_depth': [2, 3, 4]}
        ),
        'XGB': (
            xgb.XGBRegressor(n_estimators=300, random_state=SEED,
                             verbosity=0, n_jobs=-1),
            {'model__learning_rate': [0.05, 0.1], 'model__max_depth': [2, 3, 4]}
        ),
    }

# ── Preprocessing pipeline ───────────────────────────────────────────────────

def make_pipe(model):
    return Pipeline([
        ('vt',     VarianceThreshold(threshold=1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  model),
    ])

# ── Evaluation ────────────────────────────────────────────────────────────────

def metrics(y_true, y_pred):
    r,  _  = pearsonr(y_true, y_pred)
    rho, _ = spearmanr(y_true, y_pred)
    rmse   = np.sqrt(mean_squared_error(y_true, y_pred))
    r2     = r2_score(y_true, y_pred)
    return {'r': r, 'rho': rho, 'rmse': rmse, 'r2': r2}

def bootstrap_ci(y_true, y_pred, n=10000, seed=SEED):
    rng = np.random.default_rng(seed)
    n_s = len(y_true)
    rs = []
    for _ in range(n):
        idx = rng.integers(0, n_s, n_s)
        try:
            r, _ = pearsonr(y_true[idx], y_pred[idx])
            rs.append(r)
        except Exception:
            pass
    rs = np.array(rs)
    return np.percentile(rs, 2.5), np.percentile(rs, 97.5)

def nested_cv(X, y, model_name, base_model, param_grid):
    outer_cv = KFold(n_splits=N_OUTER, shuffle=True, random_state=SEED)
    inner_cv = KFold(n_splits=N_INNER, shuffle=True, random_state=SEED)

    oof_pred = np.zeros(len(y))
    fold_rs  = []

    for fold, (tr_idx, te_idx) in enumerate(outer_cv.split(X)):
        X_tr, X_te = X[tr_idx], X[te_idx]
        y_tr, y_te = y[tr_idx], y[te_idx]

        pipe = make_pipe(base_model)
        gs   = GridSearchCV(pipe, param_grid, cv=inner_cv,
                            scoring='r2', n_jobs=-1)
        gs.fit(X_tr, y_tr)
        pred = gs.best_estimator_.predict(X_te)
        oof_pred[te_idx] = pred
        fold_rs.append(pearsonr(y_te, pred)[0])

    return oof_pred, fold_rs

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    rows = []
    best_per_feat = {}   # feat_set → (model_name, oof_pred, r_value)

    for feat_name, fname in FEATURE_SETS.items():
        d = np.load(os.path.join(FEAT_DIR, fname))
        X, y = d['X'].astype(np.float64), d['y'].astype(np.float64)
        print(f"\n{'='*60}")
        print(f"Feature set: {feat_name}  shape={X.shape}")
        print(f"{'='*60}")

        best_r, best_pred, best_model_name = -99, None, None

        for mname, (base_model, param_grid) in get_models().items():
            try:
                oof_pred, fold_rs = nested_cv(X, y, mname, base_model, param_grid)
                m = metrics(y, oof_pred)
                print(f"  {mname:12s}  r={m['r']:.4f}  rho={m['rho']:.4f}  "
                      f"rmse={m['rmse']:.4f}  r2={m['r2']:.4f}  "
                      f"fold_SD={np.std(fold_rs):.3f}")
                rows.append({'feature_set': feat_name, 'model': mname, **m,
                             'fold_r_sd': np.std(fold_rs)})
                if m['r'] > best_r:
                    best_r, best_pred, best_model_name = m['r'], oof_pred, mname
            except Exception as e:
                print(f"  {mname:12s}  ERROR: {e}")

        if best_pred is not None:
            lo, hi = bootstrap_ci(y, best_pred)
            print(f"\n  >> Best: {best_model_name}  OOF r={best_r:.4f}  "
                  f"95%CI [{lo:.3f}, {hi:.3f}]")
            best_per_feat[feat_name] = (best_model_name, best_pred, best_r, lo, hi)

    # Save results table
    res_df = pd.DataFrame(rows)
    res_df.to_csv(os.path.join(RES_DIR, "step04_results.csv"), index=False)
    print(f"\nSaved: {RES_DIR}/step04_results.csv")

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY — Best model per feature set")
    print("="*60)
    for feat_name, (mname, pred, r, lo, hi) in best_per_feat.items():
        print(f"  {feat_name:15s}  {mname:12s}  r={r:.4f}  95%CI [{lo:.3f},{hi:.3f}]")

    # Save best predictions for downstream use
    y_series = None
    for feat_name, fname in FEATURE_SETS.items():
        if y_series is None:
            y_series = np.load(os.path.join(FEAT_DIR, fname))['y']
        if feat_name in best_per_feat:
            _, pred, _, _, _ = best_per_feat[feat_name]
            np.save(os.path.join(MODEL_DIR, f"oof_pred_{feat_name}.npy"), pred)

    np.save(os.path.join(MODEL_DIR, "y_true.npy"), y_series)
    print(f"\nOOF predictions saved to {MODEL_DIR}/")

if __name__ == "__main__":
    main()
