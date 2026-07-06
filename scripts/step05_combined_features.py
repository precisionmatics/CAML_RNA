"""
Step 05 — Combined feature sets + stacking

Tries:
  1. ES_CS + Chromatic concatenated (4320+5940 = 10260-dim)
  2. ES_CS + WCh concatenated (4320+5940 = 10260-dim)
  3. All four concatenated (4320+5940+4320+5940 = 20520-dim)
  4. SVR-RBF stacking: train a meta-SVR on OOF predictions from step04

Best model saved per combination.
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
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score
import xgboost as xgb

warnings.filterwarnings('ignore')

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
SEED = 42

def load(fname):
    d = np.load(os.path.join(FEAT_DIR, fname), allow_pickle=True)
    return d['X'].astype(np.float64), d['y'].astype(np.float64), d['pdbs']

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

def make_pipe(model):
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  model),
    ])

def nested_cv_svr(X, y, param_grid=None):
    if param_grid is None:
        param_grid = {'model__C': [0.1,1,10,100,500],
                      'model__gamma': ['scale','auto']}
    outer = KFold(n_splits=5, shuffle=True, random_state=SEED)
    inner = KFold(n_splits=5, shuffle=True, random_state=SEED)
    oof   = np.zeros(len(y))
    for tr, te in outer.split(X):
        pipe = make_pipe(SVR(kernel='rbf'))
        gs   = GridSearchCV(pipe, param_grid, cv=inner, scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def nested_cv_model(X, y, model, param_grid):
    outer = KFold(n_splits=5, shuffle=True, random_state=SEED)
    inner = KFold(n_splits=5, shuffle=True, random_state=SEED)
    oof   = np.zeros(len(y))
    for tr, te in outer.split(X):
        pipe = make_pipe(model)
        gs   = GridSearchCV(pipe, param_grid, cv=inner, scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def main():
    X_es, y, pdbs = load('step02_es_cs_features.npz')
    X_ch, _, _    = load('step03_chromatic.npz')
    X_w,  _, _    = load('step03_weighted.npz')
    X_wch,_, _    = load('step03_wch.npz')

    combos = {
        'ES_CS+Ch':   np.hstack([X_es, X_ch]),
        'ES_CS+WCh':  np.hstack([X_es, X_wch]),
        'Ch+WCh':     np.hstack([X_ch, X_wch]),
        'All4':       np.hstack([X_es, X_ch, X_w, X_wch]),
    }

    rows = []
    best_overall_r = 0
    best_overall_pred = None
    best_overall_name = ''

    print("=== Combined feature sets (SVR-RBF) ===")
    for name, X in combos.items():
        print(f"\n{name}  shape={X.shape}")
        oof = nested_cv_svr(X, y)
        r, rho, rmse, r2 = metrics(y, oof)
        lo, hi = bootstrap_ci(y, oof)
        print(f"  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}  "
              f"95%CI [{lo:.3f},{hi:.3f}]")
        rows.append({'combo': name, 'model':'SVR_RBF', 'r':r, 'rho':rho,
                     'rmse':rmse, 'r2':r2, 'ci_lo':lo, 'ci_hi':hi})
        np.save(os.path.join(MODEL_DIR, f"oof_pred_{name}.npy"), oof)
        if r > best_overall_r:
            best_overall_r, best_overall_pred, best_overall_name = r, oof, name

    # Also try GBR and XGB on best combo
    best_combo_X = combos[best_overall_name]
    print(f"\n=== Trying GBR + XGB on best combo ({best_overall_name}) ===")
    for mname, model, pgrid in [
        ('GBR', GradientBoostingRegressor(n_estimators=300, random_state=SEED),
                {'model__learning_rate':[0.05,0.1],'model__max_depth':[2,3,4]}),
        ('XGB', xgb.XGBRegressor(n_estimators=300, random_state=SEED, verbosity=0, n_jobs=-1),
                {'model__learning_rate':[0.05,0.1],'model__max_depth':[2,3,4]}),
        ('Ridge', Ridge(),
                {'model__alpha':[1,10,100,1000]}),
    ]:
        oof = nested_cv_model(best_combo_X, y, model, pgrid)
        r, rho, rmse, r2 = metrics(y, oof)
        lo, hi = bootstrap_ci(y, oof)
        print(f"  {mname:8s}  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  "
              f"95%CI [{lo:.3f},{hi:.3f}]")
        rows.append({'combo': best_overall_name, 'model': mname, 'r':r, 'rho':rho,
                     'rmse':rmse, 'r2':r2, 'ci_lo':lo, 'ci_hi':hi})
        if r > best_overall_r:
            best_overall_r, best_overall_pred, best_overall_name = r, oof, f"{best_overall_name}+{mname}"

    # Level-1 stacking: use step04 OOF preds as features for a meta-learner
    print("\n=== Level-1 stacking on step04 OOF predictions ===")
    oof_es  = np.load(os.path.join(MODEL_DIR, 'oof_pred_ES_CS.npy'))
    oof_ch  = np.load(os.path.join(MODEL_DIR, 'oof_pred_Chromatic.npy'))
    oof_w   = np.load(os.path.join(MODEL_DIR, 'oof_pred_Weighted.npy'))
    oof_wch = np.load(os.path.join(MODEL_DIR, 'oof_pred_WCh.npy'))
    meta_X  = np.column_stack([oof_es, oof_ch, oof_w, oof_wch])

    outer = KFold(n_splits=5, shuffle=True, random_state=SEED)
    meta_oof = np.zeros(len(y))
    for tr, te in outer.split(meta_X):
        ridge = Ridge(alpha=1.0).fit(meta_X[tr], y[tr])
        meta_oof[te] = ridge.predict(meta_X[te])
    r, rho, rmse, r2 = metrics(y, meta_oof)
    lo, hi = bootstrap_ci(y, meta_oof)
    print(f"  Meta-Ridge  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  "
          f"95%CI [{lo:.3f},{hi:.3f}]")
    rows.append({'combo':'meta_stack', 'model':'Ridge', 'r':r, 'rho':rho,
                 'rmse':rmse, 'r2':r2, 'ci_lo':lo, 'ci_hi':hi})
    if r > best_overall_r:
        best_overall_r, best_overall_pred = r, meta_oof
        best_overall_name = 'meta_stack'

    pd.DataFrame(rows).to_csv(
        os.path.join(RES_DIR, 'step05_results.csv'), index=False)

    print(f"\n{'='*55}")
    print(f"STEP 05 BEST: {best_overall_name}  r={best_overall_r:.4f}")
    print(f"{'='*55}")
    np.save(os.path.join(MODEL_DIR, 'oof_pred_step05_best.npy'), best_overall_pred)

if __name__ == '__main__':
    main()
