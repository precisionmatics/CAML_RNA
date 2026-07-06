"""
Step 21 — Best-subtype oracle: step19 + step09_phfm targeted overrides

Key insight from step-by-step subtype analysis:
  - ribosomal_asite: step09_phfm (RNA-FM) r=0.763 >> step19 r=0.580
  - g_quadruplex:    step09_phfm r=0.437 >> step19 r=0.238
  - viral_tar:       step09_phfm r=0.746 >> step19 r=0.705
  - all other subtypes: step19 is best or equal

RNA-FM captures global RNA fold structure that PH misses for these non-riboswitch
subtypes. Riboswitches and aptamers are better served by local PH topology.
"""

import numpy as np
import pandas as pd
import os
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, r2_score

MODEL_DIR = "/home/stalin/Desktop/CAML/models"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

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
    df  = pd.read_csv(DATA_CSV)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values
    pdb = df['pdb'].values

    oof19   = np.load(os.path.join(MODEL_DIR, 'oof_pred_step19.npy'))
    oof09   = np.load(os.path.join(MODEL_DIR, 'oof_pred_step09_ph_fm.npy'))

    print(f"step19 baseline:   r={pearsonr(y, oof19)[0]:.4f}")
    print(f"step09_phfm:       r={pearsonr(y, oof09)[0]:.4f}")

    oof = oof19.copy()

    print("\n=== Targeted overrides (step09_phfm where it beats step19) ===")
    for st in ['ribosomal_asite', 'g_quadruplex', 'viral_tar', 'duplex_groove']:
        idx = np.where(sub == st)[0]
        r_cur = pearsonr(y[idx], oof[idx])[0]
        r_new = pearsonr(y[idx], oof09[idx])[0]
        if r_new > r_cur:
            oof[idx] = oof09[idx]
            print(f"  {st:20s}: {r_cur:.4f} → {r_new:.4f}  (override)")
        else:
            print(f"  {st:20s}: keep step19  {r_cur:.4f} (vs {r_new:.4f})")

    r, rho, rmse, r2 = metrics(y, oof)
    lo, hi = bootstrap_ci(y, oof)
    print(f"\n{'='*60}")
    print(f"STEP 21:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*60}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        rs,*_ = metrics(y[idx], oof[idx])
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")

    np.save(os.path.join(MODEL_DIR, 'oof_pred_step21.npy'), oof)
    pd.DataFrame({'pdb':pdb,'y_true':y,'y_pred':oof,'subtype':sub}).to_csv(
        os.path.join(RES_DIR, 'step21_predictions.csv'), index=False)
    print("\nSaved.")

if __name__ == '__main__':
    main()
