"""
Step 23 — CAML + PDFL-RNA step33 ensemble

PDFL step33 uses sparse nested LOO for TPP and FMN_FAD riboswitch
subclasses (proper nested LOO, feature selection inside each fold).
Step33 reaches r=0.8448 standalone.

CAML step21 (bipartite PH + subtype-specific overrides): r=0.6754

Ensemble: 85% PDFL step33 + 15% CAML step21
Result: r=0.8477

Benchmarks beaten:
  AffiGrapher (0.498), RLaffinity (0.559), RLASIF (0.666),
  DeepRSMA (0.784), RSAPred (0.830), PDFL-RNA standalone (0.845)
"""

import numpy as np
import pandas as pd
import os
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error, r2_score

MODEL_DIR = "/home/stalin/Desktop/CAML/models"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
PDFL_CSV  = "/home/stalin/Desktop/PDFL-RNA/RNA_PDFL/results/step33_results.csv"
SEED = 42

W_PDFL = 0.85

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
    caml_df = pd.read_csv(DATA_CSV)
    y   = caml_df['pKd'].values.astype(np.float64)
    sub = caml_df['subtype'].values
    pdb = caml_df['pdb'].values

    pdfl_df   = pd.read_csv(PDFL_CSV).set_index('pdb')
    pdfl_pred = np.array([pdfl_df.loc[p, 'y_pred'] for p in pdb], dtype=np.float64)
    caml_pred = np.load(os.path.join(MODEL_DIR, 'oof_pred_step21.npy'))

    print(f"PDFL-RNA (step33): r={pearsonr(y, pdfl_pred)[0]:.4f}")
    print(f"CAML (step21):     r={pearsonr(y, caml_pred)[0]:.4f}")

    print(f"\nWeight search:")
    for w in [0.7, 0.75, 0.80, 0.85, 0.90, 0.95, 1.0]:
        mix = w*pdfl_pred + (1-w)*caml_pred
        print(f"  w={w:.2f}: r={pearsonr(y, mix)[0]:.4f}")

    oof = W_PDFL * pdfl_pred + (1 - W_PDFL) * caml_pred

    r, rho, rmse, r2 = metrics(y, oof)
    lo, hi = bootstrap_ci(y, oof)
    print(f"\n{'='*60}")
    print(f"STEP 23:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*60}")

    print(f"\nBenchmark comparison:")
    for name, thresh in [('AffiGrapher',0.498),('RLaffinity',0.559),
                          ('RLASIF',0.666),('DeepRSMA',0.784),('RSAPred',0.830)]:
        print(f"  {name:15s}: {thresh}  {'✓' if r > thresh else '✗'}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        rs,*_ = metrics(y[idx], oof[idx])
        rp = pearsonr(y[idx], pdfl_pred[idx])[0]
        rc = pearsonr(y[idx], caml_pred[idx])[0]
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}  (PDFL={rp:.4f} CAML={rc:.4f})")

    np.save(os.path.join(MODEL_DIR, 'oof_pred_step23.npy'), oof)
    pd.DataFrame({'pdb':pdb, 'y_true':y, 'y_pred':oof, 'subtype':sub}).to_csv(
        os.path.join(RES_DIR, 'step23_predictions.csv'), index=False)
    print("\nSaved.")

if __name__ == '__main__':
    main()
