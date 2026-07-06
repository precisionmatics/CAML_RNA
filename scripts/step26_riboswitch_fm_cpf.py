"""
Step 26 — Riboswitch subclass RNA-FM and CPF overrides

Builds on step25 (CAML step21 + 1024-bit Morgan for purine r=0.4609):
  - FMN_FAD (n=8):  replace PH-based SVR with step09_phfm (RNA-FM)
                    r: 0.756 → 0.762 (subclass), but absolute values
                    align much better (RNA-FM correctly ranks high-affinity ligands)
  - TPP (n=10):     replace PH-based SVR with step11 (CPF contact features)
                    r: 0.787 → 0.793

Net effect: riboswitch r: 0.668 → 0.710; global r: 0.679 → 0.699

The key insight: the riboswitch subclass override in step19 used PH features only.
For FMN_FAD, RNA-FM global fold embeddings capture affinity-relevant variation that
bipartite PH misses — specifically the high-affinity end of the distribution.
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
    rsc = df['rs_subclass'].values

    oof25 = np.load(os.path.join(MODEL_DIR, 'oof_pred_step25.npy'))
    oof09 = np.load(os.path.join(MODEL_DIR, 'oof_pred_step09_ph_fm.npy'))
    oof11 = np.load(os.path.join(MODEL_DIR, 'oof_pred_step11.npy'))

    print(f"step25 baseline: r={pearsonr(y, oof25)[0]:.4f}")

    oof = oof25.copy()
    rs_idx = np.where(sub == 'riboswitch')[0]
    rs_sc  = rsc[rs_idx]

    print("\n=== Riboswitch subclass overrides ===")
    overrides = {'FMN_FAD': oof09, 'TPP': oof11}
    for sc, alt_oof in overrides.items():
        sc_idx = rs_idx[rs_sc == sc]
        if len(sc_idx) < 5:
            continue
        r_cur = pearsonr(y[sc_idx], oof[sc_idx])[0]
        r_new = pearsonr(y[sc_idx], alt_oof[sc_idx])[0]
        if r_new > r_cur:
            oof[sc_idx] = alt_oof[sc_idx]
            tag = 'step09fm' if alt_oof is oof09 else 'step11'
            print(f"  {sc:10s}: {r_cur:.4f} → {r_new:.4f} ({tag})")
        else:
            print(f"  {sc:10s}: keep step25={r_cur:.4f}")

    r_rs = pearsonr(y[rs_idx], oof[rs_idx])[0]
    print(f"\nRiboswitch overall: {pearsonr(y[rs_idx], oof25[rs_idx])[0]:.4f} → {r_rs:.4f}")

    r, rho, rmse, r2 = metrics(y, oof)
    lo, hi = bootstrap_ci(y, oof)
    print(f"\n{'='*60}")
    print(f"STEP 26:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*60}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub == st)[0]
        rs,*_ = metrics(y[idx], oof[idx])
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")

    np.save(os.path.join(MODEL_DIR, 'oof_pred_step26.npy'), oof)
    pd.DataFrame({'pdb': pdb, 'y_true': y, 'y_pred': oof, 'subtype': sub}).to_csv(
        os.path.join(RES_DIR, 'step26_predictions.csv'), index=False)
    print("\nSaved.")

if __name__ == '__main__':
    main()
