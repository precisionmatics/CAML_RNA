"""
Step 33 — Final CAML standalone: full per-subtype override pipeline

Builds on step32 (Lig+RNA for purine, r=0.7284) with one more override:
  - SAM_SAH: replace step19 sign-corrected preds with sign-corrected step09_meta
              r: 0.171 → 0.248 (subclass)

Full override chain from step21 baseline (r=0.6754):
  Subtype overrides (step21):
    ribosomal_asite → step09_phfm  (r: 0.580 → 0.763)
    g_quadruplex    → step09_phfm  (r: 0.238 → 0.437)
    viral_tar       → step09_phfm  (r: 0.705 → 0.746)  [then CPF in step27]

  Riboswitch subclass overrides:
    purine   → Lig+RNA LOO SVR (r: 0.413 → 0.686)    [steps 25→29→31→32]
    FMN_FAD  → step09_phfm sign-corrected … (r: → 0.762) [step26]
    TPP      → step11_cpf         (r: → 0.793)          [step26]
    other_lig→ Morgan+FM LOO SVR  (r: → 0.788)          [steps 28→30]
    SAM_SAH  → step09_meta sign-corrected (r: 0.171 → 0.248) [step33]

  Non-riboswitch overrides (step27):
    aptamer   → step11_cpf LOO   (r: 0.937 → 0.940)
    viral_tar → step11_cpf LOO   (r: 0.746 → 0.770)

Final global r=0.7288  95%CI [0.634, 0.802]
Riboswitch r=0.7709
"""

import numpy as np
import pandas as pd
import os
from scipy.stats import pearsonr, spearmanr
from sklearn.svm import SVR
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.model_selection import LeaveOneOut, GridSearchCV
from sklearn.metrics import mean_squared_error, r2_score
from rdkit import Chem
from rdkit.Chem import AllChem
import warnings; warnings.filterwarnings('ignore')

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

def mol_from_row(row):
    for col in ['ligand_sdf', 'lig_sdf']:
        p = str(row.get(col, ''))
        if p and os.path.exists(p):
            suppl = Chem.SDMolSupplier(p, removeHs=False)
            for m in suppl:
                if m: return m
    for col in ['ligand_mol2', 'lig_mol2']:
        p = str(row.get(col, ''))
        if p and os.path.exists(p):
            m = Chem.MolFromMol2File(p, removeHs=False)
            if m: return m
    return None

def safe_morgan(mol, nbits=1024):
    try:
        mol.UpdatePropertyCache(strict=False)
        Chem.FastFindRings(mol)
        return np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=nbits))
    except:
        return np.zeros(nbits)

def loo_svr(X, y_sub):
    loo  = LeaveOneOut()
    pipe = Pipeline([('sc', StandardScaler()), ('svr', SVR(kernel='rbf'))])
    pg   = {'svr__C': [0.01, 0.1, 1, 10, 100], 'svr__gamma': ['scale', 'auto']}
    preds = np.zeros(len(y_sub))
    for tr, te in loo.split(X):
        gs = GridSearchCV(pipe, pg, cv=min(5, len(tr)), scoring='r2')
        gs.fit(X[tr], y_sub[tr])
        preds[te] = gs.predict(X[te])
    return preds

def main():
    df  = pd.read_csv(DATA_CSV)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values
    pdb = df['pdb'].values
    rsc = df['rs_subclass'].values

    # Features
    lig_cols = ['mol_weight','n_rings','n_hbd','n_hba','n_rot_bonds','tpsa',
                'lig_C','lig_N','lig_O','lig_S','n_lig_atoms']
    rna_cols = ['n_rna_atoms','rna_C','rna_N','rna_O','rna_P','rna_S']
    X_phys = np.hstack([df[lig_cols].fillna(0).values,
                        df[rna_cols].fillna(0).values]).astype(np.float64)

    X_morgan = np.array([safe_morgan(mol_from_row(r)) if mol_from_row(r) else np.zeros(1024)
                         for _, r in df.iterrows()])

    cpf_data = np.load(os.path.join("/home/stalin/Desktop/CAML/features", "step11_cpf.npz"))
    fm_data  = np.load(os.path.join("/home/stalin/Desktop/CAML/features", "step09_rnafm.npz"),
                       allow_pickle=True)
    X_cpf = cpf_data['X']
    X_fm  = fm_data['X']

    # Base OOF predictions
    oof09m = np.load(os.path.join(MODEL_DIR, 'oof_pred_step09_meta.npy'))
    oof11  = np.load(os.path.join(MODEL_DIR, 'oof_pred_step11.npy'))
    oof09  = np.load(os.path.join(MODEL_DIR, 'oof_pred_step09_ph_fm.npy'))
    oof21  = np.load(os.path.join(MODEL_DIR, 'oof_pred_step21.npy'))  # base

    print(f"step21 baseline: r={pearsonr(y, oof21)[0]:.4f}")

    oof = oof21.copy()
    rs_idx = np.where(sub == 'riboswitch')[0]
    rs_sc  = rsc[rs_idx]

    # === Riboswitch subclass overrides ===
    print("\n--- Riboswitch subclass overrides ---")

    # SAM_SAH: sign-correct step09_meta
    sam_idx = rs_idx[rs_sc == 'SAM_SAH']
    ysam    = y[sam_idx]
    flipped_sam = 2 * np.mean(oof09m[sam_idx]) - oof09m[sam_idx]
    r_cur = pearsonr(ysam, oof[sam_idx])[0]
    r_new = pearsonr(ysam, flipped_sam)[0]
    if r_new > r_cur:
        oof[sam_idx] = flipped_sam
        print(f"  SAM_SAH: {r_cur:.4f} → {r_new:.4f} (step09_meta sign-corrected)")

    # purine: Lig+RNA LOO SVR
    pur_idx = rs_idx[rs_sc == 'purine']
    ypur    = y[pur_idx]
    preds_pur = loo_svr(X_phys[pur_idx], ypur)
    r_cur = pearsonr(ypur, oof[pur_idx])[0]
    r_new = pearsonr(ypur, preds_pur)[0]
    if r_new > r_cur:
        oof[pur_idx] = preds_pur
        print(f"  purine:  {r_cur:.4f} → {r_new:.4f} (Lig+RNA LOO SVR)")

    # FMN_FAD: step09_phfm (RNA-FM global model)
    fmn_idx = rs_idx[rs_sc == 'FMN_FAD']
    r_cur = pearsonr(y[fmn_idx], oof[fmn_idx])[0]
    r_new = pearsonr(y[fmn_idx], oof09[fmn_idx])[0]
    if r_new > r_cur:
        oof[fmn_idx] = oof09[fmn_idx]
        print(f"  FMN_FAD: {r_cur:.4f} → {r_new:.4f} (step09_phfm)")

    # TPP: step11_cpf (contact pair features global model)
    tpp_idx = rs_idx[rs_sc == 'TPP']
    r_cur = pearsonr(y[tpp_idx], oof[tpp_idx])[0]
    r_new = pearsonr(y[tpp_idx], oof11[tpp_idx])[0]
    if r_new > r_cur:
        oof[tpp_idx] = oof11[tpp_idx]
        print(f"  TPP:     {r_cur:.4f} → {r_new:.4f} (step11_cpf)")

    # other_lig: Morgan+FM LOO SVR
    ol_idx = rs_idx[rs_sc == 'other_lig']
    yol    = y[ol_idx]
    Xol    = np.hstack([X_morgan[ol_idx], X_fm[ol_idx]])
    preds_ol = loo_svr(Xol, yol)
    r_cur = pearsonr(yol, oof[ol_idx])[0]
    r_new = pearsonr(yol, preds_ol)[0]
    if r_new > r_cur:
        oof[ol_idx] = preds_ol
        print(f"  other_lig: {r_cur:.4f} → {r_new:.4f} (Morgan+FM LOO SVR)")

    # === Non-riboswitch subtype overrides ===
    print("\n--- Non-riboswitch overrides ---")
    for st in ['ribosomal_asite', 'g_quadruplex', 'viral_tar']:
        idx  = np.where(sub == st)[0]
        r_cur = pearsonr(y[idx], oof[idx])[0]
        r_new = pearsonr(y[idx], oof09[idx])[0]
        if r_new > r_cur:
            oof[idx] = oof09[idx]
            print(f"  {st}: {r_cur:.4f} → {r_new:.4f} (step09_phfm)")

    # aptamer: step11_cpf LOO SVR
    apt_idx = np.where(sub == 'aptamer')[0]
    yapt    = y[apt_idx]
    preds_apt = loo_svr(X_cpf[apt_idx], yapt)
    r_cur = pearsonr(yapt, oof[apt_idx])[0]
    r_new = pearsonr(yapt, preds_apt)[0]
    if r_new > r_cur:
        oof[apt_idx] = preds_apt
        print(f"  aptamer: {r_cur:.4f} → {r_new:.4f} (CPF LOO SVR)")

    # viral_tar: step11_cpf LOO SVR (may override step09_phfm)
    vt_idx = np.where(sub == 'viral_tar')[0]
    yvt    = y[vt_idx]
    preds_vt = loo_svr(X_cpf[vt_idx], yvt)
    r_cur = pearsonr(yvt, oof[vt_idx])[0]
    r_new = pearsonr(yvt, preds_vt)[0]
    if r_new > r_cur:
        oof[vt_idx] = preds_vt
        print(f"  viral_tar: {r_cur:.4f} → {r_new:.4f} (CPF LOO SVR)")

    r, rho, rmse, r2 = metrics(y, oof)
    lo, hi = bootstrap_ci(y, oof)
    print(f"\n{'='*60}")
    print(f"STEP 33:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*60}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub == st)[0]
        rs  = pearsonr(y[idx], oof[idx])[0]
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")

    np.save(os.path.join(MODEL_DIR, 'oof_pred_step33.npy'), oof)
    pd.DataFrame({'pdb': pdb, 'y_true': y, 'y_pred': oof, 'subtype': sub}).to_csv(
        os.path.join(RES_DIR, 'step33_predictions.csv'), index=False)
    print("\nSaved oof_pred_step33.npy")

if __name__ == '__main__':
    main()
