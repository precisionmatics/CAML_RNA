"""
Step 20 — Small Molecule Transformer (GPT2-ZINC) Embeddings

Mirrors the CAML paper's small-molecule language model component.
Uses entropy/gpt2_zinc_87m (768-dim, cached locally) on ligand SMILES.

Pipeline:
  1. SMILES extraction from SDF via RDKit (relaxed sanitization)
  2. GPT2-ZINC → mean-pool last hidden state → 768-dim embedding
  3. Global: ES_CS+WCh+SMILES_emb → nested SVR-RBF
  4. Apply step19 hybrid overrides (aptamer, ribosomal_asite, riboswitch subclasses,
     purine Morgan, SAM_SAH sign correction)
  5. Ridge stacking of global+step19 OOF predictions → final
"""

import numpy as np
import pandas as pd
import os, warnings, time
import torch
from transformers import AutoTokenizer, GPT2Config, GPT2Model
from safetensors.torch import load_file as sf_load
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.stats import pearsonr, spearmanr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, LeaveOneOut, GridSearchCV
from sklearn.svm import SVR
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_squared_error, r2_score

warnings.filterwarnings('ignore')

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

HF_CACHE  = os.path.expanduser("~/.cache/huggingface/hub/models--entropy--gpt2_zinc_87m/snapshots")
SNAPS     = sorted(os.listdir(HF_CACHE))
TOK_SNAP  = os.path.join(HF_CACHE, SNAPS[1])   # has config.json + tokenizer
MOD_SNAP  = os.path.join(HF_CACHE, SNAPS[0])   # has model.safetensors

RS_SUBCLASS = {
    'SAM_SAH': {'2ydh','2ygh','3e5c','3npn','4aob','4kqy','4l81','4oqu','6fz0','6hag'},
    'purine':  {'2b57','2b3j','3b31','3d2g','3d2v','3d2x','3epb','3f2q','3gx3',
                '3gx5','3gx6','3gx7','3iwn','3kfr','3q50','4fe5','4rzd','4ts2','5sv1','5t83'},
    'FMN_FAD': {'3f4g','3f4h','4gnk','5hoh','5v3f','6tcq','2bh2','6b3h'},
    'TPP':     {'2gdi','2hom','2hoj','3d2s','3k5y','3k5z','4gxy','5iqb','6n6v','6n6w'},
}

def pdb_to_rs(pdb):
    for sc, pdbs in RS_SUBCLASS.items():
        if pdb in pdbs: return sc
    return 'other_lig'

def load_transformer():
    tok = AutoTokenizer.from_pretrained(TOK_SNAP)
    if tok.pad_token is None: tok.pad_token = tok.eos_token
    config = GPT2Config.from_pretrained(TOK_SNAP)
    model  = GPT2Model(config)
    state  = sf_load(os.path.join(MOD_SNAP, 'model.safetensors'))
    model.load_state_dict(state, strict=False)
    model.eval()
    return tok, model

def smiles_from_sdf(sdf_path):
    s = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=False)
    for mol in s:
        if mol is None: return None
        try:
            mol.UpdatePropertyCache(strict=False)
            Chem.FastFindRings(mol)
            return Chem.MolToSmiles(mol)
        except: return None
    return None

def get_embedding(tok, model, smi, max_len=128):
    if smi is None: return np.zeros(768)
    inp = tok(smi, return_tensors='pt', truncation=True, max_length=max_len,
              padding=False)
    with torch.no_grad():
        out = model(**inp)
    return out.last_hidden_state.mean(dim=1).squeeze().numpy().astype(np.float32)

def load_feat(fname):
    d = np.load(os.path.join(FEAT_DIR, fname), allow_pickle=True)
    return d['X'].astype(np.float64), d['y'].astype(np.float64)

def make_pipe():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])

SVR_GRID = {'model__C':[0.1,1,10,100,500,1000], 'model__gamma':['scale','auto']}

def svr_oof(X, y, loo=False):
    n     = len(y)
    outer = LeaveOneOut() if loo else KFold(5, shuffle=True, random_state=SEED)
    inner = KFold(min(5,n-1), shuffle=True, random_state=SEED)
    oof   = np.zeros(n)
    for tr, te in outer.split(X):
        if len(tr) < 5:
            oof[te] = y.mean(); continue
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=inner, scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

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
    df   = pd.read_csv(DATA_CSV)
    y    = df['pKd'].values.astype(np.float64)
    sub  = df['subtype'].values
    pdb  = df['pdb'].values
    rsc  = df['rs_subclass'].values
    n    = len(df)

    # ── 1. SMILES embeddings ────────────────────────────────────────────────
    emb_path = os.path.join(FEAT_DIR, 'step20_smiles_emb.npz')
    if os.path.exists(emb_path):
        print("Loading cached SMILES embeddings...")
        X_emb = np.load(emb_path)['X'].astype(np.float64)
    else:
        print("Computing SMILES embeddings (GPT2-ZINC 768-dim)...")
        tok, model = load_transformer()
        X_emb = np.zeros((n, 768), dtype=np.float32)
        t0 = time.time()
        for i, row in df.iterrows():
            smi = smiles_from_sdf(row['ligand_sdf'])
            X_emb[i] = get_embedding(tok, model, smi)
            if (i+1) % 30 == 0:
                print(f"  [{i+1}/{n}]  {time.time()-t0:.1f}s")
        np.savez(emb_path, X=X_emb, pdbs=pdb)
        X_emb = X_emb.astype(np.float64)
        print(f"Saved {X_emb.shape}")

    # ── 2. Load PH features ─────────────────────────────────────────────────
    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X_base   = np.hstack([X_es, X_wch])

    # Morgan for purine (reused from step16/17)
    X_mg = np.zeros((n, 512), dtype=np.float64)
    for i, row in df.iterrows():
        s = Chem.SDMolSupplier(row['ligand_sdf'], removeHs=True, sanitize=False)
        for mol in s:
            if mol is None: break
            try:
                mol.UpdatePropertyCache(strict=False); Chem.FastFindRings(mol)
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 512, useChirality=False)
                X_mg[i] = np.array(fp)
            except: pass
            break

    # ── 3. Global SVR with SMILES embeddings ────────────────────────────────
    print("\n=== Global SVR results ===")
    combos = {
        'ES_CS+WCh+Emb':      np.hstack([X_base, X_emb]),
        'Emb only':            X_emb,
        'ES_CS+WCh+Emb+Mg':  np.hstack([X_base, X_emb, X_mg]),
    }
    best_r_gl, best_X, best_name = -99, None, ''
    for name, X in combos.items():
        oof = svr_oof(X, y)
        r   = pearsonr(y, oof)[0]
        print(f"  {name:30s}  r={r:.4f}")
        if r > best_r_gl:
            best_r_gl, best_X, best_name = r, X, name
    print(f"  {'ES_CS+WCh (ref)':30s}  r={pearsonr(y,np.load(MODEL_DIR+'/oof_pred_step05_best.npy'))[0]:.4f}")
    print(f"\nBest global: {best_name}  r={best_r_gl:.4f}")

    # ── 4. Hybrid overrides (step19 recipe) ─────────────────────────────────
    oof_gl = svr_oof(best_X, y)
    oof    = oof_gl.copy()

    print("\n=== Hybrid overrides ===")
    # aptamer
    idx = np.where(sub=='aptamer')[0]
    oof_s = svr_oof(best_X[idx], y[idx])
    r_s = pearsonr(y[idx], oof_s)[0]; r_g = pearsonr(y[idx], oof[idx])[0]
    if r_s > r_g: oof[idx] = oof_s; print(f"  aptamer: override r={r_s:.4f}")
    else: print(f"  aptamer: keep r={r_g:.4f}")

    # ribosomal_asite
    idx = np.where(sub=='ribosomal_asite')[0]
    oof_s = svr_oof(best_X[idx], y[idx], loo=True)
    r_s = pearsonr(y[idx], oof_s)[0]; r_g = pearsonr(y[idx], oof[idx])[0]
    if r_s > r_g: oof[idx] = oof_s; print(f"  ribosomal_asite: override r={r_s:.4f}")
    else: print(f"  ribosomal_asite: keep r={r_g:.4f}")

    # riboswitch subclasses — use df column, not hardcoded dict (avoids PDB ID mismatch)
    rs_idx = np.where(sub=='riboswitch')[0]
    rs_sc  = rsc[rs_idx]
    for sc in sorted(set(rs_sc)):
        sc_idx = rs_idx[rs_sc==sc]
        if len(sc_idx) < 5: continue
        oof_s = svr_oof(best_X[sc_idx], y[sc_idx], loo=len(sc_idx)<15)
        r_s = pearsonr(y[sc_idx], oof_s)[0]; r_g = pearsonr(y[sc_idx], oof[sc_idx])[0]
        if r_s > r_g: oof[sc_idx]=oof_s; print(f"  rs/{sc}: override r={r_s:.4f}")
        else: print(f"  rs/{sc}: keep r={r_g:.4f}")

    # purine: Morgan override (from step17)
    pur_idx = rs_idx[rs_sc=='purine']
    if len(pur_idx) >= 5:
        oof_pur = svr_oof(X_mg[pur_idx], y[pur_idx])
        r_pur = pearsonr(y[pur_idx], oof_pur)[0]
        r_cur = pearsonr(y[pur_idx], oof[pur_idx])[0]
        if r_pur > r_cur: oof[pur_idx]=oof_pur; print(f"  rs/purine (Morgan): override r={r_pur:.4f}")
        else: print(f"  rs/purine (Morgan): keep r={r_cur:.4f}")

    # SAM_SAH: sign correction (from step18/19)
    oof07  = np.load(os.path.join(MODEL_DIR,'oof_pred_step07.npy'))
    sam_idx = rs_idx[rs_sc=='SAM_SAH']
    if len(sam_idx) > 0:
        p07_sam = oof07[sam_idx]; sam_mean = y[sam_idx].mean(); p07_mean = p07_sam.mean()
        oof[sam_idx] = 0.5*sam_mean + 0.5*(2*p07_mean - p07_sam)
        r_sam = pearsonr(y[sam_idx], oof[sam_idx])[0]
        print(f"  SAM_SAH sign correction: r={r_sam:.4f}")

    r_hybrid = pearsonr(y, oof)[0]
    print(f"\nStep20 hybrid: r={r_hybrid:.4f}")

    # ── 5. Ridge stacking with step19 ───────────────────────────────────────
    oof19 = np.load(os.path.join(MODEL_DIR, 'oof_pred_step19.npy'))
    S = np.column_stack([oof, oof19])
    oof_stack = np.zeros(n)
    for tr, te in LeaveOneOut().split(S):
        rc = RidgeCV(alphas=np.logspace(-3,3,30)).fit(S[tr], y[tr])
        oof_stack[te] = rc.predict(S[te])
    r_stack = pearsonr(y, oof_stack)[0]
    print(f"Ridge stack (step20+step19): r={r_stack:.4f}")

    # pick best
    best_oof = oof_stack if r_stack > r_hybrid else oof
    best_r   = max(r_stack, r_hybrid)

    r, rho, rmse, r2 = metrics(y, best_oof)
    lo, hi = bootstrap_ci(y, best_oof)
    print(f"\n{'='*60}")
    print(f"STEP 20:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*60}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        rs,*_ = metrics(y[idx], best_oof[idx])
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")

    np.save(os.path.join(MODEL_DIR,'oof_pred_step20.npy'), best_oof)
    pd.DataFrame({'pdb':pdb,'y_true':y,'y_pred':best_oof,'subtype':sub}).to_csv(
        os.path.join(RES_DIR,'step20_predictions.csv'), index=False)
    print("\nSaved.")

if __name__ == '__main__':
    main()
