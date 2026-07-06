"""
Step 09 — RNA-FM pocket sequence embeddings

Extracts the RNA sequence of the binding pocket (from nucleic_acid PDB),
encodes it with RNA-FM (multimolecule), and uses mean-pooled embeddings
as additional features.

Tries:
  1. RNA-FM alone (SVR-RBF)
  2. Hybrid step07 OOF + RNA-FM meta-learner (Ridge)
  3. PH (ES_CS+WCh) + RNA-FM combined (SVR-RBF)
"""

import numpy as np
import pandas as pd
import os, warnings, torch
import warnings
warnings.filterwarnings('ignore')

from multimolecule import RnaFmModel, RnaTokenizer
from Bio.PDB import PDBParser
from scipy.stats import pearsonr, spearmanr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, GridSearchCV
from sklearn.svm import SVR
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_squared_error, r2_score

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

NUCLEOTIDE_MAP = {
    'A':'A','ADE':'A','DA':'A',
    'G':'G','GUA':'G','DG':'G',
    'C':'C','CYT':'C','DC':'C',
    'U':'U','URA':'U',
    'T':'U','DT':'U',  # treat T as U for RNA-FM
}

def extract_pocket_sequence(pocket_pdb):
    """Get RNA sequence from pocket PDB (unique residues in order)."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('r', pocket_pdb)
    seen   = set()
    seq    = []
    for model in struct:
        for chain in model:
            for res in chain:
                rn  = res.get_resname().strip().upper()
                rid = (chain.id, res.get_id())
                nuc = NUCLEOTIDE_MAP.get(rn, None)
                if nuc and rid not in seen:
                    seq.append(nuc)
                    seen.add(rid)
        break
    return ''.join(seq) if seq else 'N'

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
    rs  = [pearsonr(y_true[idx:=rng.integers(0,ns,ns)],
                    y_pred[idx])[0] for _ in range(n)]
    return np.percentile(rs, 2.5), np.percentile(rs, 97.5)

def make_pipe_svr():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])

SVR_GRID = {'model__C':[0.1,1,10,100,500], 'model__gamma':['scale','auto']}

def nested_cv_svr(X, y):
    outer = KFold(5, shuffle=True, random_state=SEED)
    inner = KFold(5, shuffle=True, random_state=SEED)
    oof   = np.zeros(len(y))
    for tr, te in outer.split(X):
        gs = GridSearchCV(make_pipe_svr(), SVR_GRID, cv=inner,
                          scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def main():
    df  = pd.read_csv(DATA_CSV)
    n   = len(df)
    y   = df['pKd'].values.astype(np.float64)

    # Load RNA-FM
    device    = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Loading RNA-FM on {device}...")
    tokenizer = RnaTokenizer.from_pretrained('multimolecule/rnafm')
    model     = RnaFmModel.from_pretrained('multimolecule/rnafm').to(device)
    model.eval()

    # Extract sequences
    print("Extracting pocket sequences and encoding with RNA-FM...")
    embs = np.zeros((n, 640), dtype=np.float32)

    for i, row in df.iterrows():
        seq = extract_pocket_sequence(row['pocket_file'])
        seq = seq[:510]  # RNA-FM max ~512 tokens
        with torch.no_grad():
            inp = tokenizer(seq, return_tensors='pt',
                            padding=True, truncation=True,
                            max_length=512).to(device)
            out = model(**inp)
            # mean-pool over sequence length
            emb = out.last_hidden_state.mean(dim=1).squeeze().cpu().numpy()
        embs[i] = emb
        if (i+1) % 20 == 0:
            print(f"  [{i+1}/{n}]  seq_len={len(seq)}")

    np.savez(os.path.join(FEAT_DIR, 'step09_rnafm.npz'),
             X=embs, y=y.astype(np.float32), pdbs=df['pdb'].values)
    print(f"Saved RNA-FM embeddings: {embs.shape}")

    X_fm = embs.astype(np.float64)

    # 1. RNA-FM alone
    print("\nRNA-FM alone (SVR-RBF):")
    oof_fm = nested_cv_svr(X_fm, y)
    r, rho, rmse, _ = metrics(y, oof_fm)
    lo, hi = bootstrap_ci(y, oof_fm)
    print(f"  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  95%CI [{lo:.3f},{hi:.3f}]")
    np.save(os.path.join(MODEL_DIR, 'oof_pred_rnafm.npy'), oof_fm)

    # 2. PH + RNA-FM
    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X_ph = np.hstack([X_es, X_wch])
    X_ph_fm = np.hstack([X_ph, X_fm])
    print(f"\nPH + RNA-FM (SVR-RBF) shape={X_ph_fm.shape}:")
    oof_ph_fm = nested_cv_svr(X_ph_fm, y)
    r, rho, rmse, r2 = metrics(y, oof_ph_fm)
    lo, hi = bootstrap_ci(y, oof_ph_fm)
    print(f"  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}  "
          f"95%CI [{lo:.3f},{hi:.3f}]")
    np.save(os.path.join(MODEL_DIR, 'oof_pred_step09_ph_fm.npy'), oof_ph_fm)

    # 3. Meta-learner: hybrid step07 + RNA-FM OOF as features
    oof_hybrid = np.load(os.path.join(MODEL_DIR, 'oof_pred_step07.npy'))
    meta_X = np.column_stack([oof_hybrid, oof_fm])
    outer  = KFold(5, shuffle=True, random_state=SEED)
    meta_oof = np.zeros(n)
    for tr, te in outer.split(meta_X):
        m = Ridge(alpha=1.0).fit(meta_X[tr], y[tr])
        meta_oof[te] = m.predict(meta_X[te])
    r, rho, rmse, r2 = metrics(y, meta_oof)
    lo, hi = bootstrap_ci(y, meta_oof)
    print(f"\nMeta (hybrid07 + RNA-FM OOF → Ridge):")
    print(f"  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}  "
          f"95%CI [{lo:.3f},{hi:.3f}]")
    np.save(os.path.join(MODEL_DIR, 'oof_pred_step09_meta.npy'), meta_oof)

    print(f"\n{'='*55}")
    best_r = max(
        pearsonr(y, oof_fm)[0],
        pearsonr(y, oof_ph_fm)[0],
        pearsonr(y, meta_oof)[0],
    )
    print(f"Step 09 best r={best_r:.4f}")
    print(f"{'='*55}")

if __name__ == '__main__':
    main()
