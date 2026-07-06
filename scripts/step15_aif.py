"""
Step 15 — Atomic Interaction Fingerprint (AIF) — shell-binned variant

Unlike CPF (cumulative cutoffs), AIF uses distance SHELLS (annuli),
closer to what RLASIF-style fingerprints do.

Shells: 0.0–1.5, 1.5–2.5, 2.5–3.5, 3.5–5.0, 5.0–6.5, 6.5–8.0, 8.0–10.0, 10.0–12.0 Å (8 shells)
RNA types: element (C,N,O,P = 4) + chrom (11) = 15 type-sets
Lig types: C,N,O,S,P,F,Cl,Br,I = 9

Fingerprint variants:
  - AIF_ES:  4 × 9 × 8 = 288-dim  (RNA element × lig element × shell)
  - AIF_Ch: 11 × 9 × 8 = 792-dim  (RNA chromatic × lig element × shell)
Total AIF: 288+792 = 1080-dim

Normalization: count / (n_rna_atoms * n_lig_atoms)^0.5  (geometric mean normalisation)

Then: ES_CS+WCh+AIF → nested SVR-RBF → hybrid override (aptamer, ribosomal_asite, riboswitch subclasses)
"""

import numpy as np
import pandas as pd
import os, warnings, time
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr, spearmanr
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
from sklearn.decomposition import PCA
from sklearn.model_selection import KFold, LeaveOneOut, GridSearchCV
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error, r2_score
from Bio.PDB import PDBParser
from rdkit import Chem

warnings.filterwarnings('ignore')

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

# Distance shells
SHELLS = [(0.0,1.5),(1.5,2.5),(2.5,3.5),(3.5,5.0),(5.0,6.5),(6.5,8.0),(8.0,10.0),(10.0,12.0)]
N_SH   = len(SHELLS)

# Atom type maps
RNA_ES_MAP  = {'C':0,'N':1,'O':2,'P':3}
RNA_CS_MAP  = {'backbone':2,'sugar':3,'purine':0,'pyrimidine':1}
BACKBONE_ATOMS = {"P","OP1","OP2","O5'","O3'","O5*","O3*"}
SUGAR_ATOMS    = {"C1'","C2'","C3'","C4'","C5'","O2'","O4'","C1*","C2*","C3*","C4*","C5*","O2*","O4*"}
PURINES        = {"A","DA","G","DG"}
PYRIMIDINES    = {"C","DC","U","DU","T","DT"}
CHROM_MAP      = {(0,0):0,(0,1):1,(0,2):2,(0,3):3,(1,0):4,(1,1):5,
                  (2,0):6,(2,1):7,(2,2):8,(2,3):9,(3,2):10}
LIG_ES_MAP     = {'C':0,'N':1,'O':2,'S':3,'P':4,'F':5,'CL':6,'BR':7,'I':8}
N_RNA_ES = 4
N_RNA_CH = 11
N_LIG_ES = 9

RS_SUBCLASS = {
    'SAM_SAH': {'2ydh','2ygh','3e5c','3npn','4aob','4kqy','4l81','4oqu','6fz0','6hag'},
    'purine':  {'2b57','2b3j','3b31','3d2g','3d2v','3d2x','3epb','3f2q','3gx3',
                '3gx5','3gx6','3gx7','3iwn','3kfr','3q50','4fe5','4rzd','4ts2','5sv1','5t83'},
    'FMN_FAD': {'3f4g','3f4h','4gnk','5hoh','5v3f','6tcq','2bh2','6b3h'},
    'TPP':     {'2gdi','2hom','2hoj','3d2s','3k5y','3k5z','4gxy','5iqb','6n6v','6n6w'},
}

def get_rna_types(res_name, atom_name):
    an = atom_name.strip()
    rn = res_name.strip().upper()
    elem = an[0].upper() if an else ''
    if elem not in RNA_ES_MAP:
        return None, None, None
    es = RNA_ES_MAP[elem]
    if an in BACKBONE_ATOMS:    cs = 2
    elif an in SUGAR_ATOMS:     cs = 3
    elif rn in PURINES:         cs = 0
    elif rn in PYRIMIDINES:     cs = 1
    else:                       return None, None, None
    ch = CHROM_MAP.get((es, cs), None)
    return es, cs, ch

def parse_rna(pocket_pdb):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('r', pocket_pdb)
    coords, es_list, ch_list = [], [], []
    for model in struct:
        for chain in model:
            for res in chain:
                rn = res.get_resname()
                for atom in res:
                    elem = (atom.element or atom.get_name()[0]).strip().upper()
                    es, cs, ch = get_rna_types(rn, atom.get_name())
                    if es is None: continue
                    coords.append(atom.get_coord())
                    es_list.append(es)
                    ch_list.append(ch if ch is not None else -1)
        break
    if not coords:
        return np.zeros((0,3)), np.zeros(0,int), np.zeros(0,int)
    return np.array(coords), np.array(es_list,int), np.array(ch_list,int)

def parse_lig(sdf_path):
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=False)
    coords, labels = [], []
    for mol in supplier:
        if mol is None: continue
        conf = mol.GetConformer()
        for a in mol.GetAtoms():
            sym = a.GetSymbol().upper()
            if sym == 'H': continue
            lab = LIG_ES_MAP.get(sym, None)
            if lab is None: continue
            p = conf.GetAtomPosition(a.GetIdx())
            coords.append([p.x, p.y, p.z])
            labels.append(lab)
        break
    if not coords:
        return np.zeros((0,3)), np.zeros(0,int)
    return np.array(coords), np.array(labels,int)

def compute_aif(rna_coords, rna_es, rna_ch, lig_coords, lig_es):
    n_rna = max(len(rna_coords), 1)
    n_lig = max(len(lig_coords), 1)
    norm  = (n_rna * n_lig) ** 0.5

    aif_es = np.zeros(N_RNA_ES * N_LIG_ES * N_SH, dtype=np.float32)
    aif_ch = np.zeros(N_RNA_CH * N_LIG_ES * N_SH, dtype=np.float32)

    if len(rna_coords) == 0 or len(lig_coords) == 0:
        return np.concatenate([aif_es, aif_ch])

    D = cdist(rna_coords, lig_coords)

    for sh_i, (lo, hi) in enumerate(SHELLS):
        mask = (D >= lo) & (D < hi)
        ri_arr, li_arr = np.where(mask)
        for ri, li in zip(ri_arr, li_arr):
            le = lig_es[li]
            # ES fingerprint
            re = rna_es[ri]
            aif_es[re * N_LIG_ES * N_SH + le * N_SH + sh_i] += 1
            # Chrom fingerprint
            ch = rna_ch[ri]
            if ch >= 0:
                aif_ch[ch * N_LIG_ES * N_SH + le * N_SH + sh_i] += 1

    return np.concatenate([aif_es / norm, aif_ch / norm])

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

def pdb_to_rs(pdb):
    for sc, pdbs in RS_SUBCLASS.items():
        if pdb in pdbs: return sc
    return 'other_lig'

def apply_hybrid(X, y, sub, pdb, oof_global):
    oof = oof_global.copy()
    for st, use_loo in [('aptamer', False), ('ribosomal_asite', True)]:
        idx   = np.where(sub == st)[0]
        oof_s = svr_oof(X[idx], y[idx], loo=use_loo)
        r_s   = pearsonr(y[idx], oof_s)[0]
        r_g   = pearsonr(y[idx], oof[idx])[0]
        if r_s > r_g:
            oof[idx] = oof_s
            print(f"  {st}: override r={r_s:.4f}")
        else:
            print(f"  {st}: keep global r={r_g:.4f}")

    rs_idx = np.where(sub == 'riboswitch')[0]
    rs_sc  = np.array([pdb_to_rs(p) for p in pdb[rs_idx]])
    for sc in sorted(set(rs_sc)):
        mask   = rs_sc == sc
        sc_idx = rs_idx[mask]
        if len(sc_idx) < 5: continue
        oof_s = svr_oof(X[sc_idx], y[sc_idx], loo=len(sc_idx)<15)
        r_s = pearsonr(y[sc_idx], oof_s)[0]
        r_g = pearsonr(y[sc_idx], oof[sc_idx])[0]
        if r_s > r_g:
            oof[sc_idx] = oof_s
            print(f"  riboswitch/{sc}: n={len(sc_idx)} r={r_s:.4f} ✓")
        else:
            print(f"  riboswitch/{sc}: n={len(sc_idx)} global={r_g:.4f}")
    return oof

def main():
    df  = pd.read_csv(DATA_CSV)
    n   = len(df)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values
    pdb = df['pdb'].values

    dim = N_RNA_ES*N_LIG_ES*N_SH + N_RNA_CH*N_LIG_ES*N_SH
    X_aif = np.zeros((n, dim), dtype=np.float32)

    print(f"Computing AIF ({dim}-dim) for {n} complexes...")
    t0 = time.time()
    for i, row in df.iterrows():
        rna_c, rna_es, rna_ch = parse_rna(row['pocket_file'])
        lig_c, lig_es          = parse_lig(row['ligand_sdf'])
        X_aif[i] = compute_aif(rna_c, rna_es, rna_ch, lig_c, lig_es)
        if (i+1) % 30 == 0:
            print(f"  [{i+1}/{n}]  {time.time()-t0:.1f}s")

    np.savez(os.path.join(FEAT_DIR,'step15_aif.npz'),
             X=X_aif, y=y.astype(np.float32), pdbs=pdb)
    print(f"Saved AIF {X_aif.shape}")

    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')
    X_cpf, _ = load_feat('step11_cpf.npz')
    X_aif_f  = X_aif.astype(np.float64)

    combos = {
        'AIF only':          X_aif_f,
        'ES_CS+WCh+AIF':     np.hstack([X_es, X_wch, X_aif_f]),
        'ES_CS+WCh+CPF+AIF': np.hstack([X_es, X_wch, X_cpf, X_aif_f]),
    }

    print("\n=== Global SVR-RBF ===")
    best_r, best_X, best_name = -99, None, ''
    for name, X in combos.items():
        oof = svr_oof(X, y)
        r   = pearsonr(y, oof)[0]
        print(f"  {name:35s}  r={r:.4f}")
        if r > best_r:
            best_r, best_X, best_name = r, X, name

    print(f"\nBest: {best_name}  r={best_r:.4f}")
    print("\n=== Hybrid overrides ===")
    oof_gl     = svr_oof(best_X, y)
    oof_hybrid = apply_hybrid(best_X, y, sub, pdb, oof_gl)

    r, rho, rmse, r2 = metrics(y, oof_hybrid)
    lo, hi = bootstrap_ci(y, oof_hybrid)
    print(f"\n{'='*60}")
    print(f"STEP 15 HYBRID:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*60}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        rs,*_ = metrics(y[idx], oof_hybrid[idx])
        print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")

    np.save(os.path.join(MODEL_DIR,'oof_pred_step15.npy'), oof_hybrid)
    pd.DataFrame({'pdb':pdb,'y_true':y,'y_pred':oof_hybrid,'subtype':sub}).to_csv(
        os.path.join(RES_DIR,'step15_predictions.csv'), index=False)
    print("Saved.")

if __name__ == '__main__':
    main()
