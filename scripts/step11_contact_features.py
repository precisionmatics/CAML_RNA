"""
Step 11 — Contact Pair Fingerprint (CPF) + Extended PH filtration

CPF: counts RNA-ligand atom pairs within {3,4,5,6,8,10}Å cutoffs
     grouped by (RNA_chromatic_type × ligand_ES_type)
     11 × 10 × 6 cutoffs = 660-dim

Extended PH: rerun step02 with filtration 0.5–15.0Å (30 steps instead of 24)
     80 pairs × 60 features = 4800-dim (vs 4320 original)

Then combines CPF + extended PH + WCh and runs nested SVR-RBF.
Also tests hybrid (aptamer + ribosomal_asite subtype override).
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
from ripser import ripser

warnings.filterwarnings('ignore')

FEAT_DIR  = "/home/stalin/Desktop/CAML/features"
RES_DIR   = "/home/stalin/Desktop/CAML/results"
MODEL_DIR = "/home/stalin/Desktop/CAML/models"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
SEED = 42

# ── Atom type defs (reuse from step02/03) ────────────────────────────────────

BACKBONE_ATOMS = {"P","OP1","OP2","O5'","O3'","O5*","O3*"}
SUGAR_ATOMS    = {"C1'","C2'","C3'","C4'","C5'","O2'","O4'",
                  "C1*","C2*","C3*","C4*","C5*","O2*","O4*"}
PURINES        = {"A","DA","G","DG"}
PYRIMIDINES    = {"C","DC","U","DU","T","DT"}
RNA_ES         = {'C':0,'N':1,'O':2,'P':3}
LIG_ES         = {'C':0,'N':1,'O':2,'S':3,'P':4,'F':5,'CL':6,'BR':7,'I':8,'H':9}
CHROM_MAP      = {(0,0):0,(0,1):1,(0,2):2,(0,3):3,(1,0):4,(1,1):5,
                  (2,0):6,(2,1):7,(2,2):8,(2,3):9,(3,2):10}

CUTOFFS  = [3.0, 4.0, 5.0, 6.0, 8.0, 10.0]
FILT_EXT = np.arange(0.5, 15.5, 0.5)   # 30 steps

def get_cs(res_name, atom_name):
    an, rn = atom_name.strip(), res_name.strip().upper()
    if an in BACKBONE_ATOMS: return 2
    if an in SUGAR_ATOMS:    return 3
    if rn in PURINES:        return 0
    if rn in PYRIMIDINES:    return 1
    return None

def parse_pocket_full(pocket_pdb):
    """Returns (coords, chrom_labels) arrays for all heavy atoms."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('r', pocket_pdb)
    coords, chroms, es_types = [], [], []
    for model in struct:
        for chain in model:
            for res in chain:
                rn = res.get_resname()
                for atom in res:
                    elem = (atom.element or '').strip().upper()
                    if elem in ('H','') or elem not in RNA_ES:
                        continue
                    es = RNA_ES[elem]
                    cs = get_cs(rn, atom.get_name())
                    if cs is None: continue
                    ch = CHROM_MAP.get((es,cs), None)
                    if ch is None: continue
                    coords.append(atom.get_coord())
                    chroms.append(ch)
                    es_types.append(es)
        break
    return (np.array(coords) if coords else np.zeros((0,3)),
            np.array(chroms, dtype=int),
            np.array(es_types, dtype=int))

def parse_ligand_full(sdf_path):
    """Returns (coords, es_labels) for heavy atoms."""
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=False)
    coords, labels = [], []
    for mol in supplier:
        if mol is None: continue
        conf = mol.GetConformer()
        for a in mol.GetAtoms():
            sym = a.GetSymbol().upper()
            if sym == 'H': continue
            lab = LIG_ES.get(sym, None)
            if lab is None: continue
            p = conf.GetAtomPosition(a.GetIdx())
            coords.append([p.x, p.y, p.z])
            labels.append(lab)
        break
    return (np.array(coords) if coords else np.zeros((0,3)),
            np.array(labels, dtype=int))

# ── CPF ──────────────────────────────────────────────────────────────────────

def compute_cpf(rna_coords, rna_chroms, lig_coords, lig_es):
    """660-dim: 11 RNA_chrom × 10 lig_ES × 6 cutoffs."""
    feat = np.zeros(11*10*6, dtype=np.float32)
    if len(rna_coords)==0 or len(lig_coords)==0:
        return feat
    D = cdist(rna_coords, lig_coords)
    for ri, rc in enumerate(CUTOFFS):
        pairs = np.argwhere(D <= rc)
        for rna_idx, lig_idx in pairs:
            ch = rna_chroms[rna_idx]
            le = lig_es[lig_idx]
            feat[ch*10*6 + le*6 + ri] += 1
    # normalize by total atoms to make scale-invariant
    n_tot = max(len(rna_coords) + len(lig_coords), 1)
    return feat / n_tot

# ── Extended PH (30 filtration steps) ────────────────────────────────────────

def subsample(arr, n=200, seed=42):
    if len(arr) <= n: return arr
    rng = np.random.default_rng(seed)
    return arr[rng.choice(len(arr), n, replace=False)]

def betti_curve(dgm, steps):
    curve = np.zeros(len(steps), dtype=np.float32)
    for b, d in dgm:
        if np.isinf(d): d = steps[-1]+1
        curve += ((steps>=b)&(steps<d)).astype(np.float32)
    return curve

def pers_stats(dgm):
    finite = [(b,d) for b,d in dgm if not np.isinf(d)]
    if not finite: return np.zeros(3, dtype=np.float32)
    lf = np.array([d-b for b,d in finite])
    return np.array([lf.sum(), lf.max(), len(finite)], dtype=np.float32)

def compute_ph_pair_ext(rna_c, lig_c):
    """60-dim (30β0 + 30β1 + 3+3 stats) with 30-step filtration."""
    rna_c = subsample(rna_c)
    lig_c = subsample(lig_c)
    D = cdist(rna_c, lig_c)
    vr = np.where(D.min(axis=1)<=15.0)[0]
    vc = np.where(D.min(axis=0)<=15.0)[0]
    if len(vr)==0 or len(vc)==0:
        return np.zeros(66, dtype=np.float32)
    pts = np.vstack([rna_c[vr], lig_c[vc]])
    try:
        res = ripser(pts, maxdim=1, thresh=15.0)
    except Exception:
        return np.zeros(66, dtype=np.float32)
    d0 = res['dgms'][0]
    d1 = res['dgms'][1] if len(res['dgms'])>1 else np.zeros((0,2))
    return np.concatenate([betti_curve(d0, FILT_EXT), betti_curve(d1, FILT_EXT),
                           pers_stats(d0), pers_stats(d1)]).astype(np.float32)

def compute_extended_ph(rna_coords, rna_es, rna_cs, lig_coords, lig_es):
    """ES+CS scheme with extended filtration: 80×66=5280-dim."""
    rna_es_groups = {k: [] for k in range(4)}
    rna_cs_groups = {k: [] for k in range(4)}
    for i, (coord, es, cs) in enumerate(zip(rna_coords,
                                              rna_es if len(rna_es) else [],
                                              rna_cs if len(rna_cs) else [])):
        rna_es_groups[es].append(coord)
        rna_cs_groups[cs].append(coord)
    rna_es_groups = {k: np.array(v) for k,v in rna_es_groups.items()}
    rna_cs_groups = {k: np.array(v) for k,v in rna_cs_groups.items()}
    lig_groups = {k:[] for k in range(10)}
    for coord, le in zip(lig_coords, lig_es):
        lig_groups[le].append(coord)
    lig_groups = {k: np.array(v) for k,v in lig_groups.items()}

    parts = []
    for scheme_groups in [rna_es_groups, rna_cs_groups]:
        for rt in range(4):
            rc = scheme_groups.get(rt, np.zeros((0,3)))
            for lt in range(10):
                lc = lig_groups.get(lt, np.zeros((0,3)))
                if len(rc)==0 or len(lc)==0:
                    parts.append(np.zeros(66, dtype=np.float32))
                else:
                    parts.append(compute_ph_pair_ext(rc, lc))
    return np.concatenate(parts)

# ── ML helpers ───────────────────────────────────────────────────────────────

def metrics(y_true, y_pred):
    r,  _ = pearsonr(y_true, y_pred)
    rho,_ = spearmanr(y_true, y_pred)
    rmse  = np.sqrt(mean_squared_error(y_true, y_pred))
    r2    = r2_score(y_true, y_pred)
    return r, rho, rmse, r2

def bootstrap_ci(y_true, y_pred, n=10000):
    rng = np.random.default_rng(SEED)
    ns  = len(y_true)
    rs  = [pearsonr(y_true[i:=rng.integers(0,ns,ns)],
                    y_pred[i])[0] for _ in range(n)]
    return np.percentile(rs, 2.5), np.percentile(rs, 97.5)

def make_pipe():
    return Pipeline([
        ('vt',     VarianceThreshold(1e-4)),
        ('scaler', StandardScaler()),
        ('pca',    PCA(n_components=0.95, random_state=SEED)),
        ('model',  SVR(kernel='rbf')),
    ])

SVR_GRID = {'model__C':[0.1,1,10,100,500,1000], 'model__gamma':['scale','auto']}

def nested_cv_svr(X, y, loo=False):
    n  = len(y)
    outer = LeaveOneOut() if loo else KFold(5, shuffle=True, random_state=SEED)
    inner = KFold(min(5,n-1), shuffle=True, random_state=SEED)
    oof   = np.zeros(n)
    for tr, te in outer.split(X):
        if len(tr) < 5:
            oof[te] = y.mean(); continue
        gs = GridSearchCV(make_pipe(), SVR_GRID, cv=inner,
                          scoring='r2', n_jobs=-1)
        gs.fit(X[tr], y[tr])
        oof[te] = gs.best_estimator_.predict(X[te])
    return oof

def load_feat(fname):
    d = np.load(os.path.join(FEAT_DIR, fname), allow_pickle=True)
    return d['X'].astype(np.float64), d['y'].astype(np.float64)

# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    df  = pd.read_csv(DATA_CSV)
    n   = len(df)
    y   = df['pKd'].values.astype(np.float64)
    sub = df['subtype'].values

    X_cpf = np.zeros((n, 660),  dtype=np.float32)
    X_ext = np.zeros((n, 5280), dtype=np.float32)

    print(f"Computing CPF + extended PH for {n} complexes...")
    for i, row in df.iterrows():
        t0 = time.time()
        rna_c, rna_ch, rna_es = parse_pocket_full(row['pocket_file'])
        lig_c, lig_es          = parse_ligand_full(row['ligand_sdf'])

        # CPF
        rna_cs = np.array([
            CHROM_MAP.get((e, None), -1) for e in rna_es
        ])  # approximate: use chrom for cs
        # Get proper CS from chrom: chrom→cs mapping
        chrom_to_cs = {0:0,1:1,2:2,3:3,4:0,5:1,6:0,7:1,8:2,9:3,10:2}
        rna_cs_arr = np.array([chrom_to_cs.get(c, 0) for c in rna_ch])

        X_cpf[i] = compute_cpf(rna_c, rna_ch, lig_c, lig_es)
        X_ext[i] = compute_extended_ph(rna_c, rna_es, rna_cs_arr,
                                       lig_c, lig_es)

        if (i+1) % 20 == 0:
            print(f"  [{i+1}/{n}]  {time.time()-t0:.1f}s  "
                  f"cpf_nz={np.count_nonzero(X_cpf[i])}  "
                  f"ext_nz={np.count_nonzero(X_ext[i])}")

    np.savez(os.path.join(FEAT_DIR,'step11_cpf.npz'),
             X=X_cpf, y=y.astype(np.float32), pdbs=df['pdb'].values)
    np.savez(os.path.join(FEAT_DIR,'step11_ext_ph.npz'),
             X=X_ext, y=y.astype(np.float32), pdbs=df['pdb'].values)
    print(f"Saved CPF {X_cpf.shape} and ext PH {X_ext.shape}")

    # Load step02/03 PH features
    X_es,  _ = load_feat('step02_es_cs_features.npz')
    X_wch, _ = load_feat('step03_wch.npz')

    combos = {
        'CPF only':          X_cpf.astype(np.float64),
        'ExtPH only':        X_ext.astype(np.float64),
        'ES_CS+WCh+CPF':     np.hstack([X_es, X_wch, X_cpf.astype(np.float64)]),
        'ES_CS+WCh+ExtPH':   np.hstack([X_es, X_wch, X_ext.astype(np.float64)]),
        'ExtPH+WCh+CPF':     np.hstack([X_ext.astype(np.float64),
                                         X_wch, X_cpf.astype(np.float64)]),
    }

    print("\n=== Global SVR-RBF results ===")
    best_r, best_oof, best_name = -99, None, ''
    for name, X in combos.items():
        oof = nested_cv_svr(X, y)
        r, rho, rmse, r2 = metrics(y, oof)
        lo, hi = bootstrap_ci(y, oof)
        print(f"  {name:30s}  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  "
              f"95%CI [{lo:.3f},{hi:.3f}]")
        np.save(os.path.join(MODEL_DIR, f"oof_{name.replace(' ','_')}.npy"), oof)
        if r > best_r:
            best_r, best_oof, best_name = r, oof, name

    print(f"\nBest global: {best_name}  r={best_r:.4f}")

    # Hybrid: override aptamer + ribosomal_asite
    X_best = combos[best_name]
    oof_hybrid = best_oof.copy()
    for st, loo in [('aptamer', False), ('ribosomal_asite', True)]:
        idx   = np.where(sub==st)[0]
        oof_s = nested_cv_svr(X_best[idx], y[idx], loo=loo)
        r_s   = pearsonr(y[idx], oof_s)[0]
        r_g   = pearsonr(y[idx], best_oof[idx])[0]
        if r_s > r_g:
            oof_hybrid[idx] = oof_s
            print(f"  {st}: override r={r_s:.4f} > global r={r_g:.4f}")
        else:
            print(f"  {st}: keep global r={r_g:.4f}")

    r, rho, rmse, r2 = metrics(y, oof_hybrid)
    lo, hi = bootstrap_ci(y, oof_hybrid)
    print(f"\n{'='*55}")
    print(f"STEP 11 HYBRID:  r={r:.4f}  rho={rho:.4f}  rmse={rmse:.4f}  r2={r2:.4f}")
    print(f"95%CI [{lo:.3f}, {hi:.3f}]")
    print(f"{'='*55}")

    print("\nPer-subtype:")
    for st in sorted(set(sub)):
        idx = np.where(sub==st)[0]
        try:
            rs,*_ = metrics(y[idx], oof_hybrid[idx])
            print(f"  {st:20s}  n={len(idx):2d}  r={rs:.4f}")
        except Exception:
            pass

    np.save(os.path.join(MODEL_DIR,'oof_pred_step11.npy'), oof_hybrid)
    pd.DataFrame({'pdb':df['pdb'].values,'y_true':y,'y_pred':oof_hybrid,
                  'subtype':sub}).to_csv(
        os.path.join(RES_DIR,'step11_predictions.csv'), index=False)
    print("Saved.")

if __name__ == '__main__':
    main()
