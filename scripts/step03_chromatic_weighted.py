"""
Step 03 — Chromatic and Chemical-Weight PH features (CAML-RNA)

Chromatic (Ch): RNA atoms split into 11 functional subtypes
  (element × structural category) × 10 ligand ES types = 110 pairs → 5940-dim

Weighted (W): same 80 pairs as ES+CS but distances scaled by
  d_w(i,j) = d(i,j) / sqrt(w_i * w_j)   [cross-distances only]
  → 4320-dim (same shape as step02, different values)

WeightedChromatic (WCh): chemical weights applied within chromatic pairs
  → 5940-dim

Saves step03_chromatic.npz, step03_weighted.npz, step03_wch.npz
"""

import numpy as np
import pandas as pd
import os, time
from scipy.spatial.distance import cdist
from ripser import ripser
from Bio.PDB import PDBParser
from rdkit import Chem

# ── Reuse parsing helpers from step02 ────────────────────────────────────────

BACKBONE_ATOMS = {"P","OP1","OP2","O5'","O3'","O5*","O3*"}
SUGAR_ATOMS    = {"C1'","C2'","C3'","C4'","C5'","O2'","O4'",
                  "C1*","C2*","C3*","C4*","C5*","O2*","O4*"}
PURINES        = {"A","DA","G","DG"}
PYRIMIDINES    = {"C","DC","U","DU","T","DT"}
RNA_ES         = {'C':0,'N':1,'O':2,'P':3}
LIG_ES         = {'C':0,'N':1,'O':2,'S':3,'P':4,'F':5,'CL':6,'BR':7,'I':8,'H':9}

FILTRATION_STEPS = np.arange(0.5, 12.5, 0.5)   # 24 steps
MAX_ATOMS = 200

# ── Chemical weights (RNA atom types, 11 chromatic labels) ────────────────────
# Label encoding: (es, cs) → chromatic_id
# ES: C=0,N=1,O=2,P=3  |  CS: purine=0,pyrimidine=1,backbone=2,sugar=3
CHROM_MAP = {
    (0,0): 0,   # C_purine
    (0,1): 1,   # C_pyrimidine
    (0,2): 2,   # C_backbone
    (0,3): 3,   # C_sugar
    (1,0): 4,   # N_purine
    (1,1): 5,   # N_pyrimidine
    (2,0): 6,   # O_purine
    (2,1): 7,   # O_pyrimidine
    (2,2): 8,   # O_backbone
    (2,3): 9,   # O_sugar
    (3,2): 10,  # P_backbone
}
N_CHROM = 11

# RNA chemical weights per chromatic subtype
RNA_WEIGHTS = {
    0:  1.10,   # C_purine    – aromatic C, base stacking
    1:  1.10,   # C_pyrimidine
    2:  0.90,   # C_backbone  – aliphatic, van der Waals
    3:  0.95,   # C_sugar
    4:  1.65,   # N_purine    – strong H-bond donor/acceptor (N1,N3,N7,N9)
    5:  1.65,   # N_pyrimidine – N3, N1 Watson-Crick contacts
    6:  1.35,   # O_purine    – carbonyl O6, N7-adjacent O (weaker donor)
    7:  1.55,   # O_pyrimidine – C2=O, C4=O strong H-bond acceptors
    8:  1.45,   # O_backbone  – non-bridging phosphate O (OP1, OP2, O5', O3')
    9:  1.55,   # O_sugar     – 2'-OH: unique to RNA, donor+acceptor
    10: 1.80,   # P_backbone  – formal negative charge, strongest electrostatic
}

# Ligand ES weights
LIG_WEIGHTS = {0:1.00, 1:1.55, 2:1.40, 3:1.20, 4:1.35,
               5:0.90, 6:0.85, 7:0.80, 8:1.05, 9:0.70}

# ── PDB/SDF parsing ───────────────────────────────────────────────────────────

def get_cs(res_name, atom_name):
    an, rn = atom_name.strip(), res_name.strip().upper()
    if an in BACKBONE_ATOMS: return 2
    if an in SUGAR_ATOMS:    return 3
    if rn in PURINES:        return 0
    if rn in PYRIMIDINES:    return 1
    return None

def parse_pocket(pocket_pdb):
    """Returns list of (coord, es, cs, chrom) for each heavy atom."""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure('r', pocket_pdb)
    atoms = []
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
                    if cs is None:
                        continue
                    chrom = CHROM_MAP.get((es, cs), None)
                    if chrom is None:
                        continue
                    atoms.append((atom.get_coord(), es, cs, chrom))
        break
    return atoms

def parse_ligand(sdf_path):
    """Returns list of (coord, es) for heavy atoms."""
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=False)
    atoms = []
    for mol in supplier:
        if mol is None: continue
        conf = mol.GetConformer()
        for a in mol.GetAtoms():
            sym = a.GetSymbol().upper()
            if sym == 'H': continue
            lab = LIG_ES.get(sym, None)
            if lab is None: continue
            p = conf.GetAtomPosition(a.GetIdx())
            atoms.append((np.array([p.x, p.y, p.z]), lab))
        break
    return atoms

# ── PH helpers ───────────────────────────────────────────────────────────────

def subsample(coords, n=MAX_ATOMS, seed=42):
    if len(coords) <= n: return coords
    rng = np.random.default_rng(seed)
    return coords[rng.choice(len(coords), n, replace=False)]

def betti_curve(dgm):
    curve = np.zeros(len(FILTRATION_STEPS), dtype=np.float32)
    for b, d in dgm:
        if np.isinf(d): d = FILTRATION_STEPS[-1] + 1.0
        curve += ((FILTRATION_STEPS >= b) & (FILTRATION_STEPS < d)).astype(np.float32)
    return curve

def pers_stats(dgm):
    finite = [(b,d) for b,d in dgm if not np.isinf(d)]
    if not finite: return np.zeros(3, dtype=np.float32)
    lf = np.array([d-b for b,d in finite])
    return np.array([lf.sum(), lf.max(), len(finite)], dtype=np.float32)

def ph_pair(pts):
    """Run Ripser on pts, return 54-dim feature."""
    try:
        res = ripser(pts, maxdim=1, thresh=12.0)
    except Exception:
        return np.zeros(54, dtype=np.float32)
    d0 = res['dgms'][0]
    d1 = res['dgms'][1] if len(res['dgms'])>1 else np.zeros((0,2))
    return np.concatenate([betti_curve(d0), betti_curve(d1),
                           pers_stats(d0), pers_stats(d1)]).astype(np.float32)

def run_bipartite(rna_c, lig_c, w_rna=None, w_lig=None):
    """
    Compute 54-dim PH feature for one bipartite pair.
    If w_rna/w_lig provided, apply chemical-weight distance scaling.
    """
    rna_c = subsample(rna_c)
    lig_c = subsample(lig_c)

    D = cdist(rna_c, lig_c)
    vr = np.where(D.min(axis=1) <= 12.0)[0]
    vc = np.where(D.min(axis=0) <= 12.0)[0]
    if len(vr) == 0 or len(vc) == 0:
        return np.zeros(54, dtype=np.float32)

    rna_c = rna_c[vr]
    lig_c = lig_c[vc]

    if w_rna is not None and w_lig is not None:
        # Apply chemical-weight scaling to cross distances only
        # Build full coordinate set, then modify the distance matrix
        # Ripser accepts a precomputed distance matrix
        n_r, n_l = len(rna_c), len(lig_c)
        pts = np.vstack([rna_c, lig_c])
        n   = n_r + n_l
        D_full = cdist(pts, pts)
        # Scale cross-distances
        for i in range(n_r):
            for j in range(n_r, n):
                scaled = D_full[i, j] / np.sqrt(w_rna * w_lig)
                D_full[i, j] = scaled
                D_full[j, i] = scaled
        try:
            res = ripser(D_full, maxdim=1, thresh=12.0, distance_matrix=True)
        except Exception:
            return np.zeros(54, dtype=np.float32)
        d0 = res['dgms'][0]
        d1 = res['dgms'][1] if len(res['dgms'])>1 else np.zeros((0,2))
        return np.concatenate([betti_curve(d0), betti_curve(d1),
                               pers_stats(d0), pers_stats(d1)]).astype(np.float32)
    else:
        return ph_pair(np.vstack([rna_c, lig_c]))

# ── Per-complex feature builders ──────────────────────────────────────────────

def build_chromatic_features(rna_atoms, lig_atoms, weighted=False):
    """110 pairs × 54 = 5940-dim."""
    # Group by type
    rna_by_chrom = {k: [] for k in range(N_CHROM)}
    for (coord, es, cs, chrom) in rna_atoms:
        rna_by_chrom[chrom].append(coord)
    rna_by_chrom = {k: np.array(v) for k, v in rna_by_chrom.items()}

    lig_by_es = {k: [] for k in range(10)}
    for (coord, es) in lig_atoms:
        lig_by_es[es].append(coord)
    lig_by_es = {k: np.array(v) for k, v in lig_by_es.items()}

    parts = []
    for r_ch in range(N_CHROM):
        rna_c = rna_by_chrom.get(r_ch, np.array([]).reshape(0,3))
        w_r   = RNA_WEIGHTS[r_ch] if weighted else None
        for l_es in range(10):
            lig_c = lig_by_es.get(l_es, np.array([]).reshape(0,3))
            if len(rna_c) == 0 or len(lig_c) == 0:
                parts.append(np.zeros(54, dtype=np.float32))
            else:
                w_l = LIG_WEIGHTS[l_es] if weighted else None
                parts.append(run_bipartite(rna_c, lig_c, w_r, w_l))
    return np.concatenate(parts)

def build_weighted_esccs(rna_atoms, lig_atoms):
    """Weighted version of step02: 80 ES+CS pairs × 54 = 4320-dim."""
    # Group RNA by ES and CS
    rna_es = {k: [] for k in range(4)}
    rna_cs = {k: [] for k in range(4)}
    es_cs_w = {}  # (scheme, type) → (weight, coords_list)
    for (coord, es, cs, chrom) in rna_atoms:
        rna_es[es].append(coord)
        rna_cs[cs].append(coord)
    rna_es = {k: np.array(v) for k,v in rna_es.items()}
    rna_cs = {k: np.array(v) for k,v in rna_cs.items()}

    # Average RNA weight for ES label (average over chromatic subtypes present)
    # For ES: we use element-level weight (mean of chromatic weights for that element)
    ES_WEIGHTS = {}
    for (es, cs), chrom in CHROM_MAP.items():
        ES_WEIGHTS.setdefault(es, []).append(RNA_WEIGHTS[chrom])
    ES_WEIGHTS = {k: np.mean(v) for k, v in ES_WEIGHTS.items()}

    # CS-level weights (mean across elements in that category)
    CS_WEIGHTS = {}
    for (es, cs), chrom in CHROM_MAP.items():
        CS_WEIGHTS.setdefault(cs, []).append(RNA_WEIGHTS[chrom])
    CS_WEIGHTS = {k: np.mean(v) for k, v in CS_WEIGHTS.items()}

    lig_by_es = {k: [] for k in range(10)}
    for (coord, es) in lig_atoms:
        lig_by_es[es].append(coord)
    lig_by_es = {k: np.array(v) for k,v in lig_by_es.items()}

    parts = []
    for r_es in range(4):
        rna_c = rna_es.get(r_es, np.array([]).reshape(0,3))
        w_r = ES_WEIGHTS[r_es]
        for l_es in range(10):
            lig_c = lig_by_es.get(l_es, np.array([]).reshape(0,3))
            if len(rna_c)==0 or len(lig_c)==0:
                parts.append(np.zeros(54, dtype=np.float32))
            else:
                parts.append(run_bipartite(rna_c, lig_c, w_r, LIG_WEIGHTS[l_es]))
    for r_cs in range(4):
        rna_c = rna_cs.get(r_cs, np.array([]).reshape(0,3))
        w_r = CS_WEIGHTS[r_cs]
        for l_es in range(10):
            lig_c = lig_by_es.get(l_es, np.array([]).reshape(0,3))
            if len(rna_c)==0 or len(lig_c)==0:
                parts.append(np.zeros(54, dtype=np.float32))
            else:
                parts.append(run_bipartite(rna_c, lig_c, w_r, LIG_WEIGHTS[l_es]))
    return np.concatenate(parts)

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    df = pd.read_csv("/home/stalin/Desktop/CAML/data/dataset_clean.csv")
    n  = len(df)
    feat_dir = "/home/stalin/Desktop/CAML/features"

    X_ch  = np.zeros((n, N_CHROM*10*54), dtype=np.float32)   # 5940
    X_w   = np.zeros((n, 80*54),          dtype=np.float32)   # 4320
    X_wch = np.zeros((n, N_CHROM*10*54), dtype=np.float32)   # 5940
    failed = []

    print(f"Processing {n} complexes — Chromatic (5940) + Weighted (4320) + WCh (5940)")

    for i, row in df.iterrows():
        pdb = row['pdb']
        t0  = time.time()
        try:
            rna_atoms = parse_pocket(row['pocket_file'])
            lig_atoms = parse_ligand(row['ligand_sdf'])

            X_ch[i]  = build_chromatic_features(rna_atoms, lig_atoms, weighted=False)
            X_w[i]   = build_weighted_esccs(rna_atoms, lig_atoms)
            X_wch[i] = build_chromatic_features(rna_atoms, lig_atoms, weighted=True)

            elapsed = time.time() - t0
            if (i+1) % 10 == 0:
                print(f"  [{i+1:3d}/{n}] {pdb}  {elapsed:.1f}s  "
                      f"ch_nz={np.count_nonzero(X_ch[i])}  "
                      f"w_nz={np.count_nonzero(X_w[i])}")
        except Exception as e:
            failed.append((pdb, str(e)))
            print(f"  [{i+1:3d}/{n}] {pdb}  FAILED: {e}")

    y    = df['pKd'].values.astype(np.float32)
    pdbs = df['pdb'].values

    np.savez(os.path.join(feat_dir, "step03_chromatic.npz"),   X=X_ch,  y=y, pdbs=pdbs)
    np.savez(os.path.join(feat_dir, "step03_weighted.npz"),    X=X_w,   y=y, pdbs=pdbs)
    np.savez(os.path.join(feat_dir, "step03_wch.npz"),         X=X_wch, y=y, pdbs=pdbs)

    print(f"\nSaved chromatic {X_ch.shape}, weighted {X_w.shape}, wch {X_wch.shape}")
    if failed:
        print(f"Failed: {[f[0] for f in failed]}")

if __name__ == "__main__":
    main()
