"""
Step 02 — Element/Category-Specific PH + PSRT f-vector features (CAML-RNA)

Fixes from reviewer comments:
  1. Hydrogen removed from LIG_ES — H is not a heavy atom; RDKit strips it
     with removeHs=True, so H pairs were always zero. Now using 9 elements.
  2. PSRT f/h-vector features added — computed on the TRUE bipartite distance
     matrix (intramolecular distances = inf) per Suwayyid & Wei 2025.
  3. PH features use element-specific VR on interface atoms (same as before).

Schemes:
  ES: 4 RNA elements × 9 ligand elements = 36 pairs
  CS: 4 RNA categories × 9 ligand elements = 36 pairs
  Total: 72 pairs

Per pair features:
  PH:   24 β₀ + 24 β₁ + 3 stats_β₀ + 3 stats_β₁ = 54 features
  PSRT: 24 f₁ + 24 h₂ + 4 scalar stats            = 52 features

Output dimensions:
  X       (PH only):   72 × 54 = 3888-dim
  X_psrt  (PSRT only): 72 × 52 = 3744-dim
  X_full  (combined):  72 × 106 = 7632-dim
"""

import numpy as np
import pandas as pd
import os, time
from ripser import ripser
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser
from rdkit import Chem

# ── Atom-type dictionaries ────────────────────────────────────────────────────

RNA_ES = {'C': 0, 'N': 1, 'O': 2, 'P': 3}

# H excluded — not a heavy atom; always zero in old code, now explicitly removed
LIG_ES = {'C': 0, 'N': 1, 'O': 2, 'S': 3, 'P': 4,
           'F': 5, 'CL': 6, 'BR': 7, 'I': 8}

BACKBONE_ATOMS = {"P","OP1","OP2","O5'","O3'","O5*","O3*"}
SUGAR_ATOMS    = {"C1'","C2'","C3'","C4'","C5'","O2'","O4'",
                  "C1*","C2*","C3*","C4*","C5*","O2*","O4*"}
PURINES        = {"A","DA","G","DG"}
PYRIMIDINES    = {"C","DC","U","DU","T","DT"}

N_RNA_ES  = 4   # C N O P
N_LIG_ES  = 9   # C N O S P F Cl Br I
N_RNA_CS  = 4   # purine_base pyrimidine_base backbone sugar
N_ES_PAIRS = N_RNA_ES * N_LIG_ES   # 36
N_CS_PAIRS = N_RNA_CS * N_LIG_ES   # 36
N_PAIRS    = N_ES_PAIRS + N_CS_PAIRS  # 72

FILTRATION_STEPS = np.arange(0.5, 12.5, 0.5)  # 24 steps (same as original)
N_FILT = len(FILTRATION_STEPS)                  # 24

PH_PER_PAIR   = 2 * N_FILT + 6   # 54: β₀(24)+β₁(24)+3+3
PSRT_PER_PAIR = 2 * N_FILT + 4   # 52: f₁(24)+h₂(24)+4 scalar stats

PH_DIM   = N_PAIRS * PH_PER_PAIR    # 72 × 54 = 3888
PSRT_DIM = N_PAIRS * PSRT_PER_PAIR  # 72 × 52 = 3744
MAX_ATOMS = 200
THRESH    = 12.0

# ── Parsing ───────────────────────────────────────────────────────────────────

def get_cs_label(res_name, atom_name):
    rn, an = res_name.strip().upper(), atom_name.strip()
    if an in BACKBONE_ATOMS: return 2
    if an in SUGAR_ATOMS:    return 3
    if rn in PURINES:        return 0
    if rn in PYRIMIDINES:    return 1
    return None

def parse_rna_pocket(pocket_pdb):
    parser  = PDBParser(QUIET=True)
    struct  = parser.get_structure('rna', pocket_pdb)
    es_atoms = {k: [] for k in range(N_RNA_ES)}
    cs_atoms = {k: [] for k in range(N_RNA_CS)}
    for model in struct:
        for chain in model:
            for res in chain:
                rn = res.get_resname()
                for atom in res:
                    elem = (atom.element or '').strip().upper()
                    if elem in ('H', '') or elem not in RNA_ES:
                        continue
                    coord = atom.get_coord()
                    es_atoms[RNA_ES[elem]].append(coord)
                    cs = get_cs_label(rn, atom.get_name())
                    if cs is not None:
                        cs_atoms[cs].append(coord)
        break
    to_arr = lambda d: {k: np.array(v) if v else np.zeros((0,3)) for k,v in d.items()}
    return to_arr(es_atoms), to_arr(cs_atoms)

def parse_ligand_sdf(sdf_path):
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=False)
    lig_atoms = {k: [] for k in range(N_LIG_ES)}
    for mol in supplier:
        if mol is None: continue
        conf = mol.GetConformer()
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol().upper()
            if sym == 'H': continue
            label = LIG_ES.get(sym, None)
            if label is None: continue
            pos = conf.GetAtomPosition(atom.GetIdx())
            lig_atoms[label].append([pos.x, pos.y, pos.z])
        break
    return {k: np.array(v) if v else np.zeros((0,3)) for k,v in lig_atoms.items()}

# ── PH computation (element-specific VR on interface atoms) ───────────────────

def subsample(arr, n=MAX_ATOMS, seed=42):
    if len(arr) <= n: return arr
    return arr[np.random.default_rng(seed).choice(len(arr), n, replace=False)]

def betti_curve(dgm, steps):
    curve = np.zeros(len(steps), dtype=np.float32)
    for b, d in dgm:
        if np.isinf(d): d = steps[-1] + 1.0
        curve += ((steps >= b) & (steps < d)).astype(np.float32)
    return curve

def pers_stats(dgm):
    finite = [(b,d) for b,d in dgm if not np.isinf(d)]
    if not finite: return np.zeros(3, dtype=np.float32)
    lf = np.array([d-b for b,d in finite])
    return np.array([lf.sum(), lf.max(), float(len(finite))], dtype=np.float32)

def compute_ph_pair(rna_c, lig_c):
    rna_c, lig_c = subsample(rna_c), subsample(lig_c)
    D = cdist(rna_c, lig_c)
    vr = np.where(D.min(axis=1) <= THRESH)[0]
    vc = np.where(D.min(axis=0) <= THRESH)[0]
    if len(vr) == 0 or len(vc) == 0:
        return np.zeros(PH_PER_PAIR, dtype=np.float32)
    pts = np.vstack([rna_c[vr], lig_c[vc]])
    try:
        res = ripser(pts, maxdim=1, thresh=THRESH)
    except Exception:
        return np.zeros(PH_PER_PAIR, dtype=np.float32)
    d0 = res['dgms'][0]
    d1 = res['dgms'][1] if len(res['dgms']) > 1 else np.zeros((0,2))
    return np.concatenate([betti_curve(d0, FILTRATION_STEPS),
                           betti_curve(d1, FILTRATION_STEPS),
                           pers_stats(d0), pers_stats(d1)]).astype(np.float32)

# ── PSRT f/h-vector (TRUE bipartite: intramolecular dist = inf) ───────────────

def compute_psrt_pair(rna_c, lig_c):
    """
    Persistent Stanley-Reisner f-vector and h-vector (Suwayyid & Wei 2025).

    Uses a bipartite distance matrix where intramolecular pairs have distance
    infinity — mathematically excluding intramolecular simplices.

    For the resulting 1-dimensional bipartite complex:
      f₁(r) = |{(i,j): i∈RNA, j∈LIG, d(i,j) ≤ r}|   edge count
      f₀     = n_RNA + n_LIG                            vertex count (constant)
      h₂(r)  = 1 - f₀ + f₁(r)                          from h-polynomial

    These evolving counts under filtration are genuine PSRT invariants not
    captured by standard persistent homology Betti curves.
    """
    rna_c, lig_c = subsample(rna_c), subsample(lig_c)
    n_rna, n_lig = len(rna_c), len(lig_c)
    f0       = float(n_rna + n_lig)
    max_edge = float(n_rna * n_lig) if n_rna > 0 and n_lig > 0 else 1.0

    D_cross = cdist(rna_c, lig_c)

    f1_curve = np.array([float(np.sum(D_cross <= r))
                         for r in FILTRATION_STEPS], dtype=np.float32)
    h2_curve = (1.0 - f0 + f1_curve).astype(np.float32)

    # Scalar summaries of f₁ curve
    f1_max     = float(f1_curve[-1])
    f1_area    = float(np.trapezoid(f1_curve, FILTRATION_STEPS))
    f1_maxrate = float(np.max(np.diff(f1_curve))) if len(f1_curve) > 1 else 0.0
    f1_density = f1_max / max_edge

    stats = np.array([f1_area, f1_max, f1_maxrate, f1_density], dtype=np.float32)
    return np.concatenate([f1_curve, h2_curve, stats])

# ── Feature assembly ──────────────────────────────────────────────────────────

def compute_all_features(rna_es, rna_cs, lig_typed):
    ph_parts, psrt_parts = [], []
    # ES pairs
    for rt in range(N_RNA_ES):
        rc = rna_es.get(rt, np.zeros((0,3)))
        for lt in range(N_LIG_ES):
            lc = lig_typed.get(lt, np.zeros((0,3)))
            if len(rc) == 0 or len(lc) == 0:
                ph_parts.append(np.zeros(PH_PER_PAIR,   dtype=np.float32))
                psrt_parts.append(np.zeros(PSRT_PER_PAIR, dtype=np.float32))
            else:
                ph_parts.append(compute_ph_pair(rc, lc))
                psrt_parts.append(compute_psrt_pair(rc, lc))
    # CS pairs
    for rt in range(N_RNA_CS):
        rc = rna_cs.get(rt, np.zeros((0,3)))
        for lt in range(N_LIG_ES):
            lc = lig_typed.get(lt, np.zeros((0,3)))
            if len(rc) == 0 or len(lc) == 0:
                ph_parts.append(np.zeros(PH_PER_PAIR,   dtype=np.float32))
                psrt_parts.append(np.zeros(PSRT_PER_PAIR, dtype=np.float32))
            else:
                ph_parts.append(compute_ph_pair(rc, lc))
                psrt_parts.append(compute_psrt_pair(rc, lc))
    return np.concatenate(ph_parts), np.concatenate(psrt_parts)

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    df   = pd.read_csv("/home/stalin/Desktop/CAML/data/dataset_clean.csv")
    n    = len(df)
    y    = df['pKd'].values.astype(np.float32)
    pdbs = df['pdb'].values

    print(f"Processing {n} complexes")
    print(f"  ES: {N_ES_PAIRS} pairs  CS: {N_CS_PAIRS} pairs  Total: {N_PAIRS}")
    print(f"  PH per complex:   {PH_DIM}-dim")
    print(f"  PSRT per complex: {PSRT_DIM}-dim\n")

    X_ph   = np.zeros((n, PH_DIM),   dtype=np.float32)
    X_psrt = np.zeros((n, PSRT_DIM), dtype=np.float32)
    failed = []

    for i, row in df.iterrows():
        pdb = row['pdb']
        t0  = time.time()
        try:
            rna_es, rna_cs = parse_rna_pocket(row['pocket_file'])
            lig_typed      = parse_ligand_sdf(row['ligand_sdf'])
            X_ph[i], X_psrt[i] = compute_all_features(rna_es, rna_cs, lig_typed)
            if (i+1) % 10 == 0 or i < 3:
                print(f"  [{i+1:3d}/{n}] {pdb:6s}  {time.time()-t0:.1f}s  "
                      f"nz_ph={np.count_nonzero(X_ph[i])}  "
                      f"nz_psrt={np.count_nonzero(X_psrt[i])}")
        except Exception as e:
            failed.append((pdb, str(e)))
            print(f"  [{i+1:3d}/{n}] {pdb}  FAILED: {e}")

    X_full = np.hstack([X_ph, X_psrt])
    feat_dir = "/home/stalin/Desktop/CAML/features"
    out = os.path.join(feat_dir, "step02_es_cs_features.npz")
    np.savez(out, X=X_ph, X_psrt=X_psrt, X_full=X_full, y=y, pdbs=pdbs)

    print(f"\nSaved {out}")
    print(f"  X (PH):        {X_ph.shape}")
    print(f"  X_psrt (PSRT): {X_psrt.shape}")
    print(f"  X_full:        {X_full.shape}")
    if failed:
        print(f"Failed: {[f[0] for f in failed]}")
    print(f"Done. {n - len(failed)}/{n} succeeded.")

if __name__ == "__main__":
    main()
