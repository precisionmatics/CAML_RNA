"""
Step 02 — Bipartite Persistent Homology features (CAML-RNA)

For each RNA-ligand complex:
  - Parse pocket PDB → RNA heavy atoms with ES and CS type labels
  - Parse ligand SDF → ligand heavy atoms with ES type labels
  - For each (RNA_type, Ligand_type) pair build bipartite VR complex
  - Compute H0 and H1 Betti curves at 24 filtration steps (0.5–12.0 Å)
  - Extract 54 features per pair: 24 β0 + 24 β1 + 6 persistence stats
  - Scheme A (ES+CS): 80 pairs → 4320-dim vector
  - Save features as NPZ
"""

import numpy as np
import pandas as pd
import os, time
from ripser import ripser
from Bio.PDB import PDBParser
from rdkit import Chem

# ── Atom-type dictionaries ────────────────────────────────────────────────────

# ES label for RNA (element → int)
RNA_ES = {'C': 0, 'N': 1, 'O': 2, 'P': 3}

# Backbone atom names
BACKBONE_ATOMS = {"P","OP1","OP2","O5'","O3'","O5*","O3*"}
# Sugar atom names
SUGAR_ATOMS    = {"C1'","C2'","C3'","C4'","C5'","O2'","O4'",
                  "C1*","C2*","C3*","C4*","C5*","O2*","O4*"}
# Purine residues
PURINES        = {"A","DA","G","DG"}
# Pyrimidine residues
PYRIMIDINES    = {"C","DC","U","DU","T","DT"}

def get_cs_label(residue_name, atom_name):
    """Category-specific label: 0=purine_base, 1=pyrimidine_base,
       2=backbone, 3=sugar. Returns None if not classifiable."""
    rn = residue_name.strip().upper()
    an = atom_name.strip()
    if an in BACKBONE_ATOMS:
        return 2
    if an in SUGAR_ATOMS:
        return 3
    if rn in PURINES:
        return 0
    if rn in PYRIMIDINES:
        return 1
    return None  # modified residues, etc.

# ES label for ligand atoms
LIG_ES = {'C':0,'N':1,'O':2,'S':3,'P':4,'F':5,'CL':6,'BR':7,'I':8,'H':9}

FILTRATION_STEPS = np.arange(0.5, 12.5, 0.5)   # 24 steps
MAX_ATOMS = 200                                   # subsample cap per type

# ── PDB parsing ───────────────────────────────────────────────────────────────

def parse_rna_pocket(pocket_pdb):
    """
    Returns dict: {'ES': {0: coords_array, ...}, 'CS': {0: coords_array, ...}}
    Filters to heavy atoms only (no H).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('rna', pocket_pdb)

    es_atoms = {k: [] for k in RNA_ES.values()}   # 0,1,2,3
    cs_atoms = {k: [] for k in range(4)}           # 0,1,2,3

    for model in structure:
        for chain in model:
            for residue in chain:
                res_name = residue.get_resname().strip().upper()
                for atom in residue:
                    elem = atom.element.strip().upper() if atom.element else ''
                    if elem == 'H' or elem == '':
                        continue
                    coord = atom.get_coord()
                    # ES label
                    if elem in RNA_ES:
                        es_atoms[RNA_ES[elem]].append(coord)
                    # CS label
                    cs = get_cs_label(res_name, atom.get_name())
                    if cs is not None and elem in RNA_ES:
                        cs_atoms[cs].append(coord)
        break  # first model only

    return {
        'ES': {k: np.array(v) for k, v in es_atoms.items()},
        'CS': {k: np.array(v) for k, v in cs_atoms.items()},
    }

def parse_ligand_sdf(sdf_path):
    """Returns dict: {es_label: coords_array} for heavy atoms."""
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=True, sanitize=False)
    lig_atoms = {k: [] for k in range(10)}
    for mol in supplier:
        if mol is None:
            continue
        conf = mol.GetConformer()
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol().upper()
            if sym == 'H':
                continue
            label = LIG_ES.get(sym, None)
            if label is None:
                continue
            pos = conf.GetAtomPosition(atom.GetIdx())
            lig_atoms[label].append([pos.x, pos.y, pos.z])
        break  # first conformer
    return {k: np.array(v) for k, v in lig_atoms.items()}

# ── Bipartite PH computation ──────────────────────────────────────────────────

def subsample(arr, n=MAX_ATOMS, seed=42):
    if len(arr) <= n:
        return arr
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(arr), n, replace=False)
    return arr[idx]

def betti_curve(dgm, filtration_steps):
    """Count bars alive at each filtration step."""
    curve = np.zeros(len(filtration_steps), dtype=np.float32)
    for birth, death in dgm:
        if np.isinf(death):
            death = filtration_steps[-1] + 1.0
        alive = (filtration_steps >= birth) & (filtration_steps < death)
        curve += alive.astype(np.float32)
    return curve

def persistence_stats(dgm):
    """Returns [sum_lifetime, max_lifetime, count] for finite bars."""
    finite = [(b, d) for b, d in dgm if not np.isinf(d)]
    if not finite:
        return np.zeros(3, dtype=np.float32)
    lifetimes = np.array([d - b for b, d in finite])
    return np.array([lifetimes.sum(), lifetimes.max(), len(finite)],
                    dtype=np.float32)

def compute_bipartite_ph_pair(rna_coords, lig_coords):
    """
    Compute H0+H1 for the bipartite VR complex built from
    rna_coords ∪ lig_coords (both must be non-empty).
    Returns 54-dim feature vector.
    """
    rna_coords = subsample(rna_coords)
    lig_coords = subsample(lig_coords)

    # Keep only inter-molecular pairs within 12 Å
    from scipy.spatial.distance import cdist
    D = cdist(rna_coords, lig_coords)
    valid_rna = np.where(D.min(axis=1) <= 12.0)[0]
    valid_lig = np.where(D.min(axis=0) <= 12.0)[0]

    if len(valid_rna) == 0 or len(valid_lig) == 0:
        return np.zeros(54, dtype=np.float32)

    rna_coords = rna_coords[valid_rna]
    lig_coords = lig_coords[valid_lig]

    pts = np.vstack([rna_coords, lig_coords])

    try:
        result = ripser(pts, maxdim=1, thresh=12.0)
    except Exception:
        return np.zeros(54, dtype=np.float32)

    dgm0 = result['dgms'][0]
    dgm1 = result['dgms'][1] if len(result['dgms']) > 1 else np.zeros((0, 2))

    b0 = betti_curve(dgm0, FILTRATION_STEPS)
    b1 = betti_curve(dgm1, FILTRATION_STEPS)
    s0 = persistence_stats(dgm0)
    s1 = persistence_stats(dgm1)

    return np.concatenate([b0, b1, s0, s1]).astype(np.float32)

def compute_es_cs_features(rna_typed, lig_typed):
    """
    ES scheme: 4 RNA types × 10 lig types = 40 pairs
    CS scheme: 4 RNA types × 10 lig types = 40 pairs
    Total: 80 pairs × 54 features = 4320-dim
    """
    parts = []
    # ES pairs
    for r_type in range(4):
        rna_c = rna_typed['ES'].get(r_type, np.array([]).reshape(0,3))
        for l_type in range(10):
            lig_c = lig_typed.get(l_type, np.array([]).reshape(0,3))
            if len(rna_c) == 0 or len(lig_c) == 0:
                parts.append(np.zeros(54, dtype=np.float32))
            else:
                parts.append(compute_bipartite_ph_pair(rna_c, lig_c))
    # CS pairs
    for r_type in range(4):
        rna_c = rna_typed['CS'].get(r_type, np.array([]).reshape(0,3))
        for l_type in range(10):
            lig_c = lig_typed.get(l_type, np.array([]).reshape(0,3))
            if len(rna_c) == 0 or len(lig_c) == 0:
                parts.append(np.zeros(54, dtype=np.float32))
            else:
                parts.append(compute_bipartite_ph_pair(rna_c, lig_c))
    return np.concatenate(parts)  # 4320-dim

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    df = pd.read_csv("/home/stalin/Desktop/CAML/data/dataset_clean.csv")
    n = len(df)
    print(f"Processing {n} complexes — ES+CS bipartite PH (4320-dim each)")

    feat_dir = "/home/stalin/Desktop/CAML/features"
    os.makedirs(feat_dir, exist_ok=True)

    X = np.zeros((n, 4320), dtype=np.float32)
    failed = []

    for i, row in df.iterrows():
        pdb = row['pdb']
        t0 = time.time()
        try:
            rna_typed = parse_rna_pocket(row['pocket_file'])
            lig_typed = parse_ligand_sdf(row['ligand_sdf'])
            X[i] = compute_es_cs_features(rna_typed, lig_typed)
            elapsed = time.time() - t0
            if (i+1) % 10 == 0:
                print(f"  [{i+1:3d}/{n}] {pdb}  {elapsed:.1f}s  "
                      f"nonzero={np.count_nonzero(X[i])}")
        except Exception as e:
            failed.append((pdb, str(e)))
            print(f"  [{i+1:3d}/{n}] {pdb}  FAILED: {e}")

    y = df['pKd'].values.astype(np.float32)
    pdbs = df['pdb'].values

    out = os.path.join(feat_dir, "step02_es_cs_features.npz")
    np.savez(out, X=X, y=y, pdbs=pdbs)
    print(f"\nSaved: {out}  shape={X.shape}")
    if failed:
        print(f"Failed ({len(failed)}): {[f[0] for f in failed]}")

if __name__ == "__main__":
    main()
