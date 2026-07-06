"""
Step 01 — Dataset preparation for CAML-RNA
Reads the NA-L dataset, verifies all PDB/SDF files exist, saves clean dataset.
"""

import pandas as pd
import numpy as np
import os

DATA_CSV = "/home/stalin/Desktop/PDFL-RNA/RNA_PDFL/data/affinity/dataset.csv"
NA_L_DIR = "/home/stalin/Desktop/PDFL-RNA/NA-L"
OUT_DIR   = "/home/stalin/Desktop/CAML/data"
os.makedirs(OUT_DIR, exist_ok=True)

df = pd.read_csv(DATA_CSV)
print(f"Loaded {len(df)} entries")

# Verify files exist
missing = []
for _, row in df.iterrows():
    pdb = row['pdb']
    pocket = os.path.join(NA_L_DIR, pdb, f"{pdb}_pocket.pdb")
    ligand = os.path.join(NA_L_DIR, pdb, f"{pdb}_ligand.sdf")
    if not os.path.exists(pocket):
        missing.append((pdb, 'pocket'))
    if not os.path.exists(ligand):
        missing.append((pdb, 'ligand'))

if missing:
    print(f"WARNING: {len(missing)} missing files: {missing[:5]}")
else:
    print("All pocket + ligand files found")

# Add explicit file paths
df['pocket_file'] = df['pdb'].apply(
    lambda p: os.path.join(NA_L_DIR, p, f"{p}_pocket.pdb"))
df['ligand_sdf'] = df['pdb'].apply(
    lambda p: os.path.join(NA_L_DIR, p, f"{p}_ligand.sdf"))
df['ligand_mol2'] = df['pdb'].apply(
    lambda p: os.path.join(NA_L_DIR, p, f"{p}_ligand.mol2"))

# Print pKd stats
print(f"\npKd range: {df['pKd'].min():.2f} – {df['pKd'].max():.2f}")
print(f"pKd mean: {df['pKd'].mean():.2f}, SD: {df['pKd'].std():.2f}")

out_path = os.path.join(OUT_DIR, "dataset_clean.csv")
df.to_csv(out_path, index=False)
print(f"\nSaved: {out_path}  ({len(df)} rows)")
