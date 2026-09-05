# CAML-RNA

Reference implementation for *CAML-RNA: Commutative Algebra Machine Learning for RNA–Ligand Binding Affinity Prediction* (Arulsamy, Krishna, Kumar & Kumar).

The method combines bipartite persistent homology and persistent Stanley–Reisner theory with support-vector regression and per-subtype override models on a benchmark of 143 RNA–ligand complexes (Pearson **R = 0.7283**, LOO-CV).

## Requirements

```
python >= 3.9
numpy scipy pandas scikit-learn xgboost
ripser biopython rdkit-pypi
```

Optional (for step09): RNA-FM — https://github.com/ml4bio/RNA-FM

## Data

Curated CSV (`data/dataset_clean.csv`) and per-complex PDB/SDF files are derived from the [NA-L nucleic acid–ligand database](https://www.pdbbind.org.cn/). Contact the corresponding author for access.

## Pipeline

Scripts are run sequentially. Each writes NumPy arrays to `features/` or `models/` and CSV summaries to `results/`.

**Feature extraction**
- `step01_data_prep.py` — load 143 complexes, harmonise pKd
- `step02_bipartite_ph.py` — persistent homology + Stanley–Reisner f/h-vectors
- `step03_chromatic_weighted.py` — chromatic / weighted PH
- `step09_rnafm_features.py` — RNA-FM embeddings
- `step11_contact_features.py` — contact pair features
- `step14_ligand_desc_aug.py` — physicochemical descriptors
- `step16_morgan_riboswitch.py` — Morgan ECFP4 fingerprints
- `step20_smiles_transformer.py` — GPT2-ZINC SMILES embeddings

**Training and per-subtype models**
- `step04_ml_training.py` — baseline (R = 0.5925)
- `step05_combined_features.py` — feature concatenation + PCA
- `step06_subtype_models.py` — per-subtype SVR-RBF
- `step07_hybrid_riboswitch.py`, `step08_ligand_rna_features.py`, `step10_mkl.py`, `step12_ensemble_subtype.py`, `step13_riboswitch_boost.py`, `step15_aif.py` — feature-modality studies
- `step21_best_subtype_oracle.py` — per-subtype baseline
- `step26_riboswitch_fm_cpf.py` — RNA-FM / CPF LOO overrides
- `step33_final_overrides.py` — final result (R = 0.7283)

**Cross-validation and ablation**
- `step_10fold_cv.py` — nested 10-fold CV (R = 0.6255, Table 4)
- `step_ablation.py` — feature-modality ablation (Figure 5)

## Reproducing R = 0.7283

```bash
python scripts/step01_data_prep.py
python scripts/step02_bipartite_ph.py
python scripts/step09_rnafm_features.py
python scripts/step11_contact_features.py
python scripts/step04_ml_training.py
python scripts/step21_best_subtype_oracle.py
python scripts/step26_riboswitch_fm_cpf.py
python scripts/step33_final_overrides.py
```
