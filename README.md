# CAML-RNA: Commutative Algebra Machine Learning for RNA–Ligand Binding Affinity Prediction

Reference implementation for the manuscript *"CAML-RNA: Commutative Algebra Machine Learning for RNA–Ligand Binding Affinity Prediction"* (Arulsamy, Krishna, Kumar & Kumar). This repository contains all code used to produce the results reported in the paper (Pearson **R = 0.7283** under leave-one-out cross-validation on 143 RNA–ligand complexes).

The method combines **bipartite persistent homology (PH)** with **persistent Stanley–Reisner theory (PSRT)** — the two commutative-algebra / topological descriptors introduced in the paper — with support-vector regression and per-subtype override models.

---

## 1. Where the commutative-algebra features live

Reviewer note: the commutative-algebra (PSRT) features are implemented alongside the persistent-homology features in a single script. To make this unambiguous:

| Feature | File | Function | Lines |
|---|---|---|---|
| **Bipartite persistent homology** (β₀, β₁ Betti curves via Ripser) | `scripts/step02_bipartite_ph.py` | `compute_ph_pair` | 133–149 |
| **Persistent Stanley–Reisner f-vector and h-vector** (PSRT features, following Suwayyid & Wei 2025) | `scripts/step02_bipartite_ph.py` | `compute_psrt_pair` | 153–186 |

PSRT features are computed on the *true bipartite* distance matrix (intramolecular distances set to infinity so that only cross-part simplices are included), producing:

- `f₁(r)` — number of RNA–ligand edges with distance ≤ r, evaluated at 24 filtration radii r ∈ {0.5, 1.0, …, 12.0} Å;
- `h₂(r) = 1 − f₀ + f₁(r)` — derived from the h-polynomial of the 1-dimensional bipartite simplicial complex;
- four scalar summaries (area under f₁, final f₁, maximum first difference, density = f₁/max edges).

Repeated across the 36 element-specific and 36 category-specific atom-pair channels, this yields the 3,744-dimensional PSRT feature vector reported in the manuscript.

## 2. Pipeline overview

The pipeline is organised as sequential, self-contained scripts. Each script writes NumPy artifacts to `features/` or `models/` and CSV summaries to `results/`.

```
step01_data_prep.py         Load 143 curated NA-L complexes, harmonise pKd.
step02_bipartite_ph.py      Compute PH Betti curves + PSRT f/h-vectors (3,888 + 3,744 dims).
step03_chromatic_weighted.py  Chromatic / weighted PH variants (ablation only).
step04_ml_training.py       Baseline: 8 regressors × 4 feature sets, nested 10-fold CV.
step05_combined_features.py Feature concatenation + PCA; ES+CS+WCh SVR-RBF baseline.
step06_subtype_models.py    Per-subtype SVR-RBF training.
step07_hybrid_riboswitch.py Riboswitch-specific hybrid features.
step08_ligand_rna_features.py  Explicit ligand + RNA descriptor concatenation.
step09_rnafm_features.py    RNA-FM (640-dim) foundation-model embeddings.
step10_mkl.py               Multiple-kernel-learning hybrid.
step11_contact_features.py  Contact pair features (CPF, 660-dim).
step12_ensemble_subtype.py  Ridge stacking with one-hot subtype indicators.
step13_riboswitch_boost.py  Riboswitch-subclass OHE + MLP.
step14_ligand_desc_aug.py   Physicochemical descriptor augmentation.
step15_aif.py               Shell-binned atomic-interaction fingerprint.
step16_morgan_riboswitch.py Morgan ECFP4 fingerprint for aptamer subclass.
step20_smiles_transformer.py GPT2-ZINC SMILES embeddings (ablation only).
step21_best_subtype_oracle.py Per-subtype override baseline.
step26_riboswitch_fm_cpf.py RNA-FM/CPF LOO overrides for FMN_FAD and TPP.
step33_final_overrides.py   >>> Produces the reported R = 0.7283 <<<
step_10fold_cv.py           Nested 10-fold CV that produces the reported R₁₀ = 0.6255.
step_ablation.py            Modality ablation (Table 4, Figure 5).
```

The paper's headline result comes from `step33_final_overrides.py`; the ablation table comes from `step_ablation.py` and `step_10fold_cv.py`.

## 3. Reproducing the reported results

### 3.1 Environment

```bash
python >= 3.9
pip install numpy scipy pandas scikit-learn ripser biopython rdkit-pypi xgboost
# RNA-FM (optional, for step09): https://github.com/ml4bio/RNA-FM
```

### 3.2 Data

The 143-complex dataset is derived from the [NA-L nucleic acid–ligand database](https://www.pdbbind.org.cn/). The curated CSV (`data/dataset_clean.csv`) and per-complex PDB/SDF files are not redistributed here for licensing reasons; contact the corresponding author for access, or regenerate from the NA-L release using `step01_data_prep.py`.

### 3.3 Feature extraction

```bash
python scripts/step02_bipartite_ph.py    # writes features/step02_es_cs_features.npz
python scripts/step09_rnafm_features.py  # RNA-FM embeddings
python scripts/step11_contact_features.py  # CPF
```

### 3.4 LOO-CV result (R = 0.7283)

```bash
python scripts/step04_ml_training.py     # baseline (R = 0.5925)
python scripts/step21_best_subtype_oracle.py
python scripts/step26_riboswitch_fm_cpf.py
python scripts/step33_final_overrides.py # >>> reported R = 0.7283 <<<
```

### 3.5 Nested 10-fold ablation (Table 4, Figure 5)

```bash
python scripts/step_ablation.py
python scripts/step_10fold_cv.py
```

## 4. What is *not* in this repository

- Any persistent flag Laplacian / PDFL code. PDFL is a separate line of work by other authors and is neither used nor claimed by this paper. Earlier exploratory scripts that combined CAML-RNA out-of-fold predictions with an external PDFL model have been removed from this repository to avoid confusion; they were not part of the reported results.
- Trained model weights for third-party models (RNA-FM, GPT2-ZINC) — please obtain from the original authors.

## 5. Use of AI assistants

Consistent with journal disclosure policy: generative-AI assistants (ChatGPT, Claude) were used for grammar/style checking of the manuscript and for coding assistance during script development. All experimental design, computational results, figures, and scientific interpretations were produced and verified by the human authors. No AI-generated text or code was accepted without inspection and validation against source data.

## 6. Citation

If you use this code, please cite:

> Arulsamy S., Krishna Y., Kumar R., Kumar V. *CAML-RNA: Commutative Algebra Machine Learning for RNA–Ligand Binding Affinity Prediction.* (2026, submitted).

## 7. Contact

Corresponding author: Dr. Vanktesh Kumar (Vankteshkumar555@hotmail.com).
