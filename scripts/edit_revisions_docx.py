"""
Targeted edits to the Revisions DOCX preserving original formatting.
Opens 07_July_2026_Revised_CAML_RNA_manuscript.docx, applies all reviewer
corrections, saves as 08_August_2026_Revised_CAML_RNA_manuscript.docx.
"""

import copy
import re
from docx import Document
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.shared import Pt
import lxml.etree as etree

SRC = "/home/stalin/Desktop/CAML/Revisions/07_July_2026_Revised_CAML_RNA_manuscript.docx"
DST = "/home/stalin/Desktop/CAML/Revisions/08_August_2026_Revised_CAML_RNA_manuscript.docx"

doc = Document(SRC)
paras = doc.paragraphs


# ── helpers ──────────────────────────────────────────────────────────────────

def set_para_text(para, new_text):
    """Replace all runs in para with a single run carrying the same style as run[0]."""
    if not para.runs:
        para.add_run(new_text)
        return
    # Keep format of first run
    first_run = para.runs[0]
    font_name = first_run.font.name
    font_size = first_run.font.size
    bold      = first_run.bold
    italic    = first_run.italic
    # Clear all runs
    for run in para.runs:
        run.text = ""
    # Set text on first run
    first_run.text = new_text
    if font_name:  first_run.font.name = font_name
    if font_size:  first_run.font.size = font_size
    first_run.bold   = bold
    first_run.italic = italic


def replace_in_para(para, old, new):
    """Replace old text with new text across all runs in para."""
    full = para.text
    if old not in full:
        return False
    # Simple: clear runs, set replaced text in first run
    replaced = full.replace(old, new)
    set_para_text(para, replaced)
    return True


def global_replace(old, new):
    """Replace old→new in every paragraph."""
    count = 0
    for p in doc.paragraphs:
        if replace_in_para(p, old, new):
            count += 1
    return count


def insert_para_after(ref_para, text, bold_heading=False):
    """Insert a new paragraph immediately after ref_para."""
    new_p = OxmlElement('w:p')
    new_r = OxmlElement('w:r')
    new_t = OxmlElement('w:t')
    new_t.text = text
    new_t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
    new_r.append(new_t)
    if bold_heading:
        new_rpr = OxmlElement('w:rPr')
        new_b   = OxmlElement('w:b')
        new_rpr.append(new_b)
        new_r.insert(0, new_rpr)
    new_p.append(new_r)
    ref_para._p.addnext(new_p)
    # Return the paragraph object
    for p in doc.paragraphs:
        if p._p is new_p:
            return p
    return None


def add_table_row(table, cells):
    """Append a new row with given cell texts to table."""
    new_row = copy.deepcopy(table.rows[-1]._tr)
    tr_cells = new_row.findall(qn('w:tc'))
    for tc, text in zip(tr_cells, cells):
        for p in tc.findall('.//' + qn('w:p')):
            for r in p.findall(qn('w:r')):
                p.remove(r)
            new_r = OxmlElement('w:r')
            new_t = OxmlElement('w:t')
            new_t.text = text
            new_r.append(new_t)
            p.append(new_r)
            break
    table._tbl.append(new_row)


# ─────────────────────────────────────────────────────────────────────────────
# 1. Global R value fix: 0.7288 → 0.7283
# ─────────────────────────────────────────────────────────────────────────────
n = global_replace("0.7288", "0.7283")
print(f"  R value global fix: {n} occurrences")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Abstract [10] — full rewrite
# ─────────────────────────────────────────────────────────────────────────────
abs_text = (
    "We introduce Commutative Algebra Machine Learning for RNA (CAML-RNA), "
    "a framework that combines bipartite persistent homology (PH) and "
    "persistent Stanley–Reisner theory (PSRT) for predicting RNA–small-molecule "
    "binding affinities. RNA–ligand interactions are represented as bipartite "
    "atom-pair point clouds, enabling element-specific (ES, 36 pairs, 4 RNA × 9 ligand "
    "elements) and category-specific (CS, 36 pairs, 4 RNA structural categories × 9 "
    "ligand elements) Vietoris–Rips filtrations that extract β₀ and β₁ Betti curves "
    "(PH features, 3,888-dim combined), supplemented by persistent f-vectors and "
    "h-vectors from a true bipartite distance filtration (PSRT features, 3,744-dim, "
    "following Suwayyid and Wei). Applied to a curated benchmark of 143 RNA–ligand "
    "complexes spanning seven structurally distinct subtypes, CAML-RNA achieves "
    "Pearson R = 0.7283 (RMSE = 1.09 pKd units) under leave-one-out cross-validation "
    "and R = 0.6255 under nested 10-fold cross-validation, outperforming AffiGrapher "
    "(R = 0.498), RLaffinity (R = 0.559), and RLASIF (R = 0.666, all LOO-CV). "
    "A subtype-aware feature selection strategy achieves R = 0.940 for aptamers "
    "and R = 0.771 for the riboswitch family."
)
set_para_text(paras[10], abs_text)
print("  Abstract rewritten [10]")

# ─────────────────────────────────────────────────────────────────────────────
# 3. Introduction contributions [19] — fix element counts and framework
# ─────────────────────────────────────────────────────────────────────────────
intro_text = (
    "Here, we extend CAML to RNA–ligand binding through three RNA-specific "
    "methodological contributions: (i) bipartite persistent homology for RNA — "
    "element-specific (ES) and category-specific (CS) Vietoris–Rips filtrations "
    "on 36 RNA–ligand heavy-atom element pairs each (4 RNA types × 9 ligand "
    "element types, hydrogen excluded), capturing the multiscale topology of the "
    "binding interface via β₀/β₁ Betti curves (combined PH features, 3,888-dim); "
    "(ii) persistent f/h-vector features (PSRT) — true bipartite distance filtration "
    "(intramolecular distances set to ∞) yielding persistent f₁(r) and h₂(r) curves "
    "as genuine Stanley–Reisner invariants in the sense of Suwayyid and Wei, applied "
    "to all 72 atom-type pairs (3,744-dim); "
    "(iii) subtype-aware feature selection — a principled strategy identifying the "
    "optimal feature modality for each of seven RNA structural subtypes via nested "
    "cross-validation, evaluated under both LOO-CV and 10-fold CV."
)
set_para_text(paras[19], intro_text)
print("  Introduction contributions rewritten [19]")

# ─────────────────────────────────────────────────────────────────────────────
# 4. Benchmark text [27] — R already fixed globally; also fix surrounding text
# ─────────────────────────────────────────────────────────────────────────────
replace_in_para(paras[27], "7288", "7283")  # safety net
print("  Benchmark text checked [27]")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Hochster formula [91] — remove the "− 1"
# ─────────────────────────────────────────────────────────────────────────────
old_hochster = "βᵢ,ᵢ₊₁(k[Δ])  =  ∑ᵂ⊆V, |W|=i+1  (β̃₀(Δᵂ) − 1)"
new_hochster = "βᵢ,ᵢ₊₁(k[Δ])  =  ∑ᵂ⊆V, |W|=i+1  β̃₀(Δᵂ)"
ok = replace_in_para(paras[91], old_hochster, new_hochster)
print(f"  Hochster formula fix [91]: {ok}")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Hochster formula note [94] — clarify reduced vs unreduced
# ─────────────────────────────────────────────────────────────────────────────
old94 = "where βⱼ₋₁(Δᵂ) is the (j−1)-th Betti number of the simplicial homology of Δᵂ, and β̃₀ denotes the reduced 0-th Betti number."
new94 = ("where β̃₀(Δᵂ) = dim_k H̃₀(Δᵂ; k) is the reduced 0th Betti number "
         "(number of connected components minus one), and for j ≥ 2, "
         "β̃ⱼ₋₁(Δᵂ) = βⱼ₋₁(Δᵂ) since reduced and unreduced homology coincide in positive dimensions.")
ok = replace_in_para(paras[94], old94, new94)
print(f"  Hochster note fix [94]: {ok}")

# ─────────────────────────────────────────────────────────────────────────────
# 7. ES section [115] — remove H, fix element count
# ─────────────────────────────────────────────────────────────────────────────
old115 = ("RNA heavy atoms are typed by chemical element: Ξᴿᴺᴬ = {C, N, O, P} "
          "(four types, covering all canonical RNA heavy atoms). "
          "Ligand heavy atoms are typed by element: Ξᴸᴵᶢ = {C, N, O, S, P, F, Cl, Br, H, I} "
          "(ten types, including polar hydrogens retained in some structures). "
          "For each element pair (α, β) ∈")
new115 = ("RNA heavy atoms are typed by chemical element: Ξᴿᴺᴬ = {C, N, O, P} "
          "(four types, covering all non-hydrogen atoms in canonical RNA). "
          "Ligand heavy atoms are typed by element: Ξᴸᴵᶢ = {C, N, O, S, P, F, Cl, Br, I} "
          "(nine types; hydrogen is excluded because RDKit strips polar hydrogens and "
          "H-only channels contain no signal). "
          "For each element pair (α, β) ∈")
ok = replace_in_para(paras[115], old115, new115)
if not ok:
    # Try partial replacement
    ok = replace_in_para(paras[115], "H, I} (ten types, including polar hydrogens retained in some structures)",
                         "I} (nine types; hydrogen is excluded because RDKit strips polar hydrogens and H-only channels contain no signal)")
    ok2 = replace_in_para(paras[115], "40 element", "36 element")
    print(f"  ES section partial fix [115]: {ok}/{ok2}")
else:
    print(f"  ES section fix [115]: {ok}")

# ─────────────────────────────────────────────────────────────────────────────
# 8. ES dimensions [117] — fix filtration levels and dimensions
# ─────────────────────────────────────────────────────────────────────────────
old117 = ("Ripser [32] computes the persistent Betti numbers β0 and β1 across "
          "54 equally spaced filtration levels r ∈ [0, 12] Å. "
          "Concatenating the Betti curves over all 40 element pairs yields the "
          "ES feature vector of dimension 40 × 2 × 54 = 4,320.")
new117 = ("Ripser [32] computes β₀ and β₁ Betti curves across 24 equally spaced "
          "filtration levels r ∈ [0.5, 12.0] Å (step 0.5 Å), together with three "
          "persistence statistics per dimension (total persistence, maximum persistence, "
          "number of finite bars), giving 2 × 24 + 6 = 54 features per element pair. "
          "Concatenating over all 36 element pairs yields the ES feature vector of "
          "dimension 36 × 54 = 1,944.")
ok = replace_in_para(paras[117], old117, new117)
if not ok:
    # Partial replacements
    replace_in_para(paras[117], "54 equally spaced filtration levels r ∈ [0, 12] Å",
                    "24 equally spaced filtration levels r ∈ [0.5, 12.0] Å (step 0.5 Å)")
    replace_in_para(paras[117], "40 element pairs", "36 element pairs")
    replace_in_para(paras[117], "40 × 2 × 54 = 4,320", "36 × 54 = 1,944")
    replace_in_para(paras[117], "40 × 2 × 54 = 4320", "36 × 54 = 1,944")
print(f"  ES dimensions fix [117]: {ok}")

# ─────────────────────────────────────────────────────────────────────────────
# 9. CS section — rewrite paragraphs 119–123
# ─────────────────────────────────────────────────────────────────────────────
new_cs_intro = (
    "To complement element-level resolution, four RNA atom categories are defined "
    "based on structural role:"
)
ok = replace_in_para(paras[119], paras[119].text, new_cs_intro)
print(f"  CS intro rewrite [119]: {ok}")

# Equation paragraph [120] — keep as placeholder, update formula text
new_cs_eq = "𝒞ᴿᴺᴬ  =  { Bᵖᵘʳ,  Bᵖʸʳ,  Bᵇᵏᵇ,  Bˢᵘᵍ }	(22)"
set_para_text(paras[120], new_cs_eq)
print("  CS RNA categories updated [120]")

# Para [121] — description of RNA categories and LIG typing
new_cs_desc = (
    "where Bᵖᵘʳ comprises purine (adenine and guanine) base atoms, Bᵖʸʳ comprises "
    "pyrimidine (cytosine, uracil, thymine) base atoms, Bᵇᵏᵇ comprises phosphate "
    "backbone atoms (P, OP1, OP2, O5′, O3′), and Bˢᵘᵍ comprises ribose sugar atoms "
    "(C1′–C5′, O2′, O4′). Ligand atoms retain the same nine-element typing as in the "
    "ES scheme (Ξᴸᴵᶢ). For each of the 4 × 9 = 36 category–element pairs, a "
    "category-specific complex is constructed with the same Vietoris–Rips filtration "
    "and feature extraction as the ES scheme, yielding 36 × 54 = 1,944 CS features. "
    "The combined ES+CS PH feature vector has dimension 1,944 + 1,944 = 3,888."
)
set_para_text(paras[121], new_cs_desc)
print("  CS description rewritten [121]")

# Paras [122] and [123] — old LIG pharmacophoric categories; replace with PSRT subsection header
new_p122 = "3.4.3  Persistent f/h-Vector Features (PSRT)"
set_para_text(paras[122], new_p122)
print("  Para [122] set to PSRT subsection header")

new_p123 = (
    "In addition to the PH Betti curves, CAML-RNA computes persistent f-vector and "
    "h-vector features following Suwayyid and Wei (PSRT). These are computed using the "
    "bipartite distance filtration, which sets all intramolecular distances to ∞, "
    "ensuring that only RNA–ligand intermolecular edges form in the filtration complex. "
    "For each of the 72 atom-type pairs (36 ES + 36 CS), let nα and nβ denote the "
    "number of RNA and ligand atoms of the selected types. At filtration radius r, "
    "the number of intermolecular edges is:"
)
set_para_text(paras[123], new_p123)
print("  Para [123] set to PSRT intro text")

# Insert new paragraphs after [123] for PSRT equations and description
p123 = paras[123]
texts_to_insert = [
    ("f₁(r)  =  |{(u, v) | u ∈ Vα,  v ∈ Vβ,  d(u, v) ≤ r}|	(23)", False),
    ("The vertex count f₀ = nα + nβ is constant. The corresponding h-vector component at radius r is:", False),
    ("h₂(r)  =  1  −  f₀  +  f₁(r)	(24)", False),
    ("derived from the h-polynomial identity for dimension d = 2. "
     "Four scalar summaries of the f₁ curve (area under curve, maximum, maximum rate of increase, "
     "and fill density f₁ᵐᵃˣ / (nα nβ)) are appended. "
     "Each pair contributes 24 + 24 + 4 = 52 PSRT features (f₁ curve: 24, h₂ curve: 24, scalars: 4), "
     "giving 72 × 52 = 3,744 PSRT features total. "
     "The full CAML-RNA feature vector (PH + PSRT) has dimension 3,888 + 3,744 = 7,632.", False),
]
# Insert in reverse order so they end up in correct sequence
current = p123
for text, bold in texts_to_insert:
    new_para = insert_para_after(current, text, bold)
    if new_para:
        current = new_para
    else:
        current = p123  # fallback
print("  Inserted PSRT equation paragraphs after [123]")

# ─────────────────────────────────────────────────────────────────────────────
# 10. CPF section [126] — fix element-type pair counts and cutoffs
#     Note: paragraph indices shift by +3 after the insertions above!
#     We'll search by content instead
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "Contact Pair Features" in p.text and "CPF" in p.text and "660" in p.text:
        new_cpf = (
            "Contact Pair Features (CPF, 660-dim). RNA atoms are assigned to one of "
            "11 chromatic (element + structural-role) types (C/N/O/P in backbone, base, "
            "or sugar context, plus a catch-all), and ligand atoms to one of 10 element "
            "types (C, N, O, S, P, F, Cl, Br, I, other). For each of the 11 × 10 = 110 "
            "RNA–ligand type pairs and each of six distance cutoffs (3.5, 4.0, 4.5, 5.0, "
            "5.5, and 6.0 Å), we record the contact count, giving a 11 × 10 × 6 = "
            "660-dimensional CPF vector. CPF provides high-resolution chemical "
            "fingerprinting of the binding interface complementary to the topological "
            "Betti curve descriptors."
        )
        set_para_text(p, new_cpf)
        print("  CPF section fixed")
        break

# ─────────────────────────────────────────────────────────────────────────────
# 11. SR ideal direction fix in interpretability section [149]
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "acquires minimal generators" in p.text:
        old_sr = ("In the Stanley–Reisner terms, the β0 decay rate corresponds to the "
                  "rate at which the Stanley–Reisner ideal I(Δʳ) acquires minimal generators "
                  "as r increases, directly encoding the contact density.")
        new_sr = ("In the Stanley–Reisner terms, the β0 decay rate corresponds to the "
                  "rate at which the Stanley–Reisner ideal I(Δʳ) loses minimal generators "
                  "as r increases: each edge added to the complex eliminates one degree-2 "
                  "non-face from I(Δʳ), so rapid β₀ decay reflects a rapidly shrinking SR "
                  "ideal, directly encoding the contact density.")
        ok = replace_in_para(p, old_sr, new_sr)
        if not ok:
            replace_in_para(p, "acquires minimal generators as r increases, directly encoding",
                            "loses minimal generators as r increases: each edge added eliminates "
                            "one degree-2 non-face, directly encoding")
        print("  SR ideal direction fixed")
        break

# ─────────────────────────────────────────────────────────────────────────────
# 12. Ablation section [49] — rewrite for modality ablation
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "Figure 5 traces CAML-RNA performance across twelve key development milestones" in p.text:
        new_ablation = (
            "To assess the contribution of each feature modality, we conducted a "
            "systematic ablation study using nested 10-fold cross-validation (outer "
            "10-fold, inner 5-fold for hyperparameter tuning) with a uniform SVR-RBF "
            "pipeline on all 143 complexes. PH alone achieves R₁₀ = 0.613 and "
            "R_LOO = 0.623, confirming bipartite persistent homology as the primary "
            "topological driver. PSRT f/h-vectors alone achieve R₁₀ = 0.535, demonstrating "
            "independent signal; combining PH+PSRT does not improve over PH alone "
            "because the f₁(r) curve is highly correlated with the β₀ Betti curve "
            "(both count interface atom contacts at each radius). The strongest global "
            "combinations are PH+PSRT+RNA-FM (R₁₀ = 0.628) and PH+PSRT+CPF+RNA-FM "
            "(R₁₀ = 0.626, R_LOO = 0.632); adding all modalities drops slightly to "
            "R₁₀ = 0.609, showing that Morgan fingerprints add noise at n = 143. "
            "Subtype-specific models in step 33 push the final LOO performance to "
            "R = 0.7283. Figure 5 and Table 4 present these results."
        )
        set_para_text(p, new_ablation)
        print("  Ablation text rewritten")
        break

# ─────────────────────────────────────────────────────────────────────────────
# 13. Figure 5 caption — update
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "Figure 5" in p.text and ("ablation" in p.text.lower() or "Ablation" in p.text) and "Step" in p.text:
        new_fig5 = (
            "Figure 5. Feature modality ablation study. Pearson R (10-fold nested CV) "
            "for 13 feature combinations ranging from single modalities (PH, PSRT, CPF, "
            "RNA-FM, Morgan, Physicochemical) to multi-modality ensembles. "
            "Red bar: proposed CAML-RNA core combination (PH+PSRT+CPF+RNA-FM, R = 0.626). "
            "Blue bars: other combinations. The dashed line marks the final CAML-RNA "
            "LOO-CV performance (R = 0.7283) achieved by subtype-specific models."
        )
        set_para_text(p, new_fig5)
        print("  Figure 5 caption updated")
        break

# ─────────────────────────────────────────────────────────────────────────────
# 14. Table 1 — benchmark table updates
# ─────────────────────────────────────────────────────────────────────────────
t1 = doc.tables[0]
rows = t1.rows

def set_cell(row, col_idx, text):
    cell = row.cells[col_idx]
    for p in cell.paragraphs:
        for r in p.runs:
            r.text = ""
        if p.runs:
            p.runs[0].text = text
        else:
            p.add_run(text)
        break

# Row 1: AffiGrapher — fix n and protocol
set_cell(rows[1], 2, "143")
set_cell(rows[1], 3, "LOO-CV")
print("  Table 1 Row 1 (AffiGrapher) fixed")

# Row 2: RLaffinity — fix n and protocol
set_cell(rows[2], 2, "117")
set_cell(rows[2], 3, "10-fold CV")
print("  Table 1 Row 2 (RLaffinity) fixed")

# Row 3: RLASIF — fix n and protocol
set_cell(rows[3], 2, "143")
set_cell(rows[3], 3, "LOO-CV")
print("  Table 1 Row 3 (RLASIF) fixed")

# Row 4: CAML-RNA — R already fixed globally; fix RMSE
set_cell(rows[4], 4, "0.7283")
set_cell(rows[4], 5, "1.09")
print("  Table 1 Row 4 (CAML-RNA) fixed")

# Row 6: Replace DeepRSMA with EMMPTNet
set_cell(rows[6], 0, "EMMPTNet [new]")
set_cell(rows[6], 1, "NA-L (NL2020)")
set_cell(rows[6], 2, "143")
set_cell(rows[6], 3, "10-fold CV")
set_cell(rows[6], 4, "0.773")
set_cell(rows[6], 5, "N/A")
print("  Table 1 Row 6 replaced with EMMPTNet")

# Add DeepMIF row
add_table_row(t1, ["DeepMIF [new]", "R-SIM", "1439", "10-fold CV", "0.796", "N/A"])
print("  Table 1 DeepMIF row added")

# ─────────────────────────────────────────────────────────────────────────────
# 15. Table 2 — subtype table: Overall R value (already fixed globally)
# ─────────────────────────────────────────────────────────────────────────────
t2 = doc.tables[1]
for row in t2.rows:
    if "Overall" in row.cells[0].text or "overall" in row.cells[0].text:
        set_cell(row, 3, "0.7283")
        print("  Table 2 Overall R fixed")

# ─────────────────────────────────────────────────────────────────────────────
# 16. Add Table 4 (Ablation) — insert after the ablation text
# ─────────────────────────────────────────────────────────────────────────────
# Find the paragraph that mentions "Table 4" or add ablation table
# We'll insert a table after the conclusions section marker
# First, add a reference paragraph and then the table

# Find end of document body to append ablation table
ablation_data = [
    ["Feature Combination", "Dim", "R₁₀", "R_LOO", "95% CI"],
    ["PH (Betti curves)", "3,888", "0.613", "0.623", "[0.488, 0.709]"],
    ["PSRT (f/h-vectors)", "3,744", "0.535", "0.525", "[0.405, 0.641]"],
    ["PH+PSRT", "7,632", "0.547", "0.554", "[0.414, 0.656]"],
    ["CPF", "660", "0.300", "0.450", "[0.133, 0.444]"],
    ["RNA-FM", "640", "0.334", "0.408", "[0.168, 0.479]"],
    ["Morgan (ECFP4)", "1,024", "0.518", "0.369", "[0.391, 0.619]"],
    ["Physicochemical", "17", "0.456", "0.411", "[0.313, 0.575]"],
    ["PH+CPF", "4,548", "0.592", "0.623", "[0.467, 0.697]"],
    ["PSRT+CPF", "4,404", "0.536", "0.526", "[0.405, 0.640]"],
    ["PH+PSRT+CPF", "8,292", "0.552", "0.554", "[0.420, 0.667]"],
    ["PH+PSRT+RNA-FM", "8,272", "0.628", "0.633", "[0.515, 0.718]"],
    ["PH+PSRT+CPF+RNA-FM", "8,932", "0.626", "0.632", "[0.512, 0.716]"],
    ["All modalities", "9,973", "0.609", "0.605", "[0.494, 0.702]"],
    ["CAML-RNA (step 33, LOO)", "—", "—", "0.7283", "—"],
]

# Find ablation section paragraph to insert table caption after it
for p in doc.paragraphs:
    if "Feature modality ablation study" in p.text and "Figure 5" in p.text:
        # Insert table caption after this
        cap_para = insert_para_after(p,
            "Table 4. Ablation study: contribution of each feature modality "
            "under nested 10-fold CV. R₁₀: 10-fold CV Pearson R; "
            "R_LOO: LOO-CV Pearson R; 95% CI from 10,000 bootstrap resamples. "
            "PH: bipartite Betti curves; PSRT: persistent f/h-vectors; "
            "CPF: contact pair features (660-dim); RNA-FM: foundation model (640-dim).")
        print("  Table 4 caption inserted")
        break

# Append ablation table at end of document body
from docx.oxml import OxmlElement as _OE
def make_ablation_table(doc, data):
    from docx.shared import Inches, Pt
    table = doc.add_table(rows=len(data), cols=len(data[0]))
    table.style = 'Table Grid'
    for ri, row_data in enumerate(data):
        for ci, cell_text in enumerate(row_data):
            cell = table.rows[ri].cells[ci]
            cell.text = cell_text
            if ri == 0:
                for run in cell.paragraphs[0].runs:
                    run.bold = True
            if ri == len(data) - 1 and ci != 0:
                for run in cell.paragraphs[0].runs:
                    run.bold = True
    return table

# Add a section break paragraph then the table
doc.add_paragraph("")
doc.add_paragraph("Table 4 (Ablation Study)")
make_ablation_table(doc, ablation_data)
print("  Table 4 (Ablation) added at end of document")

# ─────────────────────────────────────────────────────────────────────────────
# 17. Save
# ─────────────────────────────────────────────────────────────────────────────
doc.save(DST)
print(f"\nSaved: {DST}")
