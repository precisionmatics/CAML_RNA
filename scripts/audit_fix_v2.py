#!/usr/bin/env python3
"""
audit_fix_v2.py
Apply all post-audit corrections from full reviewer comment review.

Fixes:
  1. Table 1 Row 4: Add CAML-RNA 10-fold CV result (R=0.6255) — Reviewer 1 Cmt 2
  2. Para body: DeepRSMA R=0.784 → R=0.796 (matches Table 1) — internal inconsistency
  3. Para [026] footnote: "RSAPred and DeepMIF" → "RSAPred and DeepRSMA" — naming fix
  4. Para [061]: Clarify bipartite VR claim (intramolecular not truly excluded for PH) — Reviewer 2 Cmt 1
  5. Para [155]: Soften β₁ = degree-3 generator (too strong) — Reviewer 2 Cmt 4
  6. Methods §3.6: Add step numbering explanation — Reviewer 1 Cmt 4
"""

from docx import Document
import re

DOCX_IN  = "Revisions/11_August_2026_Revised_CAML_RNA_manuscript.docx"
DOCX_OUT = "Revisions/12_August_2026_Revised_CAML_RNA_manuscript.docx"

doc = Document(DOCX_IN)
paras = list(doc.paragraphs)
ref_start = next(i for i, p in enumerate(paras)
                 if re.match(r'^[Rr]eferences?\s*$', p.text.strip()))
body_paras = paras[:ref_start]

def set_para_text(p, new_text):
    if not p.runs:
        return
    p.runs[0].text = new_text
    for run in p.runs[1:]:
        run.text = ''

def replace_in_para(p, old, new):
    """Try run-by-run replacement first; fallback to full-para replacement."""
    for run in p.runs:
        if old in run.text:
            run.text = run.text.replace(old, new)
            return True
    full = p.text
    if old in full:
        set_para_text(p, full.replace(old, new))
        return True
    return False

def find_para_containing(fragments, start=0, end=None):
    """Return (index, para) for first para containing ALL fragments."""
    end = end or ref_start
    for i, p in enumerate(paras[start:end], start):
        if all(f in p.text for f in fragments):
            return i, p
    return None, None

print("=== audit_fix_v2.py ===\n")

# ═══════════════════════════════════════════════════════════════════════════════
# FIX 1 — Table 1: Add CAML-RNA 10-fold CV result (Reviewer 1, Comment 2)
# ═══════════════════════════════════════════════════════════════════════════════
print("Fix 1: Table 1 — add CAML-RNA 10-fold CV result...")
tbl = doc.tables[0]
row4 = tbl.rows[4]  # CAML-RNA row

proto_cell = row4.cells[3]
r_cell     = row4.cells[4]

proto_old = proto_cell.text.strip()
r_old     = r_cell.text.strip()

if '10-fold' not in proto_old:
    for p in proto_cell.paragraphs:
        for run in p.runs:
            if 'LOO-CV' in run.text:
                run.text = run.text.replace('LOO-CV', 'LOO-CV; 10-fold CV')
    print(f"  Protocol: '{proto_old}' → '{proto_cell.text.strip()}'")

if '0.6255' not in r_old:
    for p in r_cell.paragraphs:
        for run in p.runs:
            if '0.7283' in run.text:
                run.text = run.text.replace('0.7283', '0.7283 (LOO); 0.6255 (10-fold)')
    print(f"  Pearson R: '{r_old}' → '{r_cell.text.strip()}'")

# ═══════════════════════════════════════════════════════════════════════════════
# FIX 2 — Body text: DeepRSMA R=0.784 → R=0.796 (internal inconsistency)
# ═══════════════════════════════════════════════════════════════════════════════
print("\nFix 2: Fix DeepRSMA R value 0.784 → 0.796 in body text...")
for i, p in enumerate(body_paras):
    if '0.784' in p.text:
        changed = replace_in_para(p, 'R = 0.784 on R-SIM', 'R = 0.796 on R-SIM')
        if not changed:
            changed = replace_in_para(p, '0.784', '0.796')
        if changed:
            print(f"  Para [{i}]: fixed 0.784 → 0.796 | {p.text[:80]}")

# ═══════════════════════════════════════════════════════════════════════════════
# FIX 3 — Table 1 footnote: fix "RSAPred and DeepMIF" naming
# ═══════════════════════════════════════════════════════════════════════════════
print("\nFix 3: Fix Table 1 footnote method names...")
for i, p in enumerate(body_paras):
    if 'RSAPred' in p.text and 'DeepMIF' in p.text and 'R-SIM database' in p.text:
        changed = replace_in_para(
            p,
            'RSAPred and DeepMIF were trained and evaluated on the R-SIM database; DeepMIF was evaluated on NA-L under 10-fold CV',
            'RSAPred and DeepRSMA (Huang Z et al.) were trained and evaluated on the R-SIM database; DeepMIF (Song J et al., JCIM 2026) was evaluated on NA-L under 10-fold CV'
        )
        if changed:
            print(f"  Para [{i}]: footnote method names corrected")
        else:
            # Try partial
            changed1 = replace_in_para(p, 'RSAPred and DeepMIF were trained', 'RSAPred and DeepRSMA (Huang Z et al.) were trained')
            changed2 = replace_in_para(p, '; DeepMIF was evaluated on NA-L under 10-fold CV',
                                           '; DeepMIF (Song J et al., JCIM 2026) was evaluated on NA-L under 10-fold CV')
            if changed1 or changed2:
                print(f"  Para [{i}]: footnote partially corrected")

# Also fix in table footnote cell if it exists
for ti, tbl_i in enumerate(doc.tables):
    for ri, row in enumerate(tbl_i.rows):
        for ci, cell in enumerate(row.cells):
            ct = cell.text
            if 'RSAPred and DeepMIF were trained' in ct:
                for p in cell.paragraphs:
                    replace_in_para(p,
                        'RSAPred and DeepMIF were trained and evaluated on the R-SIM database; DeepMIF was evaluated on NA-L under 10-fold CV',
                        'RSAPred and DeepRSMA (Huang Z et al.) were trained and evaluated on the R-SIM database; DeepMIF (Song J et al., JCIM 2026) was evaluated on NA-L under 10-fold CV'
                    )
                print(f"  Table {ti} Row {ri} Cell {ci}: footnote corrected")

# ═══════════════════════════════════════════════════════════════════════════════
# FIX 4 — Para [061]: Clarify bipartite VR claim (Reviewer 2, Comment 1)
# ═══════════════════════════════════════════════════════════════════════════════
print("\nFix 4: Clarify bipartite VR claim (Para ~[061])...")

NEW_BIPARTITE_SUFFIX = (
    "This construction partitions the atom set into two molecular components. "
    "For the PSRT f/h-vector features (Section 3.4.3), a strict bipartite distance "
    "filtration is applied in which all intramolecular distances are set to ∞, "
    "so that only intermolecular RNA–ligand edges enter the filtration. "
    "For the PH Betti curve features (Sections 3.4.1–3.4.2), element-specific "
    "Vietoris–Rips filtrations are computed over the union of element-selected "
    "RNA and ligand atom subsets; intramolecular edges between same-element atoms "
    "within each component may form within the VR complex, but their contribution to "
    "the Betti curves is dominated by the intermolecular signal at the filtration "
    "levels relevant to molecular recognition (r ≥ 3 Å), because the "
    "element-specific channel design separates RNA and ligand element types into "
    "dedicated pair channels that emphasize the interface geometry."
)

found_bipartite = False
for i, p in enumerate(body_paras):
    if 'excluded from the topological analysis' in p.text:
        old_text = p.text
        # Find the sentence to replace
        old_sentence = None
        for candidate in [
            'This construction naturally encodes the intermolecular character of RNA–ligand interactions: intramolecular RNA–RNA and ligand–ligand contacts are excluded from the topological analysis, ensuring that the persistent Betti curves reflect exclusively the interface geometry between the two molecular partners.',
            'intramolecular RNA–RNA and ligand–ligand contacts are excluded from the topological analysis, ensuring that the persistent Betti curves reflect exclusively the interface geometry between the two molecular partners.',
        ]:
            if candidate in old_text:
                old_sentence = candidate
                break

        if old_sentence:
            new_text = old_text.replace(old_sentence, NEW_BIPARTITE_SUFFIX)
        else:
            # Fallback: replace everything after "respectively."
            idx = old_text.find('respectively.')
            if idx >= 0:
                new_text = old_text[:idx + len('respectively.')] + ' ' + NEW_BIPARTITE_SUFFIX
            else:
                new_text = old_text + ' [BIPARTITE_CLAIM_FIX_NEEDED]'

        set_para_text(p, new_text)
        print(f"  Para [{i}]: bipartite claim updated")
        print(f"    Old: {old_text[:100]}")
        print(f"    New: {p.text[:100]}")
        found_bipartite = True
        break

if not found_bipartite:
    print("  WARNING: bipartite claim paragraph not found")

# ═══════════════════════════════════════════════════════════════════════════════
# FIX 5 — Para [155]: Soften β₁ = degree-3 generator claim (Reviewer 2, Cmt 4)
# ═══════════════════════════════════════════════════════════════════════════════
print("\nFix 5: Soften β₁ degree-3 generator claim (Para ~[155])...")

OLD_BETA1_CLAIM = "whose algebraic counterpart is a degree-3 generator of I(Δʳ) that persists over a range [rᵇᵊʳʰᵗ, rᵈᵉᵃᵗʰ] encoding the spatial extent of this enclosure."
NEW_BETA1_CLAIM = ("which, via Hochster's formula, is consistent with the presence of degree-3 "
                   "generators in I(Δʳ): the minimal vertex support of a 1-cycle in the "
                   "bipartite complex contains three vertices, so the corresponding non-face "
                   "contributes a degree-3 monomial to I(Δʳ). This is not a strict "
                   "one-to-one bijection between individual β₁ classes and cubic "
                   "generators, but provides a physically grounded algebraic interpretation of "
                   "the persistent enclosure geometry.")

found_beta1 = False
for i, p in enumerate(body_paras):
    if 'degree-3 generator' in p.text and 'algebraic counterpart' in p.text:
        old_text = p.text
        # Replace the offending sentence
        changed = replace_in_para(p,
            'whose algebraic counterpart is a degree-3 generator',
            'which, via Hochster’s formula, is consistent with the presence of degree-3 generators in I(Δʳ): the minimal vertex support of a 1-cycle in the bipartite complex contains three vertices, so the corresponding non-face contributes a degree-3 monomial'
        )
        if not changed:
            # Broader replacement
            new_text = re.sub(
                r'whose algebraic counterpart is a degree-3 generator[^.]+\.',
                ("which, via Hochster’s formula, is consistent with the presence of degree-3 "
                 "generators in I(Δʳ): the minimal vertex support of a 1-cycle in the "
                 "bipartite complex contains three vertices, so the corresponding non-face "
                 "contributes a degree-3 monomial to I(Δʳ). This is not a strict "
                 "one-to-one bijection between individual β₁ classes and cubic "
                 "generators, but provides a physically grounded algebraic interpretation of "
                 "the persistent enclosure geometry."),
                old_text
            )
            if new_text != old_text:
                set_para_text(p, new_text)
                changed = True
        if changed:
            print(f"  Para [{i}]: β₁ claim softened")
            found_beta1 = True
            break

if not found_beta1:
    # Scan more broadly
    for i, p in enumerate(body_paras):
        if 'degree-3' in p.text and ('beta1' in p.text.lower() or 'β1' in p.text or 'β' in p.text):
            print(f"  WARNING: β₁ para found at [{i}] but replacement failed: {p.text[:100]}")
            break
    if not found_beta1:
        print("  WARNING: β₁ degree-3 claim paragraph not found")

# ═══════════════════════════════════════════════════════════════════════════════
# FIX 6 — Section 3.6: Add step numbering explanation (Reviewer 1, Comment 4)
# ═══════════════════════════════════════════════════════════════════════════════
print("\nFix 6: Add step numbering explanation in ML section...")

STEP_EXPLANATION = (
    "The CAML-RNA development pipeline comprises 33 sequential computational stages "
    "(Steps 1–33). Steps 1–4 construct the base feature modalities (ES Betti "
    "curves, CS Betti curves, PSRT f/h-vectors, and the global SVR baseline). "
    "Steps 5–15 explore additional feature modalities (CPF, RNA-FM, Morgan "
    "fingerprints, physicochemical descriptors) individually and in combination. "
    "Steps 16–33 progressively apply, evaluate, and refine the subtype-aware "
    "feature override strategy, including riboswitch subclass-level refinements and "
    "sign correction for anti-correlated subgroups. The final model (Step 33) "
    "integrates all feature modalities, subtype-specific overrides, and subclass "
    "refinements. All 33 steps are available in the public code repository."
)

# Find the end of the SVR description section (around Section 3.6 ML modeling)
# Look for the 'Sign Correction' or 'Subtype-Aware' paragraph and insert after the
# paragraph that ends Section 3.6's main body (before 3.7 Evaluation Metrics).
# We'll insert after the 'Sign Correction' paragraph (Para ~[144-145])
# or just after the paragraph describing the override criterion.

insert_after_idx = None
for i, p in enumerate(body_paras):
    if 'represents the empirical ceiling' in p.text or 'SAM/SAH' in p.text and 'empirical ceiling' in p.text:
        insert_after_idx = i
        break

if insert_after_idx is None:
    # Fallback: find evaluation metrics section header and insert before it
    for i, p in enumerate(body_paras):
        if '3.7' in p.text and 'Evaluation' in p.text:
            insert_after_idx = i - 1
            break

if insert_after_idx is not None:
    # Check if step explanation already added
    already = any('33 sequential computational stages' in p.text for p in body_paras)
    if not already:
        from docx.oxml.ns import qn
        from copy import deepcopy
        import lxml.etree as etree

        # Insert a new paragraph after insert_after_idx
        ref_para = paras[insert_after_idx]
        new_para = deepcopy(ref_para._element)
        # Clear all runs text and set new content
        ref_para._element.addnext(new_para)
        # Now update the new paragraph (it was inserted after ref_para)
        # Re-read paras to find it
        doc2_paras = list(doc.paragraphs)
        for j, p2 in enumerate(doc2_paras):
            if p2._element is new_para:
                if p2.runs:
                    p2.runs[0].text = STEP_EXPLANATION
                    for run in p2.runs[1:]:
                        run.text = ''
                else:
                    # Add run
                    from docx.oxml import OxmlElement
                    r_elem = OxmlElement('w:r')
                    t_elem = OxmlElement('w:t')
                    t_elem.text = STEP_EXPLANATION
                    r_elem.append(t_elem)
                    new_para.append(r_elem)
                break
        print(f"  Step explanation inserted after Para [{insert_after_idx}]")
    else:
        print("  Step explanation already present — skipping")
else:
    print("  WARNING: insert position for step explanation not found")

# ═══════════════════════════════════════════════════════════════════════════════
# Also fix Para [027] text re: CAML-RNA 10-fold comparison
# Add sentence about 10-fold result in the body comparison paragraph
# ═══════════════════════════════════════════════════════════════════════════════
print("\nFix 7: Add 10-fold CV comparison sentence in results text...")
for i, p in enumerate(body_paras):
    if ('RLASIF' in p.text and '0.666' in p.text and 'LOO-CV' in p.text
            and 'Among the methods' in p.text):
        old_text = p.text
        # Add a sentence about 10-fold result after the RMSE CI sentence
        INSERT_AFTER = ('The 95% bootstrap confidence interval [0.634, 0.802] '
                        '(10,000 resamples, seed 42) confirms that this advantage is statistically robust.')
        NEW_SENTENCE = (' Under nested 10-fold cross-validation, CAML-RNA achieves R₁₀ = 0.6255, '
                        'which is directly comparable with DeepMIF (R = 0.773 on NA-L, 10-fold CV); '
                        'the lower 10-fold result relative to LOO-CV reflects that the subtype-specific '
                        'pipeline was optimised under LOO-CV to maximise performance on the small dataset '
                        '(n = 143), as LOO maximises available training data per prediction.')
        if INSERT_AFTER in old_text:
            changed = replace_in_para(p, INSERT_AFTER, INSERT_AFTER + NEW_SENTENCE)
            if changed:
                print(f"  Para [{i}]: 10-fold comparison sentence added")
        elif 'bootstrap confidence interval' in old_text:
            # Find and append after the CI sentence
            new_text = re.sub(
                r'(The 95\s*%?\s*bootstrap confidence interval[^.]+\.)',
                r'\1' + NEW_SENTENCE,
                old_text
            )
            if new_text != old_text:
                set_para_text(p, new_text)
                print(f"  Para [{i}]: 10-fold comparison sentence added (regex)")
        break

# ═══════════════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════════════
doc.save(DOCX_OUT)
print(f"\nSaved → {DOCX_OUT}")

# ── Verification ────────────────────────────────────────────────────────────────
print("\n=== Verification ===")
doc_v = Document(DOCX_OUT)
paras_v = list(doc_v.paragraphs)
ref_v = next(i for i, p in enumerate(paras_v) if re.match(r'^[Rr]eferences?\s*$', p.text.strip()))

# Check Table 1 CAML-RNA row
tbl_v = doc_v.tables[0]
print(f"Table 1 Row 4 Protocol: {tbl_v.rows[4].cells[3].text!r}")
print(f"Table 1 Row 4 R:        {tbl_v.rows[4].cells[4].text!r}")

# Check for 0.784 remaining
found784 = [(i, p.text[:80]) for i, p in enumerate(paras_v[:ref_v]) if '0.784' in p.text]
print(f"Remaining 0.784 occurrences: {len(found784)}")
for x in found784:
    print(f"  Para [{x[0]}]: {x[1]}")

# Check bipartite claim
found_bipartite_v = [(i, p.text[:100]) for i, p in enumerate(paras_v[:ref_v])
                     if 'excluded from the topological analysis' in p.text]
print(f"'excluded from topological analysis' still present: {len(found_bipartite_v)}")

# Check degree-3 claim
found_deg3 = [(i, p.text[:100]) for i, p in enumerate(paras_v[:ref_v])
              if 'algebraic counterpart is a degree-3 generator' in p.text]
print(f"'algebraic counterpart is a degree-3 generator' still present: {len(found_deg3)}")

# Check step explanation
found_steps = [(i, p.text[:80]) for i, p in enumerate(paras_v[:ref_v])
               if '33 sequential computational stages' in p.text]
print(f"Step explanation present: {len(found_steps) > 0}")

# Check 10-fold result in text
found_6255 = [(i, p.text[:80]) for i, p in enumerate(paras_v[:ref_v]) if '0.6255' in p.text]
print(f"0.6255 (10-fold) mentions: {len(found_6255)}")
for x in found_6255:
    print(f"  Para [{x[0]}]: {x[1]}")

# Check footnote
for i, p in enumerate(paras_v[:ref_v]):
    if 'DeepRSMA' in p.text and 'R-SIM database' in p.text:
        print(f"Footnote Para [{i}]: {p.text[:120]}")
        break
