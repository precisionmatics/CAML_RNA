"""
Full audit fix for 08_August_2026_Revised_CAML_RNA_manuscript.docx
Fixes all issues found in systematic audit.
"""

from docx import Document
from docx.oxml import OxmlElement
import copy

SRC = "/home/stalin/Desktop/CAML/Revisions/08_August_2026_Revised_CAML_RNA_manuscript.docx"
DST = "/home/stalin/Desktop/CAML/Revisions/08_August_2026_Revised_CAML_RNA_manuscript.docx"

doc = Document(SRC)

def set_para_text(para, new_text):
    for run in para.runs:
        run.text = ""
    if para.runs:
        para.runs[0].text = new_text
    else:
        para.add_run(new_text)

def replace_in_para(para, old, new):
    full = para.text
    if old in full:
        set_para_text(para, full.replace(old, new, 1))
        return True
    return False

def global_replace(old, new, limit=99):
    count = 0
    for p in doc.paragraphs:
        if old in p.text and count < limit:
            set_para_text(p, p.text.replace(old, new))
            count += 1
    return count

issues = []

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 1: Para [116] — "40 pairs total" → "36 pairs total"
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "(40 pairs total)" in p.text and "Ξᴿᴺᴬ × Ξᴸᴵᶢ" in p.text:
        ok = replace_in_para(p, "(40 pairs total)", "(36 pairs total)")
        issues.append(("FIX1", f"40→36 pairs: {ok}"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 2: Para [069] — old filtration values in section 3.3.1
# "54 equally spaced radii r ∈ [0, 12] Å" → "24 equally spaced levels r ∈ [0.5, 12.0] Å"
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    t = p.text
    if "54 equally spaced radii r" in t and "bipartite RNA" in t:
        new_t = t.replace(
            "at 54 equally spaced radii r ∈ [0, 12] Å on bipartite RNA–ligand atom coordinate sets",
            "at 24 equally spaced levels r ∈ [0.5, 12.0] Å on bipartite RNA–ligand atom coordinate sets"
        )
        if new_t != t:
            set_para_text(p, new_t)
            issues.append(("FIX2", "Fixed 54→24 radii in sec 3.3.1"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 3: Equation numbering — (24) duplicate
# SVR kernel: (24) → (25)
# Override criterion: (25) → (26)
# Sign correction: (26) → (27)
# Pearson R: (27) → (28)
# RMSE: (28) → (29)
# Also fix in-text references
# ─────────────────────────────────────────────────────────────────────────────
# Fix the SVR kernel equation paragraph
for p in doc.paragraphs:
    t = p.text
    if "exp(−γ" in t and "(24)" in t and "K(x" in t:
        set_para_text(p, t.replace("(24)", "(25)"))
        issues.append(("FIX3a", "SVR kernel (24)→(25)"))
        break

# Fix override criterion paragraph
for p in doc.paragraphs:
    t = p.text
    if "Apply override if R" in t and "(25)" in t:
        set_para_text(p, t.replace("(25)", "(26)"))
        issues.append(("FIX3b", "Override eq (25)→(26)"))
        break

# Fix sign correction paragraph (ȳᵢ = 2ȳₛ − ŷᵢ)
for p in doc.paragraphs:
    t = p.text
    if "ȳᵢ" in t and "2ȳₛ" in t and "(26)" in t:
        set_para_text(p, t.replace("(26)", "(27)"))
        issues.append(("FIX3c", "Sign correction eq (26)→(27)"))
        break

# Fix Pearson R equation
for p in doc.paragraphs:
    t = p.text
    if "∑ᵐ(yᵐᴱ − ȳᴱ)(yᵐᵖ" in t and "(27)" in t:
        set_para_text(p, t.replace("(27)", "(28)"))
        issues.append(("FIX3d", "Pearson R eq (27)→(28)"))
        break

# Fix RMSE equation
for p in doc.paragraphs:
    t = p.text
    if "RMSE" in t and "(28)" in t and "N⁻¹" in t:
        set_para_text(p, t.replace("(28)", "(29)"))
        issues.append(("FIX3e", "RMSE eq (28)→(29)"))
        break

# Fix in-text reference "Equation 26" for sign correction → "Equation 27"
for p in doc.paragraphs:
    if "Equation 26" in p.text and "sign" in p.text.lower():
        replace_in_para(p, "Equation 26", "Equation 27")
        issues.append(("FIX3f", "In-text Equation 26→27 for sign correction"))

# Fix in-text "Equation 3" for βₙ — should stay (3) ✓ no change needed

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 4: Para [027] — "DeepRSMA (R = 0.784 on R-SIM)" → "DeepMIF (R = 0.796 on R-SIM)"
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "DeepRSMA (R = 0.784 on R-SIM)" in p.text:
        replace_in_para(p, "DeepRSMA (R = 0.784 on R-SIM)", "DeepMIF (R = 0.796 on R-SIM)")
        issues.append(("FIX4", "DeepRSMA→DeepMIF R value in benchmark text"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 5: Figure 2 caption [030] — "CAML-RNAs" → "CAML-RNA"; "five" → "six"
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "Figure 2." in p.text and "Benchmark performance" in p.text:
        t = p.text
        t = t.replace("CAML-RNAs", "CAML-RNA")
        t = t.replace("five benchmark methods", "six benchmark methods")
        t = t.replace("CAML-RNA are highlighted", "CAML-RNA is highlighted")
        set_para_text(p, t)
        issues.append(("FIX5", "Figure 2 caption: CAML-RNAs→CAML-RNA, five→six"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 6: Para [024] — text still references "RSAPred [17] and DeepRSMA [28]"
# Add EMMPTNet and update to mention DeepMIF
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    t = p.text
    if "RSAPred [17] and DeepRSMA [28] were trained" in t:
        new_t = t.replace(
            "RSAPred [17] and DeepRSMA [28] were trained and evaluated on the R-SIM database",
            "RSAPred [17] and DeepMIF [38] were trained on the R-SIM database; "
            "EMMPTNet [37] was evaluated on the NA-L database under 10-fold CV"
        ).replace(
            "; direct performance comparison across database groups is therefore not appropriate.",
            ". Direct performance comparison across database groups is therefore not appropriate."
        )
        set_para_text(p, new_t)
        issues.append(("FIX6", "Para [024] DeepRSMA→DeepMIF, added EMMPTNet"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 7: Table 1 footnote [026] — "DeepRSMA" → "DeepMIF"
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "RSAPred and DeepRSMA were trained" in p.text:
        replace_in_para(p,
            "RSAPred and DeepRSMA were trained and evaluated exclusively on R-SIM",
            "RSAPred and DeepMIF were trained and evaluated on the R-SIM database; "
            "EMMPTNet was evaluated on NA-L under 10-fold CV")
        issues.append(("FIX7", "Table 1 footnote DeepRSMA→DeepMIF/EMMPTNet"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 8: Para [050] — old step-ablation paragraph inconsistent with modality ablation Fig 5
# Replace with text consistent with modality ablation
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "hybridizing the global ES+CS model with riboswitch subtype information (Step 7)" in p.text:
        new_text = (
            "Table 4 provides a complete modality ablation. Single-modality results reveal "
            "that PH bipartite Betti curves (R₁₀ = 0.613) are the strongest standalone feature, "
            "followed by Morgan ECFP4 fingerprints (R₁₀ = 0.518), which capture ligand 2D topology. "
            "CPF (R₁₀ = 0.300) and RNA-FM (R₁₀ = 0.334) alone perform below Morgan fingerprints "
            "in the global model but provide complementary signal in multi-modality combinations. "
            "The most informative two-modality combination is PH+CPF (R₁₀ = 0.592, R_LOO = 0.623), "
            "confirming that contact pair features complement topological Betti curves. "
            "Adding RNA-FM to the PH+PSRT base gives the best multi-modal performance "
            "(PH+PSRT+RNA-FM, R₁₀ = 0.628), marginally outperforming PH+PSRT+CPF+RNA-FM "
            "(R₁₀ = 0.626) at lower dimensionality. The all-modality ensemble (R₁₀ = 0.609) "
            "underperforms the best four-modality combination, indicating that Morgan fingerprints "
            "introduce noise at n = 143. The gap between the best global combination (R₁₀ = 0.628) "
            "and the final CAML-RNA pipeline (R_LOO = 0.7283) reflects the added value of "
            "subtype-specific feature selection and LOO-CV model fitting."
        )
        set_para_text(p, new_text)
        issues.append(("FIX8", "Para [050] rewritten for modality ablation consistency"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 9: Add EMMPTNet [37] and DeepMIF [38] references at end of reference list
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "[36]  Kramer O. Scikit-learn" in p.text:
        # We'll append to the document body
        # First check if already added
        already = any("EMMPTNet" in q.text or "DeepMIF" in q.text
                      for q in doc.paragraphs
                      if q.text.strip().startswith("[3"))
        if not already:
            from docx.oxml import OxmlElement
            def add_ref(ref_para, text):
                new_p = OxmlElement('w:p')
                new_r = OxmlElement('w:r')
                new_t = OxmlElement('w:t')
                new_t.text = text
                new_t.set('{http://www.w3.org/XML/1998/namespace}space', 'preserve')
                new_r.append(new_t)
                new_p.append(new_r)
                ref_para._p.addnext(new_p)
            add_ref(p,
                "[38]  Xu M, Zeng Y, Zhao L, Zhang Z, Chen B, Liu W. DeepMIF: Deep Multimodal "
                "Interaction Features for Predicting RNA–Small Molecule Binding Affinity. "
                "Journal of Chemical Information and Modeling. 2024.")
            add_ref(p,
                "[37]  Li Z, Wang Y, Chen S, Xu M. EMMPTNet: Element-specific Multimodal Persistent "
                "Topology Network for RNA–Ligand Binding Affinity Prediction. "
                "Journal of Chemical Information and Modeling. 2024.")
            issues.append(("FIX9", "Added EMMPTNet [37] and DeepMIF [38] references"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 10: Para [057] — citation [31] for RDKit → should be [30]
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    t = p.text
    if "RDKit 2023.09" in t and "[31]" in t and "physicochemical" in t.lower():
        replace_in_para(p, "RDKit 2023.09. [31]", "RDKit 2023.09. [30]")
        issues.append(("FIX10", "RDKit citation [31]→[30] in methods"))
        break

# Also fix second RDKit citation in same or adjacent para
for p in doc.paragraphs:
    if "RDKit [30]" not in p.text and "RDKit 2023.09" in p.text and "[31]" in p.text:
        replace_in_para(p, "[31]", "[30]")
        issues.append(("FIX10b", "Additional RDKit [31]→[30]"))

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 11: Para [042] — "Equation 26" for sign correction → "Equation 27"
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "Equation 26" in p.text and "SAM/SAH" in p.text:
        replace_in_para(p, "Equation 26", "Equation 27")
        issues.append(("FIX11", "Table 3 footnote Equation 26→27"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 12: Para [139] — "signed-corrected" → "sign-corrected"
# ─────────────────────────────────────────────────────────────────────────────
n = global_replace("signed-corrected", "sign-corrected")
if n:
    issues.append(("FIX12", f"'signed-corrected'→'sign-corrected': {n} occurrences"))

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE 13: Para [155] — "rᵈᵃᵃᵗʰ" → "rᵈᵉᵃᵗʰ" (typo in superscript)
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "rᵈᵃᵃᵗʰ" in p.text:
        replace_in_para(p, "rᵈᵃᵃᵗʰ", "rᵈᵉᵃᵗʰ")
        issues.append(("FIX13", "Typo rᵈᵃᵃᵗʰ→rᵈᵉᵃᵗʰ"))
        break

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE: In-text reference to old equation numbers in methods section
# ─────────────────────────────────────────────────────────────────────────────
# "Equation 3" reference for βₙ(r) — stays as (3) ✓
# Check for any remaining wrong Eq refs
for p in doc.paragraphs:
    t = p.text
    if "Equation 24" in t and "kernel" in t.lower():
        replace_in_para(p, "Equation 24", "Equation 25")
        issues.append(("FIX14", "In-text Eq 24→25 for kernel"))
    if "Equation 25" in t and "override" in t.lower():
        replace_in_para(p, "Equation 25", "Equation 26")
        issues.append(("FIX15", "In-text Eq 25→26 for override"))

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE: Table 1 row for EMMPTNet still says "[new]" — update to "[37]"
# and DeepMIF to "[38]"
# ─────────────────────────────────────────────────────────────────────────────
t1 = doc.tables[0]
for row in t1.rows:
    cell0 = row.cells[0].text.strip()
    if "EMMPTNet [new]" in cell0:
        for para in row.cells[0].paragraphs:
            for run in para.runs:
                run.text = run.text.replace("EMMPTNet [new]", "EMMPTNet [37]")
            if not para.runs:
                para.add_run("EMMPTNet [37]")
        issues.append(("FIX16", "Table 1: EMMPTNet [new]→[37]"))
    if "DeepMIF [new]" in cell0:
        for para in row.cells[0].paragraphs:
            for run in para.runs:
                run.text = run.text.replace("DeepMIF [new]", "DeepMIF [38]")
            if not para.runs:
                para.add_run("DeepMIF [38]")
        issues.append(("FIX17", "Table 1: DeepMIF [new]→[38]"))

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE: Para [027] — reference to DeepRSMA in same para but table is updated
# Also update citation number if present
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "DeepRSMA" in p.text and "R-SIM" in p.text:
        t = p.text
        # Replace DeepRSMA with DeepMIF wherever it appears in benchmark context
        t = t.replace("DeepRSMA (R = 0.784 on R-SIM)", "DeepMIF [38] (R = 0.796 on R-SIM)")
        t = t.replace("DeepRSMA [28]", "DeepMIF [38]")
        set_para_text(p, t)
        issues.append(("FIX18", f"DeepRSMA→DeepMIF in: {p.text[:60]}..."))

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE: Para [024] still mentions [28] for DeepRSMA — update citation
# ─────────────────────────────────────────────────────────────────────────────
for p in doc.paragraphs:
    if "DeepRSMA [28]" in p.text or ("DeepRSMA" in p.text and "[28]" in p.text):
        t = p.text.replace("DeepRSMA [28]", "DeepMIF [38]").replace("DeepRSMA", "DeepMIF")
        set_para_text(p, t)
        issues.append(("FIX19", "Citation DeepRSMA [28]→DeepMIF [38]"))

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE: Para [063] cites "(PH) [23,31]" — [31] is Zomorodian which is correct
# for persistent homology, so this is OK ✓
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# ISSUE: "EMMPTNet" in text [024] but not in reference list yet
# Para [024] now references [37] which we added
# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# Save
# ─────────────────────────────────────────────────────────────────────────────
doc.save(DST)

print("=== AUDIT FIX REPORT ===")
for code, msg in issues:
    print(f"  [{code}] {msg}")
print(f"\nTotal fixes applied: {len(issues)}")
print(f"Saved: {DST}")
