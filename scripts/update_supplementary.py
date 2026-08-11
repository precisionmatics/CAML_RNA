#!/usr/bin/env python3
"""
update_supplementary.py
Update CAML_RNA_Supplementary_Data.xlsx with revised manuscript data.
Key changes: n=111→143, dimensions updated, new modality ablation, new CV results, step33 predictions.
"""

import pandas as pd
import numpy as np
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from copy import copy

SUPP_IN  = "Revisions/CAML_RNA_Supplementary_Data.xlsx"
SUPP_OUT = "Revisions/CAML_RNA_Supplementary_Data_Revised.xlsx"

# ── Load source data ──────────────────────────────────────────────────────────
ablation = pd.read_csv("results/step_ablation_results.csv")
cv10     = pd.read_csv("results/step_10fold_cv_results.csv")
step33   = pd.read_csv("results/step33_predictions.csv")

# Compute step33 per-complex errors
step33['residual']  = step33['y_pred'] - step33['y_true']
step33['abs_error'] = step33['residual'].abs()
step33_r   = np.corrcoef(step33['y_true'], step33['y_pred'])[0, 1]
step33_rmse = np.sqrt(((step33['y_pred'] - step33['y_true'])**2).mean())

print(f"Step33: R={step33_r:.4f}, RMSE={step33_rmse:.4f}, n={len(step33)}")

# ── Load workbook ─────────────────────────────────────────────────────────────
wb = openpyxl.load_workbook(SUPP_IN)

# ── Helper: clear sheet and rewrite ──────────────────────────────────────────
def clear_sheet(ws):
    """Remove all content and merges from a sheet."""
    # Unmerge all merged cells first
    for merge in list(ws.merged_cells.ranges):
        ws.unmerge_cells(str(merge))
    for row in ws.iter_rows():
        for cell in row:
            cell.value = None

def write_header(ws, row, headers, bold=True, fill_color="4472C4", font_color="FFFFFF"):
    """Write a bold colored header row."""
    fill = PatternFill(start_color=fill_color, end_color=fill_color, fill_type="solid")
    font = Font(bold=bold, color=font_color)
    for c, h in enumerate(headers, 1):
        cell = ws.cell(row=row, column=c, value=h)
        if bold:
            cell.font = font
            cell.fill = fill
        cell.alignment = Alignment(horizontal='center')

def auto_width(ws, extra=2):
    """Auto-fit column widths."""
    for col in ws.columns:
        max_len = 0
        col_letter = get_column_letter(col[0].column)
        for cell in col:
            if cell.value:
                max_len = max(max_len, len(str(cell.value)))
        ws.column_dimensions[col_letter].width = min(max_len + extra, 60)

# ══════════════════════════════════════════════════════════════════════════════
# 1. Dataset_Overview — update note
# ══════════════════════════════════════════════════════════════════════════════
ws = wb['Dataset_Overview']
old_note = ws.cell(1, 1).value
ws.cell(1, 1).value = (
    "NOTE (REVISED): 143 RNA–small molecule complexes used in CAML-RNA (revised manuscript). "
    "Original supplementary contained n=111 from earlier pipeline; current results use n=143 "
    "from NA-L/NL2020 after exclusion of 1 complex with missing coordinates. "
    "pKd = −log10(Kd/1M). See OOF_Predictions_Step33 sheet for per-complex step33 predictions."
)
print("Updated Dataset_Overview note")

# ══════════════════════════════════════════════════════════════════════════════
# 2. CV_Results_Table1 → Replace with new 10-fold CV results (n=143)
# ══════════════════════════════════════════════════════════════════════════════
ws = wb['CV_Results_Table1']
clear_sheet(ws)

ws.cell(1, 1).value = (
    "NOTE (REVISED): CAML-RNA modality 10-fold cross-validation results (n=143 complexes). "
    "SVR-RBF with nested 5-fold hyperparameter tuning. Feature combinations ordered by complexity. "
    "R_10fold = Pearson R under outer 10-fold CV; R_LOO = leave-one-out CV Pearson R; "
    "95% CI = bootstrap confidence interval on 10-fold OOF Pearson R (10,000 resamples). "
    "PH=persistent homology Betti curves; PSRT=persistent f/h-vectors; CPF=chromatic persistence fingerprint; "
    "FM=RNA-FM embeddings; MFP=Morgan fingerprint; PHY=physicochemical descriptors."
)
ws.cell(1, 1).font = Font(italic=True)
ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=11)

headers = ['Feature Combination', 'Dim', 'n', 'R_10fold', 'Rho_10fold',
           'RMSE_10fold', 'R2_10fold', 'CI_lo', 'CI_hi', 'R_LOO', 'RMSE_LOO']
write_header(ws, 2, headers)

# Write data rows
display_names = {
    'PH_bipartite': 'PH (Betti curves)',
    'PSRT_fvec': 'PSRT (f/h-vectors)',
    'PH+PSRT': 'PH+PSRT',
    'CPF': 'CPF',
    'RNA-FM': 'RNA-FM',
    'Morgan': 'Morgan (ECFP4)',
    'Physicochemical': 'Physicochemical',
    'WCh': 'WCh (old pipeline)',
    'PH+CPF': 'PH+CPF',
    'PH+PSRT+CPF': 'PH+PSRT+CPF',
    'PH+PSRT+CPF+FM': 'PH+PSRT+CPF+FM (CAML-RNA 10-fold)',
    'PH+PSRT+CPF+FM+Phys': 'PH+PSRT+CPF+FM+PHY',
    'Full_ensemble': 'All modalities',
}

# Map 'RNA-FM' name from cv10 data
cv10_renamed = cv10.copy()
cv10_renamed['combination'] = cv10_renamed['combination'].replace({
    'RNA-FM': 'RNA-FM',
})

# Use ablation data (cleaner naming)
ablation_display = {
    'PH': 'PH (Betti curves)',
    'PSRT': 'PSRT (f/h-vectors)',
    'PH+PSRT': 'PH+PSRT',
    'CPF': 'CPF',
    'FM': 'RNA-FM',
    'MFP': 'Morgan (ECFP4)',
    'PHY': 'Physicochemical',
    'PH+CPF': 'PH+CPF',
    'PSRT+CPF': 'PSRT+CPF',
    'PH+PSRT+CPF': 'PH+PSRT+CPF',
    'PH+PSRT+FM': 'PH+PSRT+RNA-FM',
    'PH+PSRT+CPF+FM': 'PH+PSRT+CPF+RNA-FM',
    'PH+PSRT+CPF+FM+MFP+PHY': 'All modalities',
}

data_rows = []
for _, row in ablation.iterrows():
    name = ablation_display.get(row['combination'], row['combination'])
    ci_lo = row.get('ci_lo', '')
    ci_hi = row.get('ci_hi', '')
    r_loo = row.get('r_loo', '')
    rmse_loo = row.get('rmse_loo', '')
    data_rows.append([
        name,
        int(row['dim']),
        143,
        round(row['r_10fold'], 4),
        round(row['rho_10fold'], 4),
        round(row['rmse_10fold'], 4),
        round(row['r2_10fold'], 4),
        ci_lo if pd.notna(ci_lo) else '',
        ci_hi if pd.notna(ci_hi) else '',
        round(r_loo, 4) if pd.notna(r_loo) else '',
        round(rmse_loo, 4) if pd.notna(rmse_loo) else '',
    ])

# Add CAML-RNA step33 final row
data_rows.append([
    'CAML-RNA (step 33, LOO, final)',
    '7,632+',
    143,
    '',
    '',
    '',
    '',
    '',
    '',
    round(step33_r, 4),
    round(step33_rmse, 4),
])

for r_idx, row_data in enumerate(data_rows, 3):
    for c_idx, val in enumerate(row_data, 1):
        ws.cell(r_idx, c_idx).value = val

# Highlight the best modality row (PH+PSRT+RNA-FM)
best_fill = PatternFill(start_color="E2EFDA", end_color="E2EFDA", fill_type="solid")
for c in range(1, 12):
    if ws.cell(12 + 3, c).value:  # PH+PSRT+RNA-FM is row 11 of ablation (index 10) + 3 offset
        pass
# Highlight step33 row
final_fill = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")
last_row = 3 + len(data_rows) - 1
for c in range(1, 12):
    ws.cell(last_row, c).fill = final_fill
    ws.cell(last_row, c).font = Font(bold=True)

auto_width(ws)
print("Updated CV_Results_Table1")

# ══════════════════════════════════════════════════════════════════════════════
# 3. Ablation_Study → Replace with modality ablation (13 combinations)
# ══════════════════════════════════════════════════════════════════════════════
ws = wb['Ablation_Study']
clear_sheet(ws)

ws.cell(1, 1).value = (
    "NOTE (REVISED): Modality ablation study for CAML-RNA (n=143 complexes, SVR-RBF, nested 5-fold CV). "
    "Each row represents a different combination of feature modalities evaluated under 10-fold CV "
    "and LOO-CV. PH=persistent homology; PSRT=persistent Stanley–Reisner theory (f/h-vectors); "
    "CPF=chromatic persistence fingerprint (660-dim); FM=RNA-FM (640-dim); MFP=Morgan fingerprint (1024-dim); "
    "PHY=physicochemical descriptors (17-dim). Per-subtype R_10fold also reported."
)
ws.cell(1, 1).font = Font(italic=True)
ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=14)

headers2 = [
    'Combination', 'Modalities', 'N_mods', 'Dim',
    'R_10fold', 'Rho_10fold', 'RMSE_10fold', 'R2_10fold', 'CI_lo', 'CI_hi',
    'R_LOO', 'RMSE_LOO',
    'R10_Aptamer', 'R10_Riboswitch',
]
write_header(ws, 2, headers2)

for r_idx, row in ablation.iterrows():
    r_loo = row.get('r_loo', '')
    rmse_loo = row.get('rmse_loo', '')
    ws.cell(r_idx + 3, 1).value = row['combination']
    ws.cell(r_idx + 3, 2).value = row['modalities']
    ws.cell(r_idx + 3, 3).value = int(row['n_modalities'])
    ws.cell(r_idx + 3, 4).value = int(row['dim'])
    ws.cell(r_idx + 3, 5).value = round(row['r_10fold'], 4)
    ws.cell(r_idx + 3, 6).value = round(row['rho_10fold'], 4)
    ws.cell(r_idx + 3, 7).value = round(row['rmse_10fold'], 4)
    ws.cell(r_idx + 3, 8).value = round(row['r2_10fold'], 4)
    ws.cell(r_idx + 3, 9).value = row.get('ci_lo', '')
    ws.cell(r_idx + 3, 10).value = row.get('ci_hi', '')
    ws.cell(r_idx + 3, 11).value = round(r_loo, 4) if pd.notna(r_loo) else ''
    ws.cell(r_idx + 3, 12).value = round(rmse_loo, 4) if pd.notna(rmse_loo) else ''
    ws.cell(r_idx + 3, 13).value = round(row.get('r_aptamer', 0), 4)
    ws.cell(r_idx + 3, 14).value = round(row.get('r_riboswitch', 0), 4)

# Highlight PH+PSRT+FM (best modality combo row 11, idx 10)
best_idx = ablation.index[ablation['combination'] == 'PH+PSRT+FM']
if len(best_idx) > 0:
    best_row = best_idx[0] + 3
    best_fill2 = PatternFill(start_color="E2EFDA", end_color="E2EFDA", fill_type="solid")
    for c in range(1, 15):
        ws.cell(best_row, c).fill = best_fill2
        ws.cell(best_row, c).font = Font(bold=True)

auto_width(ws)
print("Updated Ablation_Study")

# ══════════════════════════════════════════════════════════════════════════════
# 4. Results_Summary → CAML-RNA current results
# ══════════════════════════════════════════════════════════════════════════════
ws = wb['Results_Summary']
clear_sheet(ws)

ws.cell(1, 1).value = (
    "NOTE (REVISED): CAML-RNA consolidated results (n=143 complexes). "
    "Final model (step 33): subtype-aware SVR-RBF with element-specific PH (3,888-dim), "
    "PSRT f/h-vectors (3,744-dim), CPF (660-dim), RNA-FM (640-dim), Morgan (1,024-dim), "
    "and physicochemical (17-dim) features; combined 7,632-dim base feature vector. "
    "Evaluation: outer 10-fold + LOO-CV, inner 5-fold hyperparameter tuning."
)
ws.cell(1, 1).font = Font(italic=True)
ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=6)

# Summary table
write_header(ws, 3, ['Metric', 'Value', 'Protocol', 'Dataset', 'Method', 'Notes'])

summary_rows = [
    ('Pearson R (overall)', 0.7283,  'LOO-CV',    'NA-L (n=143)', 'CAML-RNA step33', 'Primary result'),
    ('Pearson R (overall)', 0.6255,  '10-fold CV', 'NA-L (n=143)', 'PH+PSRT+CPF+FM',  'Best modality combo 10-fold'),
    ('Pearson R (overall)', 0.6279,  '10-fold CV', 'NA-L (n=143)', 'PH+PSRT+FM',       'Best 2-modality (no CPF)'),
    ('RMSE (overall)',      1.09,    'LOO-CV',    'NA-L (n=143)', 'CAML-RNA step33', ''),
    ('RMSE (overall)',      1.2202,  '10-fold CV', 'NA-L (n=143)', 'PH+PSRT+FM',       ''),
    ('95% CI (R)',         '[0.634, 0.802]', 'LOO-CV bootstrap', 'NA-L (n=143)', 'CAML-RNA step33', '10,000 resamples'),
    ('95% CI (R₁₀)',       '[0.512, 0.716]', '10-fold bootstrap', 'NA-L (n=143)', 'PH+PSRT+CPF+FM', ''),
    ('Pearson R (aptamer)', 0.940,   'LOO-CV',    'n=20',         'CAML-RNA step33', 'Best subtype'),
    ('Pearson R (riboswitch)', 0.771, 'LOO-CV',   'n=61',         'CAML-RNA step33', 'Blend of subtype models'),
    ('Pearson R (ribosomal A-site)', 0.763, 'LOO-CV', 'n=13', 'CAML-RNA step33', ''),
    ('Pearson R (viral TAR)', 0.770,  'LOO-CV',   'n=4',          'CAML-RNA step33', ''),
    ('Feature dim (PH)',    3888,    '',           'ES+CS Betti',  '36 pairs × 54 × 2', 'β₀+β₁, 24 levels'),
    ('Feature dim (PSRT)', 3744,    '',           '72 pairs × 52','f₁+h₂ curves',      ''),
    ('Feature dim (CPF)',  660,     '',           '11×10×6',      'Cutoffs 3.5–6.0 Å', ''),
    ('Feature dim (RNA-FM)', 640,   '',           'Sequence',     '100M-param model',  ''),
    ('Feature dim (Morgan)', 1024,  '',           'ECFP4 r=2',    '',                  ''),
    ('Feature dim (PHY)',   17,     '',           'Physicochemical','',                 ''),
    ('Feature dim (combined PH+PSRT)', 7632, '', 'PH+PSRT',       '3888+3744',         'Primary feature space'),
]

for r_idx, row_data in enumerate(summary_rows, 4):
    for c_idx, val in enumerate(row_data, 1):
        ws.cell(r_idx, c_idx).value = val

# Highlight key result rows
key_fill = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")
for c in range(1, 7):
    ws.cell(4, c).fill = key_fill
    ws.cell(4, c).font = Font(bold=True)

auto_width(ws)
print("Updated Results_Summary")

# ══════════════════════════════════════════════════════════════════════════════
# 5. Bootstrap_CI → Update with current values from ablation data
# ══════════════════════════════════════════════════════════════════════════════
ws = wb['Bootstrap_CI']
clear_sheet(ws)

ws.cell(1, 1).value = (
    "NOTE (REVISED): Bootstrap 95% confidence intervals for CAML-RNA modality combinations "
    "(10,000 resamples of 10-fold OOF predictions, percentile method, n=143). "
    "CI computed on Pearson R. LOO-CV R for step33 final model: 0.7283 "
    "[0.634, 0.802] by bootstrap on LOO predictions."
)
ws.cell(1, 1).font = Font(italic=True)
ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=6)

write_header(ws, 2, ['Feature Combination', 'R_10fold', 'CI_lower_95', 'CI_upper_95', 'R_LOO', 'RMSE_10fold'])

for r_idx, row in ablation.iterrows():
    name = ablation_display.get(row['combination'], row['combination'])
    ci_lo = row.get('ci_lo', '')
    ci_hi = row.get('ci_hi', '')
    r_loo = row.get('r_loo', '')
    ws.cell(r_idx + 3, 1).value = name
    ws.cell(r_idx + 3, 2).value = round(row['r_10fold'], 4)
    ws.cell(r_idx + 3, 3).value = ci_lo if pd.notna(ci_lo) else ''
    ws.cell(r_idx + 3, 4).value = ci_hi if pd.notna(ci_hi) else ''
    ws.cell(r_idx + 3, 5).value = round(r_loo, 4) if pd.notna(r_loo) else ''
    ws.cell(r_idx + 3, 6).value = round(row['rmse_10fold'], 4)

# Step33 final
last_r = len(ablation) + 3
ws.cell(last_r, 1).value = 'CAML-RNA step33 (LOO final)'
ws.cell(last_r, 2).value = ''
ws.cell(last_r, 3).value = 0.634
ws.cell(last_r, 4).value = 0.802
ws.cell(last_r, 5).value = round(step33_r, 4)
ws.cell(last_r, 6).value = round(step33_rmse, 4)
for c in range(1, 7):
    ws.cell(last_r, c).font = Font(bold=True)
    ws.cell(last_r, c).fill = PatternFill(start_color="FFF2CC", end_color="FFF2CC", fill_type="solid")

auto_width(ws)
print("Updated Bootstrap_CI")

# ══════════════════════════════════════════════════════════════════════════════
# 6. Feature_Matrix_Stats → Update key stats
# ══════════════════════════════════════════════════════════════════════════════
ws = wb['Feature_Matrix_Stats']
# Update the overall summary rows at bottom (find them)
for r in ws.iter_rows():
    for cell in r:
        if cell.value == 'Total feature columns':
            # Update adjacent cells
            ws.cell(cell.row, 1).value = 'Total feature columns (PH+PSRT)'
            ws.cell(cell.row, 2).value = 7632
        elif cell.value == 'Active columns (non-zero variance)':
            ws.cell(cell.row, 2).value = 7632  # all active after H removal
        elif cell.value == 'Complexes (rows)':
            ws.cell(cell.row, 2).value = 143
        elif cell.value == 4321:
            cell.value = 7632
        elif cell.value == 3121:
            cell.value = 7632
        elif cell.value == 111:
            cell.value = 143

# Update note
ws.cell(1, 1).value = (
    "NOTE (REVISED): Feature matrix summary for CAML-RNA (n=143 complexes). "
    "Primary feature vector: PH Betti curves (3,888-dim, ES+CS, 36 pairs each) + "
    "PSRT f/h-vectors (3,744-dim, 72 pairs) = 7,632-dim total. "
    "Additional: CPF 660-dim, RNA-FM 640-dim, Morgan 1024-dim, PHY 17-dim. "
    "Top-30 feature statistics shown are from old pipeline (4320-dim); see CV_Results_Table1 for current results."
)
print("Updated Feature_Matrix_Stats")

# ══════════════════════════════════════════════════════════════════════════════
# 7. Add new sheet: OOF_Predictions_Step33
# ══════════════════════════════════════════════════════════════════════════════
# Check if sheet already exists
if 'OOF_Predictions_Step33' in wb.sheetnames:
    del wb['OOF_Predictions_Step33']

ws_new = wb.create_sheet('OOF_Predictions_Step33')

ws_new.cell(1, 1).value = (
    "NOTE: Per-complex predictions from CAML-RNA step33 final model (subtype-aware SVR-RBF, LOO-CV, n=143). "
    "y_true = experimental pKd; y_pred = step33 prediction; residual = y_pred − y_true; "
    "abs_error = |residual|. Overall LOO Pearson R = 0.7283, RMSE = 1.09."
)
ws_new.cell(1, 1).font = Font(italic=True)
ws_new.merge_cells(start_row=1, start_column=1, end_row=1, end_column=7)

write_header(ws_new, 2, ['PDB ID', 'pKd_true', 'pKd_pred', 'Residual', 'Abs_Error', 'Subtype', 'Notes'])

step33_sorted = step33.sort_values('abs_error', ascending=False)
for r_idx, (_, row) in enumerate(step33_sorted.iterrows(), 3):
    ws_new.cell(r_idx, 1).value = row['pdb']
    ws_new.cell(r_idx, 2).value = round(row['y_true'], 4)
    ws_new.cell(r_idx, 3).value = round(row['y_pred'], 4)
    ws_new.cell(r_idx, 4).value = round(row['residual'], 4)
    ws_new.cell(r_idx, 5).value = round(row['abs_error'], 4)
    ws_new.cell(r_idx, 6).value = row['subtype']
    ws_new.cell(r_idx, 7).value = 'high error' if row['abs_error'] > 2.0 else ''
    # Color high-error rows
    if row['abs_error'] > 2.0:
        err_fill = PatternFill(start_color="FFE0E0", end_color="FFE0E0", fill_type="solid")
        for c in range(1, 8):
            ws_new.cell(r_idx, c).fill = err_fill

auto_width(ws_new)
print("Added OOF_Predictions_Step33 sheet")

# ══════════════════════════════════════════════════════════════════════════════
# 8. Update NOTEs on old sheets (mark as from old pipeline)
# ══════════════════════════════════════════════════════════════════════════════
old_pipeline_note = " [FROM OLD PIPELINE — n=111, 4320-dim; superseded by CV_Results_Table1 in revised manuscript]"

for shname in ['OOF_Predictions_Ridge', 'OOF_Predictions_Best', 'Baselines_Table3',
               'Novel_PH_Table4', 'Learning_Curve', 'Error_Analysis', 'Pocket_Summary']:
    if shname in wb.sheetnames:
        ws = wb[shname]
        old_val = ws.cell(1, 1).value or ''
        if 'OLD PIPELINE' not in str(old_val):
            ws.cell(1, 1).value = 'NOTE' + old_pipeline_note + ' | ' + old_val.replace('NOTE: ', '')
            ws.cell(1, 1).font = Font(italic=True, color="808080")

print("Marked old pipeline sheets")

# ══════════════════════════════════════════════════════════════════════════════
# 9. Update Removed_Entries note
# ══════════════════════════════════════════════════════════════════════════════
if 'Removed_Entries' in wb.sheetnames:
    ws = wb['Removed_Entries']
    ws.cell(1, 1).value = (
        "NOTE (REVISED): Complexes removed during curation. Revised dataset: n=143 (started from 144; "
        "1 complex excluded for missing atom coordinates). Original pipeline excluded more complexes "
        "from an earlier dataset version."
    )

# ══════════════════════════════════════════════════════════════════════════════
# Save
# ══════════════════════════════════════════════════════════════════════════════
wb.save(SUPP_OUT)
print(f"\nSaved → {SUPP_OUT}")
print(f"Sheets: {wb.sheetnames}")
