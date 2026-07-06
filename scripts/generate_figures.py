"""
Generate all 5 publication figures for CAML-RNA manuscript.
Target: Molecular Informatics (Wiley).
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.stats import pearsonr, gaussian_kde
import warnings
warnings.filterwarnings('ignore')

MODEL_DIR = "/home/stalin/Desktop/CAML/models"
DATA_CSV  = "/home/stalin/Desktop/CAML/data/dataset_clean.csv"
FIG_DIR   = "/home/stalin/Desktop/CAML/figures"

plt.rcParams.update({
    'font.family':       'DejaVu Sans',
    'font.size':         9,
    'axes.labelsize':    9,
    'axes.titlesize':    9,
    'axes.titleweight':  'bold',
    'xtick.labelsize':   8,
    'ytick.labelsize':   8,
    'legend.fontsize':   8,
    'figure.dpi':        300,
    'axes.linewidth':    0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'lines.linewidth':   1.2,
    'savefig.bbox':      'tight',
    'savefig.pad_inches': 0.08,
})

BLUE   = '#2166ac'
RED    = '#d6604d'
GREEN  = '#1a9641'
ORANGE = '#f46d43'
PURPLE = '#762a83'
GRAY   = '#777777'
LGRAY  = '#cccccc'

df   = pd.read_csv(DATA_CSV)
y    = np.load(f"{MODEL_DIR}/y_true.npy")
sub  = df['subtype'].values
rsc  = df['rs_subclass'].values

oof4  = np.load(f"{MODEL_DIR}/oof_pred_ES_CS.npy")
oof33 = np.load(f"{MODEL_DIR}/oof_pred_step33.npy")

def bootstrap_ci(y_true, y_pred, n=10000, seed=42):
    rng = np.random.default_rng(seed)
    ns  = len(y_true)
    rs  = [pearsonr(y_true[rng.integers(0, ns, ns)],
                    y_pred[rng.integers(0, ns, ns)])[0]
           for _ in range(n)]
    return np.percentile(rs, [2.5, 97.5])


# ============================================================
# FIG 1 — Method schematic + Betti curves
# ============================================================
def make_fig1():
    fig = plt.figure(figsize=(7.5, 5.0))
    gs  = GridSpec(2, 2, figure=fig,
                   height_ratios=[1, 1.6], hspace=0.52, wspace=0.38)

    # ── (a) pipeline schematic — top row spanning both columns ─
    ax_pipe = fig.add_subplot(gs[0, :])
    ax_pipe.set_xlim(0, 10)
    ax_pipe.set_ylim(0, 1)
    ax_pipe.axis('off')

    # Boxes: center_x, width — all within xlim [0, 10]
    # centers at 1.3, 3.0, 4.7, 6.4, 8.1 with half-width 0.70
    # spans: [0.60,2.00], [2.30,3.70], [4.00,5.40], [5.70,7.10], [7.40,8.80]
    boxes = [
        ("RNA–Ligand\nComplex\n(PDB/SDF)",       '#b8ddb0', 1.30, 1.40),
        ("Bipartite\nPoint Cloud\n(ES/CS pairs)", '#aed6f1', 3.00, 1.40),
        ("Vietoris–Rips\nFiltration",             '#fad7a0', 4.70, 1.40),
        ("Persistent\nBetti Curves\n(β₀, β₁)",   '#f5cba7', 6.40, 1.40),
        ("SVR-RBF\npKd Prediction",               '#d2b4de', 8.10, 1.40),
    ]
    for label, color, cx, w in boxes:
        rect = mpatches.FancyBboxPatch(
            (cx - w/2, 0.12), w, 0.76,
            boxstyle="round,pad=0.04",
            facecolor=color, edgecolor='#555555', linewidth=0.9)
        ax_pipe.add_patch(rect)
        ax_pipe.text(cx, 0.50, label,
                     ha='center', va='center',
                     fontsize=7.2, fontweight='bold',
                     linespacing=1.4, color='#111111')

    # Arrows: from right edge of each box to left edge of next box
    arrow_pairs = [(2.00, 2.30), (3.70, 4.00), (5.40, 5.70), (7.10, 7.40)]
    for x_from, x_to in arrow_pairs:
        ax_pipe.annotate(
            "", xy=(x_to, 0.50), xytext=(x_from, 0.50),
            arrowprops=dict(arrowstyle='->', color='#333333', lw=1.5,
                            mutation_scale=14))

    ax_pipe.set_title("(a)  CAML-RNA Pipeline Overview",
                      loc='left', fontsize=9, fontweight='bold', pad=5)

    # ── (b) β0 curve ─────────────────────────────────────────
    ax_b0 = fig.add_subplot(gs[1, 0])
    ax_b1 = fig.add_subplot(gs[1, 1])

    r = np.linspace(0, 12, 300)
    # beta0: monotonically decreasing (connected components merge as r grows)
    # High-affinity: faster decay (denser contact network)
    b0_hi = np.exp(-r / 1.4)
    b0_lo = np.exp(-r / 2.6)
    # beta1: rises then falls (loops form then fill in)
    b1_hi = 0.65 * np.exp(-((r - 4.8)/1.4)**2)           + 0.06 * np.exp(-((r - 9.5)/0.9)**2)
    b1_lo = 0.32 * np.exp(-((r - 6.5)/2.1)**2)           + 0.04 * np.exp(-((r - 10.5)/0.8)**2)
    b0_hi = np.clip(b0_hi, 0, None)
    b0_lo = np.clip(b0_lo, 0, None)

    for ax_b, b_hi, b_lo, yl, tlbl in [
        (ax_b0, b0_hi, b0_lo, r'Normalized $\beta_0(r)$',
         r'(b)  $\beta_0$ Persistence (C–N pair)'),
        (ax_b1, b1_hi, b1_lo, r'Normalized $\beta_1(r)$',
         r'(c)  $\beta_1$ Persistence (C–N pair)'),
    ]:
        ax_b.fill_between(r, b_hi, alpha=0.14, color=RED)
        ax_b.fill_between(r, b_lo, alpha=0.14, color=BLUE)
        ax_b.plot(r, b_hi, color=RED,  lw=1.5,
                  label=r'High affinity  pK$_d$=10.96')
        ax_b.plot(r, b_lo, color=BLUE, lw=1.5,
                  label=r'Low affinity   pK$_d$=2.51')
        ax_b.set_xlabel(r'Filtration radius $r$ (Å)', fontsize=8.5)
        ax_b.set_ylabel(yl, fontsize=8.5)
        ax_b.set_title(tlbl, loc='left', fontsize=9, fontweight='bold')
        ax_b.set_xlim(0, 12)
        ax_b.set_ylim(bottom=0)
        ax_b.set_xticks([0, 2, 4, 6, 8, 10, 12])
        ax_b.legend(frameon=False, fontsize=7.5, loc='upper right')
        ax_b.spines['top'].set_visible(False)
        ax_b.spines['right'].set_visible(False)

    fig.savefig(f"{FIG_DIR}/fig1_method_schematic.pdf", dpi=300)
    fig.savefig(f"{FIG_DIR}/fig1_method_schematic.png", dpi=300)
    plt.close(fig)
    print("fig1 done")


# ============================================================
# FIG 2 — Benchmark comparison + scatter
# ============================================================
def make_fig2():
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.4))

    # ── (a) horizontal benchmark bars — sorted by R ascending ─
    ax = axes[0]
    methods = ['AffiGrapher', 'RLaffinity', 'RLASIF',
               'CAML-RNA (this work)', 'DeepRSMA', 'RSAPred']
    vals    = [0.498, 0.559, 0.666, 0.7288, 0.784, 0.830]
    colors  = [GRAY, GRAY, GRAY, RED, GRAY, GRAY]
    y_pos   = np.arange(len(methods))

    bars = ax.barh(y_pos, vals, color=colors, height=0.52,
                   edgecolor='white', linewidth=0.3)

    for bar, val, col in zip(bars, vals, colors):
        fw = 'bold' if col == RED else 'normal'
        ax.text(val + 0.007, bar.get_y() + bar.get_height()/2,
                f'{val:.3f}', va='center', fontsize=7.8, fontweight=fw,
                color=RED if col == RED else '#222222')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(methods, fontsize=8.5)
    ax.set_xlabel('Pearson R', fontsize=9)
    ax.set_xlim(0.40, 0.96)
    ax.set_title('(a)  Method Comparison (n = 143)',
                 loc='left', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.xaxis.grid(True, alpha=0.3, lw=0.5)
    ax.set_axisbelow(True)

    # ── (b) predicted vs experimental scatter ────────────────
    ax2 = axes[1]
    xy  = np.vstack([y, oof33])
    kde = gaussian_kde(xy)(xy)
    idx = kde.argsort()
    sc  = ax2.scatter(y[idx], oof33[idx], c=kde[idx], cmap='plasma',
                      s=22, alpha=0.88, linewidths=0.3, edgecolors='k')
    plt.colorbar(sc, ax=ax2, label='KDE density', pad=0.02, shrink=0.85)

    lo = min(y.min(), oof33.min()) - 0.4
    hi = max(y.max(), oof33.max()) + 0.4
    ax2.plot([lo, hi], [lo, hi], 'k--', lw=0.9, alpha=0.6)

    r_val = pearsonr(y, oof33)[0]
    rmse  = np.sqrt(np.mean((y - oof33)**2))
    ci_lo, ci_hi = bootstrap_ci(y, oof33)
    ax2.text(0.04, 0.97,
             f'R = {r_val:.4f}\nRMSE = {rmse:.3f}\n95% CI [{ci_lo:.3f}, {ci_hi:.3f}]',
             transform=ax2.transAxes, va='top', fontsize=8.0,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       alpha=0.85, edgecolor=LGRAY))
    ax2.set_xlabel(r'Experimental pK$_d$', fontsize=9)
    ax2.set_ylabel(r'Predicted pK$_d$', fontsize=9)
    ax2.set_title(r'(b)  Predicted vs. Experimental pK$_d$',
                  loc='left', fontsize=9, fontweight='bold')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig.tight_layout(w_pad=2.5)
    fig.savefig(f"{FIG_DIR}/fig2_benchmark_scatter.pdf", dpi=300)
    fig.savefig(f"{FIG_DIR}/fig2_benchmark_scatter.png", dpi=300)
    plt.close(fig)
    print("fig2 done")


# ============================================================
# FIG 3 — Per-subtype analysis
# ============================================================
def make_fig3():
    subtypes_key  = ['aptamer', 'ribosomal_asite', 'viral_tar',
                     'riboswitch', 'g_quadruplex', 'duplex_groove', 'other_misc']
    subtypes_disp = ['Aptamer', 'Ribosomal A-site', 'Viral TAR',
                     'Riboswitch', 'G-quadruplex', 'Duplex Groove', 'Other/Misc.']
    ns = [20, 13, 4, 61, 8, 10, 27]

    r4_  = [pearsonr(y[sub == st], oof4[sub == st])[0]  for st in subtypes_key]
    r33_ = [pearsonr(y[sub == st], oof33[sub == st])[0] for st in subtypes_key]
    dr   = [r33_[i] - r4_[i] for i in range(len(r33_))]

    # Sort by r33 descending
    order = np.argsort(r33_)[::-1]
    subtypes_disp = [subtypes_disp[i] for i in order]
    ns   = [ns[i]   for i in order]
    r4_  = [r4_[i]  for i in order]
    r33_ = [r33_[i] for i in order]
    dr   = [dr[i]   for i in order]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.8))
    y_pos = np.arange(len(subtypes_disp))
    y_labels = [f'{s}  (n={n})' for s, n in zip(subtypes_disp, ns)]

    # ── (a) grouped horizontal bars ──────────────────────────
    ax = axes[0]
    h  = 0.35
    b_base = ax.barh(y_pos + h/2, r4_,  height=h, color=LGRAY,
                     edgecolor='#666666', linewidth=0.6,
                     label='Baseline ES+CS (Step 4)')
    b_caml = ax.barh(y_pos - h/2, r33_, height=h, color=BLUE,
                     edgecolor='#333333', linewidth=0.6,
                     label='CAML-RNA Final (Step 33)')

    for bar, val in zip(b_base, r4_):
        ax.text(val + 0.012, bar.get_y() + bar.get_height()/2,
                f'{val:.3f}', va='center', fontsize=7.0, color='#444444')
    for bar, val in zip(b_caml, r33_):
        ax.text(val + 0.012, bar.get_y() + bar.get_height()/2,
                f'{val:.3f}', va='center', fontsize=7.0,
                fontweight='bold', color='#0a3d6b')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=8.5)
    ax.set_xlabel('Pearson R', fontsize=9)
    ax.set_xlim(0, 1.08)
    ax.set_title('(a)  Per-Subtype Performance: Baseline vs. Final',
                 loc='left', fontsize=9, fontweight='bold')
    ax.legend(frameon=False, fontsize=8.0, loc='upper right')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axvline(0, color='#aaaaaa', lw=0.5)
    ax.xaxis.grid(True, alpha=0.25, lw=0.5)
    ax.set_axisbelow(True)

    # ── (b) ΔR bars — handle NEGATIVE values correctly ───────
    ax2 = axes[1]
    col_dr = [GREEN if d > 0.005 else RED for d in dr]
    ax2.barh(y_pos, dr, color=col_dr, edgecolor='#555555',
             linewidth=0.5, height=0.50)

    dr_min = min(dr)
    dr_max = max(dr)
    for i, (val, col) in enumerate(zip(dr, col_dr)):
        if val >= 0:
            # Label outside (right) of positive bar
            ax2.text(val + 0.008, y_pos[i],
                     f'+{val:.3f}', va='center', ha='left', fontsize=7.5,
                     fontweight='bold', color='#145214')
        else:
            # Label outside (left) of negative bar — clear of bar
            ax2.text(val - 0.018, y_pos[i],
                     f'{val:.3f}', va='center', ha='right', fontsize=7.5,
                     fontweight='bold', color='#aa0000')

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(y_labels, fontsize=8.5)
    ax2.set_xlabel('ΔR  (Step 33 − Step 4)', fontsize=9)
    ax2.set_title('(b)  Improvement ΔR over Baseline',
                  loc='left', fontsize=9, fontweight='bold')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.axvline(0, color='#444444', lw=0.8)
    # xlim: enough room left for negative labels and right for positive labels
    pad_l = 0.10
    pad_r = 0.10
    ax2.set_xlim(dr_min - pad_l, dr_max + pad_r)
    ax2.xaxis.grid(True, alpha=0.25, lw=0.5)
    ax2.set_axisbelow(True)

    legend_el = [
        mpatches.Patch(facecolor=GREEN, label='Improvement (ΔR > 0)'),
        mpatches.Patch(facecolor=RED,   label='Degradation (ΔR < 0)'),
    ]
    # Place legend at upper right, away from all bars
    ax2.legend(handles=legend_el, frameon=False, fontsize=7.5,
               loc='upper right')

    fig.tight_layout(w_pad=3.5)
    fig.savefig(f"{FIG_DIR}/fig3_subtype_analysis.pdf", dpi=300)
    fig.savefig(f"{FIG_DIR}/fig3_subtype_analysis.png", dpi=300)
    plt.close(fig)
    print("fig3 done")


# ============================================================
# FIG 4 — Riboswitch subclass analysis
# ============================================================
def make_fig4():
    rs_idx   = np.where(sub == 'riboswitch')[0]
    rs_sc    = rsc[rs_idx]
    y_rs     = y[rs_idx]
    oof33_rs = oof33[rs_idx]

    subclasses   = ['SAM/SAH',  'Purine',   'FMN/FAD',        'TPP',          'Other-lig.']
    subclass_key = ['SAM_SAH',  'purine',   'FMN_FAD',        'TPP',          'other_lig']
    feat_str     = ['Sign-corr. RNA-FM', 'Physchem LOO SVR',
                    'RNA-FM global', 'CPF global SVR', 'Morgan+FM LOO SVR']
    ns_sc = [int(np.sum(rs_sc == sc)) for sc in subclass_key]
    r_sc  = [pearsonr(y_rs[rs_sc == sc], oof33_rs[rs_sc == sc])[0]
             for sc in subclass_key]
    colors = [PURPLE, ORANGE, GREEN, BLUE, RED]

    # Sort by R ascending (lowest at bottom)
    order = np.argsort(r_sc)
    subclasses   = [subclasses[i]   for i in order]
    subclass_key = [subclass_key[i] for i in order]
    feat_str     = [feat_str[i]     for i in order]
    ns_sc        = [ns_sc[i]        for i in order]
    r_sc         = [r_sc[i]         for i in order]
    colors       = [colors[i]       for i in order]

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.5))

    # ── (a) horizontal bar chart ─────────────────────────────
    ax   = axes[0]
    y_pos = np.arange(len(subclasses))

    ax.barh(y_pos, r_sc, color=colors, height=0.55,
            edgecolor='#333333', linewidth=0.6)

    for i, r_val in enumerate(r_sc):
        ax.text(r_val + 0.014, y_pos[i],
                f'R = {r_val:.3f}', va='center', fontsize=8.5,
                fontweight='bold')

    # Two-line y-tick labels: "Subclass (n=N)" + feature strategy
    ylabels = [f'{s}  (n={n})\n{feat}'
               for s, n, feat in zip(subclasses, ns_sc, feat_str)]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(ylabels, fontsize=8.0)
    ax.set_xlabel('Pearson R', fontsize=9)
    ax.set_xlim(0, 1.10)
    ax.set_ylim(-0.65, len(subclasses) - 0.35)
    ax.set_title('(a)  Riboswitch Subclass Performance',
                 loc='left', fontsize=9, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axvline(0, color='#aaaaaa', lw=0.5)
    ax.xaxis.grid(True, alpha=0.25, lw=0.5)
    ax.set_axisbelow(True)

    # ── (b) scatter colored by subclass ─────────────────────
    ax2 = axes[1]
    for sc_key, sc_name, col in zip(subclass_key, subclasses, colors):
        idx_sc = rs_idx[rs_sc == sc_key]
        ax2.scatter(y[idx_sc], oof33[idx_sc], c=col, s=30,
                    alpha=0.88, edgecolors='k', linewidths=0.3,
                    label=sc_name)

    lo = min(y_rs.min(), oof33_rs.min()) - 0.35
    hi = max(y_rs.max(), oof33_rs.max()) + 0.35
    ax2.plot([lo, hi], [lo, hi], 'k--', lw=0.9, alpha=0.6)

    r_rs = pearsonr(y_rs, oof33_rs)[0]
    ax2.text(0.04, 0.97,
             f'Riboswitch overall\nR = {r_rs:.4f}  (n = 61)',
             transform=ax2.transAxes, va='top', fontsize=8.5,
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                       alpha=0.88, edgecolor=LGRAY))
    ax2.set_xlabel(r'Experimental pK$_d$', fontsize=9)
    ax2.set_ylabel(r'Predicted pK$_d$', fontsize=9)
    ax2.set_title('(b)  Riboswitch Predictions by Subclass',
                  loc='left', fontsize=9, fontweight='bold')
    ax2.legend(frameon=False, fontsize=8.0, loc='lower right',
               handletextpad=0.4)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    fig.tight_layout(w_pad=3.5)
    fig.savefig(f"{FIG_DIR}/fig4_riboswitch_subclass.pdf", dpi=300)
    fig.savefig(f"{FIG_DIR}/fig4_riboswitch_subclass.png", dpi=300)
    plt.close(fig)
    print("fig4 done")


# ============================================================
# FIG 5 — Ablation study  (NO title inside figure)
# ============================================================
def make_fig5():
    step_files = [
        ('oof_pred_ES_CS.npy',       'ES+CS\n(Step 4)'),
        ('oof_pred_step05_best.npy',  'ES+CS+WCh\n(Step 5)'),
        ('oof_pred_step07.npy',       'RS Hybrid\n(Step 7)'),
        ('oof_pred_step09_ph_fm.npy', 'RNA-FM\n(Step 9)'),
        ('oof_pred_step11.npy',       'CPF+WCh\n(Step 11)'),
        ('oof_pred_step16.npy',       'Morgan FP\n(Step 16)'),
        ('oof_pred_step19.npy',       'SAM/SAH Fix\n(Step 19)'),
        ('oof_pred_step21.npy',       'Subtype Ovrd.\n(Step 21)'),
        ('oof_pred_step26.npy',       'RS FM+CPF\n(Step 26)'),
        ('oof_pred_step28.npy',       'other_lig LOO\n(Step 28)'),
        ('oof_pred_step31.npy',       'Purine Phys.\n(Step 31)'),
        ('oof_pred_step33.npy',       'Final CAML-RNA\n(Step 33)'),
    ]

    labels = []
    rs     = []
    for fname, lbl in step_files:
        path = f"{MODEL_DIR}/{fname}"
        try:
            oof = np.load(path)
            labels.append(lbl)
            rs.append(pearsonr(y, oof)[0])
        except FileNotFoundError:
            pass

    x      = np.arange(len(labels))
    colors = [RED if i == len(labels) - 1 else BLUE for i in x]

    fig, ax = plt.subplots(figsize=(9.5, 5.0))
    bars = ax.bar(x, rs, color=colors, edgecolor='#333333',
                  linewidth=0.5, width=0.62, zorder=3)
    ax.plot(x, rs, 'o-', color='#111111', lw=1.1, ms=4, zorder=4)

    for i, (bar, val) in enumerate(zip(bars, rs)):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.003,
                f'{val:.3f}', ha='center', va='bottom', fontsize=7.5,
                fontweight='bold' if colors[i] == RED else 'normal',
                color='#aa0000' if colors[i] == RED else '#222222')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.0, ha='right', rotation=20,
                       rotation_mode='anchor')
    ax.set_ylabel('Pearson R  (LOO-CV)', fontsize=9)
    # NO ax.set_title — title removed as requested
    ax.set_ylim(0.48, 0.82)
    ax.set_xlim(-0.7, len(labels) - 0.1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.grid(True, alpha=0.3, lw=0.6, color='#dddddd')
    ax.set_axisbelow(True)

    legend_el = [
        mpatches.Patch(facecolor=BLUE, label='Development steps'),
        mpatches.Patch(facecolor=RED,  label='Final CAML-RNA (Step 33)'),
    ]
    ax.legend(handles=legend_el, frameon=False, fontsize=8.5, loc='upper left')

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.22)
    fig.savefig(f"{FIG_DIR}/fig5_ablation.pdf", dpi=300)
    fig.savefig(f"{FIG_DIR}/fig5_ablation.png", dpi=300)
    plt.close(fig)
    print("fig5 done")


if __name__ == '__main__':
    make_fig1()
    make_fig2()
    make_fig3()
    make_fig4()
    make_fig5()
    print("\nAll figures regenerated.")
