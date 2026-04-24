#!/usr/bin/env python3
"""
Generate BAP1 KURIS Figures 1 and 2 from pre-computed Excel data tables.

Input files (relative to repository root):
    data/output/BAP1_KURIS_Figure_Variants.xlsx
    data/input/BAP1_features.json

Output files:
    figures/Figure1.{png,pdf}
    figures/Figure2.{png,pdf}

Usage:
    python scripts/make_figures.py
"""

import json
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch

# ── Resolve paths relative to repo root ──
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
EXCEL_PATH = os.path.join(REPO_ROOT, 'data', 'output', 'BAP1_KURIS_Figure_Variants.xlsx')
FEATURES_PATH = os.path.join(REPO_ROOT, 'data', 'input', 'BAP1_features.json')
FIGURES_DIR = os.path.join(REPO_ROOT, 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# Constants
# ============================================================

BAP1_LENGTH = 729
SCORE_CUTOFF_ABNORMAL = -0.026891
SCORE_CUTOFF_NORMAL = -0.018367

# ── NPG figure defaults ──
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8,
    'axes.linewidth': 0.8,
    'axes.labelsize': 8,
    'axes.titlesize': 10,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 6,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'figure.dpi': 150,
    'savefig.dpi': 300,
})

# ============================================================
# Color schemes — Nature Publishing Group palette
# ============================================================

VARIANT_TYPE_COLORS = {
    'Nonsense':      '#D84315',
    'Frameshift':    '#EF6C00',
    'Splice':        '#F9A825',
    'Start lost':    '#9E9D24',
    'Stop lost':     '#827717',
    'Inframe indel': '#00695C',
    'Missense':      '#00838F',
    'Synonymous':    '#8D6E63',
    'Intronic':      '#D7CCC8',
    'UTR':           '#3E2723',
}

CLINVAR_2026_COLORS = {
    'P/LP':        '#E64B35',
    'B/LB':        '#3C5488',
    'VUS':         '#B8B8B8',
    'Conflicting': '#7B2D8E',
}

PANEL3_COLORS = {
    'KURIS/NDD':       '#E8755A',
    'TPDS P/LP':       '#722F37',
    'N229K (proband)': '#FFC107',
    'B/LB':            '#3C5488',
}

LOLLIPOP_COLORS = PANEL3_COLORS

COND_COLORS = {
    'Cancer predisposition': '#722F37',
    'Kury-Isidor syndrome':  '#E8755A',
    'Not provided':          '#B8B8B8',
}

SCORE_CLASS_COLORS = {
    'Abnormal':      '#F4A460',
    'Indeterminate': '#B0B0B0',
    'Normal':        '#80CBC4',
}

_EV_COLOR = {
    'PS3 V.Strong': '#E64B35', 'PS3 Strong': '#F08070',
    'PS3 Moderate': '#F5AFA5', 'PS3 Supptg.': '#FDDEDE',
    'BS3 V.Strong': '#3C5488', 'BS3 Strong': '#6B82B0',
    'BS3 Moderate': '#9AB0D8', 'BS3 Supptg.': '#DEE8F5',
    'Indet.': '#E0E0E0',
}

# ============================================================
# Load data
# ============================================================

def load_data():
    """Load all data from Excel and JSON files."""
    sheets = {}
    for sheet in ['Fig1a_Lollipop', 'Fig1b_ClinVar_All', 'Fig1c_PLP_Condition',
                  'Fig2a_Consequence_Hist', 'Fig2b_ClinVar_PLP_BLB',
                  'Fig2c_KURIS_Hist', 'Fig2d_Table_All', 'Fig2e_Table_KURIS']:
        sheets[sheet] = pd.read_excel(EXCEL_PATH, sheet_name=sheet)

    with open(FEATURES_PATH) as f:
        uniprot = json.load(f)

    domain_map = {
        'UCH catalytic':          ('#90CAF9', 'UCH'),
        'ULD':                    ('#CE93D8', 'ULD'),
        'HBM-like motif':        ('#FFE082', 'HBM'),
        'Interaction with BRCA1': ('#FFAB91', 'BRCA1 binding'),
    }
    domains = []
    for feat in uniprot.get('features', []):
        desc = feat.get('description', '')
        if desc in domain_map:
            color, label = domain_map[desc]
            domains.append({
                'name': label,
                'start': feat['location']['start']['value'],
                'end': feat['location']['end']['value'],
                'color': color,
            })
    domains.sort(key=lambda d: d['start'])

    return sheets, domains


# ============================================================
# LR / evidence helpers
# ============================================================

def _lr_evidence(lr):
    """Map a likelihood ratio to ACMG evidence strength label."""
    if lr >= 18.7:
        return 'PS3 Strong'
    elif lr >= 4.33:
        return 'PS3 Moderate'
    elif lr >= 2.08:
        return 'PS3 Supptg.'
    elif lr <= 1 / 18.7:
        return 'BS3 Strong'
    elif lr <= 1 / 4.33:
        return 'BS3 Moderate'
    elif lr <= 1 / 2.08:
        return 'BS3 Supptg.'
    else:
        return 'Indet.'


def compute_contingency(df_table, path_col, path_val, benign_val):
    """Compute 3-class contingency stats from a table dataframe."""
    df_path = df_table[df_table[path_col] == path_val]
    df_benign = df_table[df_table[path_col] == benign_val]
    n1, n2 = len(df_path), len(df_benign)

    rows = []
    for sc in ['Abnormal', 'Indeterminate', 'Normal']:
        a = (df_path['Score_Threshold_Class'] == sc).sum()
        b = (df_benign['Score_Threshold_Class'] == sc).sum()
        # Haldane–Anscombe correction
        a_c, b_c, n1_c, n2_c = a + 0.5, b + 0.5, n1 + 0.5, n2 + 0.5
        lr = (a_c / n1_c) / (b_c / n2_c)
        rows.append({
            'score_class': sc, 'path': a, 'benign': b,
            'lr': lr, 'corrected': (a == 0 or b == 0),
            'evidence': _lr_evidence(lr),
        })
    return {'rows': rows, 'n_path': n1, 'n_benign': n2}


# ============================================================
# Panel plotting functions
# ============================================================

def plot_lollipop(ax, df_lolli, domains):
    """Panel a: BAP1 protein lollipop plot with staggered labels."""
    backbone_y, backbone_h = 0, 0.08
    ruler_y = backbone_y - 0.09

    # Backbone
    ax.add_patch(FancyBboxPatch(
        (0, backbone_y - backbone_h / 2), BAP1_LENGTH + 1, backbone_h,
        boxstyle='round,pad=0.002', facecolor='#D9D9D9',
        edgecolor='black', linewidth=0.6, zorder=2))

    # Domains
    for d in domains:
        ax.add_patch(FancyBboxPatch(
            (d['start'], backbone_y - backbone_h / 2),
            d['end'] - d['start'], backbone_h,
            boxstyle='round,pad=0.002', facecolor=d['color'],
            edgecolor='black', linewidth=0.6, zorder=3))

    # Ruler
    tick_y = ruler_y - 0.02
    for pos in [1, 100, 200, 300, 400, 500, 600, 700, BAP1_LENGTH]:
        ax.plot([pos, pos], [tick_y - 0.015, tick_y + 0.015],
                color='black', linewidth=0.5, zorder=2)
        ax.text(pos, tick_y - 0.025, str(pos),
                ha='center', va='top', fontsize=5, color='#333')
    ax.plot([1, BAP1_LENGTH], [tick_y, tick_y], color='black',
            linewidth=0.5, zorder=2)

    # Sort and micro-offset overlapping positions
    df_s = df_lolli.sort_values('AA_Position').reset_index(drop=True)
    n_var = len(df_s)
    pos_counts = df_s['AA_Position'].value_counts()
    pos_offset_idx = {}
    x_pos = []
    for i in range(n_var):
        p = df_s.loc[i, 'AA_Position']
        n_at = pos_counts[p]
        if p not in pos_offset_idx:
            pos_offset_idx[p] = 0
        cur = pos_offset_idx[p]
        pos_offset_idx[p] += 1
        if n_at == 1:
            x_pos.append(p)
        else:
            spread = min(18, n_at * 6)
            x_pos.append(p - spread / 2 + cur * spread / max(n_at - 1, 1))

    # Tier assignment (greedy, pack labels into fewest tiers)
    h_base, h_step = 0.12, 0.055
    char_w, pad_x = 5.5, 5
    tier_ranges = {}
    assigned_tiers = []
    for i in range(n_var):
        w = len(df_s.loc[i, 'AA_Change']) * char_w + pad_x
        xl, xr = x_pos[i] - 2, x_pos[i] + w
        placed = False
        for tier in range(50):
            if tier not in tier_ranges:
                tier_ranges[tier] = []
            if not any(xl < r and xr > l for l, r in tier_ranges[tier]):
                tier_ranges[tier].append((xl, xr))
                assigned_tiers.append(tier)
                placed = True
                break
        if not placed:
            assigned_tiers.append(50)

    heights = [h_base + t * h_step for t in assigned_tiers]

    # Draw stems and labels
    for i in range(n_var):
        row = df_s.loc[i]
        x = x_pos[i]
        cat = row['Category']
        change = row['AA_Change']
        is_star = bool(row.get('Is_N229K', False))
        color = LOLLIPOP_COLORS.get(cat, '#888888')
        h = heights[i]

        line_color = 'black' if cat == 'N229K (proband)' else color
        ax.plot([x, x], [backbone_y + backbone_h / 2, h],
                color=line_color, linewidth=0.5, zorder=5)

        if is_star:
            ax.plot(x, h, '*', markersize=8, color=color,
                    markeredgecolor='black', markeredgewidth=1.0, zorder=7)
        elif cat == 'TPDS P/LP':
            ax.plot(x, h, 's', markersize=2.5, color=color,
                    markeredgecolor='black', markeredgewidth=0.3, zorder=6)
        else:
            ew = 0.8 if cat == 'N229K (proband)' else 0.3
            ax.plot(x, h, 'o', markersize=2.5, color=color,
                    markeredgecolor='black', markeredgewidth=ew, zorder=6)

        ax.text(x + 8, h, change, ha='left', va='center', fontsize=4.5,
                fontweight='bold', color=line_color, zorder=8,
                bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                          edgecolor='black', linewidth=0.3, alpha=0.95))

    actual_max = max(heights) if heights else 1.0
    ax.set_xlim(-15, BAP1_LENGTH + 55)
    ax.set_ylim(tick_y - 0.04, actual_max + h_step + 0.03)
    ax.axis('off')

    # Legend
    handles = []
    for cat, marker in [('N229K (proband)', '*'), ('KURIS/NDD', 'o'),
                         ('TPDS P/LP', 's'), ('B/LB', 'o')]:
        if cat in df_lolli['Category'].values:
            ms = 6 if cat == 'N229K (proband)' else 3.5
            ew = 1.0 if cat == 'N229K (proband)' else 0.3
            handles.append(mlines.Line2D(
                [], [], marker=marker, color=LOLLIPOP_COLORS[cat],
                linestyle='-', linewidth=0.8, markersize=ms,
                markeredgecolor='black', markeredgewidth=ew, label=cat))
    for d in domains:
        handles.append(mpatches.Patch(
            facecolor=d['color'], edgecolor='black', linewidth=0.4,
            label=d['name']))
    handles.append(mpatches.Patch(
        facecolor='#D9D9D9', edgecolor='black', linewidth=0.4, label='Other'))
    ax.legend(handles=handles, fontsize=5, loc='lower center',
              bbox_to_anchor=(0.5, 1.0), framealpha=0.95, edgecolor='#CCC',
              ncol=len(handles), handlelength=1.2, handletextpad=0.3,
              labelspacing=0.2, columnspacing=0.6)


def plot_clinvar_bars(ax, df_cv, cons_order, color_map, group_col, label_col=None):
    """Vertical stacked bar chart of variant counts by consequence."""
    if label_col is None:
        label_col = group_col
    groups = df_cv[group_col].dropna().unique()
    x_pos = np.arange(len(cons_order))
    bottom = np.zeros(len(cons_order))
    for grp in color_map:
        if grp not in groups:
            continue
        subset = df_cv[df_cv[group_col] == grp]
        counts = np.array([(subset['Consequence'] == c).sum() for c in cons_order])
        if counts.sum() > 0:
            ax.bar(x_pos, counts, bottom=bottom, width=0.65,
                   color=color_map[grp], edgecolor='white', linewidth=0.3,
                   label=f'{grp} (n={counts.sum()})')
            bottom += counts
    ax.set_xticks(x_pos)
    ax.set_xticklabels(cons_order, rotation=45, ha='right', fontsize=5)
    ax.set_ylabel('Count', fontsize=6)
    ax.legend(fontsize=4.5, loc='upper right', framealpha=0.9, edgecolor='none')
    ax.tick_params(axis='both', labelsize=5.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def plot_histogram(ax, df_hist, group_col, group_order, color_map,
                   bin_width=0.005, show_xlabel=True, show_ylabel=True,
                   show_threshold=True):
    """Stacked histogram of functional scores by category."""
    bins = np.arange(-0.30, 0.035, bin_width)
    present = [g for g in group_order if g in df_hist[group_col].unique()]

    data = [df_hist[df_hist[group_col] == g]['MAVE_Score'].dropna().values
            for g in present]
    colors = [color_map.get(g, '#888888') for g in present]

    ax.hist(data, bins=bins, stacked=True, color=colors,
            edgecolor='white', linewidth=0.3)

    if show_threshold:
        ax.axvline(x=SCORE_CUTOFF_ABNORMAL, color='black', linestyle='--', linewidth=0.8)
        ax.axvline(x=SCORE_CUTOFF_NORMAL, color='black', linestyle=':', linewidth=0.8)

    if show_xlabel:
        ax.set_xlabel('Functional Score', fontsize=6)
    if show_ylabel:
        ax.set_ylabel('Count', fontsize=6)
    ax.set_xlim(-0.30, 0.035)

    handles = []
    for g in present:
        n = (df_hist[group_col] == g).sum()
        handles.append(mpatches.Patch(color=color_map[g], label=f'{g} ({n:,})'))
    if show_threshold:
        handles.append(mlines.Line2D([], [], color='black', linestyle='--',
                                     linewidth=0.8, label=f'Abnormal ({SCORE_CUTOFF_ABNORMAL})'))
        handles.append(mlines.Line2D([], [], color='black', linestyle=':',
                                     linewidth=0.8, label=f'Normal ({SCORE_CUTOFF_NORMAL})'))
    ax.legend(handles=handles, fontsize=4.5, loc='upper left',
              framealpha=0.9, edgecolor='none')


def plot_contingency_table(ax, stats, title):
    """Contingency table with Abnormal/Indeterminate/Normal rows."""
    ax.axis('off')
    rows_data = stats['rows']
    n1, n2 = stats['n_path'], stats['n_benign']

    col_labels = ['Classification', 'P/LP', 'B/LB', 'LR', 'Evidence']
    table_data = []
    cell_colors = []

    _plp_tint, _blb_tint = '#FDCBA4', '#B2DFDB'

    for r in rows_data:
        star = '*' if r['corrected'] else ''
        lr_str = f'{r["lr"]:.2f}{star}'
        name = r['score_class']
        ev = r['evidence'] if name != 'Indeterminate' else ''
        table_data.append([name, str(r['path']), str(r['benign']), lr_str, ev])

        cls_bg = '#B0B0B0' if name == 'Indeterminate' else SCORE_CLASS_COLORS.get(name, '#FFFFFF')
        plp_bg = _plp_tint if name == 'Abnormal' else '#FFFFFF'
        blb_bg = _blb_tint if name == 'Normal' else '#FFFFFF'
        ev_bg = _EV_COLOR.get(ev, '#FFFFFF') if ev else '#FFFFFF'
        cell_colors.append([cls_bg, plp_bg, blb_bg, '#FFFFFF', ev_bg])

    table_data.append(['Total', str(n1), str(n2), '', ''])
    cell_colors.append(['#E8E8E8'] * 5)

    header_colors = ['#444444', '#E64B35', '#3C5488', '#444444', '#444444']
    tbl = ax.table(cellText=table_data, colLabels=col_labels,
                   cellColours=cell_colors, colColours=header_colors,
                   loc='lower center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(6)
    tbl.scale(1.0, 1.0)
    tbl.auto_set_column_width(col=list(range(5)))

    for j in range(5):
        tbl[0, j].get_text().set_color('white')
        tbl[0, j].get_text().set_fontsize(6)
    for i in range(1, len(table_data) + 1):
        for j in range(5):
            tbl[i, j].get_text().set_fontsize(5.5)

    # Title centered above table
    ax.text(0.5, 0.62, title, fontsize=6, fontweight='bold',
            transform=ax.transAxes, ha='center', va='bottom')


# ============================================================
# Figure 1
# ============================================================

def make_figure1(sheets, domains):
    """Figure 1: lollipop + ClinVar bars + P/LP by condition."""
    fig = plt.figure(figsize=(7.09, 4.72))
    gs = GridSpec(2, 2, figure=fig,
                  width_ratios=[1.0, 1.0],
                  height_ratios=[0.55, 0.45],
                  hspace=0.45, wspace=0.35)

    # Panel a: lollipop
    ax_a = fig.add_subplot(gs[0, :])
    plot_lollipop(ax_a, sheets['Fig1a_Lollipop'], domains)
    ax_a.text(-0.02, 1.0, 'a', fontsize=8, fontweight='bold',
              transform=ax_a.transAxes, va='bottom', ha='right')

    # Shared consequence order
    cons_order = ['Frameshift', 'Nonsense', 'Splice', 'Missense', 'Start lost',
                  'Inframe indel', 'Stop lost', 'Synonymous', 'Intronic', 'UTR']
    present = (set(sheets['Fig1b_ClinVar_All']['Consequence'].unique()) |
               set(sheets['Fig1c_PLP_Condition']['Consequence'].unique()))
    cons_order = [c for c in cons_order if c in present]

    # Panel b: ClinVar all
    ax_b = fig.add_subplot(gs[1, 0])
    plot_clinvar_bars(ax_b, sheets['Fig1b_ClinVar_All'], cons_order,
                      CLINVAR_2026_COLORS, 'ClinVar_Class')
    ax_b.text(-0.12, 1.05, 'b', fontsize=8, fontweight='bold',
              transform=ax_b.transAxes, va='bottom', ha='right')

    # Panel c: P/LP by condition
    ax_c = fig.add_subplot(gs[1, 1])
    plot_clinvar_bars(ax_c, sheets['Fig1c_PLP_Condition'], cons_order,
                      COND_COLORS, 'Condition')
    ax_c.text(-0.12, 1.05, 'c', fontsize=8, fontweight='bold',
              transform=ax_c.transAxes, va='bottom', ha='right')

    base = os.path.join(FIGURES_DIR, 'Figure1')
    plt.savefig(f'{base}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{base}.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f'Figure 1 saved: {base}.{{png,pdf}}')


# ============================================================
# Figure 2
# ============================================================

def make_figure2(sheets):
    """Figure 2: histograms (shared x) + contingency tables."""
    fig = plt.figure(figsize=(7.09, 5.91))
    gs = GridSpec(3, 2, figure=fig,
                  width_ratios=[1.3, 0.7],
                  height_ratios=[1, 1, 1],
                  hspace=0.08, wspace=0.25)

    # Variant type order for panel a
    vtype_order = ['Nonsense', 'Frameshift', 'Splice', 'Start lost', 'Stop lost',
                   'Inframe indel', 'Missense', 'Synonymous', 'Intronic', 'UTR']

    # Panel a: consequence histogram (full MAVE data)
    ax_a = fig.add_subplot(gs[0, 0])
    plot_histogram(ax_a, sheets['Fig2a_Consequence_Hist'],
                   'Variant_Type', vtype_order, VARIANT_TYPE_COLORS,
                   show_xlabel=False)
    ax_a.tick_params(axis='x', labelbottom=False)
    ax_a.text(-0.08, 1.02, 'a', fontsize=8, fontweight='bold',
              transform=ax_a.transAxes, va='top', ha='right')

    ax_empty = fig.add_subplot(gs[0, 1])
    ax_empty.axis('off')

    # Panel b: P/LP + B/LB histogram
    ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)
    cv_order = ['B/LB', 'P/LP']
    plot_histogram(ax_b, sheets['Fig2b_ClinVar_PLP_BLB'],
                   'ClinVar_Class', cv_order, CLINVAR_2026_COLORS,
                   show_xlabel=False)
    ax_b.tick_params(axis='x', labelbottom=False)
    ax_b.text(-0.08, 1.02, 'b', fontsize=8, fontweight='bold',
              transform=ax_b.transAxes, va='top', ha='right')

    # Panel d: contingency table — all variants
    ax_d = fig.add_subplot(gs[1, 1])
    stats_all = compute_contingency(sheets['Fig2d_Table_All'],
                                    'ClinVar_Class', 'P/LP', 'B/LB')
    plot_contingency_table(ax_d, stats_all, 'All Variants (ClinVar Feb 2026)')
    ax_d.text(0.0, 0.68, 'd', fontsize=8, fontweight='bold',
              transform=ax_d.transAxes, va='bottom', ha='left')

    # Panel c: KURIS/NDD histogram
    ax_c = fig.add_subplot(gs[2, 0], sharex=ax_a)
    p3_order = ['B/LB', 'TPDS P/LP', 'KURIS/NDD', 'N229K (proband)']
    plot_histogram(ax_c, sheets['Fig2c_KURIS_Hist'],
                   'Category', p3_order, PANEL3_COLORS,
                   show_xlabel=True)
    ax_c.tick_params(axis='x', labelbottom=True)
    ax_c.text(-0.08, 1.02, 'c', fontsize=8, fontweight='bold',
              transform=ax_c.transAxes, va='top', ha='right')

    # Panel e: contingency table — KURIS missense
    ax_e = fig.add_subplot(gs[2, 1])
    stats_kuris = compute_contingency(sheets['Fig2e_Table_KURIS'],
                                      'Category', 'KURIS/NDD', 'B/LB')
    plot_contingency_table(ax_e, stats_kuris, 'KURIS/NDD Missense')
    ax_e.text(0.0, 0.68, 'e', fontsize=8, fontweight='bold',
              transform=ax_e.transAxes, va='bottom', ha='left')

    fig.subplots_adjust(bottom=0.10)

    base = os.path.join(FIGURES_DIR, 'Figure2')
    plt.savefig(f'{base}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{base}.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f'Figure 2 saved: {base}.{{png,pdf}}')


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    print('Loading data...')
    sheets, domains = load_data()
    print(f'  Loaded {len(sheets)} sheets, {len(domains)} domains')

    make_figure1(sheets, domains)
    make_figure2(sheets)
    print('Done.')
