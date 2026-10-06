#!/usr/bin/env python3
"""
Figure 2: SGE functional scores + likelihood-ratio contingency tables.
    a  functional scores by molecular consequence (all scored variants)
    b  functional scores for ClinVar P/LP vs B/LB
    c  functional scores for the KURIS/NDD, TPDS P/LP, B/LB, N229K missense sets
    d  variant-reclassification schematic (N229K: VUS under ACMG evidence
       alone vs. Likely pathogenic once PS3_strong functional evidence is added)
    e  contingency (LR + 95% CI + ACMG tier): all ClinVar P/LP vs B/LB
    f  contingency: KURIS/NDD vs B/LB missense

Step 1 rebuilds this figure's own input file from the master + clinical data
tables; step 2 renders the figure from that input file.

Input file: data/figure_inputs/figure2_variants.xlsx
Output: figures/Figure2.{png,pdf}
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, palettes, stats, tables

INPUT = config.FIGURE_INPUTS / "figure2_variants.xlsx"
ABN, NRM = stats.load_cutoffs()


# ── Step 1: rebuild the per-figure input file from the source tables ──

def build_input():
    m = tables.load_master()
    m = tables.add_panel3_category(m, tables.load_clinical())
    fig2a = (m.dropna(subset=["MAVE_Score"])
             [["HGVS_p", "HGVS_c", "variant_type", "MAVE_Score", "score_threshold_class",
               "functional_classification"]]
             .rename(columns={"variant_type": "Variant_Type", "score_threshold_class": "Score_Threshold_Class",
                              "functional_classification": "Original_Paper_Class"}).sort_values("MAVE_Score"))
    cv = m[m["clinvar_2026_slim"].isin(["P/LP", "B/LB"])]
    fig2b = (cv.dropna(subset=["MAVE_Score"])
             [["HGVS_p", "HGVS_c", "clinvar_2026_slim", "MAVE_Score", "score_threshold_class", "variant_type"]]
             .rename(columns={"clinvar_2026_slim": "ClinVar_Class", "score_threshold_class": "Score_Threshold_Class",
                              "variant_type": "Variant_Type"}).sort_values(["ClinVar_Class", "MAVE_Score"]))
    fig2d = (cv.dropna(subset=["MAVE_Score"])
             [["HGVS_p", "HGVS_c", "clinvar_2026_slim", "MAVE_Score", "score_threshold_class",
               "functional_classification", "variant_type"]]
             .rename(columns={"clinvar_2026_slim": "ClinVar_Class", "score_threshold_class": "Score_Threshold_Class",
                              "functional_classification": "Original_Paper_Class", "variant_type": "Variant_Type"})
             .sort_values(["ClinVar_Class", "Score_Threshold_Class"]))
    cats = (m.dropna(subset=["panel3_category"])
            .drop_duplicates(subset=["AA_Position", "ref_aa", "alt_aa", "panel3_category"]))
    fig2c = (cats[["HGVS_p", "HGVS_c", "AA_Position", "ref_aa", "alt_aa", "panel3_category",
                   "MAVE_Score", "score_threshold_class"]]
             .rename(columns={"ref_aa": "AA_Ref", "alt_aa": "AA_Alt", "panel3_category": "Category",
                              "score_threshold_class": "Score_Threshold_Class"})
             .sort_values(["Category", "AA_Position"]))
    e = (m[m["panel3_category"].isin(["KURIS/NDD", "B/LB"])]
         .drop_duplicates(subset=["AA_Position", "ref_aa", "alt_aa", "panel3_category"])
         .dropna(subset=["score_threshold_class"]))
    fig2e = (e[["HGVS_p", "HGVS_c", "AA_Position", "panel3_category", "MAVE_Score",
                "score_threshold_class", "functional_classification"]]
             .rename(columns={"panel3_category": "Category", "score_threshold_class": "Score_Threshold_Class",
                              "functional_classification": "Original_Paper_Class"})
             .sort_values(["Category", "Score_Threshold_Class"]))
    sheets = {"Fig2a_Consequence_Hist": fig2a, "Fig2b_ClinVar_PLP_BLB": fig2b,
              "Fig2c_KURIS_Hist": fig2c, "Fig2d_Table_All": fig2d, "Fig2e_Table_KURIS": fig2e}
    with pd.ExcelWriter(INPUT, engine="openpyxl") as w:
        for name, df in sheets.items():
            df.to_excel(w, sheet_name=name, index=False)
    print(f"INFO: wrote {INPUT.relative_to(config.REPO)}")
    return {k: v for k, v in sheets.items()}


# ── Step 2: render ──

def _recompute_class(df):
    score = pd.to_numeric(df["MAVE_Score"], errors="coerce")
    cls = pd.Series("Indeterminate", index=df.index)
    cls[score <= ABN] = "Abnormal"; cls[score >= NRM] = "Normal"
    cls[score.isna()] = df["Score_Threshold_Class"][score.isna()]
    df = df.copy(); df["Score_Threshold_Class"] = cls.values
    return df


def plot_histogram(ax, df, group_col, order, color_map, show_xlabel=True):
    bins = np.arange(-0.30, 0.035, 0.005)
    present = [g for g in order if g in df[group_col].unique()]
    data = [df[df[group_col] == g]["MAVE_Score"].dropna().values for g in present]
    ax.hist(data, bins=bins, stacked=True, color=[color_map.get(g, "#888") for g in present],
            edgecolor="white", linewidth=0.3)
    ax.axvline(ABN, color="black", ls="--", lw=0.8); ax.axvline(NRM, color="black", ls=":", lw=0.8)
    if show_xlabel:
        ax.set_xlabel("Functional Score", fontsize=6)
    ax.set_ylabel("Count", fontsize=6); ax.set_xlim(-0.30, 0.035)
    handles = [mpatches.Patch(color=color_map[g], label=f"{g} ({(df[group_col] == g).sum():,})") for g in present]
    handles += [mlines.Line2D([], [], color="black", ls="--", lw=0.8, label=f"Abnormal ({ABN:.6f})"),
                mlines.Line2D([], [], color="black", ls=":", lw=0.8, label=f"Normal ({NRM:.6f})")]
    ax.legend(handles=handles, fontsize=4.5, loc="upper left", framealpha=0.9, edgecolor="none")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def contingency(df, path_col, path_val, benign_val, class_order):
    dp, db = df[df[path_col] == path_val], df[df[path_col] == benign_val]
    n1, n2 = len(dp), len(db)
    rows = []
    for sc in class_order:
        a = int((dp["Score_Threshold_Class"] == sc).sum()); b = int((db["Score_Threshold_Class"] == sc).sum())
        lr, lo, hi, corr = stats.lr_and_ci(a, b, n1, n2)
        rows.append(dict(score_class=sc, path=a, benign=b, lr=lr, lo=lo, hi=hi,
                         corrected=corr, evidence=stats.evidence_tier(lr, "compact")))
    return rows, n1, n2


def _fmt_lr(v):
    """3 decimals below 1 (so a halved value like 0.0056 reads as '0.006', not
    a '.2f'-rounded '0.01'), 2 decimals from 1 up. Matches the tables' own
    adaptive precision (bap1figs.stats.fmt_lr) in spirit, tuned for this
    panel's tighter layout."""
    return f"{v:.3f}" if v < 1 else f"{v:.2f}"


def plot_contingency(ax, rows, n1, n2, title, indet_label):
    ax.axis("off")
    plp_tint, blb_tint = "#FDCBA4", "#B2DFDB"
    data, colors = [], []
    for r in rows:
        star = "*" if r["corrected"] else ""
        name = r["score_class"]
        lr_str = "" if name == indet_label else f'{_fmt_lr(r["lr"])}{star}\n({_fmt_lr(r["lo"])}–{_fmt_lr(r["hi"])})'
        ev = r["evidence"] if name != indet_label else ""
        data.append([name, str(r["path"]), str(r["benign"]), lr_str, ev])
        cls_bg = ("#B0B0B0" if name in (indet_label, "unchanged") else
                  "#F4A460" if name in ("Abnormal", "depleted") else
                  "#80CBC4" if name in ("Normal", "enriched") else "#FFFFFF")
        colors.append([cls_bg, plp_tint if name in ("Abnormal", "depleted") else "#FFFFFF",
                       blb_tint if name in ("Normal", "enriched") else "#FFFFFF", "#FFFFFF",
                       palettes.EVIDENCE_CELL_COLORS.get(ev, "#FFFFFF") if ev else "#FFFFFF"])
    data.append(["Total", str(n1), str(n2), "", ""]); colors.append(["#E8E8E8"] * 5)
    tbl = ax.table(cellText=data, colLabels=["Classification", "P/LP", "B/LB", "LR (95% CI)", "Evidence"],
                   cellColours=colors, colColours=["#444444", "#E64B35", "#3C5488", "#444444", "#444444"],
                   loc="lower center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(6); tbl.scale(1.0, 1.3); tbl.auto_set_column_width(col=list(range(5)))
    for j in range(5):
        tbl[0, j].get_text().set_color("white"); tbl[0, j].get_text().set_fontsize(6)
    for i in range(1, len(data) + 1):
        for j in range(5):
            tbl[i, j].get_text().set_fontsize(5 if j == 3 else 5.5)
    ax.text(0.5, 0.86, title, fontsize=6, fontweight="bold", transform=ax.transAxes, ha="center", va="bottom")


def _schem_box(ax, xc, yc, w, h, title, subtitle, color):
    ax.add_patch(mpatches.FancyBboxPatch(
        (xc - w / 2, yc - h / 2), w, h, boxstyle="round,pad=0,rounding_size=0.012",
        facecolor=color, edgecolor="black", linewidth=0.6, mutation_aspect=1))
    if not subtitle:
        ax.text(xc, yc, title, ha="center", va="center", fontsize=4.8, fontweight="bold")
        return
    # Vertically center the title+subtitle block as a whole (rather than anchoring
    # to the box top) so taller boxes with more subtitle lines don't push text
    # toward the bottom edge.
    n_lines = subtitle.count("\n") + 1
    title_h, sub_line_h, gap = 0.050, 0.040, 0.014
    block_h = title_h + gap + n_lines * sub_line_h
    top_y = yc + block_h / 2
    ax.text(xc, top_y, title, ha="center", va="top", fontsize=4.8, fontweight="bold")
    ax.text(xc, top_y - title_h - gap, subtitle, ha="center", va="top",
            fontsize=3.6, style="italic", linespacing=1.25)


def _schem_arrow(ax, x, y_top, y_tip, style):
    """Dashed/dotted shaft with a short solid segment for the arrowhead --
    applying the dash/dot linestyle to the whole FancyArrowPatch (shaft +
    head) makes the head itself render broken/warped, so the head is drawn
    separately as a solid mini-arrow."""
    head_len = 0.03
    ax.plot([x, x], [y_top, y_tip + head_len], color="black", lw=0.9, linestyle=style)
    ax.annotate("", xy=(x, y_tip), xytext=(x, y_tip + head_len),
                arrowprops=dict(arrowstyle="-|>", lw=0.9, color="black", mutation_scale=7))


def plot_reclass_schematic(ax):
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ev = palettes.EVIDENCE_CELL_COLORS
    xL, xR, xPlus = 0.255, 0.765, 0.52
    ax.text(xL, 0.985, "Prior classification", ha="center", va="top", fontsize=5, fontweight="bold")
    ax.text(xR, 0.985, "Reclassification using\nfunctional evidence", ha="center", va="top",
            fontsize=5, fontweight="bold", linespacing=1.3)
    _schem_box(ax, xL, 0.81, 0.47, 0.20, "PS2_moderate",
               "(De novo, confirmed\nwith non-specific\nphenotype)", ev["PS3 Moderate"])
    ax.text(xL, 0.675, "+", ha="center", va="center", fontsize=7.5, fontweight="bold")
    _schem_box(ax, xL, 0.56, 0.47, 0.16, "PS4_supp",
               "(observed in ≥2\nadditional individuals)", ev["PS3 Supptg."])
    ax.text(xL, 0.445, "+", ha="center", va="center", fontsize=7.5, fontweight="bold")
    _schem_box(ax, xL, 0.33, 0.47, 0.16, "PM2_supp",
               "(absent in pop.\ndatabases)", ev["PS3 Supptg."])
    _schem_box(ax, xR, 0.56, 0.43, 0.16, "PS3_strong",
               "(Functionally abnormal\nscore: −0.20017)", ev["PS3 Strong"])
    ax.text(xPlus, 0.56, "+", ha="center", va="center", fontsize=7.5, fontweight="bold")
    _schem_arrow(ax, xL, 0.25, 0.155, "dashed")
    _schem_arrow(ax, xR, 0.48, 0.155, "dotted")
    _schem_box(ax, xL, 0.08, 0.30, 0.12, "VUS", None, ev["Indet."])
    _schem_box(ax, xR, 0.08, 0.40, 0.12, "Likely pathogenic", None, ev["PS3 Moderate"])


def add_proband_arrow(ax, x, y_tip, y_top):
    """Points at the N229K proband's bar (bin [-0.205, -0.200), stacked height 2)."""
    ax.annotate("", xy=(x, y_tip), xytext=(x, y_top),
                arrowprops=dict(arrowstyle="-|>", lw=1.1, color="black", mutation_scale=8))


def _label(ax, letter, va, ha, x, y):
    ax.text(x, y, letter, fontsize=10, fontweight="bold", transform=ax.transAxes, va=va, ha=ha)


def render(sheets):
    config.set_rcparams()
    sheets = {k: _recompute_class(v) if "Score_Threshold_Class" in v.columns else v for k, v in sheets.items()}
    fig = plt.figure(figsize=(7.09, 5.91))
    gs = GridSpec(3, 2, figure=fig, width_ratios=[1.3, 0.7], height_ratios=[1, 1, 1], hspace=0.08, wspace=0.25)
    ax_a = fig.add_subplot(gs[0, 0])
    plot_histogram(ax_a, sheets["Fig2a_Consequence_Hist"], "Variant_Type", palettes.VARIANT_TYPE_ORDER,
                   palettes.VARIANT_TYPE_COLORS, show_xlabel=False)
    ax_a.tick_params(axis="x", labelbottom=False); _label(ax_a, "a", "top", "right", -0.08, 1.02)
    ax_d = fig.add_subplot(gs[0, 1])
    plot_reclass_schematic(ax_d)
    _label(ax_d, "d", "top", "right", -0.08, 1.02)
    ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)
    plot_histogram(ax_b, sheets["Fig2b_ClinVar_PLP_BLB"], "ClinVar_Class", ["B/LB", "P/LP"],
                   palettes.CLINVAR_2026_COLORS, show_xlabel=False)
    ax_b.tick_params(axis="x", labelbottom=False); _label(ax_b, "b", "top", "right", -0.08, 1.02)
    ax_e = fig.add_subplot(gs[1, 1])
    plot_contingency(ax_e, *contingency(sheets["Fig2d_Table_All"], "ClinVar_Class", "P/LP", "B/LB",
                                        ["Abnormal", "Indeterminate", "Normal"]),
                     "All Variants (ClinVar Feb 2026)", "Indeterminate")
    _label(ax_e, "e", "bottom", "left", 0.0, 0.92)
    ax_c = fig.add_subplot(gs[2, 0], sharex=ax_a)
    plot_histogram(ax_c, sheets["Fig2c_KURIS_Hist"], "Category",
                   ["B/LB", "TPDS P/LP", "KURIS/NDD", "N229K (proband)"], palettes.PANEL3_COLORS, show_xlabel=True)
    add_proband_arrow(ax_c, -0.2025, 2.15, 4.4)
    _label(ax_c, "c", "top", "right", -0.08, 1.02)
    ax_f = fig.add_subplot(gs[2, 1])
    plot_contingency(ax_f, *contingency(sheets["Fig2e_Table_KURIS"], "Category", "KURIS/NDD", "B/LB",
                                        ["Abnormal", "Indeterminate", "Normal"]),
                     "KURIS/NDD Missense", "Indeterminate")
    _label(ax_f, "f", "bottom", "left", 0.0, 0.92)
    fig.subplots_adjust(bottom=0.10)
    for ext in ("png", "pdf"):
        fig.savefig(config.FIGURES / f"Figure2.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("INFO: wrote figures/Figure2.{png,pdf}")


def main():
    render(build_input())


if __name__ == "__main__":
    main()
