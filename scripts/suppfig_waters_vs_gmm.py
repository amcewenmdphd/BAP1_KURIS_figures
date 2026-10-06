#!/usr/bin/env python3
"""
Supplementary figure: Waters standardized vs GMM-cutoff classification.
    a-c  broken-axis functional-score histograms by Waters class, ClinVar class,
         and molecular consequence
    d-i  contingency (Waters class x GMM class) for the molecular / ClinVar /
         KURIS benign and pathogenic control sets

Performance metrics (sensitivity/specificity/PPV/NPV) are in the companion
table (table_waters_vs_gmm_performance.py), not on this figure.

Step 1 rebuilds this figure's input file from the master + clinical data tables;
step 2 renders it.

Input file: data/figure_inputs/waters_vs_gmm_variants.tsv
Output: figures/SuppFig_Waters_vs_GMM.{png,pdf}
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, palettes, stats, tables

INPUT = config.FIGURE_INPUTS / "waters_vs_gmm_variants.tsv"
ABN, NRM = stats.load_cutoffs()
CV_COLORS = palettes.CLINVAR_2026_COLORS


def build_input():
    m = tables.load_master()
    m = tables.add_panel3_category(m, tables.load_clinical())
    out = pd.DataFrame({
        "score": m["MAVE_Score"], "waters": m["functional_classification"],
        "gmm": m["score_threshold_class"], "cv_class": m["clinvar_2026_slim"],
        "consequence_short": m["variant_type"], "panel3_category": m["panel3_category"]})
    out.to_csv(INPUT, sep="\t", index=False)
    print(f"INFO: wrote {INPUT.relative_to(config.REPO)}")
    return out


def _matrix(sub, waters_rows, gmm_cols):
    mat = np.zeros((3, 3), dtype=int)
    for i, wc in enumerate(waters_rows):
        for j, gc in enumerate(gmm_cols):
            mat[i, j] = int(((sub["waters"] == wc) & (sub["gmm"] == gc)).sum())
    return mat


def render(df):
    config.set_rcparams()
    fig = plt.figure(figsize=(190 / 25.4, 230 / 25.4))
    gs = fig.add_gridspec(5, 6, height_ratios=[1.35, 1.35, 1.35, 1.0, 1.0], hspace=0.85, wspace=0.40)

    from bap1figs import plots
    plots.broken_hist(fig, gs[0, :], df, "waters", ["depleted", "unchanged", "enriched"],
                      palettes.WATERS_COLORS, ABN, NRM, ncol=3, show_xticklabels=False, letter="a")
    plots.broken_hist(fig, gs[1, :], df, "cv_class", ["P/LP", "VUS", "Conflicting", "B/LB"],
                      CV_COLORS, ABN, NRM, ncol=4, show_xticklabels=False, letter="b")
    _, ax_c_bot = plots.broken_hist(fig, gs[2, :], df, "consequence_short", palettes.VARIANT_TYPE_ORDER,
                                    palettes.VARIANT_TYPE_COLORS, ABN, NRM, ncol=5,
                                    show_xticklabels=True, letter="c")
    ax_c_bot.set_xlabel("Functional score")

    waters_rows = ["depleted", "unchanged", "enriched"]
    row_labels = ["Depleted (LoF)", "Unchanged", "Enriched (GoF)"]
    gmm_cols = ["Abnormal", "Indeterminate", "Normal"]
    trunc = df["consequence_short"].isin(palettes.TRUNCATING_TYPES)
    kuris_path = df[df["panel3_category"] == "KURIS/NDD"]
    kuris_ben = df[df["panel3_category"] == "B/LB"]
    panels = [
        ("d", "Synonymous (benign)", df[df["consequence_short"] == "Synonymous"], palettes.VARIANT_TYPE_COLORS["Synonymous"], 3, 0),
        ("f", "ClinVar B/LB (benign)", df[df["cv_class"] == "B/LB"], CV_COLORS["B/LB"], 3, 2),
        ("h", "KURIS B/LB miss. (benign)", kuris_ben, CV_COLORS["B/LB"], 3, 4),
        ("e", "Truncating (path)", df[trunc], palettes.VARIANT_TYPE_COLORS["Nonsense"], 4, 0),
        ("g", "ClinVar P/LP (path)", df[df["cv_class"] == "P/LP"], CV_COLORS["P/LP"], 4, 2),
        ("i", "KURIS/NDD miss. (path)", kuris_path, CV_COLORS["P/LP"], 4, 4),
    ]
    for letter, title, sub, color, row, col0 in panels:
        ax = fig.add_subplot(gs[row, col0:col0 + 2])
        mat = _matrix(sub, waters_rows, gmm_cols)
        plots.contingency_panel(ax, row_labels, gmm_cols, mat, color, letter,
                                f"{title} (n={int(mat.sum()):,})",
                                show_y=(col0 == 0), show_x=(row == 4))
    fig.subplots_adjust(left=0.075, right=0.985, top=0.94, bottom=0.03)
    plots.save(fig, "SuppFig_Waters_vs_GMM")


def main():
    render(build_input())


if __name__ == "__main__":
    main()
