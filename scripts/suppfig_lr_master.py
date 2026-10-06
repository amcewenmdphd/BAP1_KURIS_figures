#!/usr/bin/env python3
"""
Supplementary figure: master contingency forest. Single continuous LR axis
(benign <- 1 -> pathogenic) with one row per class-assignment method x control
set; pathogenic-class LR (red) to the right, benign-class LR (blue) to the left,
on ACMG evidence-strength bands.

Input: tables/Master_Contingency_Tables.tsv (built by table_master_contingency.py).
Output: figures/SuppFig_LR_master.{png,pdf}
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, plots
import matplotlib.pyplot as plt

CONTINGENCY = config.TABLES / "Master_Contingency_Tables.tsv"

GMM = "GMM score cutoffs"
WAT = "Waters: SGE FDR + directionality"
TEJ = "Tejura: SGE FDR + directionality"

# (control set, assignment method, short row label)
DISPLAY_ROWS = [
    ("ClinVar 2025 (all variants)", GMM, "GMM · ClinVar 2025"),
    ("ClinVar Feb 2026 (all variants)", GMM, "GMM · ClinVar 2026"),
    ("ClinVar Feb 2026 (all variants)", WAT, "Waters · ClinVar 2026"),
    ("ClinVar 2025 (all variants)", TEJ, "Tejura · ClinVar 2025"),
    ("ClinVar Feb 2026 (missense only)", GMM, "GMM · ClinVar 2026"),
    ("ClinVar Feb 2026 (missense only)", WAT, "Waters · ClinVar 2026"),
    ("ClinVar 2025 (missense only)", TEJ, "Tejura · ClinVar 2025"),
    ("KURIS/NDD vs B/LB missense", GMM, "GMM"),
    ("KURIS/NDD vs B/LB missense", WAT, "Waters"),
    ("KURIS/NDD vs B/LB missense", TEJ, "Tejura"),
    ("Waters ClinVar ≥1* (Table 1)", WAT, "Waters ClinVar ≥1*"),
    ("Waters Systematic (Table 1)", WAT, "Waters Systematic"),
]
DISPLAY_GROUPS = [("All variants", 4), ("Missense only", 3), ("KURIS controls", 3), ("Waters Table 1", 2)]


def _parse_ci(s):
    lo, hi = s.strip("()").split("–")
    return float(lo), float(hi)


def _lookup(cont, truth, method_label, which):
    # every classifier plots its Abnormal-class LR (pathogenic) and Normal-class LR (benign).
    cls = "Abnormal" if which == "abnormal" else "Normal"
    row = cont[(cont["Assignment method"] == method_label) & (cont["Control set"] == truth)
               & (cont["Functional class"] == cls)]
    if row.empty or "not computable" in str(row.iloc[0]["95% CI"]):
        return None
    r = row.iloc[0]
    lo, hi = _parse_ci(r["95% CI"])
    return float(str(r["Likelihood ratio"]).rstrip("*")), lo, hi


def main():
    config.set_rcparams()
    cont = pd.read_csv(CONTINGENCY, sep="\t")
    fig, ax = plt.subplots(figsize=(180 / 25.4, 165 / 25.4))
    total = len(DISPLAY_ROWS)

    los, his = [min(plots.ACMG_TICKS)], [max(plots.ACMG_TICKS)]
    for truth, ml, _ in DISPLAY_ROWS:
        for which in ("abnormal", "normal"):
            r = _lookup(cont, truth, ml, which)
            if r:
                los.append(r[1]); his.append(r[2])
    xlim = (min(los) / 1.7, max(his) * 1.7)
    plots.forest_axis(ax, xlim)

    for i, (truth, ml, short) in enumerate(DISPLAY_ROWS):
        y = total - 1 - i
        for which, col in (("abnormal", None), ("normal", None)):
            r = _lookup(cont, truth, ml, which)
            if r:
                from bap1figs import palettes
                plots.forest_marker(ax, y, *r,
                                    palettes.PATH_MARKER if which == "abnormal" else palettes.BEN_MARKER)
    ax.set_yticks(range(total))
    ax.set_yticklabels([s for _, _, s in DISPLAY_ROWS][::-1], fontsize=7)
    ax.set_ylim(-0.6, total - 0.4)

    boundary = 0
    for name, n in DISPLAY_GROUPS[:-1]:
        boundary += n
        ax.axhline(total - 0.5 - boundary, color="#BBBBBB", lw=0.6, zorder=1)
    starts = np.cumsum([0] + [g[1] for g in DISPLAY_GROUPS])
    for (name, n), s in zip(DISPLAY_GROUPS, starts):
        ymid = total - 1 - (s + (n - 1) / 2)
        ax.text(1.015, ymid, name.split(" (")[0], transform=ax.get_yaxis_transform(),
                rotation=90, va="center", ha="left", fontsize=6.5, fontweight="bold")

    plots.forest_legends(fig, "Abnormal", "Normal")
    fig.subplots_adjust(left=0.20, right=0.94, top=0.88, bottom=0.09)
    plots.save(fig, "SuppFig_LR_master")


if __name__ == "__main__":
    main()
