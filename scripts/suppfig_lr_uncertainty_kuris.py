#!/usr/bin/env python3
"""
Supplementary figure: KURIS/NDD LR stability. Single continuous LR axis with one
row per estimator, in two blocks -- the baseline LR (log-Wald / bootstrap /
Bayesian x2 / LOO) and the Tavtigian OddsPath (log-Wald / bootstrap / LOO).
Pathogenic (abnormal) LR red on the right, benign (normal) LR blue on the left;
the fraction of resamples in the full-data ACMG tier annotated above each marker.

Input: tables/Supplemental_LR_uncertainty_methods.tsv (built by table_lr_uncertainty.py).
Output: figures/SuppFig_LR_uncertainty_KURIS.{png,pdf}
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, palettes, plots
import matplotlib.pyplot as plt

STABILITY = config.TABLES / "Supplemental_LR_uncertainty_methods.tsv"
# (estimand, method, short row label). Marker shape distinguishes the three estimands.
ROWS = [
    ("LR (Haldane)", "log-Wald", "log-Wald"),
    ("LR (Haldane)", "Bootstrap", "Bootstrap"),
    ("LR (Haldane)", "LOO", "LOO"),
    ("OddsPath (Tavtigian)", "log-Wald", "log-Wald"),
    ("OddsPath (Tavtigian)", "Bootstrap", "Bootstrap"),
    ("OddsPath (Tavtigian)", "LOO", "LOO"),
    ("Bayesian", "Beta(0.5, 0.5)", "Beta(0.5,0.5)"),
    ("Bayesian", "Beta(1, 1)", "Beta(1,1)"),
]
GROUPS = [("LR (Haldane)", 3), ("OddsPath (Tavtigian)", 3), ("Bayesian", 2)]
MARKERS = {"LR (Haldane)": "o", "OddsPath (Tavtigian)": "D", "Bayesian": "s"}


def _clean(v):
    return float(str(v).rstrip("*"))


def _num(v):
    return pd.to_numeric(v, errors="coerce")


def main():
    config.set_rcparams()
    df = pd.read_csv(STABILITY, sep="\t")
    df = df[df["Truth set"] == "KURIS/NDD vs B/LB missense"]
    fig, ax = plt.subplots(figsize=(180 / 25.4, 150 / 25.4))

    # Axis spans the point estimates and the ACMG bands, not the extreme CIs; a CI
    # that runs past the axis is clamped and drawn with an arrowhead.
    ests = _num(df["Point estimate"].map(lambda v: str(v).rstrip("*")))
    xlim = (min(ests.min(), min(plots.ACMG_TICKS)) / 1.7,
            max(ests.max(), max(plots.ACMG_TICKS)) * 1.7)
    plots.forest_axis(ax, xlim)

    n = len(ROWS)
    for i, (estimand, method, _label) in enumerate(ROWS):
        y = n - i - 1
        shape = MARKERS[estimand]
        for lr_type, col in (("Abnormal", palettes.PATH_MARKER), ("Normal", palettes.BEN_MARKER)):
            row = df[(df["LR type"] == lr_type) & (df["Estimand"] == estimand)
                     & (df["Method"] == method)].iloc[0]
            est = _clean(row["Point estimate"])
            lo, hi = float(_num(row["2.5% (or min)"])), float(_num(row["97.5% (or max)"]))
            plots.forest_marker(ax, y, est, lo, hi, col, marker=shape, clip=xlim)
            pct = str(row.get("% in same tier", ""))
            if pct and "N/A" not in pct:
                ax.annotate(pct, (est, y), textcoords="offset points", xytext=(0, 7),
                            ha="center", va="bottom", fontsize=5.5, color="#222222", zorder=5)

    ax.set_yticks(range(n))
    ax.set_yticklabels([r[2] for r in ROWS][::-1], fontsize=8.5)
    ax.set_ylim(-0.6, n - 0.4)

    # separators + right-side labels between estimand blocks
    starts = np.cumsum([0] + [g[1] for g in GROUPS])
    for boundary in starts[1:-1]:
        ax.axhline(n - 0.5 - boundary, color="#BBBBBB", lw=0.6, zorder=1)
    for (name, k), s in zip(GROUPS, starts):
        ymid = n - 1 - (s + (k - 1) / 2)
        ax.text(1.015, ymid, name, transform=ax.get_yaxis_transform(), rotation=90,
                va="center", ha="left", fontsize=7, fontweight="bold")

    plots.forest_legends(fig, "Abnormal", "Normal", y_marker=1.0, y_tier=0.92)
    fig.subplots_adjust(left=0.15, right=0.93, top=0.76, bottom=0.16)
    plots.save(fig, "SuppFig_LR_uncertainty_KURIS")


if __name__ == "__main__":
    main()
