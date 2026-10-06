#!/usr/bin/env python3
"""
Supplementary Figure 1: SGE GMM calibration. Three panels from the master data
table + the fitted GMM thresholds:
    a  two-component GMM densities with the abnormal/normal score cutoffs
    b  per-consequence score strip plot
    c  cumulative score distribution by consequence, with the cutoffs

Input: master data table + GMM_thresholds.json.
Output: figures/SuppFig_GMM.{png,pdf}
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, palettes, plots, tables

XMIN, XMAX = -0.30, 0.05
ABN_COLOR, NRM_COLOR = "#F08A80", "#7B7BE8"
CATEGORY_ORDER = ["Frameshift", "Canonical Splice", "Stop Gained", "Start Lost", "Stop Lost",
                  "Missense", "Splice Region", "In-frame Indel", "Intron", "Synonymous", "UTR"]
SCORE_COL, CONS_COL = "functional_score", "display_consequence"


def plot_densities(ax, thr):
    ma, sa = thr["component_abnormal"]["mean"], thr["component_abnormal"]["sd"]
    mn, sn = thr["component_normal"]["mean"], thr["component_normal"]["sd"]
    x = np.linspace(XMIN, XMAX, 3000)
    da, dn = norm.pdf(x, ma, sa), norm.pdf(x, mn, sn)
    ax.fill_between(x, da, color=ABN_COLOR, alpha=0.85, lw=0, label="abnormal", zorder=2)
    ax.fill_between(x, dn, color=NRM_COLOR, alpha=0.85, lw=0, label="normal", zorder=2)
    ymax = max(da.max(), dn.max()) * 1.08
    # black dashed = abnormal cutoff, black dotted = normal cutoff (matches panel c)
    for t, ls, ha, dx in ((thr["abnormal_cutoff"], "--", "right", -0.004),
                          (thr["normal_cutoff"], ":", "left", 0.004)):
        ax.axvline(t, color="black", ls=ls, lw=1.2, zorder=3)
        ax.text(t + dx, ymax * 0.5, "%.4f" % t, rotation=90, va="center", ha=ha,
                fontsize=6, color="black")
    ax.set_ylim(0, ymax); ax.set_xlim(XMIN, XMAX)
    ax.set_ylabel("GMM component\ndensity")
    ax.legend(loc="upper left", frameon=False, handlelength=1.0, handleheight=1.0)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def plot_strip(ax, df):
    strip = df.dropna(subset=[SCORE_COL])
    present = [c for c in CATEGORY_ORDER if (strip[CONS_COL] == c).any()]
    offsets = list(range(len(present)))[::-1]
    ax.eventplot([strip.loc[strip[CONS_COL] == c, SCORE_COL].values for c in present],
                 orientation="horizontal", lineoffsets=offsets, linelengths=0.8,
                 linewidths=0.5, colors="black")
    ax.set_yticks(offsets); ax.set_yticklabels(present)
    ax.set_ylim(-0.7, len(present) - 0.3); ax.set_xlim(XMIN, XMAX)
    ax.set_xlabel(r"$\mathit{BAP1}$ SGE score")
    ax.grid(axis="x", color="#DDDDDD", lw=0.6, zorder=0); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def plot_cumulative(ax, df, thr):
    cdf = df.dropna(subset=[SCORE_COL, CONS_COL]).sort_values(SCORE_COL)
    gb = cdf.groupby(CONS_COL)[SCORE_COL]
    cdf = cdf.assign(cum_pct=(gb.cumcount() + 1) / gb.transform("count") * 100.0)
    present = [c for c in CATEGORY_ORDER if c in cdf[CONS_COL].unique()]
    for c in present:
        sub = cdf[cdf[CONS_COL] == c]
        ax.plot(sub[SCORE_COL], sub["cum_pct"], color=palettes.DISPLAY_CONSEQUENCE_COLORS[c],
                lw=1.6, label=c, zorder=3)
    ax.axvline(thr["abnormal_cutoff"], color="black", ls="--", lw=1.0, zorder=4)
    ax.axvline(thr["normal_cutoff"], color="black", ls=":", lw=1.0, zorder=4)
    ax.text(thr["abnormal_cutoff"], 101.5, "abnormal", rotation=90, va="bottom", ha="right",
            fontsize=5.5, color="black")
    ax.text(thr["normal_cutoff"], 101.5, "normal", rotation=90, va="bottom", ha="left",
            fontsize=5.5, color="black")
    ax.set_xlim(XMIN, XMAX); ax.set_ylim(0, 101)
    ax.set_xlabel(r"$\mathit{BAP1}$ SGE score"); ax.set_ylabel("Cumulative % of variants")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def main():
    config.set_rcparams()
    df = tables.load_master()
    df[SCORE_COL] = pd.to_numeric(df[SCORE_COL], errors="coerce")
    with open(config.GMM_THRESHOLDS) as fh:
        thr = json.load(fh)

    fig = plt.figure(figsize=(190 / 25.4, 140 / 25.4))
    gs = GridSpec(2, 2, height_ratios=[1, 2.0], width_ratios=[1, 1],
                  hspace=0.30, wspace=0.30, figure=fig)
    ax_a = fig.add_subplot(gs[0, 0]); ax_b = fig.add_subplot(gs[1, 0], sharex=ax_a)
    ax_leg = fig.add_subplot(gs[0, 1]); ax_c = fig.add_subplot(gs[1, 1])
    plot_densities(ax_a, thr); plot_strip(ax_b, df); plot_cumulative(ax_c, df, thr)

    # variant-type legend in the freed top-right cell (never overlaps the curves)
    ax_leg.axis("off")
    handles, labels = ax_c.get_legend_handles_labels()
    ax_leg.legend(handles, labels, loc="center", ncol=2, frameon=False, fontsize=6.5,
                  title="Variant type", title_fontsize=7.5, handlelength=1.6,
                  columnspacing=1.2, labelspacing=0.5)

    plt.setp(ax_a.get_xticklabels(), visible=False)
    for ax, lab in ((ax_a, "a"), (ax_b, "b"), (ax_c, "c")):
        ax.set_title(lab, loc="left", fontweight="bold", fontsize=10)
    plots.save(fig, "SuppFig_GMM")


if __name__ == "__main__":
    main()
