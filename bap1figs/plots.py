"""Shared plotting helpers: panel labels, contingency heatmaps, broken-axis
histograms, single-axis ACMG forests, and a save helper. Used by every figure so
the visual language is consistent.
"""

import numpy as np
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

from . import config, palettes


def panel_label(ax, letter, x=-0.06, y=1.04, size=10):
    """Bold lowercase panel identifier at the top-left of an axis (uniform style)."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=size,
            fontweight="bold", va="bottom", ha="right")


def save(fig, name, directory=None):
    directory = directory or config.FIGURES
    for ext in ("png", "pdf"):
        out = directory / f"{name}.{ext}"
        fig.savefig(out, bbox_inches="tight")
        print(f"INFO: wrote {out.relative_to(config.REPO)}")
    plt.close(fig)


# ── Contingency heatmap panel (rows x cols matrix, count-annotated) ──

def contingency_panel(ax, rows, cols, matrix, base_color, letter="", desc="",
                      show_y=True, show_x=True):
    cmap = LinearSegmentedColormap.from_list("cell", ["#ffffff", base_color])
    total = matrix.sum() or 1
    ax.imshow(matrix / total, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            frac = matrix[i, j] / total
            ax.text(j, i, f"{int(matrix[i, j]):,}", ha="center", va="center",
                    fontsize=7, fontweight="bold",
                    color="white" if frac > 0.35 else "black")
    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(rows)))
    ax.set_xticklabels(cols if show_x else [], rotation=30, ha="right", fontsize=6)
    ax.set_yticklabels(rows if show_y else [], fontsize=6)
    if desc:
        ax.set_title(desc, fontsize=6.5, loc="left")
    if letter:
        panel_label(ax, letter)
    ax.tick_params(length=0)


# ── Broken (split) linear-axis stacked histogram ──

def draw_stacked(ax, subset, group_col, order, color_map, edges, score_col="score"):
    bottoms = np.zeros(len(edges) - 1)
    for cat in order:
        vals = subset.loc[subset[group_col] == cat, score_col].values
        if len(vals) == 0:
            continue
        counts, _ = np.histogram(vals, bins=edges)
        ax.bar(edges[:-1], counts, width=np.diff(edges), align="edge",
               bottom=bottoms, color=color_map[cat],
               label=f"{cat} (n={len(vals):,})", linewidth=0)
        bottoms += counts
    return bottoms


def broken_hist(fig, subspec, df, group_col, order, color_map, abn_cut, nrm_cut,
                ncol=1, show_xticklabels=True, letter=None, score_col="score"):
    """Stacked histogram on a broken linear count axis: a zoomed lower panel shows
    the sparse pathogenic-range bins, a compressed upper panel carries the benign
    peak. Returns (ax_top, ax_bot)."""
    inner = subspec.subgridspec(2, 1, height_ratios=[1.0, 2.3], hspace=0.06)
    ax_top = fig.add_subplot(inner[0])
    ax_bot = fig.add_subplot(inner[1], sharex=ax_top)
    subset = df.dropna(subset=[score_col, group_col])
    edges = np.arange(-0.30, 0.021, 0.01)
    draw_stacked(ax_top, subset, group_col, order, color_map, edges, score_col)
    totals = draw_stacked(ax_bot, subset, group_col, order, color_map, edges, score_col)

    left = edges[:-1] < nrm_cut
    region_max = float(totals[left].max()) if left.any() else float(totals.max())
    brk = max(region_max * 1.25, 10)
    ax_bot.set_ylim(0, brk)
    ax_top.set_ylim(brk, float(totals.max()) * 1.08)
    ax_top.locator_params(axis="y", nbins=3)
    ax_bot.locator_params(axis="y", nbins=4)
    for ax in (ax_top, ax_bot):
        ax.axvline(abn_cut, color="black", ls="--", lw=0.8)
        ax.axvline(nrm_cut, color="black", ls=":", lw=0.8)
        ax.set_xlim(-0.30, 0.02)
        ax.spines[["top", "right"]].set_visible(False)
    ax_top.spines["bottom"].set_visible(False)
    ax_bot.spines["top"].set_visible(False)
    ax_top.tick_params(axis="x", length=0, labelbottom=False)
    dkw = dict(marker=[(-1, -0.5), (1, 0.5)], markersize=7, linestyle="none",
               color="k", mec="k", mew=1, clip_on=False)
    ax_top.plot([0], [0], transform=ax_top.transAxes, **dkw)   # break marks on the left only
    ax_bot.plot([0], [1], transform=ax_bot.transAxes, **dkw)
    ax_bot.set_ylabel("Count", y=0.7)
    if not show_xticklabels:
        ax_bot.tick_params(axis="x", labelbottom=False)
    ax_top.legend(loc="lower left", bbox_to_anchor=(0.10, 1.02), frameon=False,
                  fontsize=6, ncol=ncol, columnspacing=1.0, handlelength=1.2,
                  borderaxespad=0.0)
    if letter:
        ax_top.text(-0.065, 1.05, letter, transform=ax_top.transAxes, fontsize=10,
                    fontweight="bold", va="bottom", ha="right")
    return ax_top, ax_bot


# ── Single continuous-axis ACMG forest (benign <- 1 -> pathogenic) ──

ACMG_TICKS = [1 / 350, 1 / 18.7, 1 / 4.33, 1 / 2.08, 1, 2.08, 4.33, 18.7, 350]


def forest_axis(ax, xlim):
    """Draw the PS3 (red) / BS3 (blue) ACMG bands, the LR=1 line, and the log axis."""
    for lo, hi, color in palettes.PS3_BANDS + palettes.BS3_BANDS:
        ax.axvspan(max(lo, xlim[0]), min(hi, xlim[1]), color=color, alpha=0.40, zorder=0)
    ax.axvline(1.0, color="#555555", lw=1.0, zorder=1)
    ax.set_xscale("log")
    ax.set_xlim(xlim)
    ax.set_xticks(ACMG_TICKS)
    ax.set_xticklabels([("%g" % t if t >= 1 else "1/%g" % round(1 / t, 2)) for t in ACMG_TICKS],
                       fontsize=8, rotation=30, ha="right")
    ax.set_xlabel("Likelihood ratio", fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", length=5, width=1.1, color="#222222", direction="out")
    ax.tick_params(axis="y", length=0)


def forest_marker(ax, y, pt, lo, hi, color, marker="o", s=34, cap=0.15,
                  bar_color="black", bar_lw=1.0, alpha=1.0, clip=None):
    """One forest row. If lo/hi are None or NaN the interval is omitted (a
    point-only estimator). `marker` sets the point shape. `clip=(x0, x1)` caps a
    CI that runs past the axis and draws an arrowhead there instead of a cap."""
    import math
    has_ci = (lo is not None and hi is not None
              and not (isinstance(lo, float) and math.isnan(lo))
              and not (isinstance(hi, float) and math.isnan(hi)))
    if has_ci:
        lo_off = hi_off = False
        if clip is not None:
            if lo < clip[0]:
                lo, lo_off = clip[0], True
            if hi > clip[1]:
                hi, hi_off = clip[1], True
        ax.plot([lo, hi], [y, y], color=bar_color, lw=bar_lw, zorder=3, alpha=alpha)
        for xx, off, arrow in ((lo, lo_off, "<"), (hi, hi_off, ">")):
            if off:                                      # ran off the axis -> arrowhead
                ax.plot([xx], [y], marker=arrow, color=bar_color, ms=5, zorder=3, alpha=alpha)
            else:
                ax.plot([xx, xx], [y - cap, y + cap], color=bar_color, lw=bar_lw, zorder=3, alpha=alpha)
    ax.scatter([pt], [y], s=s, color=color, edgecolor="black", lw=0.6, marker=marker,
               zorder=4, alpha=alpha)


def forest_legends(fig, path_label, ben_label, y_marker=1.0, y_tier=0.955):
    """Marker + ACMG evidence-band legend, ordered left→right to match the x-axis:
    normal (benign, blue) on the left, abnormal (pathogenic, red) on the right, and
    the evidence bands running BS3 (benign) then PS3 (pathogenic)."""
    markers = [Line2D([0], [0], marker="o", color="w", markerfacecolor=palettes.BEN_MARKER,
                      markeredgecolor="black", markersize=7, label=ben_label),
               Line2D([0], [0], marker="o", color="w", markerfacecolor=palettes.PATH_MARKER,
                      markeredgecolor="black", markersize=7, label=path_label)]
    # BS3 reversed (Very Strong … Supporting) then PS3 (Supporting … Very Strong),
    # so the two-row legend reads in the same order as the axis (benign left of 1).
    bs3 = [mpatches.Patch(color=c, alpha=0.6, label=f"BS3 {t}")
           for (_, _, c), t in zip(palettes.BS3_BANDS, palettes.ACMG_TIERS)][::-1]
    ps3 = [mpatches.Patch(color=c, alpha=0.6, label=f"PS3 {t}")
           for (_, _, c), t in zip(palettes.PS3_BANDS, palettes.ACMG_TIERS)]
    leg1 = fig.legend(handles=markers, loc="upper center", bbox_to_anchor=(0.5, y_marker),
                      ncol=2, frameon=False, fontsize=7.5)
    fig.add_artist(leg1)
    fig.legend(handles=bs3 + ps3, loc="upper center", bbox_to_anchor=(0.5, y_tier),
               ncol=4, frameon=False, fontsize=6.3)
