#!/usr/bin/env python3
"""
Figure 1: BAP1 variant landscape.
    a  protein lollipop (KURIS/NDD, TPDS P/LP, B/LB, N229K missense on domains)
    b  ClinVar variants by molecular consequence
    c  P/LP variants by curated condition (cancer vs Kury-Isidor)

Step 1 rebuilds this figure's input file from the master + clinical data tables
(panel a) and the ClinVar-VCV intermediates (panels b, c); step 2 renders it.

Input file: data/figure_inputs/figure1_variants.xlsx
Output: figures/Figure1.{png,pdf}
"""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, palettes, stats, tables

BAP1_LENGTH = 729
INPUT = config.FIGURE_INPUTS / "figure1_variants.xlsx"
VCV_CLASSIFIED = config.SOURCE / "BAP1_ClinVar_VCV_classified.tsv"
PLP_CONDITION = config.SOURCE / "BAP1_ClinVar_PLP_by_condition.tsv"
FEATURES = config.RAW / "BAP1_features.json"
ABN, NRM = stats.load_cutoffs()


# ── Step 1: rebuild input file ──

def build_input():
    m = tables.load_master()
    m = tables.add_panel3_category(m, tables.load_clinical())
    cats = (m.dropna(subset=["panel3_category"])
            .drop_duplicates(subset=["AA_Position", "ref_aa", "alt_aa", "panel3_category"]))
    fig1a = pd.DataFrame({
        "AA_Change": cats["ref_aa"].astype(str) + cats["AA_Position"].astype(str) + cats["alt_aa"].astype(str),
        "AA_Position": cats["AA_Position"], "AA_Ref": cats["ref_aa"], "AA_Alt": cats["alt_aa"],
        "Category": cats["panel3_category"], "MAVE_Score": cats["MAVE_Score"],
        "Is_N229K": cats["panel3_category"] == "N229K (proband)"})
    fig1a["Score_Class"] = fig1a["MAVE_Score"].apply(
        lambda s: "Abnormal" if s <= ABN else ("Normal" if s >= NRM else "Indeterminate"))
    fig1a = fig1a.sort_values(["Category", "AA_Position"])

    fig1b = (pd.read_csv(VCV_CLASSIFIED, sep="\t")
             .rename(columns={"ProteinChange": "Protein", "MolecularConsequence": "Consequence_Detail",
                              "simple_consequence": "Consequence", "cv2026_class": "ClinVar_Class",
                              "AggGermlineClass": "Germline_Class"})
             [["VCV", "Protein", "Consequence_Detail", "Consequence", "ClinVar_Class", "StarRating", "Germline_Class"]]
             .sort_values(["ClinVar_Class", "Consequence"]))
    fig1c = (pd.read_csv(PLP_CONDITION, sep="\t")
             .rename(columns={"ProteinChange": "Protein", "MolecularConsequence": "Consequence_Detail",
                              "simple_consequence": "Consequence", "curated_condition": "Condition"})
             [["VCV", "Protein", "Consequence_Detail", "Consequence", "Condition", "StarRating"]]
             .sort_values(["Condition", "Consequence"]))
    with pd.ExcelWriter(INPUT, engine="openpyxl") as w:
        fig1a.to_excel(w, sheet_name="Fig1a_Lollipop", index=False)
        fig1b.to_excel(w, sheet_name="Fig1b_ClinVar_All", index=False)
        fig1c.to_excel(w, sheet_name="Fig1c_PLP_Condition", index=False)
    print(f"INFO: wrote {INPUT.relative_to(config.REPO)}")
    return {"Fig1a_Lollipop": fig1a, "Fig1b_ClinVar_All": fig1b, "Fig1c_PLP_Condition": fig1c}


def load_domains():
    with open(FEATURES) as f:
        uniprot = json.load(f)
    domains = []
    for feat in uniprot.get("features", []):
        desc = feat.get("description", "")
        if desc in palettes.DOMAIN_MAP:
            color, label = palettes.DOMAIN_MAP[desc]
            domains.append({"name": label, "start": feat["location"]["start"]["value"],
                            "end": feat["location"]["end"]["value"], "color": color})
    domains.sort(key=lambda d: d["start"])
    return domains


# ── Step 2: render ──

def plot_lollipop(ax, df_lolli, domains):
    by, bh = 0, 0.08
    ruler_y = by - 0.09
    ax.add_patch(FancyBboxPatch((0, by - bh / 2), BAP1_LENGTH + 1, bh, boxstyle="round,pad=0.002",
                                facecolor="#D9D9D9", edgecolor="black", linewidth=0.6, zorder=2))
    for d in domains:
        ax.add_patch(FancyBboxPatch((d["start"], by - bh / 2), d["end"] - d["start"], bh,
                                    boxstyle="round,pad=0.002", facecolor=d["color"], edgecolor="black",
                                    linewidth=0.6, zorder=3))
    tick_y = ruler_y - 0.02
    for pos in [1, 100, 200, 300, 400, 500, 600, 700, BAP1_LENGTH]:
        ax.plot([pos, pos], [tick_y - 0.015, tick_y + 0.015], color="black", linewidth=0.5, zorder=2)
        ax.text(pos, tick_y - 0.025, str(pos), ha="center", va="top", fontsize=5, color="#333")
    ax.plot([1, BAP1_LENGTH], [tick_y, tick_y], color="black", linewidth=0.5, zorder=2)

    df_s = df_lolli.sort_values("AA_Position").reset_index(drop=True)
    n = len(df_s)
    pos_counts = df_s["AA_Position"].value_counts()
    idx, x_pos = {}, []
    for i in range(n):
        p = df_s.loc[i, "AA_Position"]; na = pos_counts[p]; cur = idx.get(p, 0); idx[p] = cur + 1
        if na == 1:
            x_pos.append(p)
        else:
            spread = min(18, na * 6)
            x_pos.append(p - spread / 2 + cur * spread / max(na - 1, 1))
    h_base, h_step, char_w, pad_x = 0.12, 0.055, 5.5, 5
    tiers, assigned = {}, []
    for i in range(n):
        w = len(df_s.loc[i, "AA_Change"]) * char_w + pad_x
        xl, xr = x_pos[i] - 2, x_pos[i] + w
        for t in range(50):
            tiers.setdefault(t, [])
            if not any(xl < r and xr > l for l, r in tiers[t]):
                tiers[t].append((xl, xr)); assigned.append(t); break
        else:
            assigned.append(50)
    heights = [h_base + t * h_step for t in assigned]
    for i in range(n):
        row = df_s.loc[i]; x = x_pos[i]; cat = row["Category"]; change = row["AA_Change"]
        color = palettes.PANEL3_COLORS.get(cat, "#888888"); h = heights[i]
        line_color = "black" if cat == "N229K (proband)" else color
        ax.plot([x, x], [by + bh / 2, h], color=line_color, linewidth=0.5, zorder=5)
        if bool(row.get("Is_N229K", False)):
            ax.plot(x, h, "*", markersize=8, color=color, markeredgecolor="black", markeredgewidth=1.0, zorder=7)
        elif cat == "TPDS P/LP":
            ax.plot(x, h, "s", markersize=2.5, color=color, markeredgecolor="black", markeredgewidth=0.3, zorder=6)
        else:
            ew = 0.8 if cat == "N229K (proband)" else 0.3
            ax.plot(x, h, "o", markersize=2.5, color=color, markeredgecolor="black", markeredgewidth=ew, zorder=6)
        ax.text(x + 8, h, change, ha="left", va="center", fontsize=4.5, fontweight="bold", color=line_color,
                zorder=8, bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="black",
                                    linewidth=0.3, alpha=0.95))
    ax.set_xlim(-15, BAP1_LENGTH + 55)
    ax.set_ylim(tick_y - 0.04, (max(heights) if heights else 1.0) + h_step + 0.03)
    ax.axis("off")
    handles = []
    for cat, marker in [("N229K (proband)", "*"), ("KURIS/NDD", "o"), ("TPDS P/LP", "s"), ("B/LB", "o")]:
        if cat in df_lolli["Category"].values:
            ms = 6 if cat == "N229K (proband)" else 3.5
            ew = 1.0 if cat == "N229K (proband)" else 0.3
            handles.append(mlines.Line2D([], [], marker=marker, color=palettes.PANEL3_COLORS[cat], linestyle="-",
                                         linewidth=0.8, markersize=ms, markeredgecolor="black",
                                         markeredgewidth=ew, label=cat))
    for d in domains:
        handles.append(mpatches.Patch(facecolor=d["color"], edgecolor="black", linewidth=0.4, label=d["name"]))
    handles.append(mpatches.Patch(facecolor="#D9D9D9", edgecolor="black", linewidth=0.4, label="Other"))
    ax.legend(handles=handles, fontsize=5, loc="lower center", bbox_to_anchor=(0.5, 1.0), framealpha=0.95,
              edgecolor="#CCC", ncol=len(handles), handlelength=1.2, handletextpad=0.3, labelspacing=0.2,
              columnspacing=0.6)


def plot_clinvar_bars(ax, df_cv, cons_order, color_map, group_col):
    groups = df_cv[group_col].dropna().unique()
    x = np.arange(len(cons_order)); bottom = np.zeros(len(cons_order))
    for grp in color_map:
        if grp not in groups:
            continue
        sub = df_cv[df_cv[group_col] == grp]
        counts = np.array([(sub["Consequence"] == c).sum() for c in cons_order])
        if counts.sum() > 0:
            ax.bar(x, counts, bottom=bottom, width=0.65, color=color_map[grp], edgecolor="white",
                   linewidth=0.3, label=f"{grp} (n={counts.sum()})")
            bottom += counts
    ax.set_xticks(x); ax.set_xticklabels(cons_order, rotation=45, ha="right", fontsize=5)
    ax.set_ylabel("Count", fontsize=6)
    ax.legend(fontsize=4.5, loc="upper right", framealpha=0.9, edgecolor="none")
    ax.tick_params(axis="both", labelsize=5.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def render(sheets, domains):
    config.set_rcparams()
    fig = plt.figure(figsize=(7.09, 4.72))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.0, 1.0], height_ratios=[0.55, 0.45], hspace=0.45, wspace=0.35)
    ax_a = fig.add_subplot(gs[0, :])
    plot_lollipop(ax_a, sheets["Fig1a_Lollipop"], domains)
    ax_a.text(-0.02, 1.0, "a", fontsize=10, fontweight="bold", transform=ax_a.transAxes, va="bottom", ha="right")
    cons_order = ["Frameshift", "Nonsense", "Splice", "Missense", "Start lost", "Inframe indel",
                  "Stop lost", "Synonymous", "Intronic", "UTR"]
    present = (set(sheets["Fig1b_ClinVar_All"]["Consequence"].unique())
               | set(sheets["Fig1c_PLP_Condition"]["Consequence"].unique()))
    cons_order = [c for c in cons_order if c in present]
    ax_b = fig.add_subplot(gs[1, 0])
    plot_clinvar_bars(ax_b, sheets["Fig1b_ClinVar_All"], cons_order, palettes.CLINVAR_2026_COLORS, "ClinVar_Class")
    ax_b.text(-0.12, 1.05, "b", fontsize=10, fontweight="bold", transform=ax_b.transAxes, va="bottom", ha="right")
    ax_c = fig.add_subplot(gs[1, 1])
    plot_clinvar_bars(ax_c, sheets["Fig1c_PLP_Condition"], cons_order, palettes.COND_COLORS, "Condition")
    ax_c.text(-0.12, 1.05, "c", fontsize=10, fontweight="bold", transform=ax_c.transAxes, va="bottom", ha="right")
    for ext in ("png", "pdf"):
        fig.savefig(config.FIGURES / f"Figure1.{ext}", dpi=600, bbox_inches="tight")
    plt.close(fig)
    print("INFO: wrote figures/Figure1.{png,pdf}")


def main():
    render(build_input(), load_domains())


if __name__ == "__main__":
    main()
