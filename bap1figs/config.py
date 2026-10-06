"""Shared configuration: paths and matplotlib defaults for all figure/table scripts."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parent.parent
RAW = REPO / "data" / "raw"
SOURCE = REPO / "data" / "source"
FIGURE_INPUTS = REPO / "data" / "figure_inputs"
FIGURES = REPO / "figures"
TABLES = REPO / "tables"
SUPPLEMENT = REPO / "supplement"        # final numbered, captioned deliverables

# The two source tables (+ calibration) every figure input is rebuilt from.
MASTER_TABLE = SOURCE / "BAP1_master_table.tsv"
CLINICAL_TABLE = SOURCE / "BAP1_KURIS_NDD_variants.csv"
GMM_THRESHOLDS = SOURCE / "GMM_thresholds.json"

# UniProt (Q92560) feature table — source of the protein-domain annotation.
FEATURES = RAW / "BAP1_features.json"

for _d in (FIGURE_INPUTS, FIGURES, TABLES, SOURCE, SUPPLEMENT):
    _d.mkdir(parents=True, exist_ok=True)


def set_rcparams():
    """Publication defaults shared by every figure (Arial, 8 pt, vector fonts, 600 dpi)."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.linewidth": 0.8, "axes.labelsize": 8,
        "axes.titlesize": 10, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.dpi": 150, "savefig.dpi": 600,
    })
