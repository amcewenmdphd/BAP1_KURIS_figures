"""All color schemes and ACMG evidence-band palettes, shared across every figure.

One source of truth so the consequence colors on Figure 2 match the GMM-vs-Waters
supplement, the ClinVar colors match everywhere, etc.
"""

# ── Molecular consequence (Figure 2, GMM-vs-Waters histogram c) ──
VARIANT_TYPE_COLORS = {
    "Nonsense": "#D84315", "Frameshift": "#EF6C00", "Splice": "#F9A825",
    "Start lost": "#9E9D24", "Stop lost": "#827717", "Inframe indel": "#00695C",
    "Missense": "#00838F", "Synonymous": "#8D6E63", "Intronic": "#D7CCC8",
    "UTR": "#3E2723",
}
VARIANT_TYPE_ORDER = ["Nonsense", "Frameshift", "Splice", "Start lost", "Stop lost",
                      "Inframe indel", "Missense", "Synonymous", "Intronic", "UTR"]
TRUNCATING_TYPES = ["Nonsense", "Frameshift", "Splice"]

# Finer "display_consequence" categories (GMM calibration figure) mapped onto the
# same molecular-consequence colors, so the palette matches the rest of the repo.
DISPLAY_CONSEQUENCE_COLORS = {
    "Frameshift": VARIANT_TYPE_COLORS["Frameshift"],
    "Canonical Splice": VARIANT_TYPE_COLORS["Splice"],
    "Stop Gained": VARIANT_TYPE_COLORS["Nonsense"],
    "Start Lost": VARIANT_TYPE_COLORS["Start lost"],
    "Stop Lost": VARIANT_TYPE_COLORS["Stop lost"],
    "Missense": VARIANT_TYPE_COLORS["Missense"],
    "Splice Region": "#FBC02D",                    # amber, distinct from Canonical Splice
    "In-frame Indel": VARIANT_TYPE_COLORS["Inframe indel"],
    "Intron": VARIANT_TYPE_COLORS["Intronic"],
    "Synonymous": VARIANT_TYPE_COLORS["Synonymous"],
    "UTR": VARIANT_TYPE_COLORS["UTR"],
}

# ── ClinVar classification ──
CLINVAR_2026_COLORS = {"P/LP": "#E64B35", "B/LB": "#3C5488",
                       "VUS": "#B8B8B8", "Conflicting": "#7B2D8E"}

# ── Lollipop / KURIS categories (Figure 1a, Figure 2c) ──
PANEL3_COLORS = {"KURIS/NDD": "#E8755A", "TPDS P/LP": "#722F37",
                 "N229K (proband)": "#FFC107", "B/LB": "#3C5488"}

# ── Curated condition (Figure 1c) ──
COND_COLORS = {"Cancer predisposition": "#722F37", "Kury-Isidor syndrome": "#E8755A",
               "Not provided": "#B8B8B8"}

# ── Score-threshold class (Figure 2 contingency shading) ──
SCORE_CLASS_COLORS = {"Abnormal": "#F4A460", "Indeterminate": "#B0B0B0", "Normal": "#80CBC4"}

# ── Waters functional classes (GMM-vs-Waters) ──
WATERS_COLORS = {"depleted": "#D84315", "unchanged": "#9AA4AD", "enriched": "#00695C"}

# ── ACMG evidence: pathogenic (PS3) reds, benign (BS3) blues ──
# (band lo, band hi, fill color); markers are the dark endpoints.
PS3_BANDS = [(2.08, 4.33, "#FEE0D2"), (4.33, 18.7, "#FCAE91"),
             (18.7, 350.0, "#FB6A4A"), (350.0, 1e9, "#CB181D")]
BS3_BANDS = [(1 / 4.33, 1 / 2.08, "#DEEBF7"), (1 / 18.7, 1 / 4.33, "#9ECAE1"),
             (1 / 350.0, 1 / 18.7, "#4292C6"), (1e-9, 1 / 350.0, "#08519C")]
PATH_MARKER = "#99000D"      # pathogenic-class (abnormal) LR marker
BEN_MARKER = "#08306B"       # benign-class (normal) LR marker
ACMG_TIERS = ["Supporting", "Moderate", "Strong", "Very Strong"]

# Compact evidence-tier fill colors for the Figure 2 contingency table cells.
EVIDENCE_CELL_COLORS = {
    "PS3 V.Strong": "#E64B35", "PS3 Strong": "#F08070", "PS3 Moderate": "#F5AFA5",
    "PS3 Supptg.": "#FDDEDE", "BS3 V.Strong": "#3C5488", "BS3 Strong": "#6B82B0",
    "BS3 Moderate": "#9AB0D8", "BS3 Supptg.": "#DEE8F5", "Indet.": "#E0E0E0",
}

# ── BAP1 UniProt domains (Figure 1 lollipop track) ──
DOMAIN_MAP = {
    "UCH catalytic": ("#90CAF9", "UCH"), "ULD": ("#CE93D8", "ULD"),
    "HBM-like motif": ("#FFE082", "HBM"), "Interaction with BRCA1": ("#FFAB91", "BRCA1 binding"),
}
