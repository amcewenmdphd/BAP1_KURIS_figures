"""Shared table loaders and the panel3_category classifier.

Every figure input is rebuilt from the master data table + the clinical data
table through these helpers, so the two source tables are the single origin.
"""

import json

import pandas as pd

from . import config

# Proband shown as its own lollipop category: p.Asn229Lys (c.687C>A).
N229K_PROTEIN = "p.Asn229Lys"

# BAP1 functional domains annotated for Supplementary Table 1. Each entry maps a
# UniProt (Q92560) feature description to a short label; residues are read from
# data/raw/BAP1_features.json so the annotation is reproducible from the inputs.
# Listed most-specific first: a residue is labelled by the first interval it falls
# in (e.g. ULD, 670-698, wins over the broader BRCA1-interaction region it sits in).
DOMAIN_FEATURES = [
    ("HBM-like motif", "HBM"),
    ("UCH catalytic", "UCH"),
    ("ULD", "ULD"),
    ("Interaction with BRCA1", "BRCA1-binding"),
]
# Missense residue inside the protein but in none of the above.
DOMAIN_OTHER = "Other"


def load_domain_intervals(path=None):
    """[(start, end, label), ...] for the annotated domains, in priority order,
    read from the UniProt feature JSON (data/raw/BAP1_features.json)."""
    feats = json.load(open(path or config.FEATURES))["features"]
    by_desc = {}
    for f in feats:
        desc = f.get("description")
        if desc and desc not in by_desc:
            loc = f["location"]
            by_desc[desc] = (int(loc["start"]["value"]), int(loc["end"]["value"]))
    intervals = []
    for desc, label in DOMAIN_FEATURES:
        if desc in by_desc:
            s, e = by_desc[desc]
            intervals.append((s, e, label))
    return intervals


def annotate_domain(positions, path=None):
    """Map amino-acid positions to domain labels (first matching interval, else
    'Other'; <NA> position -> <NA>). `positions` is any int-like Series/iterable."""
    intervals = load_domain_intervals(path)

    def label(p):
        if pd.isna(p):
            return pd.NA
        p = int(p)
        for s, e, lab in intervals:
            if s <= p <= e:
                return lab
        return DOMAIN_OTHER

    return pd.Series(list(positions)).map(label)


def strip_prefix(series):
    """ENST...:c.x -> c.x ; ENSP...:p.x -> p.x ; '-'/blank -> <NA>."""
    s = series.astype("string").str.replace(r"^[^:]+:", "", regex=True)
    return s.where(~s.isin(["-", "", "nan"]), other=pd.NA)


def load_master(path=None):
    """Master data table with standardized derived columns used across figures."""
    m = pd.read_csv(path or config.MASTER_TABLE, sep="\t", low_memory=False)
    m["MAVE_Score"] = pd.to_numeric(m["functional_score"], errors="coerce")
    m["HGVS_p"] = strip_prefix(m["HGVSp"])
    m["HGVS_c"] = strip_prefix(m["HGVSc"])
    m["AA_Position"] = pd.to_numeric(m["protein_position"], errors="coerce").astype("Int64")
    return m


def load_clinical(path=None):
    """The curated KURIS/NDD clinical data table (18 variants incl. the proband)."""
    return pd.read_csv(path or config.CLINICAL_TABLE)


def add_panel3_category(m, clinical):
    """Classify missense variants into KURIS/NDD, TPDS P/LP, B/LB, N229K (proband),
    from the master ClinVar-2026 class and the clinical data table. Returns m with
    a new 'panel3_category' column (None for variants not in any category)."""
    clin_hgvsc = set(clinical["HGVS_c"].astype(str))
    n229k = set(clinical.loc[clinical["Protein_Change"] == N229K_PROTEIN, "HGVS_c"].astype(str))
    kuris = clin_hgvsc - n229k

    def classify(row):
        if row["variant_type"] != "Missense":
            return None
        if row["AA_Position"] is not pd.NA and row["AA_Position"] == 1:
            return None                                    # exclude M1 / start-lost
        hc = str(row["HGVS_c"])
        if hc in n229k:
            return "N229K (proband)"
        if hc in kuris:
            return "KURIS/NDD"
        if row["clinvar_2026_slim"] == "P/LP":
            return "TPDS P/LP"
        if row["clinvar_2026_slim"] == "B/LB":
            return "B/LB"
        return None

    m = m.copy()
    m["panel3_category"] = m.apply(classify, axis=1)
    return m


def kuris_controls(m):
    """KURIS control sets for the LR analyses, derived from panel3_category:
    path = KURIS/NDD missense, benign = B/LB missense. Requires add_panel3_category."""
    if "panel3_category" not in m.columns:
        raise ValueError("call add_panel3_category(m, clinical) first")
    return (m[m["panel3_category"] == "KURIS/NDD"].copy(),
            m[m["panel3_category"] == "B/LB"].copy())
