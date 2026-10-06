#!/usr/bin/env python3
"""
Supplementary Table 2 — variants shown in Figure 1.

One row per variant in the Figure 1a lollipop (KURIS/NDD, the N229K proband,
TPDS-associated ClinVar P/LP missense, and the ClinVar B/LB missense controls):
HGVS g./c./p. (with reference sequences — see Supplementary Table 1, named
again in the footnote here since Gene Symbol/HGNC ID are constant for this
single-gene table), SGE functional score/class, and ClinVar classification /
accession / review status. The "SGE evidence (this study)" column carries the
per-class ACMG code (PS3 Strong / BS3 Moderate) from the KURIS/NDD-vs-B/LB
comparison in Supplementary Table 4 (GMM cutoffs row); it is left blank for the
TPDS P/LP variants, which are not part of that comparison (see footnote).

Uses the same categorization used to build Figure 1a (bap1figs.tables.
add_panel3_category), so this table's 59 rows are exactly Figure 1's lollipop
variants. Genomic HGVS is built from the master table's GRCh38 pos/ref/alt
against NC_000003.12; HGVS c./p. substitute the RefSeq NM_004656.4 /
NP_004647.1 accessions for the Ensembl ENST/ENSP ones the master table stores
(CCDS-matched, so cDNA/protein numbering is unchanged — see
Supplementary Table 1).

Output: tables/Figure1_Variants.{xlsx,tsv,pdf} + the numbered, captioned
supplement/Supplementary_Table_2.pdf.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, tables, tablepdf, supplement

GRCH38_CHROM = "NC_000003.12"
NM_ACCESSION = "NM_004656.4"
NP_ACCESSION = "NP_004647.1"

# Per-class ACMG evidence code from the GMM-cutoffs row of the KURIS/NDD-vs-B/LB
# comparison (Supplementary Table 4 / Master_Contingency_Tables.tsv): Abnormal ->
# PS3 Strong, Normal -> BS3 Moderate. Applies only to the categories that make up
# that comparison (KURIS/NDD, the N229K proband, and the B/LB controls).
SGE_EVIDENCE_BY_CLASS = {"Abnormal": "PS3 Strong", "Normal": "BS3 Moderate"}
SGE_EVIDENCE_CATEGORIES = {"KURIS/NDD", "N229K (proband)", "B/LB"}

CATEGORY_LABEL = {
    "KURIS/NDD": "KURIS/NDD",
    "N229K (proband)": "KURIS/NDD (proband)",
    "TPDS P/LP": "Cancer predisposition (ClinVar P/LP)",
    "B/LB": "ClinVar B/LB (control)",
}
CATEGORY_ORDER = ["N229K (proband)", "KURIS/NDD", "TPDS P/LP", "B/LB"]

COLUMNS = [
    "Gene Symbol", "HGNC Gene ID", "Variant Category", "HGVS g. (GRCh38)",
    "HGVS c.", "HGVS p.", "Domain", "SGE Functional Score",
    "SGE Functional Classification", "SGE Evidence (this study)",
    "ClinVar Classification (Aggregate Germline)", "ClinVar Accession",
    "ClinVar Review Status (stars)",
]
# Gene Symbol/HGNC Gene ID are constant (single-gene table) and kept in the
# xlsx/tsv for completeness, but dropped from the printed PDF grid — where they
# are named once in the footnote instead — to leave width for the columns that
# actually vary per variant.
PDF_DROP_COLUMNS = ["Gene Symbol", "HGNC Gene ID"]


def _hgvs_g(row):
    return f"{GRCH38_CHROM}:g.{int(row['pos'])}{row['ref']}>{row['alt']}"


def _hgvs_c(hgvs_c_short):
    return f"{NM_ACCESSION}:{hgvs_c_short}"


def _hgvs_p(hgvs_p_short):
    return f"{NP_ACCESSION}:{hgvs_p_short}"


def _vcv(variation_id):
    if pd.isna(variation_id):
        return pd.NA
    return f"VCV{int(variation_id):09d}"


def build():
    m = tables.load_master()
    m = tables.add_panel3_category(m, tables.load_clinical())
    cats = (m.dropna(subset=["panel3_category"])
            .drop_duplicates(subset=["AA_Position", "ref_aa", "alt_aa", "panel3_category"])
            .copy())

    cats["Variant Category"] = cats["panel3_category"].map(CATEGORY_LABEL)
    cats["_order"] = cats["panel3_category"].map(CATEGORY_ORDER.index)

    sge_evidence = cats.apply(
        lambda r: SGE_EVIDENCE_BY_CLASS.get(r["score_threshold_class"])
        if r["panel3_category"] in SGE_EVIDENCE_CATEGORIES else pd.NA, axis=1)

    out = pd.DataFrame({
        "Gene Symbol": "BAP1",
        "HGNC Gene ID": "HGNC:950",
        "Variant Category": cats["Variant Category"].values,
        "HGVS g. (GRCh38)": cats.apply(_hgvs_g, axis=1).values,
        "HGVS c.": cats["HGVS_c"].map(_hgvs_c).values,
        "HGVS p.": cats["HGVS_p"].map(_hgvs_p).values,
        "Domain": tables.annotate_domain(cats["AA_Position"].values).values,
        "SGE Functional Score": cats["MAVE_Score"].values,
        "SGE Functional Classification": cats["score_threshold_class"].values,
        "SGE Evidence (this study)": sge_evidence.values,
        "ClinVar Classification (Aggregate Germline)": cats["clinvar_2026_classification"].values,
        "ClinVar Accession": [_vcv(v) for v in cats["clinvar_2026_variation_id"].values],
        "ClinVar Review Status (stars)": cats["clinvar_2026_stars"].values,
        "_order": cats["_order"].values,
        "_score": cats["MAVE_Score"].values,
    })
    out = out.sort_values(["_order", "_score"]).drop(columns=["_order", "_score"])
    return out.reset_index(drop=True)[COLUMNS]


def _pdf_view(df):
    d = df.drop(columns=PDF_DROP_COLUMNS).copy()
    d["SGE Functional Score"] = d["SGE Functional Score"].map(
        lambda v: "" if pd.isna(v) else f"{v:.3f}")
    d["ClinVar Review Status (stars)"] = d["ClinVar Review Status (stars)"].map(
        lambda v: "" if pd.isna(v) else str(int(v)))
    for c in ["SGE Evidence (this study)", "ClinVar Classification (Aggregate Germline)",
              "ClinVar Accession"]:
        d[c] = d[c].fillna("")
    return d.rename(columns={
        "HGVS g. (GRCh38)": "HGVS g.", "SGE Functional Score": "SGE score",
        "SGE Functional Classification": "SGE class",
        "SGE Evidence (this study)": "SGE evidence",
        "ClinVar Classification (Aggregate Germline)": "ClinVar (agg. germline)",
        "ClinVar Accession": "ClinVar VCV",
        "ClinVar Review Status (stars)": "Stars"})


def main():
    df = build()
    xlsx = config.TABLES / "Figure1_Variants.xlsx"
    tsv = config.TABLES / "Figure1_Variants.tsv"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        df.to_excel(w, sheet_name="Figure1_Variants", index=False)
    df.to_csv(tsv, sep="\t", index=False)

    kw = dict(
        aligns={"SGE score": "right", "Stars": "center", "SGE class": "center",
                "SGE evidence": "center", "Domain": "center"},
        group_cols=["HGVS c."],           # unique per row -> a rule between every variant
    )
    tablepdf.render_table_pdf(_pdf_view(df), "Figure1_Variants", **kw)
    supplement.render_supp_table("Figure1_Variants", _pdf_view(df), **kw)
    n = df["Variant Category"].value_counts()
    print(f"INFO: wrote tables/{xlsx.name} + .tsv + .pdf ({len(df)} variants: "
          f"{dict(n)})")


if __name__ == "__main__":
    main()
