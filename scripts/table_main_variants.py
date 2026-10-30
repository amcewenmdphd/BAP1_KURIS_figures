#!/usr/bin/env python3
"""
Main-article Table 1 -- every variant shown or named in the main article.

The journal requires the main article to stand alone: variants described in the
main text or figures must be fully described there, not only in the supplement.
This minimised table gives, for each of the 58 missense variants in the Figure 1a
lollipop (the UDN proband N229K, the 17 KURIS/NDD pathogenic controls, the
TPDS-associated ClinVar P/LP missense variants, and the ClinVar B/LB missense
controls): gene symbol, HGNC gene ID, HGVS g./c./p. (full reference-sequence
accessions) and the single-letter amino-acid notation used in the figures.

Rows and HGVS strings come from the same builder as Supplementary Table 2
(table_supp_figure1_variants.build), so the two tables cannot disagree; scores,
ClinVar details, and curated clinical evidence stay in Supplementary Tables 1-2.
The single-letter codes are derived from the HGVS p. and asserted equal to the
ref_aa/position/alt_aa labels Figure 1a draws.

Output: tables/Table1_Variants.{xlsx,tsv,pdf}
"""

import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from bap1figs import config, tables, tablepdf
import table_supp_figure1_variants as supp2

AA3TO1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V",
}
CATEGORY_ORDER = ["KURIS/NDD (proband)", "KURIS/NDD",
                  "Cancer predisposition (ClinVar P/LP)", "ClinVar B/LB (control)"]
GENE = "BAP1; HGNC:950"      # "BAP1" is italicized when rendered (rich_cols)

TITLE = "Table 1. BAP1 missense variants described in this study."
CAPTION = (
    "All 58 missense variants shown in Figure 1a: the UDN proband variant "
    "(p.Asn229Lys), the 17 KURIS/NDD pathogenic controls, the ClinVar "
    "pathogenic/likely-pathogenic missense variants associated with BAP1 tumor "
    "predisposition syndrome (TPDS), and the 28 ClinVar benign/likely-benign "
    "missense controls. The single-letter notation is the variant label used in "
    "Figures 1 and 2. SGE functional scores, ClinVar classifications, and the "
    "clinical evidence for each variant are in Supplementary Tables 1 and 2.")
NOTE = [
    "HGVS descriptions use the reference sequences NC_000003.12 (GRCh38), "
    "NM_004656.4 (RefSeq MANE Select transcript), and NP_004647.1 (protein) "
    "for BAP1 (HGNC:950).",
]


def _one_letter(hgvs_p):
    m = re.search(r"p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})$", hgvs_p)
    if not m:
        raise ValueError(f"unparseable HGVS p.: {hgvs_p}")
    ref, pos, alt = m.groups()
    return f"{AA3TO1[ref]}{pos}{AA3TO1[alt]}", int(pos)


def build():
    s2df = supp2.build()
    coded = s2df["HGVS p."].map(_one_letter)
    out = pd.DataFrame({
        "Gene Symbol": s2df["Gene Symbol"],
        "HGNC Gene ID": s2df["HGNC Gene ID"],
        "Variant Category": s2df["Variant Category"],
        "HGVS g. (GRCh38)": s2df["HGVS g. (GRCh38)"],
        "HGVS c.": s2df["HGVS c."],
        "HGVS p.": s2df["HGVS p."],
        "Single-letter code": [c for c, _ in coded],
        "_pos": [p for _, p in coded],
    })
    out["_order"] = out["Variant Category"].map(CATEGORY_ORDER.index)
    out = out.sort_values(["_order", "_pos"]).drop(columns=["_order", "_pos"])
    return out.reset_index(drop=True)


def _check_against_figure(df):
    """The codes must be exactly the labels Figure 1a draws (ref_aa + position + alt_aa)."""
    m = tables.load_master()
    m = tables.add_panel3_category(m, tables.load_clinical())
    cats = (m.dropna(subset=["panel3_category"])
            .drop_duplicates(subset=["AA_Position", "ref_aa", "alt_aa", "panel3_category"]))
    fig_labels = set(cats["ref_aa"].astype(str) + cats["AA_Position"].astype(str)
                     + cats["alt_aa"].astype(str))
    assert set(df["Single-letter code"]) == fig_labels, "codes differ from Figure 1a labels"
    assert df["Single-letter code"].is_unique and df["HGVS c."].is_unique


def _pdf_view(df):
    d = df.copy()
    d.insert(0, "Gene; HGNC ID", GENE)
    d = d.drop(columns=["Gene Symbol", "HGNC Gene ID"])
    return d.rename(columns={"Variant Category": "Category",
                             "HGVS g. (GRCh38)": "HGVS g.",
                             "Single-letter code": "Variant"})


def main():
    df = build()
    _check_against_figure(df)
    xlsx = config.TABLES / "Table1_Variants.xlsx"
    tsv = config.TABLES / "Table1_Variants.tsv"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        df.to_excel(w, sheet_name="Table1_Variants", index=False)
    df.to_csv(tsv, sep="\t", index=False)
    tablepdf.render_table_pdf(
        _pdf_view(df), "Table1_Variants", title=TITLE, caption=CAPTION, note=NOTE,
        group_cols=["Category"], orient="landscape", rich_cols=["Gene; HGNC ID"],
        aligns={"Variant": "center"})
    n = df["Variant Category"].value_counts()
    print(f"INFO: wrote tables/{xlsx.name} + .tsv + .pdf ({len(df)} variants: {dict(n)})")


if __name__ == "__main__":
    main()
