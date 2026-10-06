#!/usr/bin/env python3
"""
Supplementary Table 2 — phenotypic overlap between the UDN proband and the
reported KURIS/Kury-Isidor cohort.

For each clinical feature: its frequency in the Kury et al. 2022 KURIS cohort
(n=11) and whether it is present in the UDN participant (with detail). This is a
hand-curated table; its content lives in the input data file below and is simply
typeset here (Excel + TSV + PDF), matching the other supplementary tables.

Input:  data/source/SuppTable2_phenotypic_overlap.csv
Output: tables/SuppTable2_phenotypic_overlap.{xlsx,tsv,pdf}
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, tablepdf, supplement

INPUT = config.SOURCE / "SuppTable2_phenotypic_overlap.csv"


def build():
    df = pd.read_csv(INPUT)
    return df.rename(columns=lambda c: c.replace("Kury et al", "Küry et al"))


def main():
    df = build()
    df.to_csv(config.TABLES / "SuppTable2_phenotypic_overlap.tsv", sep="\t", index=False)
    df.to_excel(config.TABLES / "SuppTable2_phenotypic_overlap.xlsx", index=False)
    kw = dict(
        aligns={"KURIS Cohort (n=11) from Küry et al. 2022": "center"},
        orient="portrait",
        note="HC = head circumference; SNHL = sensorineural hearing loss.")
    tablepdf.render_table_pdf(df, "SuppTable2_phenotypic_overlap", **kw)
    supplement.render_supp_table("SuppTable2_phenotypic_overlap", df, **kw)
    print(f"INFO: wrote tables/SuppTable2_phenotypic_overlap.{{xlsx,tsv,pdf}} ({len(df)} features)")


if __name__ == "__main__":
    main()
