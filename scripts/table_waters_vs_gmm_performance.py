#!/usr/bin/env python3
"""
Waters vs GMM classification performance: TP/FP/FN/TN, indeterminate counts, and
sensitivity / specificity / PPV / NPV (Wilson 95% CIs) for each classifier
against two truth sets -- clinical (ClinVar P/LP vs B/LB) and molecular
(truncating vs synonymous). Recomputed from the master data table.

Waters' classifier folds unchanged + enriched into 'normal' (its own definition,
no indeterminate); GMM uses the score-threshold class.

Output: tables/Waters_vs_GMM_performance.{tsv,xlsx}
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, palettes, stats, tables, tablepdf, supplement

WATERS_CALL = {"depleted": "Abnormal", "unchanged": "Normal", "enriched": "Normal"}


def _perf(path, benign, call_col):
    pc, bc = path[call_col], benign[call_col]
    TP = int((pc == "Abnormal").sum()); FN = int((pc == "Normal").sum())
    IndP = int((pc == "Indeterminate").sum())
    FP = int((bc == "Abnormal").sum()); TN = int((bc == "Normal").sum())
    IndB = int((bc == "Indeterminate").sum())
    return {"TP": TP, "FP": FP, "FN": FN, "TN": TN,
            "Indeterminate (path)": IndP, "Indeterminate (benign)": IndB,
            "Sensitivity": stats.fmt_prop_ci(TP, TP + FN), "Specificity": stats.fmt_prop_ci(TN, TN + FP),
            "PPV": stats.fmt_prop_ci(TP, TP + FP), "NPV": stats.fmt_prop_ci(TN, TN + FN)}


def build():
    m = tables.load_master()
    m["waters_call"] = m["functional_classification"].map(WATERS_CALL)
    m["gmm_call"] = m["score_threshold_class"]
    clin = np.where(m["clinvar_2026_slim"] == "P/LP", "path",
                    np.where(m["clinvar_2026_slim"] == "B/LB", "benign", "na"))
    mol = np.where(m["variant_type"].isin(palettes.TRUNCATING_TYPES), "path",
                   np.where(m["variant_type"] == "Synonymous", "benign", "na"))
    rows = []
    for title, truth in (("Clinical (ClinVar P/LP vs B/LB)", clin),
                         ("Molecular (Truncating vs Synonymous)", mol)):
        t = pd.Series(truth, index=m.index)
        for method, col in (("Waters", "waters_call"), ("GMM", "gmm_call")):
            rows.append({"Truth set": title, "Method": method,
                         **_perf(m[t == "path"], m[t == "benign"], col)})
    return pd.DataFrame(rows)


def _pdf_view(perf):
    return perf.rename(columns={
        "Indeterminate (path)": "Indet. (path)", "Indeterminate (benign)": "Indet. (benign)"})


def main():
    perf = build()
    perf.to_csv(config.TABLES / "Waters_vs_GMM_performance.tsv", sep="\t", index=False)
    perf.to_excel(config.TABLES / "Waters_vs_GMM_performance.xlsx", index=False)
    kw = dict(
        aligns={c: "center" for c in ["TP", "FP", "FN", "TN", "Indet. (path)",
                "Indet. (benign)", "Sensitivity", "Specificity", "PPV", "NPV"]},
        group_cols=["Truth set"])
    tablepdf.render_table_pdf(_pdf_view(perf), "Waters_vs_GMM_performance", **kw)
    supplement.render_supp_table("Waters_vs_GMM_performance", _pdf_view(perf), **kw)
    print(f"INFO: wrote tables/Waters_vs_GMM_performance.{{tsv,xlsx,pdf}} ({len(perf)} rows)")


if __name__ == "__main__":
    main()
