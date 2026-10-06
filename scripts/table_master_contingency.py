#!/usr/bin/env python3
"""
Master contingency table: likelihood ratio (with log-Wald 95% CI and ACMG
evidence tier) for every functional-class assignment method x control set x
functional class, recomputed from the master data table + the clinical data
table.

Three assignment methods, each with its class definition made explicit:
  - GMM cutoffs        Abnormal / Indeterminate / Normal by functional-score cutoffs.
  - Waters             SGE FDR + directionality calls (depleted / unchanged / enriched);
                       Abnormal = depleted, Normal = unchanged + enriched (no indeterminate).
  - Tejura             same SGE FDR + directionality calls; Abnormal = depleted,
                       Indeterminate = enriched, Normal = unchanged.

Control sets: all ClinVar P/LP vs B/LB (2025 / 2026), missense-only,
KURIS/NDD vs B/LB missense, and the three Waters et al. 2024 Table 1 validation
sets (their published counts; LR/CI recomputed here).

Output: tables/Master_Contingency_Tables.{tsv,xlsx,pdf}
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, stats, tables, tablepdf, supplement

ABN, NRM = stats.load_cutoffs()

# Each classifier: how the call is assigned, the source column, and — per
# functional class — the raw call(s) it comprises and a plain class definition.
CLASSIFIERS = {
    # Indeterminate classes are omitted: they carry no directional evidence.
    "GMM cutoffs": dict(
        method="GMM score cutoffs",
        col="score_threshold_class",
        classes=[
            ("Abnormal", ["Abnormal"], f"score ≤ {ABN:.4f}"),
            ("Normal", ["Normal"], f"score ≥ {NRM:.4f}"),
        ]),
    "Waters": dict(
        method="Waters: SGE FDR + directionality",
        col="functional_classification",
        classes=[
            ("Abnormal", ["depleted"], "depleted"),
            ("Normal", ["unchanged", "enriched"], "unchanged + enriched"),
        ]),
    "Tejura": dict(
        method="Tejura: SGE FDR + directionality",
        col="functional_classification",
        classes=[
            ("Abnormal", ["depleted"], "depleted"),
            ("Normal", ["unchanged"], "unchanged"),
        ]),
}
# (classifier, ClinVar-slim column, control-set snapshot label)
METHOD_YEARS = [
    ("GMM cutoffs", "clinvar_2026_slim", "ClinVar Feb 2026"),
    ("GMM cutoffs", "clinvar_2025_slim", "ClinVar 2025"),
    ("Waters", "clinvar_2026_slim", "ClinVar Feb 2026"),
    ("Tejura", "clinvar_2025_slim", "ClinVar 2025"),
]
# Waters et al. 2024 Table 1: (label, N_path, N_ben, dep_path, dep_ben, ue_path, ue_ben)
WATERS_TABLE1 = [
    ("Waters ClinVar ≥2* (Table 1)", 0, 6, 0, 0, 0, 6),
    ("Waters ClinVar ≥1* (Table 1)", 7, 6, 7, 0, 0, 6),
    ("Waters Systematic (Table 1)", 2423, 138, 2419, 4, 4, 134),
]


def _extras(disp, a, b, n1, n2):
    """Tavtigian OddsPath (Abnormal/Normal classes only) with its log-Wald 95% CI
    and ACMG evidence tier. A stability check that never affects the baseline."""
    if disp not in ("Abnormal", "Normal"):
        return "—", "—", "—"
    op, lo, hi = stats.oddspath_and_ci(a, b, n1, n2, normal=(disp == "Normal"))
    if not np.isfinite(op) or op == 0:
        return "—", "—", "—"
    op_s = stats.fmt_lr(op)
    ci_s = "—" if (np.isnan(lo) or np.isnan(hi)) else f"({stats.fmt_lr(lo)}–{stats.fmt_lr(hi)})"
    return op_s, ci_s, stats.evidence_tier(op)


def _row(classifier, disp, defn, control, a, b, n1, n2, lr_s, ci_s, tier, op_s, op_ci, op_tier):
    return {
        "Assignment method": CLASSIFIERS[classifier]["method"],
        "Functional class": disp, "Class definition": defn, "Control set": control,
        "Pathogenic controls (total)": n1, "Benign controls (total)": n2,
        "Pathogenic in class": a, "Benign in class": b,
        "Likelihood ratio": lr_s, "95% CI": ci_s, "Evidence strength": tier,
        "OddsPath (Tavtigian)": op_s, "OddsPath 95% CI": op_ci, "OddsPath evidence": op_tier,
    }


def _rows_for(classifier, path_df, benign_df, control):
    spec = CLASSIFIERS[classifier]
    col = spec["col"]
    n1, n2 = len(path_df), len(benign_df)
    rows = []
    for disp, raw, defn in spec["classes"]:
        a = int(path_df[col].isin(raw).sum())
        b = int(benign_df[col].isin(raw).sum())
        lr, lo, hi, corr = stats.lr_and_ci(a, b, n1, n2)
        op_s, op_ci, op_tier = _extras(disp, a, b, n1, n2)
        rows.append(_row(classifier, disp, defn, control, a, b, n1, n2,
                         stats.fmt_lr(lr, corr), f"({stats.fmt_lr(lo)}–{stats.fmt_lr(hi)})",
                         stats.evidence_tier(lr), op_s, op_ci, op_tier))
    return rows


def _waters_table1_rows():
    rows = []
    for label, n1, n2, dp, db, up, ub in WATERS_TABLE1:
        for disp, defn, a, b in (("Abnormal", "depleted", dp, db),
                                 ("Normal", "unchanged + enriched", up, ub)):
            if n1 == 0 or n2 == 0:
                lr_s, ci_s, tier = "—", "not computable (n=0)", "—"
                op_s, op_ci, op_tier = "—", "—", "—"
            else:
                lr, lo, hi, corr = stats.lr_and_ci(a, b, n1, n2)
                lr_s = stats.fmt_lr(lr, corr)
                ci_s = f"({stats.fmt_lr(lo)}–{stats.fmt_lr(hi)})"
                tier = stats.evidence_tier(lr)
                op_s, op_ci, op_tier = _extras(disp, a, b, n1, n2)
            rows.append(_row("Waters", disp, defn, label, a, b, n1, n2, lr_s, ci_s, tier,
                             op_s, op_ci, op_tier))
    return rows


def build():
    m = tables.load_master()
    m = tables.add_panel3_category(m, tables.load_clinical())
    kuris_path, kuris_benign = tables.kuris_controls(m)
    miss = m["variant_type"] == "Missense"
    rows = []
    for classifier, year_col, snapshot in METHOD_YEARS:
        for scope, mask in (("all variants", slice(None)), ("missense only", miss)):
            sub = m if isinstance(mask, slice) else m[mask]
            rows += _rows_for(classifier, sub[sub[year_col] == "P/LP"],
                              sub[sub[year_col] == "B/LB"], f"{snapshot} ({scope})")
    for classifier in CLASSIFIERS:
        rows += _rows_for(classifier, kuris_path, kuris_benign,
                          "KURIS/NDD vs B/LB missense")
    rows += _waters_table1_rows()
    return pd.DataFrame(rows)


def _pdf_view(cont):
    d = cont.copy()
    # compress the four count columns into two "in-class / total" fractions
    d["Path. (in/N)"] = d["Pathogenic in class"].astype(str) + "/" + d["Pathogenic controls (total)"].astype(str)
    d["Ben. (in/N)"] = d["Benign in class"].astype(str) + "/" + d["Benign controls (total)"].astype(str)
    d = d.drop(columns=["Pathogenic controls (total)", "Benign controls (total)",
                        "Pathogenic in class", "Benign in class"])
    d = d.rename(columns={
        "Likelihood ratio": "LR", "Evidence strength": "Evidence",
        "OddsPath (Tavtigian)": "OddsPath", "OddsPath 95% CI": "OddsPath 95% CI",
        "OddsPath evidence": "OddsPath evid."})
    return d[["Assignment method", "Functional class", "Class definition", "Control set",
              "Path. (in/N)", "Ben. (in/N)", "LR", "95% CI", "Evidence",
              "OddsPath", "OddsPath 95% CI", "OddsPath evid."]]


def main():
    cont = build()
    cont.to_csv(config.TABLES / "Master_Contingency_Tables.tsv", sep="\t", index=False)
    cont.to_excel(config.TABLES / "Master_Contingency_Tables.xlsx", index=False)
    kw = dict(
        aligns={"Path. (in/N)": "center", "Ben. (in/N)": "center", "LR": "right",
                "95% CI": "center", "Evidence": "center", "OddsPath": "right",
                "OddsPath 95% CI": "center", "OddsPath evid.": "center"},
        group_cols=["Assignment method", "Control set"], orient="landscape",
        note=["Baseline: LR = (a/N_path)/(b/N_ben); 95% CI is log-Wald. The OddsPath columns "
              "are a stability check and do not affect the baseline.",
              "OddsPath = Tavtigian/Brnich likelihood ratio with +1 added to whichever "
              "in-class count is the error cell for that row -- benign (a false positive) "
              "for Abnormal rows: (a/(n1+1))/((b+1)/(n2+1)); pathogenic (a false negative) "
              "for Normal rows: ((a+1)/(n1+1))/(b/(n2+1)). Its 95% CI is log-Wald and "
              "'OddsPath evid.' is its ACMG tier."],
        star_note="* Haldane–Anscombe correction (+0.5 to all four cells) applied "
                  "when a functional-class cell is empty.")
    tablepdf.render_table_pdf(_pdf_view(cont), "Master_Contingency_Tables", **kw)
    supplement.render_supp_table("Master_Contingency_Tables", _pdf_view(cont), **kw)
    print(f"INFO: wrote tables/Master_Contingency_Tables.{{tsv,xlsx,pdf}} ({len(cont)} rows)")


if __name__ == "__main__":
    main()
