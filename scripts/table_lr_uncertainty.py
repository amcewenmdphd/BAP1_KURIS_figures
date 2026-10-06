#!/usr/bin/env python3
"""
LR-stability (uncertainty) table: for each control set and functional-class
direction, two estimands are each stress-tested for stability --
  - LR (Haldane) = (a/N_path)/(b/N_ben); and
  - OddsPath (Tavtigian) = (a/(N_path+1))/((b+1)/(N_ben+1)) --
each under log-Wald, bootstrap, and leave-one-out, plus a Bayesian estimand
giving the posterior-median LR and 95% credible interval under Beta(0.5, 0.5)
(Jeffreys) and Beta(1, 1) (uniform) priors, and the fraction of resamples/draws
in that estimand's full-data ACMG tier.

The baseline LR and its log-Wald CI are unchanged. Recomputed from the master +
clinical data tables; resampling is seeded for reproducibility.

Output: tables/Supplemental_LR_uncertainty_methods.{tsv,xlsx}
"""

import functools
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, stats, tables, tablepdf, supplement

RNG_SEED = 54321
N_BOOT = 10000
N_BAYES = 20000


def _lr_point(a, b, n1, n2):
    if a == 0 or b == 0:
        a, b, n1, n2 = a + 0.5, b + 0.5, n1 + 1.0, n2 + 1.0
    return (a / n1) / (b / n2)


def _wald_baseline(a, b, n1, n2):
    _, lo, hi, _ = stats.lr_and_ci(a, b, n1, n2)
    return lo, hi


def _wald_tavtigian(a, b, n1, n2, normal=False):
    _, lo, hi = stats.oddspath_and_ci(a, b, n1, n2, normal=normal)
    return lo, hi


def _resample(p_ind, b_ind, n1, n2, point_fn, wald_fn, rng):
    """log-Wald CI + bootstrap-percentile CI + leave-one-out range for one point
    estimator, plus the fraction of resamples in the full-data ACMG tier."""
    a, b = int(p_ind.sum()), int(b_ind.sum())
    point = point_fn(a, b, n1, n2)
    lo_w, hi_w = wald_fn(a, b, n1, n2)
    boot = np.array([point_fn(int(p_ind[rng.integers(0, n1, n1)].sum()),
                              int(b_ind[rng.integers(0, n2, n2)].sum()), n1, n2)
                     for _ in range(N_BOOT)])
    bootstrap = (float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5)))
    loo_vals = ([point_fn(int(p_ind[np.arange(n1) != i].sum()), b, n1 - 1, n2) for i in range(n1)]
                + [point_fn(a, int(b_ind[np.arange(n2) != j].sum()), n1, n2 - 1) for j in range(n2)])
    loo = (min(loo_vals), max(loo_vals))
    full_tier = stats.evidence_tier(point)
    pct = {"bootstrap": 100.0 * float(np.mean([stats.evidence_tier(x) == full_tier for x in boot])),
           "loo": 100.0 * float(np.mean([stats.evidence_tier(x) == full_tier for x in loo_vals]))}
    return dict(a=a, b=b, point=point, wald=(lo_w, hi_w), bootstrap=bootstrap, loo=loo,
                pct=pct, full_tier=full_tier)


def _bayes(a, b, n1, n2, prior, rng):
    """Posterior of the LR under independent Beta priors; returns (median,
    (2.5%, 97.5%) credible interval, fraction of draws in the median's ACMG tier)."""
    draws = rng.beta(a + prior[0], n1 - a + prior[1], N_BAYES) / rng.beta(b + prior[0], n2 - b + prior[1], N_BAYES)
    median = float(np.median(draws))
    ci = (float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5)))
    tier = stats.evidence_tier(median)
    pct = 100.0 * float(np.mean([stats.evidence_tier(x) == tier for x in draws]))
    return median, ci, pct


def _row(truth, lr_type, estimand, a, b, n1, n2, method, point, star, lo, hi, pct, n_iter):
    return {
        "Truth set": truth, "LR type": lr_type, "Estimand": estimand,
        "Pathogenic in class (a)": a, "Benign in class (b)": b,
        "N pathogenic total": n1, "N benign total": n2, "Method": method,
        "Point estimate": stats.fmt_lr(point, star),
        "2.5% (or min)": "" if lo is None else round(lo, 6),
        "97.5% (or max)": "" if hi is None else round(hi, 6),
        "Evidence tier": stats.evidence_tier(point), "N iterations": n_iter,
        "% in same tier": pct,
    }


def build():
    m = tables.load_master()
    m = tables.add_panel3_category(m, tables.load_clinical())
    kuris_path, kuris_benign = tables.kuris_controls(m)
    col = "score_threshold_class"
    miss = m["variant_type"] == "Missense"
    sets = {
        "Clinical - All variants": (m[m["clinvar_2026_slim"] == "P/LP"], m[m["clinvar_2026_slim"] == "B/LB"]),
        "Clinical - Missense only": (m[miss & (m["clinvar_2026_slim"] == "P/LP")], m[miss & (m["clinvar_2026_slim"] == "B/LB")]),
        "KURIS/NDD vs B/LB missense": (kuris_path, kuris_benign),
    }
    rng = np.random.default_rng(RNG_SEED)
    rows = []
    for truth, (pdf, bdf) in sets.items():
        for lr_type in ("Abnormal", "Normal"):
            pl, bl = pdf[col].values, bdf[col].values
            n1, n2 = len(pl), len(bl)
            p_ind = (pl == lr_type).astype(float); b_ind = (bl == lr_type).astype(float)
            a, b = int(p_ind.sum()), int(b_ind.sum())
            corr = (a == 0 or b == 0)

            # Estimands 1-2: LR (Haldane) and OddsPath (Tavtigian), each by log-Wald / Bootstrap / LOO.
            # OddsPath's +1 pseudo-count goes to whichever cell is the error cell for this
            # lr_type (benign for Abnormal, pathogenic for Normal) -- see stats.oddspath_and_ci.
            is_normal = (lr_type == "Normal")
            for E, point_fn, wald_fn, star in (
                    ("LR (Haldane)", _lr_point, _wald_baseline, corr),
                    ("OddsPath (Tavtigian)", functools.partial(stats.lr_tavtigian, normal=is_normal),
                     functools.partial(_wald_tavtigian, normal=is_normal), False)):
                r = _resample(p_ind, b_ind, n1, n2, point_fn, wald_fn, rng)
                rows += [
                    _row(truth, lr_type, E, a, b, n1, n2, "log-Wald",
                         r["point"], star, *r["wald"], "N/A (analytic)", ""),
                    _row(truth, lr_type, E, a, b, n1, n2, "Bootstrap",
                         r["point"], star, *r["bootstrap"], f'{r["pct"]["bootstrap"]:.1f}%', N_BOOT),
                    _row(truth, lr_type, E, a, b, n1, n2, "LOO",
                         r["point"], star, *r["loo"], f'{r["pct"]["loo"]:.1f}%', n1 + n2),
                ]
            # Estimand 3: Bayesian credible interval (posterior median point) for two priors.
            for label, prior in (("Beta(0.5, 0.5)", (0.5, 0.5)), ("Beta(1, 1)", (1, 1))):
                med, ci, pct = _bayes(a, b, n1, n2, prior, rng)
                rows.append(_row(truth, lr_type, "Bayesian", a, b, n1, n2, label,
                                 med, False, *ci, f"{pct:.1f}%", N_BAYES))
    return pd.DataFrame(rows)


def _pdf_view(stab):
    d = stab.copy()
    for c in ("2.5% (or min)", "97.5% (or max)"):
        d[c] = d[c].map(lambda v: "" if v == "" or pd.isna(v) else f"{float(v):.3g}")
    d = d.drop(columns=["N iterations"])
    return d.rename(columns={
        "Pathogenic in class (a)": "a", "Benign in class (b)": "b",
        "N pathogenic total": "N path", "N benign total": "N benign",
        "Point estimate": "Estimate", "2.5% (or min)": "Lo", "97.5% (or max)": "Hi",
        "% in same tier": "% in tier"})


def main():
    stab = build()
    stab.to_csv(config.TABLES / "Supplemental_LR_uncertainty_methods.tsv", sep="\t", index=False)
    stab.to_excel(config.TABLES / "Supplemental_LR_uncertainty_methods.xlsx", index=False)
    kw = dict(
        aligns={"a": "center", "b": "center", "N path": "center", "N benign": "center",
                "Estimate": "right", "Lo": "right", "Hi": "right",
                "Evidence tier": "center", "% in tier": "center"},
        group_cols=["Truth set", "LR type", "Estimand"], orient="landscape",
        note=[
            "Estimands: LR (Haldane) = (a/N_path)/(b/N_ben); OddsPath (Tavtigian) adds +1 to "
            "whichever in-class count is the error cell for that row -- benign (a/(N_path+1))/"
            "((b+1)/(N_ben+1)) for Abnormal rows, pathogenic ((a+1)/(N_path+1))/(b/(N_ben+1)) "
            "for Normal rows; Bayesian = posterior-median LR under Beta(0.5,0.5) (Jeffreys) "
            "and Beta(1,1) (uniform) priors. The LR (Haldane) point and its log-Wald CI are the "
            "primary baseline; the others are stability checks.",
            "Lo / Hi: analytic 95% CI (log-Wald); 95% percentile interval (Bootstrap); min and "
            "max across refits (LOO = leave-one-out, each control dropped once); 95% credible "
            "interval (Bayesian). '% in tier' is the fraction of resamples/draws in the "
            "estimand's full-data ACMG tier.",
        ],
        star_note="* Haldane–Anscombe correction (+0.5 to all four cells) applied "
                  "when a functional-class cell is empty.")
    tablepdf.render_table_pdf(_pdf_view(stab), "Supplemental_LR_uncertainty_methods", **kw)
    supplement.render_supp_table("Supplemental_LR_uncertainty_methods", _pdf_view(stab), **kw)
    print(f"INFO: wrote tables/Supplemental_LR_uncertainty_methods.{{tsv,xlsx,pdf}} ({len(stab)} rows)")


if __name__ == "__main__":
    main()
