"""Shared statistics: GMM classification, likelihood ratios, ACMG evidence tiers,
and Wilson confidence intervals. Used by every figure and table so the numbers
are identical everywhere.
"""

import json
import math

import numpy as np

from . import config


# ── GMM score-threshold classification ──

def load_cutoffs(path=None):
    with open(path or config.GMM_THRESHOLDS) as fh:
        thr = json.load(fh)
    return thr["abnormal_cutoff"], thr["normal_cutoff"]


def gmm_class(score, abn=None, nrm=None):
    if abn is None or nrm is None:
        abn, nrm = load_cutoffs()
    if score is None or (isinstance(score, float) and math.isnan(score)):
        return "Indeterminate"
    if score <= abn:
        return "Abnormal"
    if score >= nrm:
        return "Normal"
    return "Indeterminate"


# ── Likelihood ratio with log-Wald CI ──

def lr_and_ci(a, b, n1, n2, alpha=0.05):
    """LR = (a/n1)/(b/n2) with a log-Wald 95% CI. Haldane-Anscombe +0.5 is added
    to all four cells of the 2x2 (n1 -> n1+1, n2 -> n2+1) only when a class cell
    (a or b) is empty -- the standard ENIGMA/ACMG-SVI convention. Returns
    (lr, lo, hi, corrected)."""
    corrected = (a == 0 or b == 0)
    if corrected:
        a_, b_ = a + 0.5, b + 0.5
        n1_, n2_ = n1 + 1.0, n2 + 1.0
    else:
        a_, b_, n1_, n2_ = float(a), float(b), float(n1), float(n2)
    lr = (a_ / n1_) / (b_ / n2_)
    se = math.sqrt(1 / a_ - 1 / n1_ + 1 / b_ - 1 / n2_)
    z = 1.959963984540054 if alpha == 0.05 else _z(alpha)
    return lr, lr * math.exp(-z * se), lr * math.exp(z * se), corrected


def _z(alpha):
    from scipy.stats import norm
    return norm.ppf(1 - alpha / 2)


# ── Extra LR estimators for stability vetting (do NOT replace the baseline) ──

def oddspath_and_ci(a, b, n1, n2, alpha=0.05, normal=False):
    """Tavtigian/Brnich OddsPath for class membership, with a log-Wald 95% CI.

    OddsPath is the posterior-odds / prior-odds ratio, which by Bayes' theorem
    equals the likelihood ratio LR+ = P(in class | path) / P(in class | benign) =
    sensitivity / (1 - specificity). It is NOT the cross-product odds ratio
    (TP*TN)/(FP*FN). The correction adds one hypothetical misclassification to
    whichever cell is the actual error cell for the class being scored, so it stays
    conservative in both directions:
      - Abnormal class (normal=False, the default): a = pathogenic-in-class (TP),
        b = benign-in-class (a false positive) -> +1 goes to b:
            OddsPath = (a / (n1 + 1)) / ((b + 1) / (n2 + 1))
      - Normal class (normal=True): a = pathogenic-in-class (a false negative),
        b = benign-in-class (TN) -> +1 goes to a instead:
            OddsPath = ((a + 1) / (n1 + 1)) / (b / (n2 + 1))
    Adding the pseudo-count to the correct-call cell instead (e.g. b for the
    Normal class) would understate the correction and isn't used here.

    A stability check reported alongside the baseline LR, never replacing it.

    Returns (oddspath, lo, hi); lo/hi are NaN when the class with no pseudo-count
    (a for Abnormal, b for Normal) is 0 -- i.e. no data on that side at all.
    """
    if normal:
        a_, b_, n1_, n2_ = float(a + 1), float(b), float(n1 + 1), float(n2 + 1)
    else:
        a_, b_, n1_, n2_ = float(a), float(b + 1), float(n1 + 1), float(n2 + 1)
    if a_ == 0 or b_ == 0:
        return 0.0, float("nan"), float("nan")
    op = (a_ / n1_) / (b_ / n2_)
    se = math.sqrt(1 / a_ - 1 / n1_ + 1 / b_ - 1 / n2_)
    z = 1.959963984540054 if alpha == 0.05 else _z(alpha)
    return op, op * math.exp(-z * se), op * math.exp(z * se)


def lr_tavtigian(a, b, n1, n2, normal=False):
    """Tavtigian OddsPath point estimate only (see oddspath_and_ci)."""
    return oddspath_and_ci(a, b, n1, n2, normal=normal)[0]


def bayes_lr_ci(a, b, n1, n2, prior=(0.5, 0.5), draws=20000, seed=54321, alpha=0.05):
    """Bayesian credible interval for the baseline LR = (a/n1)/(b/n2) with
    independent Beta priors on each proportion. prior=(0.5,0.5) is Jeffreys;
    prior=(1,1) is uniform. Returns (lo, hi, draws_array) for the LR."""
    rng = np.random.default_rng(seed)
    pa = rng.beta(a + prior[0], n1 - a + prior[1], draws)
    pb = rng.beta(b + prior[0], n2 - b + prior[1], draws)
    ratio = pa / pb
    lo = float(np.percentile(ratio, 100 * alpha / 2))
    hi = float(np.percentile(ratio, 100 * (1 - alpha / 2)))
    return lo, hi, ratio


# ── ACMG/Tavtigian (2018) evidence tiers ──

_PS3 = [(2.08, "Supporting"), (4.33, "Moderate"), (18.7, "Strong"), (350.0, "Very Strong")]
_BS3 = [(1 / 2.08, "Supporting"), (1 / 4.33, "Moderate"), (1 / 18.7, "Strong"), (1 / 350.0, "Very Strong")]
_COMPACT = {"Supporting": "Supptg.", "Moderate": "Moderate", "Strong": "Strong",
            "Very Strong": "V.Strong"}


def evidence_tier(lr, style="full"):
    """Map an LR to its ACMG evidence tier. style='full' -> 'PS3 Strong';
    style='compact' -> 'PS3 Strong'/'PS3 Supptg.'/'Indet.' (Figure 2 labels)."""
    if not np.isfinite(lr) or lr <= 0:
        return "Indeterminate" if style == "full" else "Indet."
    if lr >= 2.08:
        name = "Supporting"
        for thr, nm in _PS3:
            if lr >= thr:
                name = nm
        return f"PS3 {name}" if style == "full" else f"PS3 {_COMPACT[name]}"
    if lr <= 1 / 2.08:
        name = "Supporting"
        for thr, nm in _BS3:
            if lr <= thr:
                name = nm
        return f"BS3 {name}" if style == "full" else f"BS3 {_COMPACT[name]}"
    return "Indeterminate" if style == "full" else "Indet."


# ── Wilson score CI for a binomial proportion ──

def wilson_ci(k, n, z=1.96):
    if n == 0:
        return float("nan"), float("nan"), float("nan")
    p = k / n
    denom = 1 + z ** 2 / n
    centre = (p + z ** 2 / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / denom
    return p, max(0.0, centre - half), min(1.0, centre + half)


def fmt_lr(lr, corrected=False):
    if lr >= 100:
        s = f"{lr:.0f}"
    elif lr >= 1:
        s = f"{lr:.2f}"
    else:
        s = f"{lr:.4f}"
    return s + ("*" if corrected else "")


def fmt_prop_ci(k, n):
    p, lo, hi = wilson_ci(k, n)
    return "—" if n == 0 else f"{p:.3f} ({lo:.3f}–{hi:.3f})"
