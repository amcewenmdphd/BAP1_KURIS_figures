#!/usr/bin/env python3
"""
Build the BAP1 KURIS master data table from the Waters et al. BAP1 SGE release.

This is the single, self-contained, traceable script that produces the annotated
master table every downstream figure and analysis is built from. It starts from
Waters Extended Data 1 (the authoritative SGE dataset), keeps **every** original
row and column unchanged, and appends derived analysis columns. The original
Waters workbook is opened read-only and is never modified.

Sources (read-only)
-------------------
    data/raw/Waters_ExtendedData1.xlsx   (sheet 'extended_data_1', 18,108 variants,
        87 columns; https://github.com/team113sanger/Waters_BAP1_SGE)
    data/raw/BAP1_cvfg_annotations.tsv   (IGVF Coding Variant Focus Group BAP1
        annotations: ClinVar 2018/2025, gnomAD v4.1, REVEL, AlphaMissense,
        SpliceAI; merged by HGVSc to add the figure/table truth sets and predictors)
    data/raw/ClinVar_Feb_2026.txt        (ClinVar BAP1 VCV export, February 2026 --
        the version cited in the paper; provides the clinvar_2026_* columns)

Outputs
-------
    data/source/BAP1_master_table.tsv          annotated master table (all input
                                               columns + the derived columns below)
    data/source/GMM_thresholds.json            fitted GMM parameters, thresholds
                                               and full provenance
    data/source/BAP1_master_table_dictionary.md  data dictionary for derived columns

Derived columns appended (originals are left untouched)
-------------------------------------------------------
    simplified_consequence   single VEP class collapsed from ``vep_consequence`` by
                             descending severity; exactly reproduces the original
                             notebook's ``simplified_consequence`` column.
    display_consequence      same classes, human-readable labels for figures.
    gmm_posterior_abnormal   GMM posterior probability of the abnormal component.
    gmm_posterior_normal     GMM posterior probability of the normal component.
    functional_consequence   posterior-based call: functionally_abnormal /
                             functionally_normal / indeterminate (>= 0.95 posterior).
    score_threshold_class    score-based call used by the figures:
                             Abnormal (score <= abnormal cutoff) /
                             Normal   (score >= normal cutoff) / Indeterminate.
    standardized_class       Waters functional_classification standardized to
                             Abnormal / Normal / Not specified (GMM-vs-Waters figure).
    clinvar_{2018,2025,2026}_{classification,slim,review_status,stars}
                             ClinVar at three timepoints (CVFG) -- the calibration
                             truth sets used by the LR / KURIS-NDD figures and tables.
    gnomad_v4_af, gnomad_v4_faf95_max, revel_score, alphamissense_pathogenicity,
    alphamissense_class, spliceai_max_delta
                             gnomAD v4.1 and in-silico predictors (CVFG).

Exclusions
----------
None. The Waters authors pre-filtered pam_flag=Y variants before releasing
Extended Data 1 (its ``pam_flag`` column is uniformly ``N``), so ED1 as delivered
is already the intended calibration set. All 18,108 variants are used; ``snvre``,
``pam_codon=Y``, and discordant variants are all retained.

GMM method (exactly as the original BAP1_GMM notebook)
------------------------------------------------------
A two-component Gaussian mixture is seeded from the empirical means of two anchor
sets and fit to their pooled scores:
  * abnormal component -> ``stop_gained`` variants inside an hg38 gene-body window
    (trims the first two exons and the extreme 3' end, limiting NMD-escape / assay
    edge effects);
  * normal component   -> ``synonymous_variant`` variants inside their central 95%.
Every scored variant then receives component posteriors; the abnormal/normal score
thresholds are the points where the posterior crosses 0.95.

Usage
-----
    python scripts/make_master_table.py
"""

import json
import os

import numpy as np
import pandas as pd
import sklearn
from sklearn.mixture import GaussianMixture

# ── Paths (resolved from repo root) ──
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
SOURCE_XLSX = os.path.join(REPO_ROOT, 'data', 'raw', 'Waters_ExtendedData1.xlsx')
SOURCE_SHEET = 'extended_data_1'
SOURCE_HEADER_ROW = 2   # 0-indexed: title row, blank row, then the header
CVFG_TSV = os.path.join(REPO_ROOT, 'data', 'raw', 'BAP1_cvfg_annotations.tsv')
CLINVAR_FEB2026_TXT = os.path.join(REPO_ROOT, 'data', 'raw', 'ClinVar_Feb_2026.txt')
OUTPUT_DIR = os.path.join(REPO_ROOT, 'data', 'source')
MASTER_TSV = os.path.join(OUTPUT_DIR, 'BAP1_master_table.tsv')
THRESHOLDS_JSON = os.path.join(OUTPUT_DIR, 'GMM_thresholds.json')
DICTIONARY_MD = os.path.join(OUTPUT_DIR, 'BAP1_master_table_dictionary.md')

# ── Source column names ──
SCORE_COL = 'functional_score'
VEP_COL = 'vep_consequence'
POSITION_COL = 'pos'

# ── GMM parameters (from the original BAP1_GMM notebook) ──
NONSENSE_POS_LOW = 52402674     # hg38 window trimming first two exons / 3' end
NONSENSE_POS_HIGH = 52409714
POSTERIOR_CUTOFF = 0.95
RANDOM_STATE = 54345


# ClinVar clinical significance -> collapsed label.
CLINVAR_SLIM = {
    'Pathogenic': 'P/LP',
    'Likely pathogenic': 'P/LP',
    'Pathogenic/Likely pathogenic': 'P/LP',
    'Benign': 'B/LB',
    'Likely benign': 'B/LB',
    'Benign/Likely benign': 'B/LB',
    'Uncertain significance': 'VUS',
    'Conflicting classifications of pathogenicity': 'Conflicting',
    'Conflicting interpretations of pathogenicity': 'Conflicting',
}

# Waters functional_classification -> standardized clinical-direction class.
STANDARDIZED_CLASS = {
    'depleted': 'Abnormal',
    'unchanged': 'Normal',
    'enriched': 'Not specified',
}

# CVFG ClinVar timepoints -> output column prefix (2018/2025 trajectory).
# The 2026 timepoint is taken from the February export (cited in the paper).
CVFG_CLINVAR_TIMEPOINTS = {
    'clinvar.201801': 'clinvar_2018',
    'clinvar.202501': 'clinvar_2025',
}

# ClinVar review status text -> review star level.
CLINVAR_REVIEW_STARS = {
    'practice guideline': 4,
    'reviewed by expert panel': 3,
    'criteria provided, multiple submitters, no conflicts': 2,
    'criteria provided, conflicting classifications': 1,
    'criteria provided, conflicting interpretations': 1,
    'criteria provided, single submitter': 1,
    'no assertion criteria provided': 0,
    'no classification provided': 0,
    'no classifications from unflagged records': 0,
}


def simplify_consequence(vep):
    """Collapse a (possibly compound) VEP consequence string to a single class by
    descending severity. This exactly reproduces the ``simplified_consequence``
    column of the original GMM notebook's input, verified against BAP1scores.tsv
    (14,165/14,165)."""
    v = str(vep)
    if 'stop_gained' in v:
        return 'stop_gained'
    if 'frameshift' in v:
        return 'frameshift_variant'
    if 'splice_acceptor' in v or 'splice_donor' in v:
        return 'splice_site_variant'
    if 'start_lost' in v:
        return 'start_lost'
    if 'stop_lost' in v:
        return 'stop_lost'
    if 'inframe' in v:
        return 'inframe_indel'
    if 'missense' in v:
        return 'missense_variant'
    if 'splice_region' in v:
        return 'splicing_variant'
    if 'synonymous' in v or 'stop_retained' in v:
        return 'synonymous_variant'
    if 'intron' in v:
        return 'intron_variant'
    if 'UTR' in v.upper() or 'prime' in v:
        return 'UTR_variant'
    return 'other'


# simplified_consequence -> human-readable label (notebook's _rename_consequences).
DISPLAY_LABELS = {
    'missense_variant': 'Missense',
    'frameshift_variant': 'Frameshift',
    'synonymous_variant': 'Synonymous',
    'intron_variant': 'Intron',
    'stop_gained': 'Stop Gained',
    'stop_lost': 'Stop Lost',
    'splice_site_variant': 'Canonical Splice',
    'splicing_variant': 'Splice Region',
    'UTR_variant': 'UTR',
    'start_lost': 'Start Lost',
    'inframe_indel': 'In-frame Indel',
}

NONSENSE_LABEL = 'stop_gained'
SYNONYMOUS_LABEL = 'synonymous_variant'


# Coarser consequence grouping used by the Figure 1/2 panels (distinct from the
# GMM's simplified_consequence -- folds splice_region into Intronic/Synonymous and
# calls stop_gained "Nonsense"). Reproduces the notebook's CONSEQUENCE_MAP.
VARIANT_TYPE_MAP = {
    'stop_gained': 'Nonsense',
    'stop_gained,splice_region_variant': 'Nonsense',
    'stop_gained,frameshift_variant': 'Nonsense',
    'stop_gained,inframe_deletion': 'Nonsense',
    'stop_gained,start_lost': 'Nonsense',
    'frameshift_variant': 'Frameshift',
    'frameshift_variant,splice_region_variant': 'Frameshift',
    'frameshift_variant,start_lost': 'Frameshift',
    'frameshift_variant,stop_lost': 'Frameshift',
    'splice_acceptor_variant': 'Splice',
    'splice_acceptor_variant,coding_sequence_variant,intron_variant': 'Splice',
    'splice_acceptor_variant,intron_variant': 'Splice',
    'splice_donor_variant': 'Splice',
    'splice_donor_variant,coding_sequence_variant,intron_variant': 'Splice',
    'splice_region_variant,intron_variant': 'Intronic',
    'splice_region_variant,synonymous_variant': 'Synonymous',
    'start_lost': 'Start lost',
    'start_lost,inframe_deletion': 'Start lost',
    'stop_lost': 'Stop lost',
    'stop_lost,inframe_deletion': 'Stop lost',
    'stop_retained_variant': 'Synonymous',
    'inframe_deletion': 'Inframe indel',
    'inframe_deletion,splice_region_variant': 'Inframe indel',
    'inframe_insertion': 'Inframe indel',
    'missense_variant': 'Missense',
    'missense_variant,splice_region_variant': 'Missense',
    'synonymous_variant': 'Synonymous',
    'intron_variant': 'Intronic',
    '3_prime_UTR_variant': 'UTR',
    '5_prime_UTR_variant': 'UTR',
}


def map_variant_type(vep):
    """Figure-panel variant type from the full VEP consequence (CONSEQUENCE_MAP,
    with the notebook's substring fallback for any unlisted compound term)."""
    v = str(vep)
    if v in VARIANT_TYPE_MAP:
        return VARIANT_TYPE_MAP[v]
    if 'frameshift' in v:
        return 'Frameshift'
    if 'stop_gained' in v or 'nonsense' in v:
        return 'Nonsense'
    if 'splice_acceptor' in v or 'splice_donor' in v:
        return 'Splice'
    if 'missense' in v:
        return 'Missense'
    if 'synonymous' in v:
        return 'Synonymous'
    if 'inframe' in v:
        return 'Inframe indel'
    if 'intron' in v or 'splice_region' in v:
        return 'Intronic'
    if 'UTR' in v:
        return 'UTR'
    if 'start_lost' in v:
        return 'Start lost'
    if 'stop_lost' in v:
        return 'Stop lost'
    return 'Other'


def load_source():
    """Read Waters Extended Data 1 read-only; return the DataFrame unchanged."""
    df = pd.read_excel(SOURCE_XLSX, sheet_name=SOURCE_SHEET, header=SOURCE_HEADER_ROW)
    print('INFO: loaded %d variants x %d columns from %s (read-only)'
          % (len(df), df.shape[1], os.path.relpath(SOURCE_XLSX, REPO_ROOT)))
    return df


def annotate_consequence(df):
    df['simplified_consequence'] = df[VEP_COL].map(simplify_consequence)
    df['display_consequence'] = df['simplified_consequence'].map(DISPLAY_LABELS)
    df['variant_type'] = df[VEP_COL].map(map_variant_type)
    print('INFO: simplified_consequence classes: %s'
          % df['simplified_consequence'].value_counts().to_dict())
    print('INFO: variant_type (figure) classes: %s'
          % df['variant_type'].value_counts().to_dict())
    return df


def annotate_standardized_class(df):
    """Standardized clinical-direction class from the Waters functional_classification
    (depleted -> Abnormal, unchanged -> Normal, enriched -> Not specified). Used by
    the GMM-vs-Waters comparison figure."""
    df['standardized_class'] = df['functional_classification'].map(STANDARDIZED_CLASS)
    return df


def merge_cvfg(df):
    """Merge the CVFG (IGVF Coding Variant Focus Group) annotations used by the
    paper figures/tables onto the master table, joined by Ensembl HGVSc with a
    genomic pos+ref+alt fallback (18,108/18,108 coverage). Adds:
      * ClinVar at three timepoints (2018/2025/2026): classification, slim label,
        review status, stars -- the calibration truth sets;
      * gnomAD v4.1 allele frequency and FAF95;
      * REVEL, AlphaMissense, SpliceAI predictors.
    (The Waters functional class itself is already present as functional_classification;
    its standardized form is added separately.)"""
    cv = pd.read_csv(CVFG_TSV, sep='\t', dtype=str)

    # Two lookup keys: Ensembl HGVSc (primary) and genomic pos+ref+alt (fallback).
    cv_hk = cv.set_index(cv['raw_hgvs_nt'].astype(str).str.strip())
    gkey = (cv['mapped_hgvs_g_start'].astype(str) + '_' + cv['mapped_hgvs_g_ref'].astype(str)
            + '_' + cv['mapped_hgvs_g_alt'].astype(str))
    cv_gk = cv.set_index(gkey)

    df_hk = df['HGVSc'].astype(str).str.strip()
    df_gk = (pd.to_numeric(df[POSITION_COL], errors='coerce').astype('Int64').astype(str)
             + '_' + df['ref'].astype(str) + '_' + df['alt'].astype(str))

    def pull(src_col):
        """Map a CVFG column onto df by HGVSc, filling gaps via the genomic key."""
        by_h = df_hk.map(cv_hk[src_col].groupby(level=0).first())
        by_g = df_gk.map(cv_gk[src_col].groupby(level=0).first())
        return by_h.fillna(by_g)

    # ClinVar timepoints.
    for src_prefix, out_prefix in CVFG_CLINVAR_TIMEPOINTS.items():
        cls = pull('%s.clinical_significance' % src_prefix)
        df['%s_classification' % out_prefix] = cls
        df['%s_slim' % out_prefix] = cls.map(CLINVAR_SLIM)
        df['%s_review_status' % out_prefix] = pull('%s.review_status' % src_prefix)
        df['%s_stars' % out_prefix] = pull('%s.stars' % src_prefix)

    # gnomAD v4.1 and in-silico predictors.
    df['gnomad_v4_af'] = pull('gnomad.v4_1.allele_frequency')
    df['gnomad_v4_faf95_max'] = pull('gnomad.v4_1.faf95_max')
    df['revel_score'] = pull('revel.score')
    df['alphamissense_pathogenicity'] = pull('alphamissense.pathogenicity')
    df['alphamissense_class'] = pull('alphamissense.class')
    df['spliceai_max_delta'] = pull('spliceai.max_delta_score')

    print('INFO: CVFG merged; ClinVar populated: 2018=%d, 2025=%d'
          % (df['clinvar_2018_classification'].notna().sum(),
             df['clinvar_2025_classification'].notna().sum()))
    print('INFO: standardized_class: %s'
          % df['standardized_class'].value_counts().to_dict())
    return df


def merge_clinvar_feb2026(df):
    """Merge the ClinVar February 2026 VCV export (cited in the paper) as the 2026
    timepoint, by genomic pos+ref+alt via Canonical SPDI. Adds clinvar_2026_*."""
    cv = pd.read_csv(CLINVAR_FEB2026_TXT, sep='\t', dtype=str)

    def spdi_key(spdi):
        try:
            _, pos0, ref, alt = str(spdi).split(':')
            return '%d_%s_%s' % (int(pos0) + 1, ref, alt)   # SPDI is 0-based
        except (ValueError, AttributeError):
            return None

    cv['_key'] = cv['Canonical SPDI'].map(spdi_key)
    cv = cv[cv['_key'].notna()].drop_duplicates('_key', keep='first')
    pull = {
        'Germline classification': 'clinvar_2026_classification',
        'Germline review status': 'clinvar_2026_review_status',
        'Germline date last evaluated': 'clinvar_2026_date',
        'Condition(s)': 'clinvar_2026_condition',
        'Oncogenicity classification': 'clinvar_2026_oncogenicity',
        'VariationID': 'clinvar_2026_variation_id',
    }
    cv_small = cv[['_key'] + list(pull)].rename(columns=pull).set_index('_key')

    key = (pd.to_numeric(df[POSITION_COL], errors='coerce').astype('Int64').astype(str)
           + '_' + df['ref'].astype(str) + '_' + df['alt'].astype(str))
    for col in pull.values():
        df[col] = key.map(cv_small[col])
    df['clinvar_2026_slim'] = df['clinvar_2026_classification'].map(CLINVAR_SLIM)
    df['clinvar_2026_stars'] = df['clinvar_2026_review_status'].map(CLINVAR_REVIEW_STARS)

    print('INFO: ClinVar Feb 2026 merged onto %d variants; slim: %s'
          % (int(df['clinvar_2026_classification'].notna().sum()),
             df['clinvar_2026_slim'].value_counts().to_dict()))
    return df


# ClinVar VCVs whose aggregate "germline" classification is not a true germline
# call: every germline-track SCV submission for these is flagged Origin=somatic
# in the raw SCV export (data/raw/BAP1_ClinVar_SCV.tsv) -- i.e. each was found
# by tumor-only sequencing, not inherited/germline testing, despite ClinVar's
# aggregate germline field carrying a P/LP classification for it:
#   VCV000800323  p.Ser596Thr / c.1787G>C  "Likely pathogenic", from multiple
#                 myeloma tumor sequencing (pos=52403241, ref=C, alt=G, hg38)
# Same exclusion used in build_clinical_intermediates.py's `false_germline` set
# (Figure 1b/1c); applied here too so every clinvar_* truth set used for
# LR/ACMG-evidence calibration (Tables 3-5, Figure 2e) excludes it as well.
FALSE_GERMLINE_LOCI = {(52403241, 'C', 'G')}


def exclude_false_germline_clinvar(df):
    """Null every clinvar_2018/2025/2026_* column for the loci in
    FALSE_GERMLINE_LOCI, so they read as ClinVar-unclassified (like any other
    variant absent from ClinVar) rather than as a germline P/LP control."""
    pos = pd.to_numeric(df[POSITION_COL], errors='coerce')
    mask = pd.Series(False, index=df.index)
    for p, ref, alt in FALSE_GERMLINE_LOCI:
        mask |= (pos == p) & (df['ref'] == ref) & (df['alt'] == alt)
    cols = [c for c in df.columns if c.startswith(('clinvar_2018', 'clinvar_2025', 'clinvar_2026'))]
    n = int(mask.sum())
    if n:
        df.loc[mask, cols] = np.nan
        print('INFO: excluded %d false-germline (somatic-only-origin) ClinVar record(s) '
              'from clinvar_* columns (HGVSc: %s)' % (n, df.loc[mask, 'HGVSc'].tolist()))
    else:
        print('WARNING: FALSE_GERMLINE_LOCI matched 0 rows -- check POSITION_COL/ref/alt encoding')
    return df


def fit_gmm(df):
    """Fit the 2-component GMM on the stop_gained and synonymous anchor sets,
    exactly as the original notebook did. No variants are excluded: the authors
    pre-filtered pam_flag=Y upstream, so ED1 as delivered is already the intended
    calibration set."""
    score = pd.to_numeric(df[SCORE_COL], errors='coerce')
    pos = pd.to_numeric(df[POSITION_COL], errors='coerce')

    nonsense = ((df['simplified_consequence'] == NONSENSE_LABEL)
                & (pos > NONSENSE_POS_LOW) & (pos <= NONSENSE_POS_HIGH))

    syn_all = df['simplified_consequence'] == SYNONYMOUS_LABEL
    syn_scores = score[syn_all].dropna()
    syn_low, syn_high = syn_scores.quantile(0.025), syn_scores.quantile(0.975)
    synonymous = syn_all & (score >= syn_low) & (score <= syn_high)

    mean_abn = score[nonsense].mean()
    mean_nrm = score[synonymous].mean()
    print('INFO: abnormal anchor (stop_gained, in-window): n=%d mean=%.6f'
          % (int(nonsense.sum()), mean_abn))
    print('INFO: normal anchor (synonymous central 95%%): n=%d mean=%.6f'
          % (int(synonymous.sum()), mean_nrm))

    fit_scores = pd.concat([score[nonsense], score[synonymous]]).dropna().values.reshape(-1, 1)
    gm = GaussianMixture(
        n_components=2,
        means_init=np.array([mean_abn, mean_nrm]).reshape(-1, 1),
        random_state=RANDOM_STATE,
    ).fit(fit_scores)
    return gm, score


def classify(df, gm, score):
    """Add posteriors, functional_consequence and score_threshold_class for ALL
    variants (excluded ones included -- flagged, not dropped)."""
    scored_mask = score.notna()
    proba = gm.predict_proba(score[scored_mask].values.reshape(-1, 1))

    df['gmm_posterior_abnormal'] = np.nan
    df['gmm_posterior_normal'] = np.nan
    df.loc[scored_mask, 'gmm_posterior_abnormal'] = proba[:, 0]
    df.loc[scored_mask, 'gmm_posterior_normal'] = proba[:, 1]

    fc = pd.Series('indeterminate', index=df.index)
    fc[df['gmm_posterior_abnormal'] >= POSTERIOR_CUTOFF] = 'functionally_abnormal'
    fc[df['gmm_posterior_normal'] >= POSTERIOR_CUTOFF] = 'functionally_normal'
    # Failsafe: nothing at/above the highest normal-called score is abnormal.
    normal_max = score[(fc == 'functionally_normal')].max()
    fc[score >= normal_max] = 'functionally_normal'
    fc[~scored_mask] = ''
    df['functional_consequence'] = fc.values
    return df


def find_threshold(gm, score, component):
    """Score where the posterior for ``component`` (0=abnormal, 1=normal) crosses
    the cutoff, refined on a dense grid."""
    scored = score.dropna().values
    post = gm.predict_proba(scored.reshape(-1, 1))[:, component]
    guess = scored[np.abs(post - POSTERIOR_CUTOFF).argmin()]
    grid = np.linspace(guess - 0.01, guess + 0.01, 10000)
    gpost = gm.predict_proba(grid.reshape(-1, 1))[:, component]
    return float(grid[np.abs(gpost - POSTERIOR_CUTOFF).argmin()])


def add_threshold_class(df, score, thr_abnormal, thr_normal):
    cls = pd.Series('Indeterminate', index=df.index)
    cls[score <= thr_abnormal] = 'Abnormal'
    cls[score >= thr_normal] = 'Normal'
    cls[score.isna()] = ''
    df['score_threshold_class'] = cls.values
    return df


_CLINVAR_DERIVED = []
for _tp in ('clinvar_2018', 'clinvar_2025', 'clinvar_2026'):
    _CLINVAR_DERIVED += ['%s_classification' % _tp, '%s_slim' % _tp,
                         '%s_review_status' % _tp, '%s_stars' % _tp]
# February 2026 export carries extra fields beyond the CVFG timepoints.
_CLINVAR_DERIVED += ['clinvar_2026_date', 'clinvar_2026_condition',
                     'clinvar_2026_oncogenicity', 'clinvar_2026_variation_id']

DERIVED_COLUMNS = [
    'simplified_consequence', 'display_consequence', 'variant_type',
    'gmm_posterior_abnormal', 'gmm_posterior_normal', 'functional_consequence',
    'score_threshold_class', 'standardized_class',
] + _CLINVAR_DERIVED + [
    'gnomad_v4_af', 'gnomad_v4_faf95_max',
    'revel_score', 'alphamissense_pathogenicity', 'alphamissense_class',
    'spliceai_max_delta',
]


def write_dictionary(df):
    rows = [
        ('simplified_consequence', 'Single VEP class collapsed from vep_consequence by severity; exactly reproduces the notebook column.'),
        ('display_consequence', 'Human-readable consequence label (same GMM classes) for the GMM figure.'),
        ('variant_type', 'Coarser Figure 1/2 consequence group (Nonsense/Splice/Intronic/...); folds splice_region into Intronic/Synonymous.'),
        ('gmm_posterior_abnormal', 'GMM posterior probability of the abnormal component (blank if no score).'),
        ('gmm_posterior_normal', 'GMM posterior probability of the normal component (blank if no score).'),
        ('functional_consequence', 'Posterior call: functionally_abnormal / functionally_normal / indeterminate (>=0.95).'),
        ('score_threshold_class', 'Score-threshold call used by figures: Abnormal / Normal / Indeterminate.'),
        ('standardized_class', 'Waters functional_classification standardized: depleted->Abnormal, unchanged->Normal, enriched->Not specified (GMM-vs-Waters figure).'),
        ('clinvar_2018_classification', 'ClinVar clinical significance, Jan 2018 (CVFG), joined by HGVSc.'),
        ('clinvar_2018_slim', 'Collapsed 2018 class: P/LP, B/LB, VUS, Conflicting.'),
        ('clinvar_2018_review_status', 'ClinVar 2018 review status.'),
        ('clinvar_2018_stars', 'ClinVar 2018 review star level.'),
        ('clinvar_2025_classification', 'ClinVar clinical significance, Jan 2025 (CVFG) -- calibration truth set.'),
        ('clinvar_2025_slim', 'Collapsed 2025 class: P/LP, B/LB, VUS, Conflicting.'),
        ('clinvar_2025_review_status', 'ClinVar 2025 review status.'),
        ('clinvar_2025_stars', 'ClinVar 2025 review star level.'),
        ('clinvar_2026_classification', 'ClinVar germline classification, February 2026 export (cited in the paper) -- primary calibration truth set.'),
        ('clinvar_2026_slim', 'Collapsed Feb 2026 class: P/LP, B/LB, VUS, Conflicting.'),
        ('clinvar_2026_review_status', 'ClinVar Feb 2026 germline review status.'),
        ('clinvar_2026_stars', 'ClinVar Feb 2026 review star level (mapped from review status).'),
        ('clinvar_2026_date', 'ClinVar Feb 2026 germline date last evaluated.'),
        ('clinvar_2026_condition', 'ClinVar Feb 2026 reported condition(s).'),
        ('clinvar_2026_oncogenicity', 'ClinVar Feb 2026 oncogenicity classification, if any.'),
        ('clinvar_2026_variation_id', 'ClinVar VariationID (Feb 2026).'),
        ('gnomad_v4_af', 'gnomAD v4.1 allele frequency (CVFG).'),
        ('gnomad_v4_faf95_max', 'gnomAD v4.1 max filtering allele frequency, 95% (CVFG).'),
        ('revel_score', 'REVEL missense pathogenicity score (CVFG).'),
        ('alphamissense_pathogenicity', 'AlphaMissense pathogenicity score (CVFG).'),
        ('alphamissense_class', 'AlphaMissense class: likely_pathogenic / likely_benign / ambiguous (CVFG).'),
        ('spliceai_max_delta', 'SpliceAI maximum delta score (CVFG).'),
    ]
    with open(DICTIONARY_MD, 'w') as fh:
        fh.write('# BAP1 master table -- derived columns\n\n')
        fh.write('Built by `scripts/make_master_table.py` from '
                 '`data/raw/Waters_ExtendedData1.xlsx` (read-only). All original '
                 'Waters columns are preserved unchanged; the columns below are appended.\n\n')
        fh.write('| column | description |\n|---|---|\n')
        for name, desc in rows:
            fh.write('| `%s` | %s |\n' % (name, desc))
    print('INFO: wrote %s' % os.path.relpath(DICTIONARY_MD, REPO_ROOT))


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    df = load_source()
    original_columns = list(df.columns)

    df = annotate_consequence(df)

    gm, score = fit_gmm(df)
    means = gm.means_.ravel()
    sds = np.sqrt(gm.covariances_.ravel())
    weights = gm.weights_.ravel()

    df = classify(df, gm, score)
    thr_abnormal = find_threshold(gm, score, component=0)
    thr_normal = find_threshold(gm, score, component=1)
    df = add_threshold_class(df, score, thr_abnormal, thr_normal)
    print('INFO: abnormal cutoff (score <= %.6f); normal cutoff (score >= %.6f)'
          % (thr_abnormal, thr_normal))
    print('INFO: score_threshold_class counts: %s'
          % df['score_threshold_class'].value_counts().to_dict())

    df = annotate_standardized_class(df)
    df = merge_cvfg(df)
    df = merge_clinvar_feb2026(df)
    df = exclude_false_germline_clinvar(df)

    # ── Write master table: originals first (unchanged), then derived columns ──
    ordered = original_columns + [c for c in DERIVED_COLUMNS if c in df.columns]
    df[ordered].to_csv(MASTER_TSV, sep='\t', index=False)
    print('INFO: wrote %s (%d rows x %d cols)'
          % (os.path.relpath(MASTER_TSV, REPO_ROOT), len(df), len(ordered)))

    # ── Provenance / thresholds ──
    provenance = {
        'source_file': os.path.relpath(SOURCE_XLSX, REPO_ROOT),
        'source_sheet': SOURCE_SHEET,
        'cvfg_annotation_source': os.path.relpath(CVFG_TSV, REPO_ROOT),
        'cvfg_merge_key': 'Ensembl HGVSc, genomic pos+ref+alt fallback',
        'clinvar_2026_source': os.path.relpath(CLINVAR_FEB2026_TXT, REPO_ROOT),
        'clinvar_timepoints': ['2018-01 (CVFG)', '2025-01 (CVFG)', '2026-02 (paper)'],
        'n_variants': int(len(df)),
        'exclusion_rule': ('none; authors pre-filtered pam_flag=Y upstream, so '
                           'ED1 as delivered is the calibration set'),
        'random_state': RANDOM_STATE,
        'posterior_cutoff': POSTERIOR_CUTOFF,
        'sklearn_version': sklearn.__version__,
        'gmm_covariance_type': gm.covariance_type,
        'gmm_tol': gm.tol,
        'gmm_max_iter': gm.max_iter,
        'gmm_reg_covar': gm.reg_covar,
        'gmm_n_iter_converged': int(gm.n_iter_),
        'nonsense_window_hg38': [NONSENSE_POS_LOW, NONSENSE_POS_HIGH],
        'abnormal_cutoff': thr_abnormal,
        'normal_cutoff': thr_normal,
        'component_abnormal': {'mean': float(means[0]), 'sd': float(sds[0]),
                               'weight': float(weights[0])},
        'component_normal': {'mean': float(means[1]), 'sd': float(sds[1]),
                             'weight': float(weights[1])},
    }
    with open(THRESHOLDS_JSON, 'w') as fh:
        json.dump(provenance, fh, indent=2)
    print('INFO: wrote %s' % os.path.relpath(THRESHOLDS_JSON, REPO_ROOT))

    write_dictionary(df)


if __name__ == '__main__':
    main()
