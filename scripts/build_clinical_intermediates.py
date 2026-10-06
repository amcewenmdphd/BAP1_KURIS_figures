#!/usr/bin/env python3
"""
Rebuild the clinical / ClinVar-VCV intermediate tables from the raw input files,
so the clinical figure panels are reproducible end to end (raw ClinVar +
NDD cohort sources -> intermediates -> figure workbook -> figures).

This is a faithful port of the clinical-curation cells of the original
BAP1_KURIS_ScoreThreshold_Figures notebook (cells 7, 9, 11, 13, 15, 17, 19).

Raw inputs (data/raw/)
    BAP1_ClinVar_VCV.tsv        ClinVar variant (VCV) records, XML-flattened
    BAP1_ClinVar_SCV.tsv        ClinVar submission (SCV) records
    DECIPHER_BAP1_NDD.csv       DECIPHER NDD patients
    Kaplanis_et_al_2020_BAP1.tsv  Kaplanis 2020 de novo NDD variants
    Kury_et_al_2022_Table1.csv  Kury 2022 Kury-Isidor cohort (Table 1)

Outputs (data/source/)
    BAP1_ClinVar_SCV_curated.tsv    germline SCVs with study-assigned condition
    BAP1_KURIS_NDD_merged.tsv       NDD/KURIS variants merged across all sources
    BAP1_Cancer_P_LP_merged.tsv     cancer-predisposition P/LP variants
    BAP1_ClinVar_VCV_classified.tsv all classified ClinVar variants (Fig 1b source)
    BAP1_ClinVar_PLP_by_condition.tsv  P/LP germline variants by condition (Fig 1c source)

Usage
    python scripts/make_clinical_tables.py
"""

import os
import re

import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
DATA_IN = os.path.join(REPO_ROOT, 'data', 'raw')
DATA_OUT = os.path.join(REPO_ROOT, 'data', 'source')

PLP_CLASSES = {'Pathogenic', 'Likely pathogenic', 'Pathogenic/Likely pathogenic'}
BLB_CLASSES = {'Benign', 'Likely benign', 'Benign/Likely benign'}


def _s(v):
    """Safe string: '' for NaN/None."""
    return '' if pd.isna(v) else str(v).strip()


# ── ClinVar VCV parse (notebook cell 7) ──

def build_vcv_bap1(df_vcv_all):
    def parse_spdi(spdi):
        if pd.isna(spdi):
            return None, None, None
        parts = str(spdi).split(':')
        if len(parts) == 4:
            return int(parts[1]) + 1, parts[2], parts[3]   # 0-based -> 1-based
        return None, None, None

    v = df_vcv_all[
        (df_vcv_all['GeneSymbol'] == 'BAP1')
        & df_vcv_all['AggGermlineClass'].notna()
        & (df_vcv_all['AggGermlineClass'] != '')
        & (df_vcv_all['AggGermlineClass'] != 'not provided')
    ].copy()

    spdi = v['CanonicalSPDI'].apply(parse_spdi)
    v['spdi_pos'] = pd.to_numeric([x[0] for x in spdi], errors='coerce')
    v['spdi_ref'] = [x[1] for x in spdi]
    v['spdi_alt'] = [x[2] for x in spdi]

    v['_vlen'] = pd.to_numeric(v['GRCh38_variantLength'], errors='coerce')
    v = v[(v['_vlen'] < 50) | v['_vlen'].isna()].copy()
    v = v[v['spdi_pos'].notna() | v['MolecularConsequence'].notna()].copy()

    def map_germline_class(sig):
        if pd.isna(sig) or sig in ('', 'not provided'):
            return None
        sig = str(sig).strip()
        if sig in PLP_CLASSES:
            return 'P/LP'
        if sig in BLB_CLASSES:
            return 'B/LB'
        if sig == 'Uncertain significance':
            return 'VUS'
        if 'onflicting' in sig:
            return 'Conflicting'
        return None

    v['cv2026_class'] = v['AggGermlineClass'].apply(map_germline_class)
    return v


def simplify_consequence_cv(mc):
    if pd.isna(mc):
        return 'Other'
    mc = str(mc).lower()
    if 'frameshift' in mc: return 'Frameshift'
    if 'nonsense' in mc: return 'Nonsense'
    if 'splice acceptor' in mc: return 'Splice'
    if 'splice donor' in mc: return 'Splice'
    if 'missense' in mc and 'initiator' in mc: return 'Start lost'
    if 'missense' in mc: return 'Missense'
    if 'inframe' in mc: return 'Inframe indel'
    if 'intron' in mc: return 'Intronic'
    if 'synonymous' in mc: return 'Synonymous'
    if 'stop lost' in mc: return 'Stop lost'
    if 'utr' in mc: return 'UTR'
    return 'Other'


# ── VCV-level condition curation for Fig 1c (notebook cell 9) ──

# ClinVar ConditionXRefs used for condition assignment. Every cross-reference seen
# on a germline BAP1 submission in the February 2026 ClinVar release is catalogued
# here. This mapping was curated and MANUALLY VERIFIED: each code's definition was
# checked against its official OMIM / MedGen / MONDO / Orphanet / HPO entry, and its
# cancer / Kury-Isidor / not-provided assignment was reviewed by hand. The key is
# the lowercased fragment that identifies the code in the (lowercased) ConditionXRefs
# string, the value is the official ontology term. Matching is substring, so tokens
# are kept distinct. The CANCER_KW / NDD_KW keyword lists below were verified the
# same way against the distinct BAP1 condition strings in this release.
XREF_KURIS = {
    '619762':   'OMIM #619762: Kury-Isidor syndrome; KURIS',
    'cn306819': 'MedGen CN306819: Kury-Isidor syndrome; KURIS (retired ID; replaced by C5676925 / MONDO:0859230)',
    'c5676925': 'MedGen C5676925: Kury-Isidor syndrome; KURIS (current ID; not in the Feb 2026 export)',
}
XREF_CANCER = {
    '614327':       'OMIM #614327: Tumor predisposition syndrome 1; TPDS1',
    '606661':       'OMIM #606661: Melanoma, uveal, susceptibility to, 2; UVM2',
    'c3280492':     'MedGen C3280492: Tumor susceptibility linked to germline BAP1 mutations',
    'c0027672':     'MedGen C0027672: Hereditary cancer-predisposing syndrome',
    'c4707290':     'MedGen C4707290: BAP1 tumor predisposition syndrome',
    'c1333600':     'MedGen C1333600: Hereditary cancer',
    'c0553580':     'MedGen C0553580: Ewing sarcoma',
    'cn204945':     'MedGen CN204945: Uveal melanoma',
    'c0025149':     'MedGen C0025149: Medulloblastoma',
    'c0206716':     'MedGen C0206716: Ganglioglioma',
    '0013692':      'MONDO:0013692: BAP1-related tumor predisposition syndrome',
    'orpha289539':  'Orphanet ORPHA:289539: BAP1-related tumor predisposition syndrome',
    'orphanet:618': 'Orphanet ORPHA:618: Familial melanoma',
    '0007716':      'HP:0007716: Intraocular melanoma',
    '0009592':      'HP:0009592: Astrocytoma',
    'orpha29073':   'Orphanet ORPHA:29073: Multiple myeloma',
}
# Catalogued but deliberately NOT mapped (placeholders / non-BAP1) -> Not provided:
XREF_EXCLUDED = {
    'cn517202': 'MedGen CN517202: not provided (ClinVar placeholder)',
    'c3661900': 'MedGen C3661900: not provided (ClinVar placeholder)',
    'cn169374': 'MedGen CN169374: not specified (ClinVar placeholder)',
    'cn235283': 'MedGen CN235283: ClinVar concept (nonspecific)',
    'c2677590': 'MedGen C2677590: RFT1-congenital disorder of glycosylation (non-BAP1)',
}
KURIS_XREFS = set(XREF_KURIS)
CANCER_XREFS = set(XREF_CANCER)
STRIP_CONDITIONS = {'Neoplasm', 'Neoplasms'}
LIST_EVERYTHING = {
    'Melanoma, uveal, susceptibility to, 2',
    'Kury-Isidor syndrome',
    'BAP1-related tumor predisposition syndrome',
}
CONDITION_PRIORITY = ['Kury-Isidor syndrome', 'Cancer predisposition', 'Not provided']

# Unified keyword lists, shared by the VCV-level (Fig 1c) and SCV-level curation so
# the two systems cannot diverge. 'bap1' is deliberately NOT a cancer keyword: a bare
# "BAP1-related condition/disorder" carries no directional signal and is left as
# Not provided. Anything that is neither NDD nor cancer folds into Not provided.
CANCER_KW = [
    'tumor', 'tumour', 'cancer', 'carcinoma', 'melanoma', 'mesothelioma',
    'neoplasm', 'neoplastic', 'predispos', 'astrocytoma', 'medulloblastoma',
    'ganglioglioma', 'ewing', 'renal', 'myeloma', 'uveal', 'malignant',
    'oncogenic', 'basal cell', 'leukemia', 'leukaemia', 'sarcoma', 'thymoma',
    'lymphoma', 'blastoma', 'glioma',
]
NDD_KW = [
    'neurodevelopmental', 'developmental delay', 'intellectual disab',
    'kury-isidor', 'kuris', 'kury', 'autism', 'seizure', 'epilep',
    'microcephaly', 'hypotonia',
]
# Generic BAP1 condition names that carry no directional signal: they match neither
# NDD nor cancer ('bap1' is deliberately not a cancer keyword), so they fall through
# to Not provided. Catalogued here as an explicit double-check.
AMBIGUOUS_CONDITIONS = [
    'BAP1-related condition',
    'BAP1-related disorder',
]


def build_cv_classified_and_plp(df_vcv_bap1, df_scv_all):
    onco_vcvs = set(df_scv_all[df_scv_all['InterpretationType'] == 'oncogenicity']['VCV'].unique())

    # False germline: all germline SCVs marked somatic origin
    germ_per_vcv = df_scv_all[df_scv_all['InterpretationType'] == 'germline'].groupby('VCV')['SCV'].count()
    som_per_vcv = df_scv_all[(df_scv_all['InterpretationType'] == 'germline')
                             & (df_scv_all['Origin'] == 'somatic')].groupby('VCV')['SCV'].count()
    false_germline = {vid for vid in som_per_vcv.index
                      if som_per_vcv[vid] == germ_per_vcv.get(vid, 0)}

    def classify_xref(xref_str):
        if pd.isna(xref_str) or xref_str == '':
            return None
        xr = xref_str.lower()
        has_kuris = any(x in xr for x in KURIS_XREFS)
        has_cancer = any(x in xr for x in CANCER_XREFS)
        if has_kuris and has_cancer:
            return None
        if has_kuris:
            return 'Kury-Isidor syndrome'
        if has_cancer:
            return 'Cancer predisposition'
        return None

    def vcv_condition_from_scvs(vcv_id):
        germ = df_scv_all[(df_scv_all['VCV'] == vcv_id)
                          & (df_scv_all['InterpretationType'] == 'germline')]
        if germ.empty:
            return 'Not provided'
        conditions = set()
        for _, s in germ.iterrows():
            cond_text = str(s.get('Condition', '')).lower()
            xref_cls = classify_xref(str(s.get('ConditionXRefs', '')))
            if xref_cls:
                conditions.add(xref_cls)
                continue
            if any(kw in cond_text for kw in NDD_KW):
                conditions.add('Kury-Isidor syndrome')
            elif any(kw in cond_text for kw in CANCER_KW):
                conditions.add('Cancer predisposition')
        if 'Kury-Isidor syndrome' in conditions:
            return 'Kury-Isidor syndrome'
        if 'Cancer predisposition' in conditions:
            return 'Cancer predisposition'
        return 'Not provided'

    def simplify_single_condition(cond):
        c = cond.lower().strip()
        if any(k in c for k in NDD_KW):
            return 'Kury-Isidor syndrome'
        if any(k in c for k in CANCER_KW):
            return 'Cancer predisposition'
        return 'Not provided'          # ambiguous / blank folds into Not provided

    def assign_condition(row):
        acc = row['VCV']
        if acc in false_germline:
            return 'False germline (somatic)'
        scv_cond = vcv_condition_from_scvs(acc)
        if scv_cond != 'Not provided':
            return scv_cond
        raw_parts = [p.strip() for p in str(row['AggConditions']).split('|')]
        if acc in onco_vcvs:
            raw_parts = [p for p in raw_parts if p not in STRIP_CONDITIONS]
        has_list_everything = LIST_EVERYTHING.issubset(set(raw_parts))
        cleaned = [p for p in raw_parts if not (has_list_everything and p in LIST_EVERYTHING)]
        if not cleaned:
            return 'Not provided'
        simplified = {simplify_single_condition(c) for c in cleaned}
        simplified.discard('Not provided')
        if not simplified:
            return 'Not provided'
        for cat in CONDITION_PRIORITY:
            if cat in simplified:
                return cat
        return 'Not provided'

    df_plp_cv = df_vcv_bap1[df_vcv_bap1['AggGermlineClass'].isin(PLP_CLASSES)].copy()
    df_plp_cv['is_false_germline'] = df_plp_cv['VCV'].isin(false_germline)
    df_plp_cv['simple_consequence'] = df_plp_cv['MolecularConsequence'].apply(simplify_consequence_cv)
    df_plp_cv['curated_condition'] = df_plp_cv.apply(assign_condition, axis=1)
    df_plp_germline = df_plp_cv[df_plp_cv['curated_condition'] != 'False germline (somatic)'].copy()

    # Exclude the same false-germline (somatic-only-origin) VCVs used for Fig 1c's
    # condition breakdown, so Fig 1b's per-class counts agree with it -- these are
    # variants ClinVar's aggregate germline field calls P/LP, but every germline-
    # track SCV for them is flagged Origin=somatic (found in tumor sequencing, e.g.
    # VCV000800323/multiple myeloma, VCV000633598/uveal melanoma), so they are not
    # true germline pathogenic calls.
    df_cv_classified = df_vcv_bap1[df_vcv_bap1['cv2026_class'].notna()
                                    & ~df_vcv_bap1['VCV'].isin(false_germline)].copy()
    df_cv_classified['simple_consequence'] = df_cv_classified['MolecularConsequence'].apply(simplify_consequence_cv)
    return df_cv_classified, df_plp_germline


# ── SCV-level condition curation (notebook cell 11) ──

# ── SCVs flagged for manual review ──────────────────────────────────────────
# The single source of truth for every hand-reviewed germline SCV. For each we
# hard-code (a) the FLAG -- why the automated rule cannot be trusted for this
# record -- and (b) the OUTCOME of the manual review: the StudyCondition assigned
# after reading the submitter's coded fields, reported clinical features
# (ObservedHPO), and free-text comment by hand.
#
# These are the ONLY records where a human overrides the automated assignment.
# A systematic scan of all germline SCVs (Feb 2026 ClinVar release) confirms the
# set is complete: every comment-contradicts-xref conflict is flagged here, and no
# unflagged P/LP record carries recoverable ObservedHPO. Flag categories:
#   'comment contradicts xref' -- coded xref points one way, comment the other
#                                 (these are the only flags that flip the rule's output)
#   'no coded signal'          -- condition/xref uninformative; classified from ObservedHPO
#   'confirmatory'             -- coded signal agrees with the review; changes nothing
_MANUAL_REVIEW = {
    'SCV004801702': dict(
        flag='no coded signal',
        outcome='NDD / Kury-Isidor syndrome',
        note="condition 'see cases', no xref; ObservedHPO (global developmental delay, "
             "microcephaly, poor speech, infantile hypotonia) compatible with KURIS",
    ),
    'SCV002523020': dict(
        flag='no coded signal',
        outcome='Cancer predisposition',
        note="condition 'see cases', no xref; ObservedHPO (papillary thyroid carcinoma, "
             "neoplasm, hypothyroidism) compatible with TPDS",
    ),
    'SCV005399745': dict(
        flag='confirmatory',
        outcome='Cancer predisposition',
        note="truncating variant; comment (NMD/loss of protein, TPDS1 MIM#614327) "
             "consistent with cancer; agrees with xref OMIM:614327",
    ),
    'SCV005399747': dict(
        flag='confirmatory',
        outcome='NDD / Kury-Isidor syndrome',
        note="xref OMIM:619762 (KURIS); agrees with rule",
    ),
    'SCV004911795': dict(
        flag='comment contradicts xref',
        outcome='NDD / Kury-Isidor syndrome',
        note="entry written for KURIS not TPDS despite cancer xref (MedGen:C0027672): "
             "'...for BAP1-related neurodevelopmental disorder; however, its clinical "
             "significance for BAP1-related tumor predisposition syndrome is uncertain'",
    ),
    'SCV004053060': dict(
        flag='comment contradicts xref',
        outcome='NDD / Kury-Isidor syndrome',
        note="comment describes de novo BAP1-related neurodevelopmental disorder despite "
             "cancer xref (MedGen:C0027672)",
    ),
    'SCV004403353': dict(
        flag='comment contradicts xref',
        outcome='NDD / Kury-Isidor syndrome',
        note="cites BAP1 neurodevelopmental functional paper (PMID:35051358) despite "
             "cancer xref (MedGen:C3280492)",
    ),
}
# Assignment applied by curate() (SCV -> reviewed StudyCondition).
_MANUAL_OVERRIDES = {scv: rec['outcome'] for scv, rec in _MANUAL_REVIEW.items()}
def curate_scv_conditions(df_vcv_all, df_scv_all):
    vcv_germ = df_vcv_all[
        (df_vcv_all['GeneSymbol'] == 'BAP1')
        & df_vcv_all['AggGermlineClass'].notna()
        & (df_vcv_all['AggGermlineClass'] != '')
        & (df_vcv_all['AggGermlineClass'] != 'not provided')]
    df_scv_germ = df_scv_all[
        df_scv_all['VCV'].isin(set(vcv_germ['VCV']))
        & (df_scv_all['InterpretationType'] == 'germline')].copy()

    def is_multi_phenotype(xrefs):
        return '619762' in xrefs and ('614327' in xrefs or '606661' in xrefs)

    def curate(row):
        scv_id = str(row['SCV'])
        cond = str(row['Condition']).lower()
        xrefs = str(row['ConditionXRefs']).lower()
        comment = str(row['Comment']).lower()
        if scv_id in _MANUAL_OVERRIDES:
            return _MANUAL_OVERRIDES[scv_id], 'manual review'
        if is_multi_phenotype(xrefs):
            if any(kw in comment for kw in NDD_KW):
                return 'NDD / Kury-Isidor syndrome', 'comment keyword'
            if any(kw in comment for kw in CANCER_KW):
                return 'Cancer predisposition', 'comment keyword'
            return 'Not provided', 'multi-phenotype (nonspecific)'
        if '619762' in xrefs or any(kw in cond for kw in NDD_KW):
            return 'NDD / Kury-Isidor syndrome', 'condition text/xref'
        if any(kw in cond for kw in CANCER_KW):
            return 'Cancer predisposition', 'condition text/xref'
        if any(x in xrefs for x in CANCER_XREFS):
            return 'Cancer predisposition', 'condition text/xref'
        # No directional signal in the condition (e.g. a bare "BAP1-related condition")
        # or a blank condition: fall back to the submission comment, else Not provided.
        if any(kw in comment for kw in NDD_KW):
            return 'NDD / Kury-Isidor syndrome', 'comment keyword'
        if any(kw in comment for kw in CANCER_KW):
            return 'Cancer predisposition', 'comment keyword'
        return 'Not provided', 'no directional signal'

    res = df_scv_germ.apply(curate, axis=1, result_type='expand')
    df_scv_germ['StudyCondition'] = res[0]
    df_scv_germ['ConditionMethod'] = res[1]
    # Annotate every SCV with its manual-review flag and outcome. Flagged records
    # carry the reason they need review and the reviewed result; all others are 'OK'.
    scv_str = df_scv_germ['SCV'].astype(str)
    df_scv_germ['ManualReviewFlag'] = scv_str.map(
        {s: r['flag'] for s, r in _MANUAL_REVIEW.items()}).fillna('OK')
    df_scv_germ['ManualReviewOutcome'] = scv_str.map(
        {s: '%s -- %s' % (r['outcome'], r['note']) for s, r in _MANUAL_REVIEW.items()}).fillna('OK')

    def simplify_class(c):
        c = str(c).lower()
        if c in ('pathogenic', 'likely pathogenic'):
            return 'P/LP'
        if c in ('benign', 'likely benign'):
            return 'B/LB'
        if c == 'uncertain significance':
            return 'VUS'
        return 'Other'

    df_scv_germ['SimpleClass'] = df_scv_germ['Classification'].apply(simplify_class)
    return df_scv_germ


# ── Submitter lab lookup + SCV line formatting (notebook cell 13) ──

_ORG_NAMES = {
    '3': 'OMIM', '1006': 'BCM: BMGL', '1238': 'U of Chicago Genetics',
    '25969': 'ARUP Laboratories', '26957': 'GeneDx', '61756': 'Ambry Genetics',
    '167595': 'Revvity Omics', '239772': 'PreventionGenetics', '274978': 'GDL-UMCU',
    '500031': 'Invitae', '500035': 'Mendelics', '500057': 'Sema4', '500060': 'EGL',
    '500068': 'Mayo Clinic', '500104': 'VCGS', '500105': 'Fulgent Genetics',
    '500110': 'Quest Diagnostics', '504895': 'Illumina', '505150': 'IoHGE',
    '505336': 'IHG-UKRWTH', '505613': 'MedGenet', '505691': 'GM4113, RH',
    '505849': 'Color', '505870': 'CHGT', '505925': 'Radboudumc_MUMC+',
    '505999': 'UDN', '506086': 'UKL', '506152': 'PLM-SHS', '506185': 'GenomeConnect',
    '506213': 'Lab (506213)', '506273': 'True Health Diagnostics', '506382': 'DL-UMCG',
    '506385': 'IMGAG', '506395': 'LDGA', '506453': 'VUMC GD', '506497': 'CGDx-EMC',
    '506672': 'SJCRH-MP', '507240': 'Myriad Genetics', '507439': 'IHG',
    '507558': 'NYGC', '509060': 'Lab (509060)', '509089': 'KCCC/NGS Lab Kuwait',
    '509169': 'Breakthrough Genomics',
}
_CLS_SHORT = {'Pathogenic': 'P', 'Likely pathogenic': 'LP', 'Uncertain significance': 'VUS',
              'Likely benign': 'LB', 'Benign': 'B', 'not provided': 'not provided',
              'Likely Benign': 'LB', 'likely benign': 'LB', 'benign': 'B'}
_COND_SHORT = {'NDD / Kury-Isidor syndrome': 'NDD/KURIS',
               'Cancer predisposition': 'Cancer', 'Not provided': 'Not specified'}
_COND_ORDER = {'NDD / Kury-Isidor syndrome': 0, 'Cancer predisposition': 1, 'Not provided': 2}


def _org_names_with_data(df_scv_germ):
    names = dict(_ORG_NAMES)
    for _, r in df_scv_germ[['SubmitterOrgID', 'OrgAbbreviation']].dropna(subset=['OrgAbbreviation']).iterrows():
        oid, abbr = str(r['SubmitterOrgID']), str(r['OrgAbbreviation'])
        if oid not in names and abbr and abbr != 'nan':
            names[oid] = abbr
    return names


def _make_get_lab(names):
    def get_lab(row):
        abbr = str(row.get('OrgAbbreviation', ''))
        if abbr and abbr != 'nan':
            return abbr
        return names.get(str(row.get('SubmitterOrgID', '')), f"Org#{row.get('SubmitterOrgID', '')}")
    return get_lab


# ── NDD/KURIS merged table (notebook cells 13/15/17) ──

def build_ndd_merged(df_vcv_all, df_scv_germ, df_decipher, df_kaplanis, get_lab):
    kury_raw = pd.read_csv(os.path.join(DATA_IN, 'Kury_et_al_2022_Table1.csv'))
    kury_hgvs_c = list(kury_raw.loc[kury_raw['Feature'] == 'HGVS_c'].iloc[0][1:].dropna().values)
    kury_hgvs_p = list(kury_raw.loc[kury_raw['Feature'] == 'HGVS_p'].iloc[0][1:].dropna().values)
    kury_families = list(kury_raw.columns[1:])
    kury_variants = {}
    for i, hc in enumerate(kury_hgvs_c):
        fam = kury_families[i] if i < len(kury_families) else f'F{i+1}'
        hp = kury_hgvs_p[i] if i < len(kury_hgvs_p) else ''
        kury_variants.setdefault(hc, {'hgvs_c': hc, 'hgvs_p': hp, 'families': []})['families'].append(fam)

    dec_variants = {}
    for _, r in df_decipher.iterrows():
        hc = _s(r.get('cDNA_change', r.get('HGVS_c', '')))
        if not hc:
            continue
        detail = f"DECIPHER {_s(r.get('Patient', ''))}"
        if _s(r.get('Pathogenicity', '')):
            detail += f": {_s(r.get('Pathogenicity', ''))}"
        if _s(r.get('Phenotypes', '')):
            detail += f" ({_s(r.get('Phenotypes', ''))})"
        dec_variants.setdefault(hc, {'hgvs_c': hc,
                                     'hgvs_p': _s(r.get('Protein_change', r.get('HGVS_p', ''))),
                                     'patients': []})['patients'].append(detail)

    kap_variants = {}
    for _, r in df_kaplanis.iterrows():
        hc = _s(r.get('HGVS_c', ''))
        if not hc:
            continue
        kap_variants.setdefault(hc, {'hgvs_c': hc, 'hgvs_p': _s(r.get('HGVS_p', '')),
                                     'ids': []})['ids'].append(f"{_s(r.get('id', ''))} ({_s(r.get('study', ''))})")

    cv_ndd = df_scv_germ[df_scv_germ['StudyCondition'] == 'NDD / Kury-Isidor syndrome'].copy()
    cv_ndd_vcvs = sorted(cv_ndd['VCV'].unique())
    all_for_ndd = df_scv_germ[df_scv_germ['VCV'].isin(set(cv_ndd_vcvs))].copy()

    vcv_hgvs = {}
    for _, r in df_vcv_all[df_vcv_all['VCV'].isin(set(cv_ndd_vcvs))].iterrows():
        m = re.search(r'(c\.\S+)', _s(r.get('VariantName', '')))
        if m:
            vcv_hgvs[r['VCV']] = m.group(1)

    all_hgvs = set(kury_variants) | set(dec_variants) | set(kap_variants) | set(vcv_hgvs.values())
    all_hgvs.discard('')

    rows = []
    for hc in sorted(all_hgvs):
        row = {'HGVS_c': hc}
        hp = (kury_variants.get(hc, {}).get('hgvs_p', '')
              or dec_variants.get(hc, {}).get('hgvs_p', '')
              or kap_variants.get(hc, {}).get('hgvs_p', ''))
        row['Protein'] = hp
        vcv_id = next((v for v, h in vcv_hgvs.items() if h == hc), '')
        row['VCV'] = vcv_id
        if vcv_id:
            vr = df_vcv_all[df_vcv_all['VCV'] == vcv_id].iloc[0]
            cons = _s(vr.get('MolecularConsequence', ''))
            row['Consequence'] = cons.split(';')[0].split('(')[0].strip() if cons else ''
            _chr, _start = _s(vr.get('GRCh38_Chr', '')), _s(vr.get('GRCh38_start', ''))
            row['GRCh38_Position'] = f'chr{_chr}:{_start}' if _chr and _start else ''
            row['dbSNP'] = _s(vr.get('dbSNP_rsID', ''))
            row['AggGermlineClass'] = _s(vr.get('AggGermlineClass', ''))
            row['StarRating'] = _s(vr.get('StarRating', ''))
        else:
            row.update({'Consequence': '', 'GRCh38_Position': '', 'dbSNP': '',
                        'AggGermlineClass': '', 'StarRating': ''})

        kuris_parts = []
        if hc in kury_variants:
            fams = ', '.join(kury_variants[hc]['families'])
            details = []
            for fam in kury_variants[hc]['families']:
                if fam in kury_raw.columns:
                    for _, kr in kury_raw.iterrows():
                        feat, val = _s(kr['Feature']), _s(kr[fam])
                        if feat in ('Sex', 'Age_at_last_investigation', 'Developmental_delay_or_ID',
                                    'Seizures', 'Behavioral_anomalies', 'Cardiac_anomalies',
                                    'Eye_anomalies', 'Hands', 'Feet', 'Facial_dysmorphism',
                                    'Other_recurrent_signs', 'Growth_failure') \
                                and val and val not in ('', '-', 'N/A', 'nan', 'not applicable'):
                            details.append(f'{feat}: {val}')
            detail = f'Kury et al. 2022 ({fams})'
            if details:
                detail += ' [' + '; '.join(details) + ']'
            kuris_parts.append(detail)
        if hc in dec_variants:
            kuris_parts.extend(dec_variants[hc]['patients'])
        if hc in kap_variants:
            kuris_parts.append(f"Kaplanis et al. 2020 ({', '.join(kap_variants[hc]['ids'])})")
        if vcv_id:
            for _, sr in cv_ndd[cv_ndd['VCV'] == vcv_id].iterrows():
                cls = _CLS_SHORT.get(_s(sr['Classification']), _s(sr['Classification']))
                line = f"{_s(sr['SCV'])}: {cls} (NDD/KURIS), {get_lab(sr)}"
                if _s(sr.get('DateLastEvaluated', '')):
                    line += f", {_s(sr.get('DateLastEvaluated', ''))}"
                comment = _s(sr.get('Comment', ''))
                if comment:
                    line += f" [{comment[:147] + '...' if len(comment) > 150 else comment}]"
                kuris_parts.append(line)
        row['KURIS_Sources'] = '\n'.join(kuris_parts)

        other_parts = []
        if vcv_id:
            for _, sr in all_for_ndd[(all_for_ndd['VCV'] == vcv_id)
                                     & (all_for_ndd['StudyCondition'] != 'NDD / Kury-Isidor syndrome')].iterrows():
                cls = _CLS_SHORT.get(_s(sr['Classification']), _s(sr['Classification']))
                cond = _COND_SHORT.get(_s(sr['StudyCondition']), _s(sr['StudyCondition']))
                other_parts.append(f"{_s(sr['SCV'])}: {cls} ({cond}), {get_lab(sr)}")
        row['Other_ClinVar_SCVs'] = '\n'.join(other_parts)

        agg = row['AggGermlineClass'].lower()
        has_ndd_plp = vcv_id and any(
            _s(r2['Classification']).lower() in ('pathogenic', 'likely pathogenic')
            for _, r2 in cv_ndd[cv_ndd['VCV'] == vcv_id].iterrows())
        has_dec_plp = hc in dec_variants and any(
            'Pathogenic' in p or 'Likely pathogenic' in p for p in dec_variants[hc]['patients'])
        in_kap = hc in kap_variants
        in_dec = hc in dec_variants
        in_cv_ndd = vcv_id != '' and len(cv_ndd[cv_ndd['VCV'] == vcv_id]) > 0
        if 'benign' in agg or 'likely benign' in agg:
            row['Include'], row['InclusionNotes'] = 'FALSE', 'Benign/LB aggregate classification in ClinVar'
        elif hc in kury_variants:
            row['Include'], row['InclusionNotes'] = 'TRUE', 'Present in Kury et al. 2022'
        elif has_ndd_plp:
            row['Include'], row['InclusionNotes'] = 'TRUE', 'P/LP for NDD/KURIS in ClinVar'
        elif has_dec_plp:
            row['Include'], row['InclusionNotes'] = 'TRUE', 'P/LP in DECIPHER'
        elif in_kap and not in_dec and not in_cv_ndd:
            row['Include'], row['InclusionNotes'] = 'FALSE', 'Kaplanis only (no other NDD source)'
        else:
            row['Include'], row['InclusionNotes'] = 'FALSE', 'No P/LP assertion for NDD/KURIS'
        rows.append(row)

    return pd.DataFrame(rows)


# ── Cancer P/LP merged table (notebook cell 19) ──

def build_cancer_merged(df_vcv_all, df_scv_germ, get_lab):
    vcv_plp = df_vcv_all[(df_vcv_all['GeneSymbol'] == 'BAP1')
                         & (df_vcv_all['AggGermlineClass'].isin(PLP_CLASSES))].copy()
    rows = []
    for _, vr in vcv_plp.iterrows():
        vcv_id = vr['VCV']
        m = re.search(r'(c\.\S+)', _s(vr.get('VariantName', '')))
        cons = _s(vr.get('MolecularConsequence', ''))
        _chr, _start = _s(vr.get('GRCh38_Chr', '')), _s(vr.get('GRCh38_start', ''))
        scvs = df_scv_germ[df_scv_germ['VCV'] == vcv_id]
        cancer_lines = []
        for _, sr in scvs[scvs['StudyCondition'] == 'Cancer predisposition'].iterrows():
            cls = _CLS_SHORT.get(_s(sr['Classification']), _s(sr['Classification']))
            line = f"{_s(sr['SCV'])}: {cls} (Cancer), {get_lab(sr)}"
            if _s(sr.get('DateLastEvaluated', '')):
                line += f", {_s(sr.get('DateLastEvaluated', ''))}"
            cancer_lines.append(line)
        other_lines = []
        for _, sr in scvs[scvs['StudyCondition'] != 'Cancer predisposition'].iterrows():
            cls = _CLS_SHORT.get(_s(sr['Classification']), _s(sr['Classification']))
            cond = _COND_SHORT.get(_s(sr['StudyCondition']), _s(sr['StudyCondition']))
            other_lines.append(f"{_s(sr['SCV'])}: {cls} ({cond}), {get_lab(sr)}")
        rows.append({
            'VCV': vcv_id, 'HGVS_c': m.group(1) if m else '',
            'Protein': _s(vr.get('ProteinChange', '')),
            'Consequence': cons.split(';')[0].split('(')[0].strip() if cons else '',
            'GRCh38_Position': f'chr{_chr}:{_start}' if _chr and _start else '',
            'dbSNP': _s(vr.get('dbSNP_rsID', '')),
            'AggGermlineClass': _s(vr.get('AggGermlineClass', '')),
            'StarRating': _s(vr.get('StarRating', '')),
            'Cancer_SCVs': '\n'.join(cancer_lines),
            'Other_ClinVar_SCVs': '\n'.join(other_lines),
        })
    return pd.DataFrame(rows)


def main():
    df_vcv_all = pd.read_csv(os.path.join(DATA_IN, 'BAP1_ClinVar_VCV.tsv'), sep='\t', dtype=str)
    df_scv_all = pd.read_csv(os.path.join(DATA_IN, 'BAP1_ClinVar_SCV.tsv'), sep='\t', dtype=str)
    df_decipher = pd.read_csv(os.path.join(DATA_IN, 'DECIPHER_BAP1_NDD.csv'))
    df_kaplanis = pd.read_csv(os.path.join(DATA_IN, 'Kaplanis_et_al_2020_BAP1.tsv'), sep='\t')
    print('INFO: VCV %d rows, SCV %d rows' % (len(df_vcv_all), len(df_scv_all)))

    df_vcv_bap1 = build_vcv_bap1(df_vcv_all)
    df_cv_classified, df_plp_germline = build_cv_classified_and_plp(df_vcv_bap1, df_scv_all)
    df_scv_germ = curate_scv_conditions(df_vcv_all, df_scv_all)
    get_lab = _make_get_lab(_org_names_with_data(df_scv_germ))

    df_ndd_merged = build_ndd_merged(df_vcv_all, df_scv_germ, df_decipher, df_kaplanis, get_lab)
    df_cancer = build_cancer_merged(df_vcv_all, df_scv_germ, get_lab)

    # Persist intermediates (drop private helper columns). The curated SCV file is
    # restricted to missense SCVs to keep it short — the manuscript's variant-level
    # analysis is missense-focused, and the NDD/Cancer merges above already consumed
    # the full in-memory table, so this filter does not affect any downstream output.
    vcv_cons = df_vcv_all.set_index('VCV')['MolecularConsequence']
    df_scv_out = df_scv_germ.copy()
    df_scv_out.insert(df_scv_out.columns.get_loc('VariationID') + 1,
                      'MolecularConsequence',
                      df_scv_out['VCV'].map(vcv_cons).fillna(''))
    df_scv_out = df_scv_out[
        df_scv_out['MolecularConsequence'].str.contains('missense', case=False, na=False)]
    scv_cols = [c for c in df_scv_out.columns if not c.startswith('_')]
    _write(df_scv_out[scv_cols], 'BAP1_ClinVar_SCV_curated.tsv')
    _write(df_ndd_merged, 'BAP1_KURIS_NDD_merged.tsv')
    _write(df_cancer, 'BAP1_Cancer_P_LP_merged.tsv')

    # Fig 1b / 1c sources (VCV-level).
    fig1b = df_cv_classified[['VCV', 'ProteinChange', 'MolecularConsequence',
                              'simple_consequence', 'cv2026_class', 'StarRating',
                              'AggGermlineClass']]
    _write(fig1b, 'BAP1_ClinVar_VCV_classified.tsv')
    fig1c = df_plp_germline[['VCV', 'ProteinChange', 'MolecularConsequence',
                             'simple_consequence', 'curated_condition', 'StarRating']]
    _write(fig1c, 'BAP1_ClinVar_PLP_by_condition.tsv')


def _write(df, name):
    path = os.path.join(DATA_OUT, name)
    df.to_csv(path, sep='\t', index=False)
    print('INFO: wrote %s (%d rows)' % (os.path.relpath(path, REPO_ROOT), len(df)))


if __name__ == '__main__':
    main()
