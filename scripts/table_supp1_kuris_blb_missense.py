#!/usr/bin/env python3
"""
Supplementary Table 1 — KURIS/NDD and B/LB missense variants.

Every KURIS/NDD variant (from the clinical data table, excl. the UDN proband,
N229K — shown in Figure 1 and Supplementary Table 2 instead) and every ClinVar
B/LB missense variant (from the master data table), each labelled
with the gene (BAP1; HGNC:950) and annotated with its SGE functional score, GMM
functional class, BAP1 protein domain, and ClinVar aggregate germline
classification / accession / star rating. KURIS rows additionally carry the
curated ClinVar-lab / Kury / DECIPHER evidence from the clinical data table.
HGVS g./c./p. are given in full against the reference sequences (GRCh38
NC_000003.12; RefSeq NM_004656.4 / NP_004647.1 — MANE Select, CCDS2853.1-
matched to the Ensembl-canonical ENST00000460680.6 / ENSP00000417132.1 the
master table's HGVS_c/HGVS_p are computed against, so cDNA/protein numbering
is unchanged).

Fully recreatable from the source tables:
  - SGE score + GMM class  <- master data table (by HGVS c)
  - Domain                 <- UniProt feature table (data/raw/BAP1_features.json)
  - KURIS ClinVar / Kury / DECIPHER columns  <- clinical data table
  - B/LB ClinVar columns   <- master data table (Feb 2026 ClinVar snapshot)
  - HGVS g. (GRCh38)       <- master data table, by HGVS c (both KURIS and B/LB;
                              the clinical table's own hg38_pos is 0-based/SPDI,
                              not the 1-based VCF/HGVS position ClinVar reports,
                              so it is not used here)

Output: tables/SuppTable1_KURIS_NDD_and_BLB_Missense.{xlsx,tsv}
"""

import html
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bap1figs import config, tables, tablepdf, supplement

GRCH38_CHROM = "NC_000003.12"
NM_ACCESSION = "NM_004656.4"
NP_ACCESSION = "NP_004647.1"

GENE_HGNC = "BAP1; HGNC:950"      # constant; "BAP1" is italicized when rendered (rich_cols)

COLUMNS = [
    "Gene; HGNC ID", "Category", "HGVS g", "HGVS p", "HGVS c",
    "SGE Functional Score", "Domain", "Functional Classification",
    "ClinVar Classification (Aggregate Germline)", "ClinVar Accession",
    "ClinVar Stars", "ClinVar Lab Entries Related to KURIS", "Küry et al",
    "DECIPHER",
]


def _hgvs_g(pos, ref, alt):
    if pd.isna(pos) or pd.isna(ref) or pd.isna(alt):
        return pd.NA
    return f"{GRCH38_CHROM}:g.{int(pos)}{ref}>{alt}"


def _hgvs_c(short):
    return pd.NA if pd.isna(short) else f"{NM_ACCESSION}:{short}"


def _hgvs_p(short):
    return pd.NA if pd.isna(short) else f"{NP_ACCESSION}:{short}"


def _vcv(variation_id):
    """ClinVar VariationID -> zero-padded VCV accession (matches clinical table)."""
    if pd.isna(variation_id):
        return pd.NA
    return f"VCV{int(variation_id):09d}"


def _stars(x):
    if pd.isna(x):
        return pd.NA
    return int(x)


def _truncated_scv_ids():
    """SCV ids whose ClinVar Comment field is at (or past) ClinVar's hard
    1000-character export cap — i.e. genuinely cut off mid-word server-side,
    not something our formatting broke. Checked against the raw SCV table
    (data/raw/BAP1_ClinVar_SCV.tsv) so the marker below is exact, not a
    punctuation guess (which would misfire on legitimate comma-separated
    lists like DECIPHER phenotype terms that also end mid-word)."""
    raw = pd.read_csv(config.RAW / "BAP1_ClinVar_SCV.tsv", sep="\t", low_memory=False)
    lens = raw["Comment"].astype("string").str.len()
    return set(raw.loc[lens >= 1000, "SCV"])


# SCVs whose ClinVar_KURIS_SCVs entry in the source CSV has been hand-abridged
# into a complete, deliberately shortened summary (not the verbatim, server-cut
# comment) -- so they must NOT be flagged by _mark_truncated even though their
# original raw ClinVar comment was long enough to hit the 1000-char export cap.
MANUALLY_ABRIDGED_SCV = {"SCV004911795", "SCV005399747"}

_TRUNCATED_SCV = _truncated_scv_ids() - MANUALLY_ABRIDGED_SCV


def _mark_truncated(s):
    """If the last SCV entry quoted in s is one of ClinVar's own
    server-truncated comments, replace the ']' our own '[...]' wrapper adds
    with '…]' so the cutoff reads as a cutoff instead of a silent, confusing
    stop mid-word."""
    if not isinstance(s, str) or not s.endswith("]"):
        return s
    ids = re.findall(r"SCV\d+", s)
    if ids and ids[-1] in _TRUNCATED_SCV:
        return s[:-1].rstrip() + "…]"
    return s


def _clean(series):
    """Tidy ClinVar free text: decode HTML entities (&uuml; -> ü, &amp; -> &, …),
    drop the encoded carriage returns (_x000D_), remove any dangling incomplete
    entity left by ClinVar's 1000-char comment truncation (e.g. '(K&u'), mark
    the truncation itself with an ellipsis (see _mark_truncated), and collapse
    whitespace."""
    return (series.astype("string")
            .str.replace("_x000D_", " ", regex=False)
            .map(lambda s: html.unescape(s) if isinstance(s, str) else s)
            .str.replace(r"&[A-Za-z]{1,7};?(?=[)\]\s]|$)", "", regex=True)  # leftover broken entity
            .str.replace(r"\s+", " ", regex=True).str.strip()
            .map(_mark_truncated))


def build_kuris(m):
    """KURIS/NDD rows from the clinical data table, excluding the UDN proband
    (N229K — it has its own row in Figure 1 and Supplementary Table 2; this
    table is the calibration truth set, not a variant inventory). Scores + GMM
    class + domain joined from the master table by HGVS c (so they match
    Figure 2)."""
    clin = tables.load_clinical()
    clin = clin[clin["Protein_Change"] != tables.N229K_PROTEIN].reset_index(drop=True)
    mm = m.drop_duplicates(subset="HGVS_c").set_index("HGVS_c")
    sub = mm.reindex(clin["HGVS_c"])
    # Genomic position from the master table (1-based, matches ClinVar's GRCh38
    # VCF coordinates), not the clinical table's own hg38_pos (0-based/SPDI).
    hgvs_g = [_hgvs_g(p, r, a) for p, r, a in
              zip(sub["pos"].values, sub["ref"].values, sub["alt"].values)]
    out = pd.DataFrame({
        "Gene; HGNC ID": GENE_HGNC,
        "Category": "KURIS/NDD",
        "HGVS g": hgvs_g,
        "HGVS p": [_hgvs_p(s) for s in clin["Protein_Change"].values],
        "HGVS c": [_hgvs_c(s) for s in clin["HGVS_c"].values],
        "SGE Functional Score": sub["MAVE_Score"].values,
        "Domain": tables.annotate_domain(sub["AA_Position"].values).values,
        "Functional Classification": sub["score_threshold_class"].values,
        "ClinVar Classification (Aggregate Germline)": clin["ClinVar_Classification"].values,
        "ClinVar Accession": clin["VCV"].values,
        "ClinVar Stars": clin["Star_Rating"].values,
        "ClinVar Lab Entries Related to KURIS": _clean(clin["KURIS_ClinVar_SCVs"]).values,
        "Küry et al": _clean(clin["Kury_et_al"]).values,
        "DECIPHER": _clean(clin["DECIPHER"]).values,
    })
    return out.sort_values("SGE Functional Score").reset_index(drop=True)


def build_blb(m):
    """ClinVar B/LB missense rows from the master data table (Feb 2026 snapshot)."""
    blb = (m[m["panel3_category"] == "B/LB"]
           .drop_duplicates(subset="HGVS_c")
           .sort_values("MAVE_Score"))
    hgvs_g = [_hgvs_g(p, r, a) for p, r, a in
              zip(blb["pos"], blb["ref"], blb["alt"])]
    return pd.DataFrame({
        "Gene; HGNC ID": GENE_HGNC,
        "Category": "B/LB",
        "HGVS g": hgvs_g,
        "HGVS p": [_hgvs_p(s) for s in blb["HGVS_p"].values],
        "HGVS c": [_hgvs_c(s) for s in blb["HGVS_c"].values],
        "SGE Functional Score": blb["MAVE_Score"].values,
        "Domain": tables.annotate_domain(blb["AA_Position"].values).values,
        "Functional Classification": blb["score_threshold_class"].values,
        "ClinVar Classification (Aggregate Germline)": blb["clinvar_2026_classification"].values,
        "ClinVar Accession": [_vcv(v) for v in blb["clinvar_2026_variation_id"].values],
        "ClinVar Stars": [_stars(s) for s in blb["clinvar_2026_stars"].values],
        "ClinVar Lab Entries Related to KURIS": pd.NA,
        "Küry et al": pd.NA,
        "DECIPHER": pd.NA,
    }).reset_index(drop=True)


def build():
    m = tables.load_master()
    m = tables.add_panel3_category(m, tables.load_clinical())
    combined = pd.concat([build_kuris(m), build_blb(m)], ignore_index=True)
    return combined[COLUMNS]


def _pdf_view(df):
    d = df.copy()
    d["SGE Functional Score"] = d["SGE Functional Score"].map(
        lambda v: "" if pd.isna(v) else f"{v:.3f}")
    d["ClinVar Stars"] = d["ClinVar Stars"].map(
        lambda v: "" if pd.isna(v) else str(int(v)))
    for c in ["ClinVar Accession", "ClinVar Classification (Aggregate Germline)",
              "ClinVar Lab Entries Related to KURIS", "Küry et al", "DECIPHER"]:
        d[c] = d[c].fillna("")
    return d.rename(columns={
        "HGVS g": "HGVS g.", "HGVS p": "HGVS p.", "HGVS c": "HGVS c.",
        "SGE Functional Score": "SGE score", "Functional Classification": "GMM class",
        "ClinVar Classification (Aggregate Germline)": "ClinVar (agg. germline)",
        "ClinVar Accession": "ClinVar VCV", "ClinVar Stars": "Stars",
        "ClinVar Lab Entries Related to KURIS": "ClinVar KURIS lab entry (SCV)",
        "Küry et al": "Küry et al. 2022", "DECIPHER": "DECIPHER"})


def main():
    df = build()
    xlsx = config.TABLES / "SuppTable1_KURIS_NDD_and_BLB_Missense.xlsx"
    tsv = config.TABLES / "SuppTable1_KURIS_NDD_and_BLB_Missense.tsv"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        df.to_excel(w, sheet_name="KURIS_NDD_and_BLB_Missense", index=False)
    df.to_csv(tsv, sep="\t", index=False)
    kw = dict(
        aligns={"SGE score": "right", "Stars": "center", "GMM class": "center",
                "Domain": "center", "ClinVar VCV": "left"},
        group_cols=["HGVS c."],
        detail_cols=["ClinVar KURIS lab entry (SCV)", "Küry et al. 2022", "DECIPHER"],
        rich_cols=["Gene; HGNC ID"])
    tablepdf.render_table_pdf(_pdf_view(df), "SuppTable1_KURIS_NDD_and_BLB_Missense", **kw)
    supplement.render_supp_table("SuppTable1_KURIS_NDD_and_BLB_Missense", _pdf_view(df), **kw)
    n_kuris = int((df["Category"] == "KURIS/NDD").sum())
    n_blb = int((df["Category"] == "B/LB").sum())
    print(f"INFO: wrote tables/{xlsx.name} + .tsv + .pdf "
          f"({len(df)} variants: {n_kuris} KURIS/NDD + {n_blb} B/LB)")


if __name__ == "__main__":
    main()
