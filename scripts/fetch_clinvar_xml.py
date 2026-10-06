#!/usr/bin/env python3
"""
Download the full ClinVar XML release and extract BAP1 variant data.

This is the UPSTREAM PROVENANCE for the two raw ClinVar inputs
(data/raw/BAP1_ClinVar_VCV.tsv and data/raw/BAP1_ClinVar_SCV.tsv). It is NOT run
by scripts/run_all.py: it needs a one-time ~3 GB download and the `lxml` package,
and the parsed TSVs are committed as raw inputs so the rest of the pipeline never
has to re-run it. Kept here so the provenance of those TSVs is self-documenting.

The committed TSVs are from the February 2026 ClinVar release (records last
updated through 2026-02-01), the same release as data/raw/ClinVar_Feb_2026.txt.
Because the download URL points at "…_00-latest", re-running this will fetch
whatever the current monthly release is; pin the release date when citing it.

Downloads ClinVarVCVRelease_00-latest.xml.gz from the NCBI FTP site
(single ~3 GB download), then stream-parses it to extract all BAP1 variants
with complete VCV- and SCV-level details.

Usage:
    python scripts/fetch_clinvar_xml.py
    python scripts/fetch_clinvar_xml.py --gene BAP1 --keep-full-xml
    python scripts/fetch_clinvar_xml.py --xml-path /path/to/already/downloaded.xml.gz

Output files (in data/raw/):
    BAP1_ClinVar_variants.xml    BAP1-only XML (extracted from full release)
    BAP1_ClinVar_VCV.tsv         Variant-level data (one row per VCV)
    BAP1_ClinVar_SCV.tsv         Submission-level data (one row per SCV)

The full release XML (~3 GB) is deleted after parsing unless --keep-full-xml
is specified. If you already have the file, pass --xml-path to skip download.

Requirements:
    pip install lxml pandas
"""
import argparse
import gzip
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from lxml import etree

CLINVAR_URLS = [
    # Current filename: ClinVarVCVRelease (renamed from ClinVarVariationRelease)
    # Monthly releases at https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/"
    "ClinVarVCVRelease_00-latest.xml.gz",
]
# Minimum expected size for the full release (currently ~3 GB)
MIN_DOWNLOAD_SIZE = 500 * 1024 * 1024  # 500 MB
GENE = "BAP1"
DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "raw"

# ACMG/AMP evidence codes pattern
ACMG_RE = re.compile(
    r"\b(PVS1|PS[1-4]|PM[1-6]|PP[1-5]|BA1|BS[1-4]|BP[1-7])"
    r"(?:_(Strong|Moderate|Supporting|Very_?Strong|Stand_?Alone))?\b",
    re.IGNORECASE,
)

# ClinVar review-status → star rating
STAR_RATING = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, single submitter": 1,
    "no assertion for the individual variant": 0,
    "no assertion criteria provided": 0,
    "no classification provided": 0,
    "no classification for the single variant": 0,
}


# ── Download ────────────────────────────────────────────────────────
def _has_command(name):
    """Check if a command is available on PATH."""
    import shutil
    return shutil.which(name) is not None


def _validate_gzip(path):
    """Check that a file is actually a gzip file (not an HTML error page)."""
    if not path.exists():
        return False, "File does not exist"
    size = path.stat().st_size
    if size < MIN_DOWNLOAD_SIZE:
        # Read first bytes to diagnose
        with open(path, "rb") as f:
            header = f.read(64)
        if header[:2] != b'\x1f\x8b':
            snippet = header[:200].decode("utf-8", errors="replace")
            return False, (f"Not a gzip file (got {size:,} bytes starting with: "
                           f"{snippet!r}). Server likely returned an error page.")
        return False, f"File too small ({size:,} bytes, expected >500 MB)"
    return True, "OK"


def download_clinvar_xml(dest_path):
    """Download full ClinVar XML release via curl or wget with resume support."""
    print(f"Downloading ClinVar XML release to {dest_path}")
    print(f"  This is ~3 GB — it may take 10-30 min depending on connection.")
    print()

    use_curl = _has_command("curl")
    use_wget = _has_command("wget")

    if not use_curl and not use_wget:
        print("ERROR: Neither curl nor wget found. Download manually:")
        print(f"  curl -L -C - -o {dest_path} {CLINVAR_URLS[0]}")
        sys.exit(1)

    for url in CLINVAR_URLS:
        print(f"  Trying: {url}")

        if use_curl:
            tool = "curl"
            cmd = [
                "curl",
                "--location",           # Follow redirects
                "--continue-at", "-",   # Resume partial downloads
                "--retry", "5",
                "--retry-delay", "5",
                "--connect-timeout", "60",
                "--fail",               # Fail on HTTP errors (no error pages)
                "--user-agent", "ClinVar-downloader/1.0",
                "--progress-bar",
                "--output", str(dest_path),
                url,
            ]
        else:
            tool = "wget"
            cmd = [
                "wget",
                "--continue",
                "--show-progress",
                "--timeout=60",
                "--tries=5",
                "-O", str(dest_path),
                url,
            ]

        print(f"  Using: {tool}")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"  Download failed with {tool} (exit code {e.returncode})")
            # Remove failed partial file before trying next URL
            if dest_path.exists() and dest_path.stat().st_size < MIN_DOWNLOAD_SIZE:
                dest_path.unlink()
            continue

        # Validate the download
        ok, msg = _validate_gzip(dest_path)
        if ok:
            size_gb = dest_path.stat().st_size / (1024**3)
            print(f"\nDownload complete: {size_gb:.2f} GB")
            return
        else:
            print(f"  Validation failed: {msg}")
            dest_path.unlink(missing_ok=True)
            continue

    # All URLs failed
    print(f"\nERROR: All download attempts failed.")
    print(f"  Please download manually and re-run with --xml-path:")
    print(f"")
    print(f"  curl -L -o {dest_path} {CLINVAR_URLS[0]}")
    print(f"  python {sys.argv[0]} --xml-path {dest_path}")
    sys.exit(1)


# ── Stream parse ────────────────────────────────────────────────────
def is_gene_match(va_element, gene):
    """Check if a VariationArchive element belongs to the target gene."""
    # Check VariationName attribute (fast path)
    name = va_element.get("VariationName", "")
    if f"({gene})" in name or f"{gene}:" in name:
        return True
    # Check Gene elements inside the XML (thorough path)
    for g in va_element.iter("Gene"):
        if g.get("Symbol", "") == gene:
            return True
    return False


def _text(element, tag, default=""):
    """Safely get text of a child element."""
    el = element.find(tag)
    return (el.text or default) if el is not None else default


def _get_seq_location(parent, assembly):
    """Extract SequenceLocation attributes for a given assembly."""
    for sl in parent.iter("SequenceLocation"):
        if sl.get("Assembly") == assembly:
            return {
                "Chr": sl.get("Chr", ""),
                "Accession": sl.get("Accession", ""),
                "start": sl.get("start", ""),
                "stop": sl.get("stop", ""),
                "display_start": sl.get("display_start", ""),
                "display_stop": sl.get("display_stop", ""),
                "variantLength": sl.get("variantLength", ""),
                "referenceAlleleVCF": sl.get("referenceAlleleVCF", ""),
                "alternateAlleleVCF": sl.get("alternateAlleleVCF", ""),
                "positionVCF": sl.get("positionVCF", ""),
                "forDisplay": sl.get("forDisplay", ""),
            }
    return {}


def _collect_xrefs(parent, target_db=None):
    """Collect XRef elements, optionally filtered by DB name."""
    refs = []
    for xref in parent.iter("XRef"):
        db = xref.get("DB", "")
        xid = xref.get("ID", "")
        xtype = xref.get("Type", "")
        if target_db is None or db == target_db:
            refs.append({"DB": db, "ID": xid, "Type": xtype})
    return refs


def _extract_acmg_codes(text_sources):
    """Extract ACMG/AMP codes from a list of text strings."""
    codes = set()
    for text in text_sources:
        for match in ACMG_RE.finditer(text):
            code = match.group(1).upper()
            modifier = match.group(2)
            if modifier:
                code += f"_{modifier}"
            codes.add(code)
    return ";".join(sorted(codes)) if codes else ""


def _extract_pmids(parent):
    """Extract all PubMed IDs from Citation elements under parent."""
    pmids = set()
    for citation in parent.iter("Citation"):
        for cid in citation.iter("ID"):
            if cid.get("Source", "") == "PubMed" and cid.text:
                pmids.add(cid.text)
    return ";".join(sorted(pmids)) if pmids else ""


def parse_variation_archive(va):
    """
    Extract comprehensive VCV- and SCV-level records from one
    VariationArchive element.

    Returns (vcv_record, scv_records):
        vcv_record  — dict with one row of variant-level data
        scv_records — list of dicts, one per SCV submission
    """
    # ================================================================
    # VCV-level (variant) data
    # ================================================================
    vcv = va.get("Accession", "")
    variation_id = va.get("VariationID", "")
    variant_name = va.get("VariationName", "")
    variation_type_attr = va.get("VariationType", "")
    date_created_va = va.get("DateCreated", "")
    date_updated_va = va.get("DateLastUpdated", "")
    num_submitters = va.get("NumberOfSubmitters", "")
    num_submissions = va.get("NumberOfSubmissions", "")
    record_type = va.get("RecordType", "")

    record_status = _text(va, "RecordStatus")
    species = _text(va, "Species")

    # ── SimpleAllele fields ──
    allele_id = ""
    variant_type = ""
    variant_length = ""
    canonical_spdi = ""
    protein_changes = []
    mol_consequences = []
    hgvs_list = []
    functional_consequence = ""
    functional_consequence_xref = ""
    cytogenetic_location = ""
    dbsnp_rsid = ""
    all_xrefs = []
    grch38 = {}
    grch37 = {}
    global_maf_value = ""
    global_maf_source = ""
    global_maf_minor_allele = ""

    # Gene-level info
    gene_symbol = ""
    gene_id = ""
    hgnc_id = ""
    gene_full_name = ""

    sa_el = va.find(".//SimpleAllele")
    if sa_el is not None:
        allele_id = sa_el.get("AlleleID", "")

        # VariantType
        vt = sa_el.find("VariantType")
        if vt is not None:
            variant_type = vt.text or ""

        # VariantLength
        vl = sa_el.find("VariantLength")
        if vl is not None:
            variant_length = vl.text or ""

        # CanonicalSPDI
        spdi = sa_el.find("CanonicalSPDI")
        if spdi is not None:
            canonical_spdi = spdi.text or ""

        # ProteinChange(s)
        for pc in sa_el.iter("ProteinChange"):
            if pc.text:
                protein_changes.append(pc.text)

        # MolecularConsequence(s)
        for mc in sa_el.iter("MolecularConsequence"):
            mc_type = mc.get("Type", "")
            mc_db = mc.get("DB", "")
            mc_id = mc.get("ID", "")
            mol_consequences.append(
                f"{mc_type}" + (f" ({mc_db}:{mc_id})" if mc_id else ""))

        # HGVS list
        for hgvs_el in sa_el.iter("HGVS"):
            expr = hgvs_el.text or ""
            htype = hgvs_el.get("Type", "")
            assembly = hgvs_el.get("Assembly", "")
            change = hgvs_el.get("Change", "")
            hgvs_list.append(
                f"{expr}" + (f" [{htype}]" if htype else "")
                + (f" [{assembly}]" if assembly else ""))

        # FunctionalConsequence
        fc = sa_el.find("FunctionalConsequence")
        if fc is not None:
            functional_consequence = fc.get("Value", "")
            fc_xrefs = _collect_xrefs(fc)
            if fc_xrefs:
                functional_consequence_xref = ";".join(
                    f"{x['DB']}:{x['ID']}" for x in fc_xrefs)

        # Location (allele-level)
        loc = sa_el.find("Location")
        if loc is not None:
            cyto = loc.find("CytogeneticLocation")
            if cyto is not None:
                cytogenetic_location = cyto.text or ""
            grch38 = _get_seq_location(loc, "GRCh38")
            grch37 = _get_seq_location(loc, "GRCh37")

        # XRefs (dbSNP, ClinGen, etc.)
        xref_list = sa_el.find("XRefList")
        if xref_list is not None:
            for xref in xref_list:
                db = xref.get("DB", "")
                xid = xref.get("ID", "")
                xtype = xref.get("Type", "")
                all_xrefs.append(f"{db}:{xid}")
                if db == "dbSNP":
                    dbsnp_rsid = xid

        # GlobalMinorAlleleFrequency
        gmaf = sa_el.find("GlobalMinorAlleleFrequency")
        if gmaf is not None:
            global_maf_value = gmaf.get("Value", "")
            global_maf_source = gmaf.get("Source", "")
            global_maf_minor_allele = gmaf.get("MinorAllele", "")

        # Gene info
        for gene_el in sa_el.iter("Gene"):
            gene_symbol = gene_el.get("Symbol", "")
            gene_id = gene_el.get("GeneID", "")
            hgnc_id = gene_el.get("HGNC_ID", "")
            gene_full_name = gene_el.get("FullName", "")
            break

    # ── ClassifiedRecord / Classifications ──
    agg_germline = ""
    agg_germline_review = ""
    agg_germline_date = ""
    agg_germline_date_created = ""
    agg_somatic = ""
    agg_somatic_review = ""
    agg_onco = ""
    agg_onco_review = ""

    cr = va.find(".//ClassifiedRecord")
    if cr is not None:
        cls = cr.find("Classifications")
        if cls is not None:
            gc = cls.find("GermlineClassification")
            if gc is not None:
                agg_germline = _text(gc, "Description")
                agg_germline_review = _text(gc, "ReviewStatus")
                agg_germline_date = gc.get("DateLastEvaluated", "")
                agg_germline_date_created = gc.get("DateCreated", "")
            sc = cls.find("SomaticClinicalImpact")
            if sc is not None:
                agg_somatic = _text(sc, "Description")
                agg_somatic_review = _text(sc, "ReviewStatus")
            oc = cls.find("OncogenicityClassification")
            if oc is not None:
                agg_onco = _text(oc, "Description")
                agg_onco_review = _text(oc, "ReviewStatus")

    star_rating = STAR_RATING.get(agg_germline_review.lower(), "")
    if star_rating == "" and agg_somatic_review:
        star_rating = STAR_RATING.get(agg_somatic_review.lower(), "")

    # ── Aggregate condition (from Classifications > ConditionList) ──
    agg_conditions = []
    if cr is not None:
        for cls_parent in cr.iter("Classifications"):
            for cond_list in cls_parent.iter("ConditionList"):
                for ts in cond_list.iter("TraitSet"):
                    for trait in ts.iter("Trait"):
                        name_el = trait.find(".//Name/ElementValue")
                        if name_el is not None and name_el.text:
                            agg_conditions.append(name_el.text)
            break
    agg_conditions_str = ";".join(agg_conditions)

    # ── RCV accessions ──
    rcv_list = []
    if cr is not None:
        for rcv in cr.iter("RCVAccession"):
            racc = rcv.get("Accession", "")
            rtitle = rcv.get("Title", "")
            rcv_list.append(f"{racc}" + (f" ({rtitle})" if rtitle else ""))
    rcv_str = ";".join(rcv_list) if rcv_list else ""

    vcv_record = {
        "VCV": vcv,
        "VariationID": variation_id,
        "VariantName": variant_name,
        "RecordStatus": record_status,
        "RecordType": record_type,
        "Species": species,
        "AlleleID": allele_id,
        "VariationType": variation_type_attr or variant_type,
        "VariantLength": variant_length,
        "DateCreated": date_created_va,
        "DateLastUpdated": date_updated_va,
        "NumSubmitters": num_submitters,
        "NumSubmissions": num_submissions,
        # Gene
        "GeneSymbol": gene_symbol,
        "GeneID": gene_id,
        "HGNC_ID": hgnc_id,
        "GeneFullName": gene_full_name,
        # Allele details
        "CanonicalSPDI": canonical_spdi,
        "ProteinChange": ";".join(protein_changes),
        "MolecularConsequence": ";".join(mol_consequences),
        "HGVSexpressions": " | ".join(hgvs_list),
        "FunctionalConsequence": functional_consequence,
        "FunctionalConsequenceXRef": functional_consequence_xref,
        "CytogeneticLocation": cytogenetic_location,
        "dbSNP_rsID": dbsnp_rsid,
        "OtherXRefs": ";".join(all_xrefs),
        # GRCh38 coordinates
        "GRCh38_Chr": grch38.get("Chr", ""),
        "GRCh38_Accession": grch38.get("Accession", ""),
        "GRCh38_start": grch38.get("start", ""),
        "GRCh38_stop": grch38.get("stop", ""),
        "GRCh38_referenceAlleleVCF": grch38.get("referenceAlleleVCF", ""),
        "GRCh38_alternateAlleleVCF": grch38.get("alternateAlleleVCF", ""),
        "GRCh38_positionVCF": grch38.get("positionVCF", ""),
        "GRCh38_variantLength": grch38.get("variantLength", ""),
        # GRCh37 coordinates
        "GRCh37_Chr": grch37.get("Chr", ""),
        "GRCh37_Accession": grch37.get("Accession", ""),
        "GRCh37_start": grch37.get("start", ""),
        "GRCh37_stop": grch37.get("stop", ""),
        "GRCh37_referenceAlleleVCF": grch37.get("referenceAlleleVCF", ""),
        "GRCh37_alternateAlleleVCF": grch37.get("alternateAlleleVCF", ""),
        "GRCh37_positionVCF": grch37.get("positionVCF", ""),
        "GRCh37_variantLength": grch37.get("variantLength", ""),
        # MAF
        "GlobalMAF": global_maf_value,
        "GlobalMAF_Source": global_maf_source,
        "GlobalMAF_MinorAllele": global_maf_minor_allele,
        # Aggregate classifications
        "AggGermlineClass": agg_germline,
        "AggGermlineReviewStatus": agg_germline_review,
        "AggGermlineDateLastEval": agg_germline_date,
        "AggGermlineDateCreated": agg_germline_date_created,
        "AggSomaticClass": agg_somatic,
        "AggSomaticReviewStatus": agg_somatic_review,
        "AggOncoClass": agg_onco,
        "AggOncoReviewStatus": agg_onco_review,
        "StarRating": star_rating,
        "AggConditions": agg_conditions_str,
        # RCV
        "RCVaccessions": rcv_str,
    }

    # ================================================================
    # SCV-level (submission) data
    # ================================================================
    scv_records = []

    for ca in va.iter("ClinicalAssertion"):
        scv_acc = ca.get("Accession", "")
        scv_ver = ca.get("Version", "")
        scv_id = ca.get("ID", "")
        scv_date_created = ca.get("DateCreated", "")
        scv_date_updated = ca.get("DateLastUpdated", "")
        submission_date = ca.get("SubmissionDate", "")

        scv_record_status = _text(ca, "RecordStatus")

        # ── ClinVarAccession element ──
        acc_el = ca.find("ClinVarAccession")
        org_id = ""
        org_type = ""
        org_category = ""
        org_abbrev = ""
        if acc_el is not None:
            org_id = acc_el.get("OrgID", "")
            org_type = acc_el.get("OrgType", "")
            org_category = acc_el.get("OrganizationCategory", "")
            org_abbrev = acc_el.get("OrgAbbreviation", "")
            if not scv_acc:
                scv_acc = acc_el.get("Accession", "")
                scv_ver = acc_el.get("Version", "")

        # ── ClinVarSubmissionID ──
        submitter = ""
        local_key = ""
        submitter_date = ""
        sub_el = ca.find("ClinVarSubmissionID")
        if sub_el is not None:
            submitter = sub_el.get("submitter", "")
            local_key = sub_el.get("localKey", "")
            submitter_date = sub_el.get("submitterDate", "")

        # ── Assertion ──
        assertion_type = _text(ca, "Assertion")

        # ── Classification ──
        interp_type = ""
        classification = ""
        date_evaluated = ""
        review_status = ""

        cls_el = ca.find("Classification")
        if cls_el is not None:
            for tag, itype in [
                ("GermlineClassification", "germline"),
                ("SomaticClinicalImpact", "somatic"),
                ("OncogenicityClassification", "oncogenicity"),
            ]:
                child = cls_el.find(tag)
                if child is not None and child.text:
                    interp_type = itype
                    classification = child.text
                    break
            if not classification:
                desc = cls_el.find("Description")
                if desc is not None:
                    classification = desc.text or ""
            rs = cls_el.find("ReviewStatus")
            if rs is not None:
                review_status = rs.text or ""
            de = cls_el.get("DateLastEvaluated", "")
            if de:
                date_evaluated = de

        # ── Assertion method and its citation ──
        assertion_method = ""
        assertion_method_citation = ""
        for attr_set in ca.iter("AttributeSet"):
            for attr in attr_set.iter("Attribute"):
                if attr.get("Type") == "AssertionMethod":
                    assertion_method = attr.text or ""
                    # Look for citation in same AttributeSet
                    for cit in attr_set.iter("Citation"):
                        for cid in cit.iter("ID"):
                            src = cid.get("Source", "")
                            if src and cid.text:
                                assertion_method_citation = (
                                    f"{src}:{cid.text}")
                                break
                    break

        # ── Condition / Trait (asserted) ──
        conditions = []
        condition_xrefs = []
        # Find the TraitSet that is a direct child or under the assertion
        # (not under ObservedIn)
        for ts in ca.findall("TraitSet"):
            ts_type = ts.get("Type", "")
            for trait in ts.iter("Trait"):
                trait_type = trait.get("Type", "")
                name = ""
                name_el = trait.find(".//Name/ElementValue")
                if name_el is not None:
                    name = name_el.text or ""
                xrefs = []
                for xref in trait.iter("XRef"):
                    db = xref.get("DB", "")
                    xid = xref.get("ID", "")
                    if db and xid:
                        xrefs.append(f"{db}:{xid}")
                conditions.append(name)
                condition_xrefs.extend(xrefs)
        condition_str = ";".join(conditions) if conditions else ""
        condition_xref_str = ";".join(condition_xrefs) if condition_xrefs else ""

        # ── ObservedIn data ──
        origins = []
        affected_statuses = []
        num_families = ""
        num_families_with_variant = ""
        num_individuals = ""
        num_tested = ""
        genders = []
        observed_phenotypes = []
        observed_hpo_ids = []
        obs_descriptions = []
        obs_pmids = set()

        for obs in ca.iter("ObservedIn"):
            for sample in obs.iter("Sample"):
                for o in sample.iter("Origin"):
                    if o.text:
                        origins.append(o.text)
                for af in sample.iter("AffectedStatus"):
                    if af.text:
                        affected_statuses.append(af.text)
                nt = sample.find("NumberTested")
                if nt is not None and nt.text:
                    num_tested = nt.text
                for gend in sample.iter("Gender"):
                    if gend.text:
                        genders.append(gend.text)
                fd = sample.find("FamilyData")
                if fd is not None:
                    num_families = fd.get("NumFamilies", "")
                    num_families_with_variant = fd.get(
                        "NumFamiliesWithVariant", "")
                nf = sample.find("NumberOfFamilies")
                if nf is not None and nf.text:
                    num_families = nf.text
                ni = sample.find("NumberOfIndividuals")
                if ni is not None and ni.text:
                    num_individuals = ni.text

            # Observed phenotypes (TraitSet under ObservedIn)
            for ts in obs.iter("TraitSet"):
                for trait in ts.iter("Trait"):
                    name_el = trait.find(".//Name/ElementValue")
                    if name_el is not None and name_el.text:
                        observed_phenotypes.append(name_el.text)
                    for xref in trait.iter("XRef"):
                        db = xref.get("DB", "")
                        xid = xref.get("ID", "")
                        if db in ("HP", "HPO") and xid:
                            observed_hpo_ids.append(xid)

            # ObservedData descriptions and citations
            for od in obs.iter("ObservedData"):
                for attr in od.iter("Attribute"):
                    if attr.get("Type") == "Description" and attr.text:
                        obs_descriptions.append(attr.text)
                for cit in od.iter("Citation"):
                    for cid in cit.iter("ID"):
                        if cid.get("Source", "") == "PubMed" and cid.text:
                            obs_pmids.add(cid.text)

        origin_str = ";".join(origins) if origins else ""
        affected_str = ";".join(affected_statuses) if affected_statuses else ""
        gender_str = ";".join(genders) if genders else ""
        obs_pheno_str = ";".join(observed_phenotypes) if observed_phenotypes else ""
        obs_hpo_str = ";".join(observed_hpo_ids) if observed_hpo_ids else ""
        obs_desc_str = " | ".join(
            d[:500] for d in obs_descriptions) if obs_descriptions else ""

        # ── Collection method ──
        methods = []
        for m in ca.iter("Method"):
            mt = m.find("MethodType")
            if mt is not None and mt.text:
                methods.append(mt.text)
        method_str = ";".join(methods) if methods else ""

        # ── Mode of inheritance ──
        moi = ""
        for attr_set in ca.iter("AttributeSet"):
            for attr in attr_set.iter("Attribute"):
                if attr.get("Type") == "ModeOfInheritance" and attr.text:
                    moi = attr.text
                    break

        # ── All AttributeSet key-value pairs (beyond assertion method / MOI) ──
        extra_attributes = []
        for attr_set in ca.iter("AttributeSet"):
            for attr in attr_set.iter("Attribute"):
                atype = attr.get("Type", "")
                if atype in ("AssertionMethod", "ModeOfInheritance"):
                    continue
                if attr.text:
                    extra_attributes.append(f"{atype}={attr.text}")

        # ── ACMG codes ──
        acmg_text_sources = []
        for attr_set in ca.iter("AttributeSet"):
            for attr in attr_set.iter("Attribute"):
                if attr.get("Type") in ("AssertionMethod", "ModeOfInheritance"):
                    continue
                if attr.text:
                    acmg_text_sources.append(attr.text)
        for comment in ca.iter("Comment"):
            if comment.text:
                acmg_text_sources.append(comment.text)
        desc_el = ca.find(".//Description")
        if desc_el is not None and desc_el.text:
            acmg_text_sources.append(desc_el.text)
        acmg_str = _extract_acmg_codes(acmg_text_sources)

        # ── All citations/PMIDs ──
        pmid_str = _extract_pmids(ca)

        # ── Comments ──
        comments = []
        for c in ca.iter("Comment"):
            if c.text:
                comments.append(c.text[:1000])
        comment_str = " | ".join(comments) if comments else ""

        # ── StudyName / SubmissionName ──
        study_name = _text(ca, "StudyName")
        submission_names = []
        for sn_list in ca.iter("SubmissionNameList"):
            for sn in sn_list.iter("SubmissionName"):
                if sn.text:
                    submission_names.append(sn.text)
        submission_name_str = ";".join(submission_names) if submission_names else ""

        scv_records.append({
            "VCV": vcv,
            "VariationID": variation_id,
            "SCV": scv_acc,
            "SCV_Version": scv_ver,
            "SCV_ID": scv_id,
            "RecordStatus": scv_record_status,
            "DateCreated": scv_date_created,
            "DateLastUpdated": scv_date_updated,
            "SubmissionDate": submission_date,
            # Submitter
            "Submitter": submitter,
            "SubmitterOrgID": org_id,
            "OrgType": org_type,
            "OrganizationCategory": org_category,
            "OrgAbbreviation": org_abbrev,
            "LocalKey": local_key,
            "SubmitterDate": submitter_date,
            # Classification
            "InterpretationType": interp_type,
            "Classification": classification,
            "DateLastEvaluated": date_evaluated,
            "ReviewStatus": review_status,
            "AssertionType": assertion_type,
            "AssertionMethod": assertion_method,
            "AssertionMethodCitation": assertion_method_citation,
            # Condition
            "Condition": condition_str,
            "ConditionXRefs": condition_xref_str,
            # Observed data
            "Origin": origin_str,
            "AffectedStatus": affected_str,
            "Gender": gender_str,
            "NumFamilies": num_families,
            "NumFamiliesWithVariant": num_families_with_variant,
            "NumIndividuals": num_individuals,
            "NumTested": num_tested,
            "CollectionMethod": method_str,
            "ModeOfInheritance": moi,
            # Observed phenotypes
            "ObservedPhenotypes": obs_pheno_str,
            "ObservedHPO": obs_hpo_str,
            "ObservedDataDescriptions": obs_desc_str,
            # Evidence
            "ACMG_Codes": acmg_str,
            "PMIDs": pmid_str,
            "ObservedDataPMIDs": ";".join(sorted(obs_pmids)),
            "ExtraAttributes": ";".join(extra_attributes),
            # Metadata
            "Comment": comment_str,
            "StudyName": study_name,
            "SubmissionNames": submission_name_str,
        })

    return vcv_record, scv_records


def stream_parse_clinvar(xml_path, gene):
    """
    Stream-parse the ClinVar XML release, extracting only variants for
    the target gene. Memory-efficient: processes one VariationArchive
    at a time, never loading the full file into memory.

    Returns (vcv_records, scv_records, gene_xml_chunks) where:
        vcv_records     — list of dicts (one per variant)
        scv_records     — list of dicts (one per submission)
        gene_xml_chunks — list of XML bytestrings for matched elements
    """
    vcv_records = []
    scv_records = []
    gene_xml_chunks = []

    # Determine how to open the file
    if str(xml_path).endswith(".gz"):
        source = gzip.open(xml_path, "rb")
    else:
        source = open(xml_path, "rb")

    n_total = 0
    n_matched = 0

    print(f"Stream-parsing {xml_path.name} for gene={gene}...")
    print(f"  (This may take 10-20 min for the full release)")

    try:
        context = etree.iterparse(
            source,
            events=("end",),
            tag="VariationArchive",
        )

        for event, elem in context:
            n_total += 1

            if n_total % 50000 == 0:
                print(f"  Processed {n_total:,} variants, "
                      f"found {n_matched} {gene} matches...",
                      flush=True)

            if is_gene_match(elem, gene):
                n_matched += 1
                # Save raw XML for this variant
                gene_xml_chunks.append(etree.tostring(elem, pretty_print=True))
                # Parse VCV + SCV records
                vcv_rec, scv_recs = parse_variation_archive(elem)
                vcv_records.append(vcv_rec)
                scv_records.extend(scv_recs)

            # Free memory: clear the element and its predecessors
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]

    finally:
        source.close()

    print(f"  Done: {n_total:,} total variants, {n_matched} {gene} matches")
    print(f"    VCV records: {len(vcv_records)}")
    print(f"    SCV records: {len(scv_records)}")
    return vcv_records, scv_records, gene_xml_chunks


# ── Main ────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Download ClinVar XML release and extract gene-specific "
                    "SCV-level data")
    parser.add_argument(
        "--gene", default=GENE,
        help="Gene symbol to extract (default: BAP1)")
    parser.add_argument(
        "--xml-path", type=Path, default=None,
        help="Path to already-downloaded ClinVarVariationRelease XML "
             "(skip download)")
    parser.add_argument(
        "--keep-full-xml", action="store_true",
        help="Keep the full ~3 GB release XML after parsing "
             "(default: delete it)")
    args = parser.parse_args()

    os.makedirs(DATA_DIR, exist_ok=True)
    gene = args.gene

    # Step 1: Get the XML file
    if args.xml_path:
        xml_path = args.xml_path
        if not xml_path.exists():
            print(f"ERROR: {xml_path} does not exist")
            sys.exit(1)
        print(f"Using existing XML: {xml_path}")
    else:
        xml_path = DATA_DIR / "ClinVarVariationRelease_00-latest.xml.gz"
        if xml_path.exists():
            size_gb = xml_path.stat().st_size / (1024**3)
            print(f"Found existing download: {xml_path} ({size_gb:.2f} GB)")
            print(f"  Delete it and re-run to force a fresh download.")
        else:
            download_clinvar_xml(xml_path)

    # Step 2: Stream-parse for target gene
    vcv_records, scv_records, gene_xml_chunks = stream_parse_clinvar(
        xml_path, gene)

    if not vcv_records:
        print(f"WARNING: No {gene} variants found in {xml_path.name}")
        print(f"  Check that the gene symbol is correct.")
        sys.exit(1)

    # Step 3: Save gene-specific XML
    gene_xml_path = DATA_DIR / f"{gene}_ClinVar_variants.xml"
    with open(gene_xml_path, "wb") as f:
        f.write(b'<?xml version="1.0" encoding="UTF-8"?>\n')
        f.write(f'<ClinVarVCVRelease gene="{gene}">\n'.encode())
        for chunk in gene_xml_chunks:
            f.write(chunk)
        f.write(b"</ClinVarVCVRelease>\n")
    size_mb = gene_xml_path.stat().st_size / (1024**2)
    print(f"\n{gene} XML saved: {gene_xml_path} ({size_mb:.1f} MB)")

    # Step 4: Save VCV-level TSV (one row per variant)
    df_vcv = pd.DataFrame(vcv_records)
    vcv_tsv_path = DATA_DIR / f"{gene}_ClinVar_VCV.tsv"
    df_vcv.to_csv(vcv_tsv_path, sep="\t", index=False)
    print(f"VCV TSV saved: {vcv_tsv_path} ({len(df_vcv)} variants)")

    # Step 5: Save SCV-level TSV (one row per submission)
    df_scv = pd.DataFrame(scv_records)
    scv_tsv_path = DATA_DIR / f"{gene}_ClinVar_SCV.tsv"
    df_scv.to_csv(scv_tsv_path, sep="\t", index=False)
    print(f"SCV TSV saved: {scv_tsv_path} ({len(df_scv)} submissions)")

    # Step 6: Clean up full XML if not keeping
    if not args.keep_full_xml and not args.xml_path:
        if xml_path.exists() and ("ClinVarVCVRelease" in xml_path.name
                                  or "ClinVarVariationRelease" in xml_path.name):
            print(f"\nRemoving full release XML: {xml_path.name}")
            xml_path.unlink()
            print(f"  (Use --keep-full-xml to retain it)")

    # Summary
    print(f"\n{'=' * 60}")
    print(f"Summary for {gene}:")
    print(f"  Variants (VCV):           {len(df_vcv)}")
    print(f"  Submissions (SCV):        {len(df_scv)}")
    print(f"")
    print(f"  VCV-level:")
    print(f"    With GRCh38 coords:     "
          f"{(df_vcv['GRCh38_Chr'] != '').sum()}")
    print(f"    With CanonicalSPDI:     "
          f"{(df_vcv['CanonicalSPDI'] != '').sum()}")
    print(f"    With dbSNP rsID:        "
          f"{(df_vcv['dbSNP_rsID'] != '').sum()}")
    print(f"    With GlobalMAF:         "
          f"{(df_vcv['GlobalMAF'] != '').sum()}")
    print(f"    With star rating:       "
          f"{(df_vcv['StarRating'] != '').sum()}")
    print(f"    With HGVS expressions:  "
          f"{(df_vcv['HGVSexpressions'] != '').sum()}")
    print(f"    With RCV accessions:    "
          f"{(df_vcv['RCVaccessions'] != '').sum()}")
    print(f"")
    print(f"  SCV-level:")
    print(f"    With ACMG codes:        "
          f"{(df_scv['ACMG_Codes'] != '').sum()}")
    print(f"    With PMIDs:             "
          f"{(df_scv['PMIDs'] != '').sum()}")
    print(f"    With origin info:       "
          f"{(df_scv['Origin'] != '').sum()}")
    print(f"    With observed pheno:    "
          f"{(df_scv['ObservedPhenotypes'] != '').sum()}")
    print(f"    With HPO terms:         "
          f"{(df_scv['ObservedHPO'] != '').sum()}")
    print(f"    With affected status:   "
          f"{(df_scv['AffectedStatus'] != '').sum()}")
    print(f"    Current records:        "
          f"{(df_scv['RecordStatus'] == 'current').sum()}")
    print(f"    Superseded records:     "
          f"{(df_scv['RecordStatus'] != 'current').sum()}")
    print(f"")
    print(f"  Aggregate germline classifications:")
    print(f"    {df_vcv['AggGermlineClass'].value_counts().to_string()}")
    print(f"")
    print(f"  SCV interpretation types:")
    print(f"    {df_scv['InterpretationType'].value_counts().to_string()}")
    print(f"")
    print(f"  Top SCV classifications:")
    print(f"    {df_scv['Classification'].value_counts().head(10).to_string()}")
    print(f"")
    print(f"  VCV columns ({len(df_vcv.columns)}):")
    for c in df_vcv.columns:
        print(f"    {c}")
    print(f"")
    print(f"  SCV columns ({len(df_scv.columns)}):")
    for c in df_scv.columns:
        print(f"    {c}")
    print(f"\nOutput files:")
    print(f"  {gene_xml_path}")
    print(f"  {vcv_tsv_path}")
    print(f"  {scv_tsv_path}")


if __name__ == "__main__":
    main()
