# BAP1 master table -- derived columns

Built by `scripts/make_master_table.py` from `data/raw/Waters_ExtendedData1.xlsx` (read-only). All original Waters columns are preserved unchanged; the columns below are appended.

| column | description |
|---|---|
| `simplified_consequence` | Single VEP class collapsed from vep_consequence by severity; exactly reproduces the notebook column. |
| `display_consequence` | Human-readable consequence label (same GMM classes) for the GMM figure. |
| `variant_type` | Coarser Figure 1/2 consequence group (Nonsense/Splice/Intronic/...); folds splice_region into Intronic/Synonymous. |
| `gmm_posterior_abnormal` | GMM posterior probability of the abnormal component (blank if no score). |
| `gmm_posterior_normal` | GMM posterior probability of the normal component (blank if no score). |
| `functional_consequence` | Posterior call: functionally_abnormal / functionally_normal / indeterminate (>=0.95). |
| `score_threshold_class` | Score-threshold call used by figures: Abnormal / Normal / Indeterminate. |
| `standardized_class` | Waters functional_classification standardized: depleted->Abnormal, unchanged->Normal, enriched->Not specified (GMM-vs-Waters figure). |
| `clinvar_2018_classification` | ClinVar clinical significance, Jan 2018 (CVFG), joined by HGVSc. |
| `clinvar_2018_slim` | Collapsed 2018 class: P/LP, B/LB, VUS, Conflicting. |
| `clinvar_2018_review_status` | ClinVar 2018 review status. |
| `clinvar_2018_stars` | ClinVar 2018 review star level. |
| `clinvar_2025_classification` | ClinVar clinical significance, Jan 2025 (CVFG) -- calibration truth set. |
| `clinvar_2025_slim` | Collapsed 2025 class: P/LP, B/LB, VUS, Conflicting. |
| `clinvar_2025_review_status` | ClinVar 2025 review status. |
| `clinvar_2025_stars` | ClinVar 2025 review star level. |
| `clinvar_2026_classification` | ClinVar germline classification, February 2026 export (cited in the paper) -- primary calibration truth set. |
| `clinvar_2026_slim` | Collapsed Feb 2026 class: P/LP, B/LB, VUS, Conflicting. |
| `clinvar_2026_review_status` | ClinVar Feb 2026 germline review status. |
| `clinvar_2026_stars` | ClinVar Feb 2026 review star level (mapped from review status). |
| `clinvar_2026_date` | ClinVar Feb 2026 germline date last evaluated. |
| `clinvar_2026_condition` | ClinVar Feb 2026 reported condition(s). |
| `clinvar_2026_oncogenicity` | ClinVar Feb 2026 oncogenicity classification, if any. |
| `clinvar_2026_variation_id` | ClinVar VariationID (Feb 2026). |
| `gnomad_v4_af` | gnomAD v4.1 allele frequency (CVFG). |
| `gnomad_v4_faf95_max` | gnomAD v4.1 max filtering allele frequency, 95% (CVFG). |
| `revel_score` | REVEL missense pathogenicity score (CVFG). |
| `alphamissense_pathogenicity` | AlphaMissense pathogenicity score (CVFG). |
| `alphamissense_class` | AlphaMissense class: likely_pathogenic / likely_benign / ambiguous (CVFG). |
| `spliceai_max_delta` | SpliceAI maximum delta score (CVFG). |
