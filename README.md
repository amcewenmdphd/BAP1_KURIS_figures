# BAP1 KURIS — figures & tables

Reproduces every figure and table in the BAP1 KURIS manuscript from the raw
inputs. All shared code lives in the `bap1figs` package.

## Reproducibility chain

```
data/raw/                         (Waters ED1, CVFG, ClinVar VCV/SCV, DECIPHER, Kaplanis, Kury, ClinVar Feb 2026)
   │  build_master_table.py
   │  build_clinical_intermediates.py
   ▼
data/source/                      (master data table + clinical data table + GMM thresholds + ClinVar intermediates)
   │  each figure/table script rebuilds its own input file …
   ▼
data/figure_inputs/               (one input file per figure)
   │  … then renders
   ▼
figures/   tables/                (working outputs)
   │  build_supplement.py         (numbers, captions, composes)
   ▼
supplement/                       (final Supplementary_Figure_N.pdf / Supplementary_Table_N.pdf + guide)
```

## Run everything

```bash
pip install -r requirements.txt
python scripts/run_all.py
```

Or run any single figure/table on its own (each is self-contained):

```bash
python scripts/figure2.py            # rebuilds data/figure_inputs/figure2_variants.xlsx, then Figure 2
python scripts/table_master_contingency.py
```
## Scripts

Upstream (raw → source tables):

| Script | Builds |
|---|---|
| `build_master_table.py` | `data/source/BAP1_master_table.tsv`, `GMM_thresholds.json` |
| `build_clinical_intermediates.py` | ClinVar SCV/VCV curation + KURIS/Cancer merges |
| `fetch_clinvar_xml.py` | provenance only — `data/raw/BAP1_ClinVar_{VCV,SCV}.tsv` (not in `run_all`) |

Figures:

| Script | Figure |
|---|---|
| `figure1.py` | Figure 1 — lollipop + ClinVar panels |
| `figure2.py` | Figure 2 — functional scores + LR contingency |
| `suppfig1_gmm.py` | Supp Fig 1 — GMM calibration |
| `suppfig_waters_vs_gmm.py` | Waters vs GMM classification |
| `suppfig_lr_master.py` | master contingency forest |
| `suppfig_lr_uncertainty_kuris.py` | KURIS/NDD LR stability |

Tables:

| Script | Table |
|---|---|
| `table_master_contingency.py` | LR by method × control set × class (+ CI + ACMG tier) |
| `table_lr_uncertainty.py` | LR stability (log-Wald / bootstrap / Bayesian / LOO) |
| `table_waters_vs_gmm_performance.py` | sensitivity / specificity / PPV / NPV |
| `table_supp1_kuris_blb_missense.py` | Supp Table 1 — KURIS/NDD + B/LB missense (HGVS g./c./p., score, domain, ClinVar; reference sequences in the footnote) |
| `table_supp_figure1_variants.py` | Supp Table 2 — all Figure 1 variants (HGVS g./c./p., SGE + ClinVar) |
| `table_supp2_phenotypic_overlap.py` | Supp Table 6 — proband vs KURIS-cohort phenotype overlap |

Each table script writes an Excel workbook, a tab-separated file, and a typeset
`.pdf` (via `bap1figs.tablepdf`) into `tables/`, plus a numbered, captioned copy
into `supplement/`.

## Supplement

`scripts/build_supplement.py` assembles the numbered, captioned supplement into
`supplement/`, in citation order:

| # | Figure | # | Table |
|---|---|---|---|
| Fig 1 | GMM calibration | Table 1 | KURIS/NDD + B/LB missense variants |
| Fig 2 | Waters vs GMM concordance | Table 2 | all variants shown in Figure 1 |
| Fig 3 | KURIS/NDD LR stability | Table 3 | classifier performance |
| Fig 4 | method × control comparison | Table 4 | master contingency table |
|  |  | Table 5 | LR uncertainty by estimator |
|  |  | Table 6 | phenotypic overlap |
