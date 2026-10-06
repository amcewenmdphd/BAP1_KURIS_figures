# BAP1 KURIS — figures & tables

Reproduces every figure and table in the BAP1 KURIS manuscript from the raw
inputs. **One script per figure/table**; all shared code (color schemes,
likelihood-ratio math, plotting helpers) lives in the `bap1figs` package.

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

Every figure/table script has two steps: **(1)** rebuild its own input file from
the two source tables (the master data table + the clinical data table), then
**(2)** render the figure/table from that input file. So a user with only the two
source tables can recreate every input file, figure, and table.

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

## The `bap1figs` shared library

| Module | Contents |
|---|---|
| `config`   | paths + matplotlib defaults (Arial, 8 pt, vector fonts, 600 dpi) |
| `palettes` | all color schemes + ACMG PS3 (red) / BS3 (blue) evidence bands |
| `stats`    | GMM classification, likelihood ratio + log-Wald CI (Haldane–Anscombe), ACMG evidence tiers, Wilson CIs |
| `tables`   | master/clinical loaders, `panel3_category` (KURIS/NDD, TPDS P/LP, B/LB, N229K), domain annotation (UniProt Q92560) |
| `plots`    | panel labels, contingency heatmaps, broken-axis histograms, single-axis ACMG forests |
| `tablepdf` | typeset PDF renderer shared by the supplementary-table scripts (title/caption aware) |
| `supplement` | canonical figure/table order, titles, and Genetics-in-Medicine-style captions |

## Scripts

Upstream (raw → source tables):

| Script | Builds |
|---|---|
| `build_master_table.py` | `data/source/BAP1_master_table.tsv`, `GMM_thresholds.json` |
| `build_clinical_intermediates.py` | ClinVar SCV/VCV curation + KURIS/Cancer merges |
| `fetch_clinvar_xml.py` | provenance only — `data/raw/BAP1_ClinVar_{VCV,SCV}.tsv` (not in `run_all`) |

### ClinVar ingestion (three paths)

ClinVar enters the pipeline three ways, all committed as raw inputs so nothing
needs re-fetching:

1. **Full XML release** — `fetch_clinvar_xml.py` downloads the ~3 GB
   `ClinVarVCVRelease_00-latest.xml.gz` and `lxml`-streams it to
   `data/raw/BAP1_ClinVar_VCV.tsv` (one row/VCV) and `…_SCV.tsv` (one row/SCV).
   These feed the clinical curation (conditions, SCV assertions, KURIS/NDD set,
   Figure 1c). Snapshot: **February 2026** release (records through 2026-02-01).
   This script is heavyweight and is **not** part of `run_all.py`.
2. **ClinVar website tab export** — `data/raw/ClinVar_Feb_2026.txt`, the
   February 2026 germline-classification timepoint used in the master table.
3. **IGVF CVFG annotations** — `data/raw/BAP1_cvfg_annotations.tsv`, supplying the
   2018 and 2025 ClinVar timepoints.

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

## Supplement (final deliverables)

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

Supplementary Table 1 (KURIS/NDD + B/LB missense variants) also carries the
required reference sequences for every HGVS g./c./p. description in the
manuscript (GRCh38 NC_000003.12; RefSeq NM_004656.4/NP_004647.1), in its
footnote — see `table_supp1_kuris_blb_missense.py`. Supplementary Table 2
(`table_supp_figure1_variants.py`) lists every variant plotted in Figure 1,
cited from the main text rather than reproduced there in full.

Each `Supplementary_Figure_N.pdf` is a US-Letter page with the figure at
double-column width and its caption below (on a following page when the figure is
too tall), following Genetics in Medicine conventions. Each
`Supplementary_Table_N.pdf` carries its title and caption above the table. The
titles, captions, and ordering live in `bap1figs/supplement.py`;
`supplement/SUPPLEMENT_GUIDE.md` is the reader's guide that lists every item.

## Notes

- GMM score cutoffs come from one file (`data/source/GMM_thresholds.json`) that
  every script reads, so calibration can never drift between figures.
- Likelihood ratios use the standard Haldane–Anscombe correction (+0.5 to all
  four cells of the 2×2 when a class cell is empty; marked `*`), matching the
  ACMG SVI likelihood-ratio calculator.
- The `SuppFig_LR_master` and `SuppFig_LR_uncertainty_KURIS` forests read the
  contingency / stability **tables** as their input files.
