# Supplementary material — reader's guide

BAP1 saturation genome editing functional evidence for KURIS/NDD variant interpretation. Figures and tables are numbered in the order they are cited; each file is a self-contained US-Letter PDF (figures carry their caption, tables carry their title and caption). Working copies of every table are also provided as `.tsv` and `.xlsx` in `tables/`.

## Supplementary figures

**Supplementary Figure 1. Fitted Gaussians for variant functional classification.**

(a) “Abnormal” (red) and “Normal” (blue) Gaussians used to draw thresholds for classifying variants as functionally abnormal or functionally normal. Estimated density from GMM modeling is on the Y-axis, and SGE score is on the X-axis. Vertical dashed line at X = −0.0286 represents the upper threshold for functionally abnormal variants. Vertical dotted line at X = −0.0197 represents the lower threshold for functionally normal variants. (b) Strip plot of molecular consequence (Y-axis) vs. SGE score (X-axis). (c) Cumulative distribution of functional scores by molecular consequence, with the abnormal (dashed) and normal (dotted) cutoffs marked.

_File:_ `supplement/Supplementary_Figure_1.pdf`

**Supplementary Figure 2. Comparing Waters classifications based on directionality and q-values to the GMM score-cutoff-derived classification.**

(a–c) Functional-score histograms on a split count axis, stratified by (a) Waters functional class, (b) ClinVar classification, and (c) molecular consequence; the abnormal (dashed) and normal (dotted) score cutoffs are overlaid. (d–i) Cross-tabulation of Waters functional class against GMM class for benign control sets (d, synonymous; f, ClinVar benign/likely-benign; h, KURIS benign/likely-benign missense) and pathogenic control sets (e, truncating; g, ClinVar pathogenic/likely-pathogenic; i, KURIS/NDD missense). In (d–i), the count is printed in each cell and the cell shading is proportional to that cell's fraction of the control set.

_File:_ `supplement/Supplementary_Figure_2.pdf`

**Supplementary Figure 3. Likelihood ratios across functional-classification methods and control sets.**

Abnormal-class (red) and normal-class (blue) likelihood ratios with log-Wald 95% confidence intervals for each functional-classification method (GMM score cutoffs, Waters, Tejura) applied to each control set: ClinVar pathogenic/likely-pathogenic vs benign/likely-benign (2025 and February 2026 releases; all variants and missense only), the KURIS/NDD vs benign/likely-benign missense set, and the Waters et al. validation sets. Shaded bands denote ACMG evidence-strength tiers.

_File:_ `supplement/Supplementary_Figure_3.pdf`

**Supplementary Figure 4. Stability of the KURIS/NDD missense likelihood-ratio estimate across estimators.**

Point estimates and intervals for the abnormal-class (red) and normal-class (blue) likelihood ratios computed under three estimands: the Haldane-corrected likelihood ratio (circles), the Tavtigian OddsPath (diamonds), and Bayesian posterior estimates under Beta(0.5, 0.5) and Beta(1, 1) priors (squares). For the likelihood ratio and OddsPath, the intervals shown are the log-Wald 95% confidence interval, the percentile bootstrap range, and the leave-one-out range; for the Bayesian estimand, the 95% credible interval. Shaded bands denote ACMG evidence-strength tiers, and the fraction of resamples falling in the full-data tier is annotated above each marker. Intervals extending beyond the axis are drawn with an arrowhead.

_File:_ `supplement/Supplementary_Figure_4.pdf`

## Supplementary tables

**Supplementary Table 1. Phenotype-specific KURIS/NDD truth set used for assay calibration.**

Missense variants identified in individuals with KURIS/NDD and benign/likely-benign ClinVar missense variants, annotated with the affected protein domain, the SGE functional score and GMM class, ClinVar classification and review status, and associated clinical information. Each variant's clinical occurrences are listed on a labelled second row.

_File:_ `supplement/Supplementary_Table_1.pdf` (data: `tables/SuppTable1_KURIS_NDD_and_BLB_Missense.tsv`, `tables/SuppTable1_KURIS_NDD_and_BLB_Missense.xlsx`)

**Supplementary Table 2: All BAP1 missense variants associated with KURIS/NDD and TPDS, and their ClinVar B/LB controls**

KURIS/NDD missense variants, TPDS-associated ClinVar pathogenic/likely-pathogenic missense variants, and ClinVar benign/likely-benign controls, with HGVS nomenclature, domain, SGE functional scores and class, and ClinVar classifications. The SGE evidence column gives the ACMG code (PS3 Strong / BS3 Moderate) from the KURIS/NDD-vs-B/LB calibration (Supplementary Table 4) for the KURIS/NDD and B/LB rows; it is blank for the TPDS P/LP rows, which are not part of that comparison.

_File:_ `supplement/Supplementary_Table_2.pdf` (data: `tables/Figure1_Variants.tsv`, `tables/Figure1_Variants.xlsx`)

**Supplementary Table 3. Performance of the Waters and GMM functional classifiers.**

Sensitivity, specificity, positive predictive value, and negative predictive value (each with a Wilson 95% confidence interval), together with the true- and false-positive and -negative counts (TP, FP, FN, TN), for the Waters classification (assigned from SGE q-values and directionality) and the GMM score-cutoff classification, against a clinical truth set (ClinVar pathogenic/likely-pathogenic vs benign/likely-benign) and a molecular truth set (truncating vs synonymous variants). All four metrics are computed only over variants each classifier called Abnormal or Normal; GMM's indeterminate calls (Indet. columns) are excluded from every denominator — e.g. GMM's specificity against the clinical truth set is 1,133/1,158, not 1,133/1,171.

_File:_ `supplement/Supplementary_Table_3.pdf` (data: `tables/Waters_vs_GMM_performance.tsv`, `tables/Waters_vs_GMM_performance.xlsx`)

**Supplementary Table 4. Contingency table of functional-classification evidence across all truth sets and class-assignment methods.**

Likelihood ratios and Tavtigian OddsPath values, each with a log-Wald 95% confidence interval and the corresponding ACMG evidence-strength tier, for every functional-classification method, control set, and functional class. Class definitions and assignment methods are given in the footnotes.

_File:_ `supplement/Supplementary_Table_4.pdf` (data: `tables/Master_Contingency_Tables.tsv`, `tables/Master_Contingency_Tables.xlsx`)

**Supplementary Table 5. Likelihood-ratio uncertainty of the functional evidence across estimators.**

Point estimates and intervals for the abnormal- and normal-class likelihood ratios under three estimands — the Haldane-corrected likelihood ratio, the Tavtigian OddsPath, and Bayesian posterior estimates — for each control set. For the likelihood ratio and OddsPath, the bounds are the log-Wald 95% confidence interval, the percentile bootstrap range, and the leave-one-out range; the Bayesian rows give the posterior median and 95% credible interval under Beta(0.5, 0.5) and Beta(1, 1) priors.

_File:_ `supplement/Supplementary_Table_5.pdf` (data: `tables/Supplemental_LR_uncertainty_methods.tsv`, `tables/Supplemental_LR_uncertainty_methods.xlsx`)

**Supplementary Table 6. Phenotypic overlap between the UDN individual and the published features of KURIS described by Küry et al. 2022.**

For each clinical feature reported in the Küry et al. 2022 KURIS cohort (n=11), its frequency in that cohort and whether it was observed in the UDN proband, showing the phenotypic overlap that supports the proband's KURIS diagnosis.

_File:_ `supplement/Supplementary_Table_6.pdf` (data: `tables/SuppTable2_phenotypic_overlap.tsv`, `tables/SuppTable2_phenotypic_overlap.xlsx`)
