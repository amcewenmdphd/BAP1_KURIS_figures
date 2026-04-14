# BAP1 KURIS Figure Generation

Standalone script to reproduce Figures 1 and 2 from pre-computed data tables.

## Requirements

- Python 3.8+
- pandas
- numpy
- matplotlib
- openpyxl

Install dependencies:

```bash
pip install pandas numpy matplotlib openpyxl
```

## Input files

| File | Description |
|------|-------------|
| `data/output/BAP1_KURIS_Figure_Variants.xlsx` | Source data for all figure panels (8 sheets) |
| `data/input/BAP1_features.json` | UniProt BAP1 domain annotations (Q92560) |

## Usage

From the repository root:

```bash
python scripts/make_figures.py
```

## Output

| File | Description |
|------|-------------|
| `figures/Figure1.png` | Figure 1 (300 dpi PNG) |
| `figures/Figure1.pdf` | Figure 1 (vector PDF) |
| `figures/Figure2.png` | Figure 2 (300 dpi PNG) |
| `figures/Figure2.pdf` | Figure 2 (vector PDF) |

## Figure panels

**Figure 1** (180 x 120 mm):
- **a** BAP1 protein lollipop plot showing missense variant positions by category
- **b** ClinVar classified variants by molecular consequence
- **c** P/LP variants by curated condition (Cancer predisposition, Kury-Isidor syndrome)

**Figure 2** (180 x 150 mm):
- **a** Functional scores (SGE) by molecular consequence (all 18,108 variants)
- **b** Functional scores for ClinVar P/LP and B/LB variants
- **c** Functional scores for KURIS/NDD and TPDS missense variants
- **d** Contingency table: all ClinVar P/LP vs B/LB variants
- **e** Contingency table: KURIS/NDD vs B/LB missense variants

Score classification thresholds: Abnormal (score <= -0.026891), Normal (score >= -0.018367).
Likelihood ratios use Haldane-Anscombe correction (+ 0.5 to all cells); corrected values marked with *.
