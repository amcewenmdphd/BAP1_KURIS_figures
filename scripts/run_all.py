#!/usr/bin/env python3
"""Run the whole pipeline in dependency order: raw inputs -> source tables ->
per-figure inputs -> figures + tables."""

import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

# (raw -> source), then tables (some feed figures), then figures.
ORDER = [
    "build_master_table.py",
    "build_clinical_intermediates.py",
    "table_master_contingency.py",
    "table_lr_uncertainty.py",
    "table_waters_vs_gmm_performance.py",
    "table_supp1_kuris_blb_missense.py",
    "table_supp_figure1_variants.py",
    "table_supp2_phenotypic_overlap.py",
    "figure1.py",
    "figure2.py",
    "suppfig1_gmm.py",
    "suppfig_lr_master.py",
    "suppfig_lr_uncertainty_kuris.py",
    "suppfig_waters_vs_gmm.py",
    "build_supplement.py",          # numbered, captioned deliverables + guide
]

for script in ORDER:
    print(f"\n=== {script} ===")
    runpy.run_path(str(HERE / script), run_name="__main__")
print("\nAll figures and tables regenerated.")
