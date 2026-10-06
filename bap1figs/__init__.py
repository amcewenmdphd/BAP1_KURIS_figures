"""bap1figs — shared library for the BAP1 KURIS figure/table scripts.

Modules:
    config    paths + matplotlib defaults
    palettes  color schemes + ACMG evidence bands
    stats     GMM classification, likelihood ratios, ACMG tiers, Wilson CIs
    tables    master/clinical loaders + panel3_category classifier
    plots     panel labels, contingency heatmaps, broken-axis histograms, ACMG forests
    tablepdf  typeset PDF renderer shared by the supplementary-table scripts
"""

from . import config, palettes, stats, tables, plots, tablepdf  # noqa: F401
