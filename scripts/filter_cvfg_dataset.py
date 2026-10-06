#!/usr/bin/env python3
"""
Filter the IGVF/CVFG integrated variant-effect dataset down to a single gene.

Source: "integrated_variant_effect_dataset.tsv.gz" (Supplementary Data 1),
        https://zenodo.org/records/18637474
Output: data/raw/<GENE>_integrated_variant_effect.tsv   (BAP1 by default)

The integrated dataset spans many genes and is large, so this streams the
gzip line by line and keeps only rows for the requested gene; it never loads the
whole file into memory. The Zenodo host is not reachable from the analysis
environment, so run this locally on the downloaded file, then commit the small
per-gene output:

    python scripts/filter_cvfg_dataset.py \\
        --input ~/Downloads/integrated_variant_effect_dataset.tsv.gz \\
        --gene BAP1

Options
    --input PATH      path to the downloaded .tsv.gz (or plain .tsv)   [required]
    --gene SYMBOL     gene symbol to keep                               [BAP1]
    --gene-col NAME   gene column name                                  [auto-detect]
    --output PATH     output path  [data/raw/<GENE>_integrated_variant_effect.tsv]
    --list-columns    just print the header and exit (use to find the gene column)
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

# Common names for the gene column across IGVF/MaveDB-style releases.
GENE_COL_CANDIDATES = [
    "gene", "gene_symbol", "genesymbol", "symbol", "hgnc_symbol", "hgnc",
    "target", "target_gene", "gene_name", "targetgene",
]


def _open(path):
    path = str(path)
    return gzip.open(path, "rt", newline="") if path.endswith(".gz") else open(path, "rt", newline="")


def detect_gene_col(header):
    lower = {h.lower(): h for h in header}
    for cand in GENE_COL_CANDIDATES:
        if cand in lower:
            return lower[cand]
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--gene", default="BAP1")
    ap.add_argument("--gene-col", default=None)
    ap.add_argument("--output", default=None)
    ap.add_argument("--list-columns", action="store_true")
    args = ap.parse_args()

    inp = Path(args.input).expanduser()
    if not inp.exists():
        sys.exit(f"ERROR: input not found: {inp}")

    repo = Path(__file__).resolve().parent.parent
    out = (Path(args.output).expanduser() if args.output
           else repo / "data" / "raw" / f"{args.gene}_integrated_variant_effect.tsv")

    with _open(inp) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)

        if args.list_columns:
            for i, h in enumerate(header, 1):
                print(f"{i:3d}  {h}")
            return

        gene_col = args.gene_col or detect_gene_col(header)
        if gene_col is None:
            sys.exit("ERROR: could not auto-detect a gene column.\n"
                     f"  Columns: {header}\n"
                     "  Re-run with --gene-col <name> (or --list-columns to inspect).")
        if gene_col not in header:
            sys.exit(f"ERROR: gene column {gene_col!r} not in header: {header}")
        gi = header.index(gene_col)

        out.parent.mkdir(parents=True, exist_ok=True)
        kept = scanned = 0
        with open(out, "w", newline="") as ofh:
            w = csv.writer(ofh, delimiter="\t")
            w.writerow(header)
            for row in reader:
                scanned += 1
                if len(row) > gi and row[gi] == args.gene:
                    w.writerow(row)
                    kept += 1

    print(f"gene column : {gene_col}")
    print(f"scanned     : {scanned:,} rows")
    print(f"kept        : {kept:,} rows for gene {args.gene!r}")
    print(f"written     : {out}")
    if kept == 0:
        print(f"\nWARNING: 0 rows matched {args.gene!r} in column {gene_col!r}. "
              "Check the exact gene symbol used in the file "
              "(run with --list-columns, and inspect a few values of that column).")


if __name__ == "__main__":
    main()
