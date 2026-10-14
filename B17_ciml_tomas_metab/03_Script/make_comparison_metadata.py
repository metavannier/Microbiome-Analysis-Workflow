#!/usr/bin/env python3
# #################################################################
#
#   Write the sample metadata of one comparison: the rows of
#   sample-metadata.tsv matching the comparison subset (e.g.
#   "Body-site=Feces;Sexe=Male"), minus the excluded samples.
#   The file keeps the QIIME 2 metadata format (same header, optional
#   '#q2:types' row) and is used to subset the feature tables and as
#   --m-metadata-file for every analysis of the comparison.
#
# #################################################################

import argparse
import sys

import pandas as pd


def parse_subset(text):
    conditions = {}
    for cond in filter(None, (c.strip() for c in text.split(";"))):
        if "=" not in cond:
            sys.exit(f"ERROR: malformed subset condition '{cond}'")
        col, values = cond.split("=", 1)
        conditions[col.strip()] = [v.strip() for v in values.split(",")]
    return conditions


def main(argv=None):
    p = argparse.ArgumentParser(description="Write a comparison metadata.")
    p.add_argument("--metadata", required=True)
    p.add_argument("--subset", default="",
                   help="'COL=VAL[,VAL];COL2=VAL' (empty = all samples)")
    p.add_argument("--group-columns", required=True,
                   help="comma-separated columns compared (checked)")
    p.add_argument("--exclude-samples", default="")
    p.add_argument("--output", required=True)
    args = p.parse_args(argv)

    md = pd.read_csv(args.metadata, sep="\t", dtype=str,
                     keep_default_na=False, skip_blank_lines=True)
    id_col = md.columns[0]
    md[id_col] = md[id_col].str.strip()
    directives = md[md[id_col].str.startswith("#")]
    samples = md[(~md[id_col].str.startswith("#")) & (md[id_col] != "")]

    keep = pd.Series(True, index=samples.index)
    for col, values in parse_subset(args.subset).items():
        if col not in samples.columns:
            sys.exit(f"ERROR: subset column '{col}' not in {args.metadata}")
        keep &= samples[col].str.strip().isin(values)
    excluded = {s.strip() for s in args.exclude_samples.split(",")
                if s.strip()}
    keep &= ~samples[id_col].isin(excluded)
    subset = samples[keep]
    if subset.empty:
        sys.exit(f"ERROR: no sample matches '{args.subset}'")

    for col in filter(None, (c.strip() for c in
                             args.group_columns.split(","))):
        if col not in subset.columns:
            sys.exit(f"ERROR: group column '{col}' not in {args.metadata}")
        counts = subset[col].str.strip().replace("", pd.NA).dropna() \
            .value_counts().sort_index()
        if len(counts) < 2:
            sys.exit(f"ERROR: column '{col}' has fewer than 2 groups in "
                     f"the subset '{args.subset}': {dict(counts)}")
        print(f"{col}: " + ", ".join(f"{g}={n}" for g, n in counts.items()))

    pd.concat([directives, subset]).to_csv(args.output, sep="\t",
                                           index=False)
    print(f"{args.output}: {len(subset)} samples")


if __name__ == "__main__":
    main()
