#!/usr/bin/env python3
# #################################################################
#
#   Table of the abundance after taxa collapse and feature filtering
#   per group, formatted as LEfSe input (port of the former lefse.R):
#
#       index   M6969   M6970   ...      <- sample IDs  (LEfSe -u 1)
#       Diet    CTL     SAM     ...      <- class       (LEfSe -c 2)
#       d__Bacteria|p__...|s__...  12.0  0.0  ...
#
#   Input: the TSV written by 'biom convert --to-tsv' from the filtered
#   collapsed table ('# Constructed from biom file' + '#OTU ID' header).
#
# #################################################################

import argparse
import sys

import pandas as pd


def main(argv=None):
    p = argparse.ArgumentParser(description="Format a table for LEfSe.")
    p.add_argument("--table-tsv", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--group-column", required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args(argv)

    table = pd.read_csv(args.table_tsv, sep="\t", skiprows=1, dtype=str)
    table = table.rename(columns={table.columns[0]: "index"})
    # Drop the empty observation-metadata column added by --header-key.
    table = table.drop(columns=[c for c in table.columns
                                if c.lower() == "taxonomy"])
    samples = list(table.columns[1:])

    md = pd.read_csv(args.metadata, sep="\t", dtype=str,
                     keep_default_na=False)
    md = md[~md[md.columns[0]].str.startswith("#")]
    md = md.set_index(md.columns[0])
    if args.group_column not in md.columns:
        sys.exit(f"ERROR: column '{args.group_column}' not in metadata")
    missing = [s for s in samples if s not in md.index]
    if missing:
        sys.exit("ERROR: samples missing from metadata: " + ", ".join(missing))

    class_row = pd.DataFrame(
        [[args.group_column] + md.loc[samples, args.group_column].tolist()],
        columns=table.columns)
    table["index"] = table["index"].str.replace('"', "") \
        .str.replace(";", "|")
    header = pd.DataFrame([table.columns.tolist()], columns=table.columns)
    lefse = pd.concat([header, class_row, table], ignore_index=True)
    lefse.to_csv(args.output, sep="\t", header=False, index=False)
    print(f"{args.output}: {len(samples)} samples, {len(table)} taxa")


if __name__ == "__main__":
    main()
