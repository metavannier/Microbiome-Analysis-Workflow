#!/usr/bin/env python3
# #################################################################
#
#   Format a QIIME 2 taxa bar plot table (level-N.csv, e.g. the
#   "Download CSV" of taxonomic level 7 in QIIME 2 View) into a
#   MetaboAnalyst-ready table:
#
#       Sample,Group,<taxon 1>,<taxon 2>,...
#       M6995,CTL,0.0706973,0.0414641,...
#
#   - one row per sample, one column per taxon (samples in rows);
#   - the Group of each sample is read from sample-metadata.tsv
#     (--metadata + --group-column), or, for backward compatibility,
#     from the second column of an existing MetaboAnalyst table
#     (--groups-from, e.g. "Ileum-cleaned-F vs M.csv");
#   - counts are converted to per-sample relative abundances
#     (the denominator is the total of ALL taxa of the sample);
#   - taxa absent from every selected sample are removed.
#
#   The script only needs pandas: it runs in the QIIME 2 container
#   (workflow) or in any Python 3 environment (manual use).
#
# #################################################################

import argparse
import csv
import sys

import pandas as pd


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Format a QIIME 2 taxa bar plot CSV (level-N.csv) "
                    "for MetaboAnalyst.")
    p.add_argument("--taxa-csv", required=True,
                   help="level-N.csv downloaded from QIIME 2 View (or "
                        "extracted from the taxa bar plot .qzv).")
    p.add_argument("--output", required=True,
                   help="Output MetaboAnalyst table (.csv).")
    grp = p.add_argument_group("group assignment (choose one)")
    grp.add_argument("--metadata",
                     help="QIIME 2 sample metadata (sample-metadata.tsv).")
    grp.add_argument("--group-column",
                     help="Metadata column giving the group of each sample "
                          "(e.g. Diet, Sexe, Name-experiment).")
    grp.add_argument("--groups-from",
                     help="Existing MetaboAnalyst table whose first two "
                          "columns are Sample and Group (old procedure).")
    sel = p.add_argument_group("sample selection")
    sel.add_argument("--samples",
                     help="TSV/CSV whose first column lists the sample IDs to "
                          "keep (e.g. the comparison metadata written by the "
                          "workflow).")
    sel.add_argument("--subset", action="append", default=[],
                     metavar="COLUMN=VALUE[,VALUE...]",
                     help="Keep samples whose metadata COLUMN equals one of "
                          "the VALUEs. Can be repeated; several conditions "
                          "may also be separated by ';'.")
    sel.add_argument("--exclude-samples", default="",
                     help="Comma-separated sample IDs to drop.")
    out = p.add_argument_group("output options")
    out.add_argument("--counts", action="store_true",
                     help="Keep raw counts instead of relative abundances.")
    out.add_argument("--keep-zero-taxa", action="store_true",
                     help="Keep taxa absent from all the selected samples.")
    out.add_argument("--min-group-prevalence", type=int, default=0,
                     help="Keep only taxa present (>0) in at least N samples "
                          "of at least one group (0 = no filter).")
    out.add_argument("--sep", default=",",
                     help="Field separator of the output (default ','). "
                          "MetaboAnalyst expects comma-separated files.")
    out.add_argument("--decimal", default=".",
                     help="Decimal mark of the output (default '.').")
    out.add_argument("--float-format", default="%.10g",
                     help="printf-style format of the values.")
    out.add_argument("--sample-label", default="Sample",
                     help="Header of the sample column (default 'Sample').")
    out.add_argument("--group-label", default="Group",
                     help="Header of the group column (default 'Group').")
    args = p.parse_args(argv)

    if args.groups_from is None and (args.metadata is None
                                     or args.group_column is None):
        p.error("give either --metadata and --group-column, "
                "or --groups-from")
    if args.subset and args.metadata is None:
        p.error("--subset requires --metadata")
    if args.sep == args.decimal:
        p.error("--sep and --decimal must be different")
    return args


def read_metadata(path):
    """Read a QIIME 2 metadata TSV as strings, indexed by sample ID."""
    md = pd.read_csv(path, sep="\t", dtype=str, comment=None,
                     skip_blank_lines=True, keep_default_na=False)
    # Drop the optional QIIME 2 '#q2:types' directive row.
    first = md.columns[0]
    md = md[~md[first].str.startswith("#")]
    md = md[md[first].str.strip() != ""]
    md[first] = md[first].str.strip()
    return md.set_index(first)


def parse_subset(items):
    """['Body-site=Ileal;Sexe=Male', 'Diet=CTL,SAM'] -> {col: [values]}"""
    conditions = {}
    for item in items:
        for cond in filter(None, (c.strip() for c in item.split(";"))):
            if "=" not in cond:
                sys.exit(f"ERROR: malformed subset condition '{cond}' "
                         "(expected COLUMN=VALUE[,VALUE])")
            col, values = cond.split("=", 1)
            conditions[col.strip()] = [v.strip() for v in values.split(",")]
    return conditions


def read_groups_from(path):
    """Sample -> Group from the first two columns of an old table."""
    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
    dialect = csv.Sniffer().sniff(sample, delimiters=",;\t")
    old = pd.read_csv(path, sep=dialect.delimiter, dtype=str,
                      usecols=[0, 1], encoding="utf-8-sig",
                      keep_default_na=False)
    old.columns = ["Sample", "Group"]
    old = old.apply(lambda col: col.str.strip())
    if old["Sample"].duplicated().any():
        sys.exit(f"ERROR: duplicated samples in {path}: "
                 + ", ".join(old["Sample"][old["Sample"].duplicated()]))
    empty = old["Group"] == ""
    if empty.any():
        print("WARNING: no group in " + path + " for samples (dropped): "
              + ", ".join(old["Sample"][empty]), file=sys.stderr)
    return old[~empty].set_index("Sample")["Group"]


def read_sample_list(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        sample = fh.read(4096)
    sep = "\t" if "\t" in sample.splitlines()[0] else ","
    ids = pd.read_csv(path, sep=sep, dtype=str, usecols=[0],
                      keep_default_na=False).iloc[:, 0].str.strip()
    return set(i for i in ids if i and not i.startswith("#"))


def taxa_columns(table, metadata_cols):
    """Taxa columns of a barplot CSV (QIIME 2 appends the metadata columns
    after the taxa). At level N every taxon has exactly N ranks, so the
    taxa are the numeric columns, not in the metadata, with the same
    number of ';' as the first taxon column."""
    first = table.columns[0]
    n_sep = first.count(";")
    cols, skipped = [], []
    for col in table.columns:
        if col in metadata_cols:
            continue
        if col.count(";") != n_sep or (
                n_sep == 0 and not col.startswith(("Unassigned", "d__", "k__"))):
            skipped.append(col)
            continue
        values = pd.to_numeric(table[col], errors="coerce")
        if not values.notna().all():
            sys.exit(f"ERROR: taxon column '{col}' contains non-numeric "
                     "values")
        cols.append(col)
    if skipped:
        print("Metadata columns of the taxa table (not taxa): "
              + ", ".join(skipped))
    return cols


def main(argv=None):
    args = parse_args(argv)

    table = pd.read_csv(args.taxa_csv, dtype=str, keep_default_na=False)
    id_col = table.columns[0]
    table[id_col] = table[id_col].str.strip()
    table = table.set_index(id_col)
    if table.index.duplicated().any():
        sys.exit("ERROR: duplicated sample IDs in " + args.taxa_csv)

    metadata = read_metadata(args.metadata) if args.metadata else None
    metadata_cols = set(metadata.columns) if metadata is not None else set()
    taxa = taxa_columns(table, metadata_cols)
    if not taxa:
        sys.exit("ERROR: no taxon column found in " + args.taxa_csv)
    counts = table[taxa].apply(pd.to_numeric)

    # ---- group of each sample ----------------------------------------
    if args.groups_from:
        groups = read_groups_from(args.groups_from)
    else:
        if args.group_column not in metadata.columns:
            sys.exit(f"ERROR: column '{args.group_column}' not found in "
                     f"{args.metadata}. Available: "
                     f"{', '.join(metadata.columns)}")
        groups = metadata[args.group_column].str.strip()
        groups = groups[groups != ""]

    # ---- sample selection ----------------------------------------------
    keep = pd.Series(True, index=counts.index)
    if metadata is not None:
        # Samples absent from the metadata are outside the analysis.
        keep &= counts.index.isin(metadata.index)
        absent = [s for s in metadata.index if s not in counts.index]
        if absent:
            print(f"{len(absent)} samples of {args.metadata} are not in the "
                  "taxa table (e.g. below the sampling depth): "
                  + ", ".join(absent))
    for col, values in parse_subset(args.subset).items():
        if col not in metadata.columns:
            sys.exit(f"ERROR: subset column '{col}' not in metadata")
        allowed = set(metadata.index[metadata[col].str.strip().isin(values)])
        keep &= counts.index.isin(allowed)
    if args.samples:
        keep &= counts.index.isin(read_sample_list(args.samples))
    excluded = {s.strip() for s in args.exclude_samples.split(",")
                if s.strip()}
    keep &= ~counts.index.isin(excluded)
    no_group = keep & ~counts.index.isin(groups.index)
    if no_group.any() and args.groups_from:
        # Old procedure: the --groups-from table also defines the samples.
        print(f"{int(no_group.sum())} samples absent from "
              f"{args.groups_from} were left out")
    elif no_group.any():
        print(f"WARNING: no '{args.group_column}' value for samples "
              "(dropped): " + ", ".join(counts.index[no_group]),
              file=sys.stderr)
    keep &= counts.index.isin(groups.index)
    if not keep.any():
        sys.exit("ERROR: no sample left after selection")

    totals = counts.sum(axis=1)
    empty = keep & (totals <= 0)
    if empty.any():
        print("WARNING: samples without any read (dropped): "
              + ", ".join(counts.index[empty]), file=sys.stderr)
        keep &= ~empty

    # Relative abundance over ALL taxa of each sample, computed before
    # any taxon is removed (identical to the previous MetaboAnalyst files).
    values = counts if args.counts else counts.div(totals, axis=0)
    values = values.loc[keep]
    sample_groups = groups.reindex(values.index)

    # ---- taxa filtering --------------------------------------------------
    n_taxa = values.shape[1]
    if not args.keep_zero_taxa:
        values = values.loc[:, (values > 0).any(axis=0)]
    if args.min_group_prevalence > 0:
        present = (values > 0).groupby(sample_groups).sum()
        values = values.loc[:, (present >= args.min_group_prevalence)
                            .any(axis=0)]

    # ---- write -------------------------------------------------------------
    out = values.copy()
    out.insert(0, args.group_label, sample_groups)
    out.index.name = args.sample_label
    out.to_csv(args.output, sep=args.sep, decimal=args.decimal,
               float_format=args.float_format, quoting=csv.QUOTE_MINIMAL)

    print(f"{args.output}: {out.shape[0]} samples x {values.shape[1]} taxa "
          f"({n_taxa - values.shape[1]} taxa removed)")
    for group, n in sample_groups.value_counts().sort_index().items():
        print(f"  {group}: {n} samples")
    if sample_groups.nunique() < 2:
        print("WARNING: fewer than 2 groups - MetaboAnalyst needs at least "
              "two groups", file=sys.stderr)


if __name__ == "__main__":
    main()
