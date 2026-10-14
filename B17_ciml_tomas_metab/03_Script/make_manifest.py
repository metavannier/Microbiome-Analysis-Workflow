#!/usr/bin/env python3
# #################################################################
#
#   Build manifest.tsv (QIIME 2 PairedEndFastqManifestPhred33V2) from
#   the paired FASTQ files of 00_RowData, e.g.
#
#     M6969_S2_L001_R1_001.fastq.gz / M6969_S2_L001_R2_001.fastq.gz
#
#   gives
#
#     sample-id  forward-absolute-filepath                                   reverse-absolute-filepath
#     M6969      $PWD/05_Output/00_trimmed/M6969_S2_L001_R1_001.fastq.gz     $PWD/05_Output/00_trimmed/M6969_S2_L001_R2_001.fastq.gz
#
#   The paths point to the reads trimmed by Trimmomatic; '$PWD' is
#   expanded by QIIME 2 at import time (the workflow runs from the
#   project directory). The sample ID is the file name part before the
#   first '_' (M6969).
#
#   Usage (from the project directory):
#       python3 03_Script/make_manifest.py --metadata sample-metadata.tsv
#
#   With --expand MANIFEST, the script instead copies an existing
#   manifest with '$PWD' (and any environment variable) replaced by its
#   value (used by the workflow before the QIIME 2 import).
#
# #################################################################

import argparse
import glob
import os
import sys


def expand(manifest, output):
    with open(manifest, encoding="utf-8-sig") as src, \
            open(output, "w") as dst:
        for i, line in enumerate(src):
            fields = line.rstrip("\r\n").split("\t")
            if i > 0 and len(fields) >= 3 and not line.startswith("#"):
                fields[1:3] = [os.path.abspath(os.path.expandvars(f.strip()))
                               for f in fields[1:3]]
            dst.write("\t".join(fields) + "\n")
    print(f"{manifest} -> {output}")


def main(argv=None):
    p = argparse.ArgumentParser(description="Build the QIIME 2 manifest.")
    p.add_argument("--expand", metavar="MANIFEST",
                   help="expand the variables of an existing manifest")
    p.add_argument("--raw-dir", default="00_RowData")
    p.add_argument("--trimmed-dir", default="$PWD/05_Output/00_trimmed",
                   help="directory written in the manifest "
                        "(default: $PWD/05_Output/00_trimmed)")
    p.add_argument("--r1-suffix", default="_R1_001.fastq.gz")
    p.add_argument("--r2-suffix", default="_R2_001.fastq.gz")
    p.add_argument("--metadata",
                   help="sample-metadata.tsv: warn about samples missing on "
                        "either side and keep only samples it lists")
    p.add_argument("--output", default="manifest.tsv")
    args = p.parse_args(argv)

    if args.expand:
        expand(args.expand, args.output)
        return

    r1_files = sorted(glob.glob(os.path.join(args.raw_dir,
                                             "*" + args.r1_suffix)))
    if not r1_files:
        sys.exit(f"ERROR: no '*{args.r1_suffix}' file in {args.raw_dir}")

    rows = {}
    for r1 in r1_files:
        name1 = os.path.basename(r1)
        name2 = name1[: -len(args.r1_suffix)] + args.r2_suffix
        if not os.path.exists(os.path.join(args.raw_dir, name2)):
            sys.exit(f"ERROR: mate file {name2} of {name1} is missing")
        sample = name1.split("_")[0]
        if sample in rows:
            sys.exit(f"ERROR: sample {sample} has several R1 files "
                     f"({rows[sample][0]} and {name1})")
        rows[sample] = (name1, name2)

    if args.metadata:
        with open(args.metadata, encoding="utf-8-sig") as fh:
            lines = [line.rstrip("\n").split("\t") for line in fh]
        known = [f[0].strip() for f in lines[1:]
                 if f and f[0].strip() and not f[0].startswith("#")]
        no_fastq = [s for s in known if s not in rows]
        no_meta = sorted(set(rows) - set(known))
        if no_fastq:
            print("WARNING: in metadata but without FASTQ: "
                  + ", ".join(no_fastq), file=sys.stderr)
        if no_meta:
            print("WARNING: FASTQ without metadata (not written): "
                  + ", ".join(no_meta), file=sys.stderr)
        rows = {s: rows[s] for s in known if s in rows}

    with open(args.output, "w") as out:
        out.write("sample-id\tforward-absolute-filepath\t"
                  "reverse-absolute-filepath\n")
        for sample, (name1, name2) in rows.items():
            out.write(f"{sample}\t{args.trimmed_dir}/{name1}\t"
                      f"{args.trimmed_dir}/{name2}\n")
    print(f"{args.output}: {len(rows)} samples")


if __name__ == "__main__":
    main()
