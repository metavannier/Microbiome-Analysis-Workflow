#!/usr/bin/env python3
# #################################################################
#
#   Summary of the taxonomic levels of a QIIME 2 taxonomy (taxonomy.tsv
#   exported from classify-sklearn): for each level, the rank it holds
#   (d__ domain, k__ kingdom, p__ phylum, ..., s__ species) and the
#   number of ASVs assigned at that level. The level numbers are those
#   of the taxa bar plot and of 'qiime taxa collapse --p-level'.
#
# #################################################################

import argparse
import collections
import csv

RANKS = {"d": "domain", "k": "kingdom", "p": "phylum", "c": "class",
         "o": "order", "f": "family", "g": "genus", "s": "species"}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--taxonomy", required=True, help="taxonomy.tsv")
    p.add_argument("--output", required=True)
    args = p.parse_args(argv)

    prefixes = collections.defaultdict(collections.Counter)
    assigned = collections.Counter()
    n_asv = 0
    with open(args.taxonomy, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 2 or row[0].startswith("#"):
                continue
            n_asv += 1
            for level, rank in enumerate(row[1].split(";"), 1):
                rank = rank.strip()
                prefix, _, name = rank.partition("__")
                if name:
                    prefixes[level][prefix] += 1
                    assigned[level] += 1

    with open(args.output, "w") as out:
        out.write("level\trank\tprefix\tassigned_ASVs\tpercent_ASVs\n")
        for level in sorted(prefixes):
            prefix = prefixes[level].most_common(1)[0][0]
            out.write(f"{level}\t{RANKS.get(prefix, prefix)}\t{prefix}__\t"
                      f"{assigned[level]}\t{100 * assigned[level] / n_asv:.1f}\n")
    print(f"{args.output}: {n_asv} ASVs, {len(prefixes)} levels")


if __name__ == "__main__":
    main()
