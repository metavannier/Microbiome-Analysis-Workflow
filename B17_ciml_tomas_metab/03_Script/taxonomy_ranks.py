#!/usr/bin/env python3
# #################################################################
#
#   Summary of the taxonomic levels of a QIIME 2 taxonomy (taxonomy.tsv
#   exported from classify-sklearn): for each level, the rank it holds
#   (d__ domain, k__ kingdom, p__ phylum, ..., s__ species), the
#   number of ASVs assigned at that level and, among them, the ASVs with
#   an informative name: placeholders (uncultured, metagenome,
#   unidentified, Incertae_Sedis, "... bacterium", "... sp.") and names
#   copied from the upper rank (rank propagation) are not counted. The
#   level numbers are those of the taxa bar plot and of
#   'qiime taxa collapse --p-level'.
#
# #################################################################

import argparse
import collections
import csv
import re

RANKS = {"d": "domain", "k": "kingdom", "p": "phylum", "c": "class",
         "o": "order", "f": "family", "g": "genus", "s": "species"}

PLACEHOLDER = re.compile(
    r"uncultured|metagenome|unidentified|unknown|unclassified|incertae_sedis"
    r"|environmental|synthetic|organism|^bacterium|_bacterium$|_sp\.?$"
    r"|_archaeon$", re.IGNORECASE)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--taxonomy", required=True, help="taxonomy.tsv")
    p.add_argument("--output", required=True)
    args = p.parse_args(argv)

    prefixes = collections.defaultdict(collections.Counter)
    assigned = collections.Counter()
    named = collections.Counter()
    n_asv = 0
    with open(args.taxonomy, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) < 2 or row[0].startswith("#"):
                continue
            n_asv += 1
            parent = None
            for level, rank in enumerate(row[1].split(";"), 1):
                rank = rank.strip()
                prefix, _, name = rank.partition("__")
                if name:
                    prefixes[level][prefix] += 1
                    assigned[level] += 1
                    if name != parent and not PLACEHOLDER.search(name):
                        named[level] += 1
                parent = name

    with open(args.output, "w") as out:
        out.write("level\trank\tprefix\tassigned_ASVs\tpercent_ASVs"
                  "\tnamed_ASVs\tpercent_named\n")
        for level in sorted(prefixes):
            prefix = prefixes[level].most_common(1)[0][0]
            out.write(f"{level}\t{RANKS.get(prefix, prefix)}\t{prefix}__\t"
                      f"{assigned[level]}\t{100 * assigned[level] / n_asv:.1f}"
                      f"\t{named[level]}\t{100 * named[level] / n_asv:.1f}\n")
    print(f"{args.output}: {n_asv} ASVs, {len(prefixes)} levels")


if __name__ == "__main__":
    main()
