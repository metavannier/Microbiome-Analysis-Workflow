#!/usr/bin/env python3
# #################################################################
#
#   Cross-check of two taxonomies of the same ASVs (QIIME 2
#   taxonomy.tsv exported from classify-sklearn), e.g. the SILVA
#   weighted (habitat-specific) classifier against the SILVA uniform
#   classifier, as recommended by SILVA for the weighted classifiers.
#
#   For each ASV, the two lineages are compared rank by rank down to
#   the deepest rank assigned by both (the uniform classifiers stop at
#   the genus). An ASV whose lineages differ is flagged "different"
#   with the first rank that differs, "one_unassigned" when only one
#   classifier assigns it: check its sequence (BLAST, exact match)
#   before interpreting it.
#
# #################################################################

import argparse
import collections
import csv
import sys

RANKS = {"d": "domain", "k": "kingdom", "p": "phylum", "c": "class",
         "o": "order", "f": "family", "g": "genus", "s": "species"}


def read_taxonomy(path):
    """Feature ID -> (taxon, confidence)"""
    taxa = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        if len(header) < 2:
            sys.exit(f"ERROR: {path} is not a QIIME 2 taxonomy table")
        for row in reader:
            if not row or row[0].startswith("#"):
                continue
            conf = row[2] if len(row) > 2 else ""
            taxa[row[0]] = (row[1].strip(), conf)
    return taxa


def lineage(taxon):
    """'d__Bacteria;k__Bacillati;...' -> [('d', 'Bacteria'), ...]"""
    ranks = []
    for item in taxon.split(";"):
        prefix, sep, name = item.strip().partition("__")
        if not sep or not name:
            break
        ranks.append((prefix, name))
    return ranks


def compare(taxon, check):
    a, b = lineage(taxon), lineage(check)
    n = min(len(a), len(b))
    if n == 0:
        return ("unassigned" if not a and not b else "one_unassigned"), "", ""
    for i in range(n):
        if a[i] != b[i]:
            return "different", RANKS.get(a[i][0], a[i][0]), \
                RANKS.get(a[n - 1][0], a[n - 1][0])
    return "same", "", RANKS.get(a[n - 1][0], a[n - 1][0])


def rank_name(taxon, prefix):
    for p, name in lineage(taxon):
        if p == prefix:
            return name
    return ""


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--taxonomy", required=True,
                   help="taxonomy.tsv of the classifier used in the analysis")
    p.add_argument("--check", required=True,
                   help="taxonomy.tsv of the control classifier")
    p.add_argument("--output", required=True)
    args = p.parse_args(argv)

    taxonomy = read_taxonomy(args.taxonomy)
    check = read_taxonomy(args.check)

    rows = []
    for fid, (taxon, conf) in taxonomy.items():
        if fid in check:
            ctaxon, cconf = check[fid]
            status, first_diff, deepest = compare(taxon, ctaxon)
        else:
            ctaxon, cconf, status, first_diff, deepest = "", "", "missing", "", ""
        rows.append([fid, status, first_diff, deepest,
                     rank_name(taxon, "g"), rank_name(ctaxon, "g"),
                     taxon, conf, ctaxon, cconf])

    order = {"different": 0, "one_unassigned": 1, "missing": 2,
             "unassigned": 3, "same": 4}
    rows.sort(key=lambda r: (order[r[1]], r[0]))
    with open(args.output, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(["Feature ID", "agreement", "first_different_rank",
                         "deepest_compared_rank", "genus", "check_genus",
                         "Taxon", "Confidence", "check_Taxon",
                         "check_Confidence"])
        writer.writerows(rows)

    counts = collections.Counter(r[1] for r in rows)
    print(f"{args.output}: {len(rows)} ASVs - "
          + ", ".join(f"{k}: {counts[k]}" for k in order if counts[k]))
    pairs = collections.Counter((r[4], r[5]) for r in rows
                                if r[1] == "different" and r[2] in
                                ("genus", "family", "order", "class",
                                 "phylum", "kingdom", "domain"))
    for (g1, g2), n in pairs.most_common(20):
        print(f"  {n} ASVs: '{g1 or '-'}' vs '{g2 or '-'}' (check)")


if __name__ == "__main__":
    main()
