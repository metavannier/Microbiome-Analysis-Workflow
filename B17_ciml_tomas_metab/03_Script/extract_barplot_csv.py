#!/usr/bin/env python3
# #################################################################
#
#   Extract the level-N.csv table stored inside a QIIME 2 taxa bar
#   plot visualization (.qzv). This is exactly the file obtained with
#   "Download CSV" for that taxonomic level in QIIME 2 View.
#
# #################################################################

import argparse
import re
import shutil
import sys
import zipfile


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--qzv", required=True, help="taxa bar plot .qzv")
    p.add_argument("--level", type=int, default=7,
                   help="taxonomic level (default 7 = deepest rank: genus with "
                        "SILVA 144, species with SILVA 138; see "
                        "taxonomy_ranks.tsv)")
    p.add_argument("--output", required=True, help="output .csv")
    args = p.parse_args(argv)

    pattern = re.compile(r"^[^/]+/data/level-(\d+)\.csv$")
    with zipfile.ZipFile(args.qzv) as qzv:
        members = {int(m.group(1)): name for name in qzv.namelist()
                   if (m := pattern.match(name))}
        if args.level not in members:
            sys.exit(f"ERROR: level-{args.level}.csv not found in "
                     f"{args.qzv}; available levels: "
                     f"{sorted(members) or 'none'}")
        with qzv.open(members[args.level]) as src, \
                open(args.output, "wb") as dst:
            shutil.copyfileobj(src, dst)
    print(f"{members[args.level]} -> {args.output}")


if __name__ == "__main__":
    main()
