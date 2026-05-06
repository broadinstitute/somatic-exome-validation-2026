#!/usr/bin/env python
# Junko Tsuji <jtsuji@broadinstitute.org>

import sys
import os.path
import fileinput

from argparse import ArgumentParser

def main(args):
    prev = []
    for x in fileinput.input(args.bed):
        x = x.rstrip().split("\t")
        if prev == []:
            prev = x[:]
            continue
        if prev[3] == x[3] and int(prev[2]) >= int(x[1]):
            prev[2] = str(max(int(x[2]), int(prev[2])))
        else:
            print("\t".join(prev))
            prev = x[:]
    print("\t".join(prev))


if __name__ == "__main__":
    parser = ArgumentParser(description="merge overlapping entries of a given gene")
    parser.add_argument("bed", help="input BED")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass

