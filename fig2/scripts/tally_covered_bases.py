#!/usr/bin/env python
# author: Junko Tsuji <jtsuji@broadinstitute.org>

import sys
import os.path

import fileinput
from argparse import ArgumentParser

def main(args):
    LCS_CUTOFF = 20

    targets = {}
    for line in fileinput.input(args.input):
        line = line.rstrip().split("\t")[:-1]

        locus = line[-5]+":"+line[-4]+"-"+line[-3]
        cov = map(int, line[3:-5])

        # ignore all sex chromosomes
        if line[0] == "Y" or line[0] == "X":
            continue

        targets.setdefault(locus, [0 for i in range(len(cov))])
        for i in range(len(cov)):
            if cov[i] >= LCS_CUTOFF:
                targets[locus][i] += 1
    for key in targets.keys():
        print("\t".join([key] + map(str, targets[key])))


if __name__ == "__main__":
    parser = ArgumentParser(description="Count covered bases per-sample")
    parser.add_argument("input", help="intersected per-base coverage")
    args = parser.parse_args()
    try:
        main(args)
    except KeyboardInterrupt: pass
