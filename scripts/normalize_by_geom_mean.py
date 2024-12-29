#!/usr/bin/env python

import argparse
from scipy.stats.mstats import gmean
from rna_folding.utils import combinatorically_complete_genotypes


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file with scores to be "
                        "normalized, space-seperated file with phenotypes in "
                        "first column and scores in second colum", type=str)
    parser.add_argument("-o", "--output", help="Output file", type=str)

    args = parser.parse_args()

    phenotypes = []
    scores = []

    with open(args.input, "r") as f:
        for raw_line in f:
            line = raw_line.strip().split(" ")

            phenotypes.append(line[0])
            scores.append(float(line[1]))

    gm = gmean(scores)
    print(gm, scores)
    scores_norm = [score/gm for score in scores]

    with open(args.output, "w") as f:
        for p, s in zip(phenotypes, scores_norm):
            f.write(p + " " + str(s) + "\n")
