#!/usr/bin/env python

"""Takes a gp_map file, i.e. a space-delimited file, where the first column
contains genotypes and the second phenotypes and returns a sorted csv file
where the ith row contains the phenotype count of the ith phenotype.

e.g:
input: 
genotype    phenotype (numeric)
AUGCGGG     2
GACGACG     3
AGAGAGU     2
Note: The input file should not have a header, this is just for clarity

output:
0
1
2

"""

import argparse
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map"
                        "file", required=True)
    parser.add_argument("-o", "--output", help="Output file"
                        "file", required=True)

    args = parser.parse_args()

    gp_map = np.genfromtxt(args.file, delimiter=" ")  # read data

    phenotypes = gp_map[:, 1].astype(int)  # get phenotype column

    values, counts = np.unique(phenotypes, return_counts=True)  # count
    stack = np.stack((values, counts), axis=1)

    np.savetxt(args.output, stack.astype(int), delimiter=" ", fmt="%i")  # write to csv
