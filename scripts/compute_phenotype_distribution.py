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
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map in "
                        "ph-to-gt format, i.e. assumes a csv file where each "
                        "phenotype is followed by all the genotypes that map "
                        "to it, separated by spaces. Assumes one-to-one "
                        "mapping. "
                        "e.g.:   ()() 10 4 2 1 "
                        "        .... 3 5 6 "
                        "        (..) 7 8 9 ")
    parser.add_argument("-o", "--output", help="Output file"
                        "file", required=True)

    args = parser.parse_args()


    with open(args.file, "r") as file:
        with open(args.output, "w") as outfile:
            for line_ in file:
                line = line_.split(" ")
                # count number of genotypes that map to this phenotype
                # -1 to exclude the phenotype at the start of the line
                count = len(line) - 1
                outfile.write(line[0] + " " + str(count) + "\n")
    file.close()
    outfile.close()
        