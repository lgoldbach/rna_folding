#!/usr/bin/env python

import argparse

from rna_folding.utils import dotbracket_to_genotype


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--phenotypes", required=True, type=str, 
                        help="File with phenotypes in dot-bracket notation. "
                        "One per line")
    parser.add_argument("-a", "--alphabet", type=str,
                        help="Which alphabet to convert to, max. 2 letters, "
                        "e.g. 'GC'", default="GC")
    parser.add_argument("-d", "--deterministic", default=True, type=bool,
                        help="Do not assign bases randomly.")
    parser.add_argument("-o", "--output", required=True, type=str, 
                        help="Output file")


    
    print("Parse files ...")
    args = parser.parse_args()

    with open(args.phenotype, "r") as file:
        for line in file:
            gt = dotbracket_to_genotype(line.strip())
            