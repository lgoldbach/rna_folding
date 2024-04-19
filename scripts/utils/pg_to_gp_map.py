#!/usr/bin/env python

import argparse

from rna_folding.parsing import gpmap_to_dict


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, 
                        help="GP map file in g-p format")
    parser.add_argument("-g", "--genotype_file", required=True, type=str,
                        help="List of genotype files in order")
    parser.add_argument("-o", "--out", required=False, type=str,
                        help="Name of output file")

    args = parser.parse_args()

    pg_dict = gpmap_to_dict(args.input, args.genotype_file)

    with open(args.out, "w") as file:
        for gt in sorted(pg_dict.keys()):
            file.write(gt + " " + " ".join(pg_dict[gt])+ "\n")


    