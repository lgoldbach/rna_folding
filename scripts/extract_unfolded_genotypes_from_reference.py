#!/usr/bin/env python

import argparse
import numpy as np
from rna_folding.utils import dict_to_gpmap


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-i", "--input", help="Input gp_map", type=str)
    parser.add_argument("-g", "--genotypes", help="Genotype file", type=str)
    parser.add_argument("-u", "--unfolded", help="The unfolded genotype", type=str, required=False)

    args = parser.parse_args()

    with open(args.genotypes, "r") as gt_file:
            genotypes = np.array([line.strip() for line in gt_file])

    with open(args.input, "r") as gp_map_file:
        for line in gp_map_file:
            line_split = line.strip().split(" ")
            ph = line_split[0]
            if ph == args.unfolded:
                genotype_ids = [gt for gt in line_split[1:]]
                # unfolded_genotypes = genotypes[genotype_ids]
    
    with open(args.output, "w") as file:
        print(genotype_ids[:10], flush=True)
        for genotype in genotype_ids:
            file.write(genotype + "\n")
    