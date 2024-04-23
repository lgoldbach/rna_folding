#!/usr/bin/env python

import argparse
import numpy as np
from rna_folding.utils import dict_to_gpmap


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-i", "--input", help="input gp_map", type=str)


    args = parser.parse_args()
    
    gp_map = {}

    with open(args.input, "r") as file:
        for ph_i, line in enumerate(file):
            line_ = line.strip().split(" ")
            ph = line_[0]

            for gt in line_[1:]:
                if gt not in gp_map:  # O(1) operation because dict
                    gp_map[gt] = [ph_i]
                else:
                    gp_map[gt].append(ph_i) 
            
    unique_suboptimal_sets = set()
    for i in gp_map:
        unique_suboptimal_sets.add(tuple(gp_map[i]))

    with open(args.output, "w") as outfile:
        for s in unique_suboptimal_sets:
            line = " ".join([str(i) for i in s])
            line += "\n"
            outfile.write(line)
