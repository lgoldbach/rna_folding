#!/usr/bin/env python

import argparse
import numpy as np
from rna_folding.utils import dict_to_gpmap


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-i", "--input", help="input gp_map", type=str)
    parser.add_argument("-r", "--ranking", help="A file containing phenotypes that defines the ranking", type=str)

    args = parser.parse_args()
    
    flat_gp_map = {}

    phenotypes = [ph.strip() for ph in open(args.ranking, "r")]
    ph_to_rank = dict([(ph, i) for i, ph in enumerate(phenotypes)])

    gp_map = open(args.input, "r")
    for line in gp_map:
        line_ = line.strip().split(" ")
        ph = line_[0]
        rank = ph_to_rank[ph]

        for gt in line_[1:]:
            if gt not in flat_gp_map:  # O(1) operation because dict
                flat_gp_map[gt] = rank  # assign phenotype
            elif rank < flat_gp_map[gt]:
                flat_gp_map[gt] = rank  # update to better ranked ph
    
    # turn into ph_to_gt dict
    ph_to_gt = dict([(ph, []) for ph in phenotypes])
    for gt in flat_gp_map:
        ph = phenotypes[flat_gp_map[gt]]  # change rank to ph
        ph_to_gt[ph].append(gt)

    # save to file
    dict_to_gpmap(ph_to_gt, args.output)
