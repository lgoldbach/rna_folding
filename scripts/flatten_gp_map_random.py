#!/usr/bin/env python

import argparse
import numpy as np
from rna_folding.utils import dict_to_gpmap
from rna_folding.parsing import gpmap_to_dict


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-i", "--input", help="input gp_map", type=str)
    parser.add_argument("-g", "--genotypes", help="genotype file", type=str)

    args = parser.parse_args()

    gp_map = gpmap_to_dict(gpmap_file=args.input)

    # flatten
    for gt in gp_map:
        gp_map[gt] = np.random.choice(gp_map[gt])
    
    # turn into ph_to_gt dict
    ph_to_gt = {}
    for gt in gp_map:
        ph = gp_map[gt]
        if ph in ph_to_gt:
            ph_to_gt[ph].append(gt)
        else:
            ph_to_gt[ph] = [gt]

    # save to file
    dict_to_gpmap(ph_to_gt, args.output)
