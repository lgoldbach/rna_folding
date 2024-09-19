#!/usr/bin/env python

import argparse
import random

from rna_folding.parsing import gpmap_to_lists, dict_to_gpmap


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="g-p map text file ", required=True)
    parser.add_argument("-o", "--output", help="output g-p map text file",
                        required=True)
    
    args = parser.parse_args()

    gt, ph = gpmap_to_lists(args.input)

    random.shuffle(gt)

    # turn into dictionary keyed by phenotype (non-redundant)
    D = {}
    for p, g in zip(ph, gt):
        if p not in D:
            D[p] = [g]
        else:
            D[p].append(g)

    dict_to_gpmap(ph_to_gt=D, file=args.output)  # save to file
