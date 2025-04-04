#!/usr/bin/env python

import argparse
import RNA
from rna_folding.parsing import dict_to_gpmap


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-g", "--genotypes", help="List of genotypes", type=str)
    parser.add_argument("-i", "--input", help="input gp_map", type=str)

    args = parser.parse_args()
    
    flat_gp_map = {} 
    flat_ge_map = {}

    genotypes = []
    with open(args.genotypes, "r") as g:
        for line in g:
            genotypes.append(line.strip())

    ## read in genotypes
    with open(args.input, "r") as gp_map: 
        for line in gp_map:
            line_ = line.strip().split(" ")
            ph = line_[0]
            for gt_id in line_[1:]:
                gt = genotypes[int(gt_id)]

                fc = RNA.fold_compound(gt)  # generate fold_compound object
                fe = fc.eval_structure(ph)  # compute free energy
                if gt not in flat_gp_map:  # O(1) operation because dict
                    flat_gp_map[gt] = ph  # assign phenotype
                    flat_ge_map[gt] = fe
                elif fe < flat_ge_map[gt]:
                    flat_gp_map[gt] = ph  # update to lower free energy ph
                    flat_ge_map[gt] = fe
                elif fe == flat_ge_map[gt]:
                    print("Same fe:", gt, ph, fe, flat_gp_map[gt], flat_ge_map[gt], flush=True)

    # turn into ph_to_gt dict
    ph_to_gt = dict([(ph, []) for ph in flat_gp_map.values()])
    for gt in flat_gp_map:
        ph_to_gt[ph].append(gt)

    # save to file
    dict_to_gpmap(ph_to_gt, args.output)
