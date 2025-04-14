#!/usr/bin/env python

import argparse
import RNA
from rna_folding.parsing import dict_to_gpmap
from rna_folding.utils import is_compatible
from rna_folding.base_pairing import BasePairing
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-g", "--genotypes", help="List of genotypes", type=str)
    parser.add_argument("-u", "--unfolded", help="Unfolded phenotype", type=str)
    parser.add_argument("-i", "--input", help="input gp_map", type=str)

    args = parser.parse_args()
    
    flat_gp_map = {} 
    flat_ge_map = {}

    bp = BasePairing(bases="UCAG", graph_path="/home/lgold/phd/research/projects/connectivity/rna_folding/data/graphs/", id=7)

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
                if fe >= 0:  # ignore free en. >= 0, because those are unfold.
                    continue
                if gt not in flat_gp_map:  # O(1) operation because dict
                    flat_gp_map[gt] = ph  # assign phenotype
                    flat_ge_map[gt] = fe
                elif fe < flat_ge_map[gt]:
                    flat_gp_map[gt] = ph  # update to lower free energy ph
                    flat_ge_map[gt] = fe
    # c = 0
    # c_bad = 0
    # for gt in flat_gp_map:
    #     ph = RNA.fold(gt)[0]
    #     if ph != flat_gp_map[gt]:
    #         c_bad += 1
    #     else:
    #         c += 1
    
    # print(c, c_bad)

    # turn into ph_to_gt dict
    ph_to_gt = {}
    for gt in genotypes:
        if gt in flat_gp_map:
            ph = flat_gp_map[gt]
        else:  # if there is no entry, then it must be unfolded
            ph = args.unfolded
        if ph in ph_to_gt:
            ph_to_gt[ph].append(gt)
        else:
            ph_to_gt[ph] = [gt]

    # save to file
    dict_to_gpmap(ph_to_gt, args.output)
