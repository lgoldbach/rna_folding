#!/usr/bin/env python

import argparse
import numpy as np
from rna_folding.parsing import dict_to_gpmap


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-i", "--input", help="input gp_map", type=str)
    parser.add_argument("-u", "--unfolded", help="The unfolded genotype", type=str, required=True)
    parser.add_argument("-g", "--unfolded_genotypes", help="Text file containing genotypes that are supposed to be unfolded", type=str, required=False)
    parser.add_argument("-r", "--ranking", help="A file containing phenotypes that defines the ranking. IMPORTANT: The last phenotype has to be the unfolded", type=str)

    args = parser.parse_args()

    phenotypes = [ph.strip() for ph in open(args.ranking, "r")]

    ph_to_rank = dict([(ph, i) for i, ph in enumerate(phenotypes)])

    r_unf = len(phenotypes)-1  # last rank to unfolded phenotype


    flat_gp_map = {}
    with open(args.input, "r") as gp_map:
        for line in gp_map:
            line_ = line.strip().split(" ")
            ph = line_[0]

            try:
                rank = ph_to_rank[ph]
            except KeyError:  # the phenotype is not part of the ranking
                rank = r_unf
                ph_to_rank[ph] = rank  # assign it the lowest rank (unfolded ph)
                

            for gt in line_[1:]:
                if gt not in flat_gp_map:  # O(1) operation because dict
                    flat_gp_map[gt] = rank  # assign phenotype
                elif rank < flat_gp_map[gt]:
                    flat_gp_map[gt] = rank  # update to better ranked ph
    
    aa = np.unique(list(flat_gp_map.values()), return_counts=True)
    print(aa, flush=True)
    # delete random genotypes from map
    count = 0
    gt_list = []
    if args.unfolded_genotypes:
        r_unf = ph_to_rank[args.unfolded]  # get rank of unfolded ph
        with open(args.unfolded_genotypes, "r") as unf_gts:
            for line in unf_gts:
                count += 1
                gt = line.strip()
                gt_list.append(gt)
                flat_gp_map[gt] = r_unf
    print(count, flush=True)
    print(np.unique(gt_list, return_counts=True), flush=True)
    # turn into ph_to_gt dict
    ph_to_gt = dict([(ph, []) for ph in phenotypes])
    for gt in flat_gp_map:
        ph = phenotypes[flat_gp_map[gt]]  # translate rank to ph
        ph_to_gt[ph].append(gt)  # append genotype to the ph it maps to


    # save to file
    dict_to_gpmap(ph_to_gt, args.output)
