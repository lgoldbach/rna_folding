#!/usr/bin/env python

import argparse
import numpy as np
from rna_folding.parsing import dict_to_gpmap


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-i", "--input", help="input gp_map", type=str)
    parser.add_argument("-d", "--dropout", help="How many nodes to drop out", type=int, required=False)
    parser.add_argument("-u", "--unfolded", help="The unfolded genotype", type=str, required=False)
    parser.add_argument("-g", "--unfolded_genotypes", help="Text file containing genotypes that are supposed to be unfolded", type=str, required=False)
    parser.add_argument("-r", "--ranking", help="A file containing phenotypes that defines the ranking", type=str)

    args = parser.parse_args()
    
    flat_gp_map = {}

    phenotypes = [ph.strip() for ph in open(args.ranking, "r")]

    # remove unfolded phenotype from ranking and add to the end
    if args.unfolded in phenotypes:
        phenotypes.remove(args.unfolded)

    ph_to_rank = dict([(ph, i) for i, ph in enumerate(phenotypes)])

    gp_map = open(args.input, "r").readlines()

    
    for line in gp_map:
        line_ = line.strip().split(" ")
        ph = line_[0]

        if ph not in ph_to_rank:  # add phenotype if it's not part of ranking
            rank = max(list(ph_to_rank.values())) + 1
            ph_to_rank[ph] = rank
            print("Added to ranking:", ph, flush=True)
            continue

        rank = ph_to_rank[ph]

        for gt in line_[1:]:
            if gt not in flat_gp_map:  # O(1) operation because dict
                flat_gp_map[gt] = rank  # assign phenotype
            elif rank < flat_gp_map[gt]:
                flat_gp_map[gt] = rank  # update to better ranked ph

    # delete random genotypes from map
    if args.dropout:
        phenotypes.append(args.unfolded)  # add unfolded phenotype to end of ranking for reference

        dropout_ids = np.random.choice(np.arange(len(flat_gp_map)), size=args.dropout, replace=False)
        dropout_keys = np.array(list(flat_gp_map.keys()))[dropout_ids]
        for key in dropout_keys:
            flat_gp_map[key] = -1  # -1 references the last unfolded phenotype

    # set all genotypes from the unfolded_genotypes file to unfolded
    if args.unfolded_genotypes:
        phenotypes.append(args.unfolded)  # add unfolded phenotype to end of ranking for reference
        with open(args.unfolded_genotypes, "r") as unf_gts:
            for line in unf_gts:
                gt = line.strip()
                flat_gp_map[int(gt)] = -1  # -1 references the last unfolded phenotype

    # turn into ph_to_gt dict
    ph_to_gt = dict([(ph, []) for ph in phenotypes])
    for gt in flat_gp_map:
        ph = phenotypes[flat_gp_map[gt]]  # change rank to ph
        ph_to_gt[ph].append(gt)
    # save to file
    dict_to_gpmap(ph_to_gt, args.output)
