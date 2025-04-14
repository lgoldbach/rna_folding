#!/usr/bin/env python

"""Generate a uniform genotype phenotype maps. If #gt mod #ph > 0 then there it
won't be perfectly uniform.

"""

import argparse

import numpy as np
from rna_folding.parsing import dict_to_gpmap


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file with genotypes", required=True, type=str)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--phenotypes", help="List of phenotypes (either this or --n_ph is required", nargs="+")
    group.add_argument("-n", "--n_ph", help="Number of phenotypes (either this or '--phenotypes' is required", type=int)
    parser.add_argument("-o", "--output", help="File output for phenotypes", required=True, type=str)

    args = parser.parse_args()
    
    # if no phenotypes given simple use a range of numbers
    if not args.phenotypes:
        phenotypes = [str(i) for i in range(1, args.n_ph+1)]
        args.n_ph = len(phenotypes)
    
    with open(args.input, "r") as f:
        genotypes_txt = np.loadtxt(f, dtype=str)
        genotypes = np.array([str(i) for i in range(len(genotypes_txt))])  # to numeric

    div = len(genotypes) // args.n_ph  # how many genotypes per phenotype
    rem = len(genotypes) % args.n_ph  # the remainder
    
    gt_range = np.arange(len(genotypes))  # get a range of same size as genotypes

    rem_gt_idx = np.random.choice(gt_range, size=rem)  # pick random indices

    rem_gt = genotypes[rem_gt_idx]  # pick genotype for the remainder

    genotypes = np.delete(genotypes, rem_gt_idx)
    gt_range = np.arange(len(genotypes))  # get a range of same size as genotypes

    # Assign the same number of random genotypes to each phenotypes
    pg_map = {}
    for ph in phenotypes:
        random_gt_idx = np.random.choice(gt_range, size=div)
        pg_map[ph] = list(genotypes[random_gt_idx])

    # Assign the remaining genotypes as evenly as possible
    already_chosen = []  # keep track of phenotypes already chosen
    for gt in rem_gt:
        while len(already_chosen) < rem:
            ph = np.random.choice(phenotypes)
            if ph in already_chosen:
                continue
            else:
                pg_map[ph].append(gt)
                already_chosen.append(ph)

    dict_to_gpmap(pg_map, args.output)  # save as txt file
