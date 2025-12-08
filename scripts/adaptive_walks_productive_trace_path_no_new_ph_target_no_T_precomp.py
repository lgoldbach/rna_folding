#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import numpy as np
import datetime
import time

from rna_folding.adaptive_walks import kimura_fixation_from_fitness, pairwise_transition_prob_dict, kimura_fixation, productive_adaptive_walk


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--gp_map", help="Input gp map.", required=True)
    parser.add_argument("-f", "--fl", help="Fitness landscape "
                        "file", required=True)
    parser.add_argument("-s", "--sample_size_walks", help="How many random walks to sample ",
                        type=int, required=False)
    parser.add_argument("-n", "--population_size", help="Populaton size", 
                        type=int, required=False)
    parser.add_argument("-m", "--max_steps", help="Maximum number of steps "
                        "per walk", type=int, required=False)
    parser.add_argument("-r", "--seed", help="random seed", type=int,
                        required=False)
    parser.add_argument("-a", "--avoid", help="Phenotypes to avoid as starting "
                        "nodes", required=False)
    parser.add_argument("-l", "--walk_lengths", help="file for output data",
                        required=True)
    parser.add_argument("-p", "--paths", help="File where all paths are saved", required=True)
    
    args = parser.parse_args()

    if args.seed:
        rng = np.random.default_rng(seed=args.seed)
    else:
        rng = np.random.default_rng(seed=1996)

    G = pickle.load(open(args.gp_map, "rb"))

    print("Load complete", flush=True)
    phenotypes = G.phenotype_set
    print("Ph. set complete", flush=True)
    
    # load fl. There will be no new target phenotype selected. All phenotype have fitness
    # in the interval [0, 1)
    ph_to_fitness = {}
    max_fit = 0
    with open(args.fl, "r") as f:
        for line in f:
            data = line.strip().split(" ")
            phenotype = data[0]
            fitness = float(data[1])
            ph_to_fitness[phenotype] = fitness
            if fitness > max_fit:
                max_fit = fitness
    print("f loop complete", flush=True)

    all_nodes = set(G.genotypes)
    print("Num. of genotypes:", len(all_nodes), flush=True)
    print("Genotype sets complete", flush=True)
    # phenotypes to avoid as starting nodes, e.g. lethal ones.
    if args.avoid:
        # get all lethal nodes
        lethal_nodes = [g for i, g in enumerate(G.genotypes) if G.phenotypes[i]==args.avoid]
    else:
        lethal_nodes = []
    # extract nodes that are not target (or lethal if applicable)
    potential_starting_nodes = list(all_nodes.difference(set(lethal_nodes)))
    
    start_gt = rng.choice(potential_starting_nodes, size=min(args.sample_size_walks, len(potential_starting_nodes)), replace=False)

    adaptive_walk_lengths = []  # store adaptive walk lenghts
    paths = []  # store whole paths of genotypes 

    print("Start adaptive walks", flush=True)
    for g in start_gt:
        # store adaptive walks by target phenotype
        path = productive_adaptive_walk(G, g,
                                fitness_function=ph_to_fitness, 
                                fixation_function=kimura_fixation,
                                population_size=args.population_size,
                                max_steps=args.max_steps,
                                rng=rng, max_fit=max_fit)
        
        if ph_to_fitness[G.map(path[-1])] == max_fit:  # walk reached target
            adaptive_walk_lengths.append(len(path))
        else:
            adaptive_walk_lengths.append(-1)  # walk didn't reach target

        paths.append(path)  # save path
    print("Finished adaptive walks", flush=True)

    # Write adaptive walk path each paths of genotypes into a single line each 
    with open(args.paths, "a") as file:
        for path in paths:
            file.write(" ".join(path) + "\n")
    print("Finish writing paths", flush=True)
    # record the adaptive walk lenghts
    with open(args.walk_lengths, "a") as file:
        for lengths in adaptive_walk_lengths:
            file.write(str(lengths) + "\n")

    