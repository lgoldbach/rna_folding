#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import numpy as np
import datetime
import time

from rna_folding.adaptive_walks import kimura_fixation_from_fitness, pairwise_transition_prob_dict, kimura_fixation, productive_adaptive_walk_w_T


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

    # read in fitness landscape
    ph_to_fitness = {}
    with open(args.fl, "r") as f:
        for line in f:
            data = line.strip().split(" ")
            phenotype = data[0]
            fitness = float(data[1])
            ph_to_fitness[phenotype] = fitness
            
            if fitness == 1:  # assign target phenotype
                target_ph = phenotype


    G = pickle.load(open(args.gp_map, "rb"))

    phenotypes = sorted(list(set(nx.get_node_attributes(G, "phenotype").values())))

    # fixartion probability
    fix_prob = lambda x, y: kimura_fixation_from_fitness(x, y, N=args.population_size)
    
    # precompute transition probabilities between any pairs of phenotypes
    T = pairwise_transition_prob_dict(f_map=ph_to_fitness, func=fix_prob)
        
    all_nodes = set(G.nodes)
    # get list of target nodes. Do not start walks from there (would be redundant)
    non_starting_nodes = [x for x,y in G.nodes(data=True) if y['phenotype']==target_ph]
    target_node_s = len(non_starting_nodes)

    # phenotypes to avoid as starting nodes, e.g. lethal ones.
    if args.avoid:
        # get all lethal nodes
        lethal_nodes = [x for x,y in G.nodes(data=True) if y['phenotype']==args.avoid]
        non_starting_nodes += lethal_nodes  # also disallow lethal nodes as starting nodes
    # extract nodes that are not target (or lethal if applicable)
    potential_starting_nodes = list(all_nodes.difference(set(non_starting_nodes)))
    
    start_gt = rng.choice(potential_starting_nodes, size=min(args.sample_size_walks, len(potential_starting_nodes)), replace=False)
    # print(f"Start walks", datetime.datetime.now().hour, datetime.datetime.now().minute, flush=True)

    adaptive_walk_lengths = []  # store adaptive walk lenghts for each phenotype
    paths = []  # store whole paths of genotypes 
    for g in start_gt:
        # store adaptive walks by target phenotype
        path = productive_adaptive_walk_w_T(G, g,
                                fitness_function=ph_to_fitness, 
                                T=T,
                                max_steps=args.max_steps,
                                rng=rng)
        
        if ph_to_fitness[G.nodes[path[-1]]["phenotype"]] == 1:  # walk reached target
            adaptive_walk_lengths.append(len(path))
        else:
            adaptive_walk_lengths.append(-1)  # walk didn't reach target
        # print("Path")
        # for i in path:
        #     print(ph_to_fitness[G.nodes[i]["phenotype"]])     
        paths.append(path)  # save path

    # Write adaptive walk path each paths of genotypes into a single line each 
    with open(args.paths, "w") as file:
        for path in paths:
            file.write(" ".join(path) + "\n")
    
    # record the adaptive walk lenghts
    with open(args.walk_lengths, "w") as file:
        for lengths in adaptive_walk_lengths:
            file.write(str(lengths) + "\n")

        