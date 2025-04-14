#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import numpy as np
import datetime
import time

from rna_folding.adaptive_walks import greedy_adaptive_walk


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="gp map "
                        "file", required=True)
    parser.add_argument("-l", "--sample_size_landscapes", help="How many random fitness  "
                        "landscapes to sample", type=int, required=False)
    parser.add_argument("-s", "--sample_size_walks", help="How many random walks to sample ",
                        type=int, required=False)
    parser.add_argument("-m", "--max_steps", help="Maximum number of steps "
                        "per walk", type=int, required=False)
    parser.add_argument("-r", "--seed", help="random seed", type=int,
                        required=False)
    parser.add_argument("-u", "--lethal_phenotype", help="Define a lethal phenotype whose fitness will be set to 0", type=str, required=False)
    parser.add_argument("-o", "--output", help="file for output data",
                        required=True)
    parser.add_argument("-d", "--seldif", help="The maximum difference in fitness, i.e. the maximum selection coefficient", type=float, required=True)
    
    args = parser.parse_args()

    if args.seed:
        rng = np.random.default_rng(seed=args.seed)
    else:
        rng = np.random.default_rng(seed=1996)

    # print(f"Start loading", datetime.datetime.now().hour, datetime.datetime.now().minute, flush=True)
    G = pickle.load(open(args.input, "rb"))
    # print(f"Done", datetime.datetime.now().hour, datetime.datetime.now().minute, flush=True)

    phenotypes = sorted(list(set(nx.get_node_attributes(G, "phenotype").values())))
    # track on how many random fitness landscapes phenotype j was reachable by 
    # phenotype i, so that navigability[ph_j][ph_i] will be a number between 0 
    # and sample_size
    navigability = {}

    adaptive_walk_lengths = {}  # store adaptive walk lenghts for each phenotype
    for target_ph in phenotypes:  # loop over target phenotypes
        adaptive_walk_lengths[target_ph] = []

        # print(f"Start {target_ph}", datetime.datetime.now().hour, datetime.datetime.now().minute, flush=True)
        for i in range(args.sample_size_landscapes):
            # assign random fitness to every phenotype
            ph_to_fitness = {}
            for ph in phenotypes:
                f = rng.uniform(0, args.seldif)  # in [0, 1) interval
                ph_to_fitness[ph] = f
            
            if args.lethal_phenotype:
                ph_to_fitness[args.lethal_phenotype] = 0

            ph_to_fitness[target_ph] = 1  # target phenotype gets 1

            # print(f"Start getting genotypes", datetime.datetime.now().hour, datetime.datetime.now().minute, flush=True)
            # start_gt = []
            # while len(start_gt) < args.sample_size_walks:
            #     candidate_gt = np.random.choice(G.nodes)
            #     if G.nodes[candidate_gt]["phenotype"] != target_ph:
            #         start_gt.append(candidate_gt)

            all_nodes = set(G.nodes)
            # get list of target nodes. Do not start walks from there (would be redundant)
            non_starting_nodes = [x for x,y in G.nodes(data=True) if y['phenotype']==target_ph]
            target_node_s = len(non_starting_nodes)
            if args.lethal_phenotype:
                # get all lethal nodes
                lethal_nodes = [x for x,y in G.nodes(data=True) if y['phenotype']==args.lethal_phenotype]
                non_starting_nodes += lethal_nodes  # also disallow lethal nodes as starting nodes
            # extract nodes that are not target (or lethal if applicable)
            potential_starting_nodes = list(all_nodes.difference(set(non_starting_nodes)))
            
            
            start_gt = rng.choice(potential_starting_nodes, size=min(args.sample_size_walks, len(potential_starting_nodes)), replace=False)
            # print(f"Start walks", datetime.datetime.now().hour, datetime.datetime.now().minute, flush=True)
            for g in start_gt:
                # store adaptive walks by target phenotype
                path = greedy_adaptive_walk(G, g,
                                     fitness_function=ph_to_fitness, 
                                     max_steps=args.max_steps,
                                     rng=rng)
                if ph_to_fitness[G.nodes[path[-1]]["phenotype"]] == 1:  # walk reached target
                    adaptive_walk_lengths[target_ph].append(len(path))
                else:
                    adaptive_walk_lengths[target_ph].append(-1)  # walk didn't reach target
            
    with open(args.output, "w") as file:
        for target_ph in adaptive_walk_lengths:
            file.write(f"{target_ph}")
            for path_len in adaptive_walk_lengths[target_ph]:
                file.write(" " + str(path_len))
            file.write("\n")
        