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
    parser.add_argument("-i", "--input", help="gp map "
                        "file", required=True)
    parser.add_argument("-l", "--sample_size_landscapes", help="How many random fitness  "
                        "landscapes to sample", type=int, required=False)
    parser.add_argument("-s", "--sample_size_walks", help="How many random walks to sample ",
                        type=int, required=False)
    parser.add_argument("-n", "--population_size", help="Populaton size", 
                        type=int, required=False)
    parser.add_argument("-m", "--max_steps", help="Maximum number of steps "
                        "per walk", type=int, required=False)
    parser.add_argument("-r", "--seed", help="random seed", type=int,
                        required=False)
    parser.add_argument("-u", "--lethal_phenotype", help="Define a lethal phenotype whose fitness will be set to 0", type=str, required=False)
    parser.add_argument("-o", "--walk_success", help="file for output data",
                        required=True)
    parser.add_argument("-p", "--paths", help="File where all paths are saved", required=True)
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
    paths={}
    for target_ph in phenotypes:  # loop over target phenotypes
        paths[target_ph] = []
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

            # pre compute fixation probability for all phenotype pairs
            fix_prob = lambda x, y: kimura_fixation_from_fitness(x, y, N=args.population_size)
            
            T = pairwise_transition_prob_dict(f_map=ph_to_fitness, func=fix_prob)

            for pair in T:
                if pair[0] == pair[1]:
                    print(pair, ph_to_fitness[pair[0]], ph_to_fitness[pair[1]], T[pair], flush=True)
                    
            # print(ph_to_fitness, flush=True)
            # for i in ph_to_fitness:
            #     for j in ph_to_fitness:
            #         print(i, j, kimura_fixation(ph_to_fitness[j]-ph_to_fitness[i], N=args.population_size), flush=True)
            
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
                path = productive_adaptive_walk_w_T(G, g,
                                     fitness_function=ph_to_fitness, 
                                     T=T,
                                     max_steps=args.max_steps,
                                     rng=rng)
                if ph_to_fitness[G.nodes[path[-1]]["phenotype"]] == 1:  # walk reached target
                    adaptive_walk_lengths[target_ph].append(len(path))
                    paths[target_ph].append(path)  # save successful path
                else:
                    adaptive_walk_lengths[target_ph].append(-1)  # walk didn't reach target
                    paths[target_ph].append(path)  # save unsuccessful path

    with open(args.paths, "w") as file:
            for ph in paths:
                file.write(ph + "\n")
                for path in paths[ph]:
                    file.write(" ".join(path) + "\n")
                    
    with open(args.walk_success, "w") as file:
        for target_ph in adaptive_walk_lengths:
            file.write(f"{target_ph}")
            for path_len in adaptive_walk_lengths[target_ph]:
                file.write(" " + str(path_len))
            file.write("\n")

    

        