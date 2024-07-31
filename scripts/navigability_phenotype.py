#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import copy

from rna_folding.utils import random_fitness_landscape_from_nx_graph, remove_nonadaptive_edges


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Neutral component graph pickle "
                        "file", required=True)
    parser.add_argument("-n", "--sample_size", help="How many random fitness  "
                        "landscapes to sample", type=int, required=False)
    parser.add_argument("-o", "--output", help="file for output data",
                        required=True)
    

    args = parser.parse_args()

    G_ = pickle.load(open(args.input, "rb"))

    phenotypes = set(nx.get_node_attributes(G_, "phenotype").values())
    
    # track on how many random fitness landscapes phenotype j was reachable by 
    # phenotype i, so that navigability[ph_j][ph_i] will be a number between 0 
    # and sample_size
    navigability = {}

    for ph_j in phenotypes:  # loop over target phenotypes
        navigability[ph_j] = {}  # create a sub-dict for each of them

        for i in range(args.sample_size):
            G = copy.deepcopy(G_)  # deep_copy gp_map
            G = random_fitness_landscape_from_nx_graph(G, peak_phenotype=ph_j)
            G = remove_nonadaptive_edges(G)

            # get all neutral components with target phenotype
            ncs = [g for g, attr in G.nodes(data=True) 
                     if attr['phenotype']==ph_j]
            
            # collect shortest paths, one for each neutral component of the
            # target phenotype
            shortest_paths_all = []  
            for nc in ncs:  # for each neutral component
                # get one shortest paths to target from all neutral_components
                # dict keyed by source (docs are wrong)
                paths = dict(nx.single_target_shortest_path(G, target=nc))

                shortest_paths_all.append(paths)
    
            # check whether at least on of the target neutral components is 
            # reachable, if so navig_temp[ph_i] = 1, else, none
            # If reachable by multiple neutral components this will still be 1
            # we only care whether it is reachable not how well for now
            navig_temp = {}
            for nc, paths in zip(ncs, shortest_paths_all):
                for node in paths:  # if node in paths it has a shortest path
                    ph_i = G.nodes[node]["phenotype"]  # get node phenotype
                    if ph_i not in navig_temp:
                        navig_temp[ph_i] = 1
            
            for ph in navig_temp:
                if ph not in navigability[ph_j]:
                    navigability[ph_j][ph] = 1
                else:
                    navigability[ph_j][ph] += 1

    with open(args.output, "w") as file:
        for ph_j in navigability:
            for ph_i in navigability[ph_j]:
                nav = navigability[ph_j][ph_i] / args.sample_size # normalize
                file.write(f"{ph_j} {ph_i} {nav}\n")
            