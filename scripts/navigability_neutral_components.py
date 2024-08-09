#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import copy
import datetime

from rna_folding.utils import random_fitness_landscape_from_nx_graph, remove_nonadaptive_edges


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Neutral component graph pickle "
                        "file", required=True)
    parser.add_argument("-n", "--sample_size", help="How many random fitness  "
                        "landscapes to sample", type=int, required=False)
    parser.add_argument("-o", "--navigability", help="file for output data",
                        required=True)
    parser.add_argument("-s", "--shortest_path_lengths", help="All shortest path lengths",
                        required=True)


    args = parser.parse_args()

    G_ = pickle.load(open(args.input, "rb"))

    phenotypes = set(nx.get_node_attributes(G_, "phenotype").values())
    
    # track on how many random fitness landscapes phenotype j was reachable by 
    # phenotype i, so that navigability[ph_j][ph_i] will be a number between 0 
    # and sample_size
    navigability = {}

    all_shortest_path_lengths = []
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
            # reachable, if so navig_temp[source_node] = 1, else, none
            # If multiple target neutral components (i.e. same phenotype) are 
            # reachable by the same  source neutral component this will
            #  still be 1. We only care whether it is reachable, not how well.
            navig_temp = {}
            for paths in shortest_paths_all: 
                for source_node in paths:  # if node in paths it has a shortest path
                    if source_node not in navig_temp:  # only count once
                        navig_temp[source_node] = 1

                    all_shortest_path_lengths.append(len(paths[source_node]))
            
            # Add to the total count of how many time source_node reached ph_j
            # across samples
            for source_node in navig_temp:
                if source_node not in navigability[ph_j]:
                    navigability[ph_j][source_node] = 1
                else:
                    navigability[ph_j][source_node] += 1

    with open(args.navigability, "w") as file:
        for ph_j in navigability:
            for source_node in navigability[ph_j]:
                nav = navigability[ph_j][source_node] / args.sample_size # normalize
                file.write(f"{ph_j} {source_node} {nav}\n")

    with open(args.shortest_path_lengths, "w") as file:
        file.write(" ".join([str(l) for l in all_shortest_path_lengths]))


            