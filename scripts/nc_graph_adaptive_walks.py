#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import numpy as np
import datetime
import time

from rna_folding.adaptive_walks import nc_graph_to_directed_graph, nc_uniform_adaptive_walk
from rna_folding.adaptive_walks import load_fl_file_to_dict, nc_graph_to_directed_graph


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nc_graph", help="Input gp map.", required=True)
    parser.add_argument("-f", "--fl", help="Fitness landscape "
                        "file", required=True)
    parser.add_argument("-s", "--sample_size_walks", help="How many random walks to sample ",
                        type=int, required=False)
    parser.add_argument("-m", "--max_steps", help="maxmum number of steps", required=True, type=int)
    parser.add_argument("-r", "--seed", help="random seed", type=int,
                        required=False)
    parser.add_argument("-p", "--paths", help="File where all paths are saved", required=True)
    
    args = parser.parse_args()

    if args.seed:
        rng = np.random.default_rng(seed=args.seed)
    else:
        rng = np.random.default_rng(seed=1996)

    nc_graph = pickle.load(open(args.nc_graph, "rb"))

    # load fl. All phenotype have fitness in the interval [0, 1)
    ph_to_f = load_fl_file_to_dict(args.fl)

    nc_to_f = {}
    for nc in nc_graph:
        ph = nc_graph.nodes[nc]["phenotype"]
        f = ph_to_f[ph]
        nc_to_f[nc] = f
    
    max_fit = max(ph_to_f.values())

    # directed graph that only allows uphill steps
    di_graph = nc_graph_to_directed_graph(nc_graph, ph_to_f)

    phenotypes = [nc_graph.nodes[nc]["phenotype"] for nc in nc_graph]
    

    start_nc = rng.choice(nc_graph.nodes, size=args.sample_size_walks, replace=True)

    paths = []  # store whole paths of ncs

    for nc in start_nc:
        # store adaptive walks by target phenotype
        path = nc_uniform_adaptive_walk(di_graph,
                                        nc, 
                                        max_steps=args.max_steps,
                                        rng=rng) 

        paths.append(path)  # save path

    # Write adaptive walk path each paths of genotypes into a single line each 
    with open(args.paths, "a") as file:
        for path in paths:
            file.write(" ".join([str(p) for p in path]) + "\n")
    