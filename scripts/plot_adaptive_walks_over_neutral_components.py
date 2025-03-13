#!/usr/bin/env python

import argparse
import pickle
import RNA
from rna_folding.adaptive_walks import kimura_fixation

import matplotlib.pyplot as plt


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--paths", help="Adaptive walk paths ", required=True)
    parser.add_argument("-n", "--nc_graph", help="neutral component graph in "
                        ".pickle format ", type=str, required=True)
    parser.add_argument("-m", "--nc_to_gt", help="File that maps neutral"
                        "components to genotypes", required=True)
    parser.add_argument("-f", "--fl", help="fitness landscape, i.e. phenotype "
                        "to fitness map", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    

    args = parser.parse_args()
    
    # read in paths
    paths = []
    with open(args.paths, "r") as f:
        for line in f:
            path = line.strip().split(" ")
            paths.append(path)

    # read file into dictionary
    gt_to_nc = {}
    with open(args.nc_to_gt, "r") as f:
        for line in f:
            data = line.strip().split(" ")
            for gt in data[1:]:
                gt_to_nc[gt] = int(data[0])

    # read fitness landscape
    ph_to_f = {}
    with open(args.fl, "r") as f:
        for line in f:
            data = line.strip().split(" ")
            ph_to_f[data[0]] = float(data[1])

    # load neutral component nx.Graph object  
    nc_graph = pickle.load(open(args.nc_graph, "rb"))

    # translate genotype path to neutral component fitness path
    nc_paths = []
    for path in paths:
        ph = nc_graph.nodes[gt_to_nc[path[0]]]["phenotype"]
        fitness = ph_to_f[ph]
        nc_path = [fitness]
        for i, gt in enumerate(path[1:]):
            try:
                nc = gt_to_nc[gt]
            except KeyError:
            ph = nc_graph.nodes[nc]["phenotype"]
            fitness = ph_to_f[ph]
            if fitness != nc_path[-1]:  # ignore neutral steps
                nc_path.append(fitness)
        nc_paths.append(nc_path)
    fig, ax = plt.subplots()

    success = 0
    for path in nc_paths:
        if path[-1] == 1: 
            success += 1
            ax.plot(range(len(path)), path, color="blue")
        else:
            ax.plot(range(len(path)), path, color="red")

    print(success/len(nc_paths), flush=True)
    plt.savefig(args.output, format="pdf", dpi=30)
