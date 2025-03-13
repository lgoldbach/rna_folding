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
        for i, line in enumerate(f):
            if i < 1000000:
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
            nc = gt_to_nc[gt]
            ph = nc_graph.nodes[nc]["phenotype"]
            fitness = ph_to_f[ph]
            if fitness != nc_path[-1]:  # ignore neutral steps
                nc_path.append(fitness)
        nc_paths.append(nc_path)

    fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, sharey=True, figsize=(10, 5))

    peaks_nc = []
    peaks_f = []
    for nc in nc_graph.nodes:
        neighbors = nc_graph.neighbors(nc)
        nc_ph = nc_graph.nodes[nc]["phenotype"]
        nc_f = ph_to_f[nc_ph]
        peak = True
        for ne in neighbors:
            ne_ph = nc_graph.nodes[ne]["phenotype"]
            ne_f = ph_to_f[ne_ph]
            if ne_f >= nc_f:
                peak = False
        if peak:
            peaks_nc.append(nc)
            peaks_f.append(nc_f)

    success = 0

    l1 = False
    l2 = False
    l3 = False

    plateau_end = 0
    peak_end = 0
    for path in nc_paths:
        if path[-1] == 1: 
            success += 1
            if not l1:
                l1 = True
                ax1.plot(range(len(path)), path, color="blue", alpha=0.2, label="Successful")
            else:
                ax1.plot(range(len(path)), path, color="blue", alpha=0.2)
            ax1.scatter(len(path)-1, path[-1], color="darkblue", marker="^", zorder=10)
        else:
            if path[-1] in peaks_f:
                peak_end += 1
                if not l2:
                    l2 = True
                    ax2.plot(range(len(path)), path, color="red", alpha=0.2, label="Unsuccessful peak-ending")
                else:
                    ax2.plot(range(len(path)), path, color="red", alpha=0.2)

                ax2.scatter(len(path)-1, path[-1], color="darkred", marker="^", zorder=10)
            else:
                plateau_end += 1
                if not l3:
                    l3 = True
                    ax2.plot(range(len(path)), path, color="green", alpha=0.2, label="Unsuccessful plateau-ending")
                else:
                    ax2.plot(range(len(path)), path, color="green", alpha=0.2)
                
                ax2.scatter(len(path)-1, path[-1], color="darkgreen", zorder=10)
    
    ax1.set_xlabel("NC step")
    ax2.set_xlabel("NC step")
    ax1.set_ylabel("Fitness")

    ax1.set_ylim(0, 1.1)
    ax1.yaxis.grid()
    ax2.yaxis.grid()
    leg = ax1.legend(loc="lower right")
    for lh in leg.legendHandles: 
        lh.set_alpha(1)
    leg = ax2.legend(loc="lower right")
    for lh in leg.legendHandles: 
        lh.set_alpha(1)

    fig.suptitle(f"Number of NC: {len(nc_graph.nodes)}, Number of peaks: {len(peaks_nc)}, Success rate: {success/len(nc_paths)}")

    sum_ends = plateau_end+peak_end
    ax2.set_title(f"Plateau end: {int((plateau_end/sum_ends)*100)}%, peak end: {int((peak_end/sum_ends)*100)}%")
    plt.tight_layout()
    
    plt.savefig(args.output, format="pdf", dpi=30)
