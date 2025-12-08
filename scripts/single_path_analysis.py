#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.stats import pearsonr
import copy
import networkx as nx

from rna_folding.adaptive_walks import load_fl_file_to_dict, nc_graph_to_directed_graph
from rna_folding.analysis import get_peaks


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--nc_paths", help="File with adaptive walks along neutral components", required=True)
    parser.add_argument("-n", "--nc_graph", help="Neutral component graph", required=True)
    parser.add_argument("-f", "--fitness_landscapes", help="Fitness landscape files", required=True)
    parser.add_argument("-g", "--gp_map", help="whole gp map", required=True)
    parser.add_argument("-i", "--ignore", help="Phenotype to ignore, e.g. the unviable", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()

    nc_graph = pickle.load(open(args.nc_graph, "rb"))
    ph_to_nc = {nc_graph.nodes[nc]["phenotype"]: nc for nc in nc_graph.nodes}
    ph_to_f = load_fl_file_to_dict(args.fitness_landscapes)

    nc_to_f = {}
    for nc in nc_graph:
        ph = nc_graph.nodes[nc]["phenotype"]
        f = ph_to_f[ph]
        nc_to_f[nc] = f

    viable_gt_sum = sum([nc_graph.nodes[nc]["size"] for nc in nc_graph.nodes if nc_graph.nodes[nc]["phenotype"] != args.ignore])

    paths = []
    with open(args.nc_paths, "r") as f:
        for c, line in enumerate(f):
            path_raw = line.strip().split(" ")
            path = [int(path_raw[0])]
            for i, nc in enumerate(path_raw[1:], start=1):
                if nc != path_raw[i-1]:
                    path.append(int(nc))
            paths.append(path)
    
    # get dict of peak id to peak fitness
    peaks_to_f = dict(zip(*get_peaks(nc_graph=nc_graph, ph_to_f=ph_to_f)))

    global_peak_f = max(peaks_to_f.values())

    local_peaks = []
    global_peaks = []
    for peak in peaks_to_f:
        if peaks_to_f[peak] == global_peak_f:  # get global peaks
            global_peaks.append(peak)
        else:  # all non-global peaks are local peaks
            local_peaks.append(peak)

    global_peak_combined_size = sum([nc_graph.nodes[nc]["size"] for nc in global_peaks])
    local_peak_combined_size = sum([nc_graph.nodes[nc]["size"] for nc in local_peaks])
    # Turn nc graph into directed graph to easily check accessibility
    di_graph = nc_graph_to_directed_graph(nc_graph, ph_to_f)

    N = 20
    
    paths_sample_id = np.random.choice(range(len(paths)), N)
    paths_sample = [paths[path_id] for path_id in paths_sample_id]

    fig, axes = plt.subplots(ncols=2, nrows=1, figsize=(10, 5), sharey=True)

    global_ending = 0
    local_ending = 0
    plateau = 0

    xmins = []
    for path in paths_sample:
        global_acc_nc_along_path = []
        local_acc_nc_along_path = []
        
        x_f = []  # x values will be fitness of ncs along path
        for nc in path:  # go along path nc by nc
            x_f.append(nc_to_f[nc])

            global_acc_peaks = []
            for peak in global_peaks:  # go through all global peaks
                if nx.has_path(di_graph, source=nc, target=peak):  # check which are accessible
                    global_acc_peaks.append(peak)
            global_acc_nc_along_path.append(global_acc_peaks)  # save all accessible peaks per step

            # repeat for local peaks
            local_acc_peaks = []
            for peak in local_peaks:  # go through all peaks
                if nx.has_path(di_graph, source=nc, target=peak):  # check which are accessible
                    local_acc_peaks.append(peak)
            local_acc_nc_along_path.append(local_acc_peaks)  # save all accessible peaks per step
    
        y_global = []
        for peaks in global_acc_nc_along_path:
            combined_peak_size = sum([nc_graph.nodes[p]["size"] for p in peaks])
            y_global.append(combined_peak_size)

        y_local = []
        for peaks in local_acc_nc_along_path:
            combined_peak_size = sum([nc_graph.nodes[p]["size"] for p in peaks])
            y_local.append(combined_peak_size)
        
        y = [yg/(yl+yg) for yl, yg in zip(y_local, y_global)]
        
        # step size counting backwards from end of path (=0)
        # e.g. path of length 3 will have x = [-2, -1, 0]
        x = [-len(path)+i+1 for i in range(len(path))]  
        xmins.append(min(x))

        x2 = [nc_to_f[nc] for nc in path]  
        
        if path[-1] in global_peaks:
            axes[0].plot(x, y, marker="x", color="green")
            global_ending += 1
        elif nc in local_peaks:
            axes[0].plot(x, y, marker="x", color="red")
            local_ending += 1
        else:
            axes[0].plot(x2, y, marker="x", color="blue")
            plateau += 1
        
        if path[-1] in global_peaks:
            axes[1].plot(x2, y, marker="x", color="green")
            global_ending += 1
        elif nc in local_peaks:
            axes[1].plot(x2, y, marker="x", color="red")
            local_ending += 1
        else:
            axes[1].plot(x2, y, marker="x", color="blue")
            plateau += 1
    
    all_paths = global_ending+local_ending+plateau
    succ_paths = global_ending/all_paths
    unsucc_paths = local_ending/all_paths
    plateau_paths = plateau/all_paths

    axes[0].bar(1, succ_paths, color="green")
    axes[0].bar(1, unsucc_paths, bottom=succ_paths, color="red")
    axes[0].bar(1, plateau, bottom=succ_paths+unsucc_paths, color="blue")

    axes[0].hlines(global_peak_combined_size/(global_peak_combined_size+local_peak_combined_size), 
                   xmin=min(xmins)-0.3, xmax=1, color="grey", linestyle="--", 
                   zorder=-10, label="|global peaks|/|all peaks|")
    axes[0].set_xlabel("Steps away from end")
    axes[0].set_ylabel("Fraction of accessible peaks that are global peaks")
        # ax.plot(x, y_global, label="global", marker="x")
        # ax.plot(x, y_local, label="local", marker="x")
    
    axes[1].set_xlabel("Fitness")
    axes[1].set_ylabel("Fraction of accessible peaks that are global peaks")

    labels = [item.get_text() for item in axes[0].get_xticklabels()]
    labels[-1] = "Navigability"

    axes[0].legend(loc="lower left")
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=10)
