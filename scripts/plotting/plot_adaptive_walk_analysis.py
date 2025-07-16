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
    parser.add_argument("-p", "--paths", help="File with adaptive walks", nargs="+", required=True)
    parser.add_argument("-n", "--nc_graph", help="Neutral component graph", required=True)
    parser.add_argument("-f", "--fitness_landscapes", help="Fitness landscape files", nargs="+", required=True)
    parser.add_argument("-g", "--nc_to_gt", help="Map from neutral components to genotypes", required=True)
    parser.add_argument("-m", "--gp_map", help="whole gp map", required=True)
    parser.add_argument("-i", "--ignore", help="Phenotype to ignore", type=str, required=False)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()

    fig, axes = plt.subplots(nrows=len(args.fitness_landscapes), 
                             ncols=4, 
                             figsize=(20, 5*len(args.fitness_landscapes)),
                             sharey='col', sharex='col')
    for row in axes:
        for axis in row:
            axis.tick_params(
            axis='x',           # changes apply to the x-axis
            which='both',       # both major and minor ticks are affected
            bottom=True,
            top=False,
            labelbottom=True)    # labels along the bottom edge are on

    nc_graph = pickle.load(open(args.nc_graph, "rb"))
    gp_map = pickle.load(open(args.gp_map, "rb"))
    number_of_gt = 2**25

    ##### need to replace this with gt to nc
    gt_to_nc = {}
    with open(args.nc_to_gt, "r") as f:
        for line_ in f:
            line = line_.strip().split(" ")
            nc = int(line[0])
            for gt in line[1:]:
                gt_to_nc[gt] = nc

    pickle.dump(gt_to_nc, open("gt_to_nc.pkl", "wb"))

    ph_to_nc = {nc_graph.nodes[nc]["phenotype"]: nc for nc in nc_graph.nodes}

    for i, fl in enumerate(args.fitness_landscapes):
        ph_to_f = load_fl_file_to_dict(fl)
    
        nc_reached = {}
        accessible_from = {}

        start_ncs = {nc: [] for nc in nc_graph.nodes}

        peaks, peak_f = get_peaks(nc_graph=nc_graph, ph_to_f=ph_to_f)
        # global_peak = peaks[np.argmax(peak_f)]
        global_peaks = list(np.array(peaks)[np.argwhere(peak_f == np.amax(peak_f))].flatten())  # get all global peaks

        peak_size_sum = sum([nc_graph.nodes[p]["size"] for p in peaks])


        di_graph = nc_graph_to_directed_graph(nc_graph, ph_to_f)

        plateau_steps = [0]*1000
        global_peak_steps = [0]*1000
        local_peak_steps = [0]*1000

        # paths to NC paths (already have this somewhere
        with open(args.paths[i], "r") as f:
            for c, line in enumerate(f):
                path = line.strip().split(" ")


                final_nc = gt_to_nc[path[-1]]  # map to nc

                start_nc = gt_to_nc[path[0]]  # phenotype reached

                if start_nc == final_nc:  # ignore paths that didnt move
                    continue

                for step, gt in enumerate(path):
                    nc = gt_to_nc[gt]
                    if nc in global_peaks:
                        global_peak_steps[step] += 1
                        
                    elif nc in peaks:
                        local_peak_steps[step] += 1
                    else:
                        plateau_steps[step] += 1                

                start_ncs[final_nc].append(start_nc)

                if final_nc not in nc_reached:
                    nc_reached[final_nc] = 1  # add to dir and count 1
                else:
                    nc_reached[final_nc] += 1  # increase count

                if final_nc not in accessible_from:  # if accessibility not yet assessed
                    accessible_from[final_nc] = 0

                    for source in di_graph.nodes:  # for all possible sources
                        if nx.has_path(di_graph, source=source, target=final_nc):  # check if accessible
                            accessible_from[final_nc] += di_graph.nodes[source]["size"]  # size of source nc
        
        per_nc_navigability_peak = []
        per_nc_accessibility_peak = []
        nc_size_peak = []
        nc_f_peak = []

        per_nc_navigability_plat = []
        per_nc_accessibility_plat = []
        nc_size_plat = []
        nc_f_plat = []
        all_paths = sum(list(nc_reached.values()))

        accesibility_global_peak = []
        navigability_global_peak = []
        peak_size_global = []
        peak_f_global = []

        peak_count = 0
        for nc in nc_reached:
            if nc in global_peaks:
                navigability_global_peak.append(nc_reached[nc]/all_paths)
                peak_size_global.append(nc_graph.nodes[nc]["size"]/peak_size_sum)
                accesibility_global_peak.append(accessible_from[nc]/number_of_gt)
                peak_f_global.append(ph_to_f[nc_graph.nodes[nc]["phenotype"]])

                peak_count += nc_reached[nc]
            # note: global peaks are contained in peaks which is 
            # why I need to elif and check global peaks first
            elif nc in peaks:
                per_nc_navigability_peak.append(nc_reached[nc]/all_paths)
                nc_size_peak.append(nc_graph.nodes[nc]["size"]/peak_size_sum)
                per_nc_accessibility_peak.append(accessible_from[nc]/number_of_gt)
                nc_f_peak.append(ph_to_f[nc_graph.nodes[nc]["phenotype"]])

                peak_count += nc_reached[nc]
            else:
                per_nc_navigability_plat.append(nc_reached[nc]/all_paths)
                nc_size_plat.append(nc_graph.nodes[nc]["size"]/peak_size_sum)
                per_nc_accessibility_plat.append(accessible_from[nc]/number_of_gt)
                nc_f_plat.append(ph_to_f[nc_graph.nodes[nc]["phenotype"]])
                

            # if accessible_from[nc] < 10 and nc_reached[nc]/all_paths > 0.1:
            #     print(nc, 
            #           accessible_from[nc],
            #           nc_graph.nodes[nc]["phenotype"], 
            #           peak_size[-1],
            #           peak_f[-1],
            #           len(start_ncs[nc]),
            #           flush=True)
            #     for nc_s in start_ncs[nc]:
            #         print(nc_graph.nodes[nc_s]["size"], flush=True)

        fraction_peaks_over_plateau = np.round(peak_count/all_paths, 2)
        sum_of_paths = sum([global_peak_steps[0], local_peak_steps[0], plateau_steps[0]])
        # pickle.dump(per_peak_navigability, open("per_peak_navigability.pkl", "wb"))
        # pickle.dump(peak_size, open("peak_size.pkl", "wb"))
        # pickle.dump(per_peak_accessibility, open("per_peak_accessibility.pkl", "wb"))
        # pickle.dump(peak_f, open("peak_f.pkl", "wb"))

        axes[i][1].plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=0.5)
        
        axes[i][0].scatter(per_nc_accessibility_peak, per_nc_navigability_peak, color="red", alpha=0.5, label="local peak")
        axes[i][1].scatter(nc_size_peak, per_nc_navigability_peak, color="red", alpha=0.5, label="local peak")
        axes[i][2].scatter(nc_f_peak, per_nc_navigability_peak, color="red", alpha=0.5, label="local peak")

        axes[i][0].scatter(per_nc_accessibility_plat, per_nc_navigability_plat, color="blue", alpha=0.5, label="plateau")
        axes[i][1].scatter(nc_size_plat, per_nc_navigability_plat, color="blue", alpha=0.5, label="plateau")
        axes[i][2].scatter(nc_f_plat, per_nc_navigability_plat, color="blue", alpha=0.5, label="plateau")
        # for global
        axes[i][0].scatter(accesibility_global_peak, navigability_global_peak, marker="*", s=50, color="green", label="global peak")
        axes[i][1].scatter(peak_size_global, navigability_global_peak, marker="*", s=50, color="green", label="global peak")
        axes[i][2].scatter(peak_f_global, navigability_global_peak, marker="*", s=50, color="green", label="global peak")

        axes[i][3].plot(range(len(plateau_steps)), plateau_steps, color="blue", label="plateau")
        axes[i][3].plot(range(len(local_peak_steps)), local_peak_steps, color="red", label="local peak")
        axes[i][3].plot(range(len(global_peak_steps)), global_peak_steps, color="green", label="global peak")

        axes[i][0].set_ylabel("Navigability")
        axes[i][1].set_ylabel("Navigability")
        axes[i][2].set_ylabel("Navigability")
        axes[i][0].set_xlabel("Fraction of genotypes from which NC is accessible")
        axes[i][1].set_xlabel("|NC|/ruggedness") 
        axes[i][2].set_xlabel("NC fitness")
        axes[i][3].set_xlabel("Steps")
        axes[i][3].set_ylabel("Number of paths\n" + "Peak ending fraction: " + str(fraction_peaks_over_plateau) + "\n" + "Number of paths: " + str(sum_of_paths))

        axes[i][0].legend()
        axes[i][1].legend()
        axes[i][2].legend()
        axes[i][3].legend()

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
