#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF


def plot_length_over_counts(path_length_file, 
                            ph_counts_file, 
                            ax, 
                            all_path_lengths,
                            **kwargs):
    path_lengths = {}
    with open(path_length_file, "r") as paths_lengths_f:
        for line in paths_lengths_f:
            line = line.strip().split()
            path_lengths[line[0]] = [int(path_len) for path_len in line[1:]]
            
            all_path_lengths += path_lengths[line[0]]  # add all lengths to a single list
    
    ph_counts = {}
    with open(ph_counts_file, "r") as ph_counts_f:
        for line in ph_counts_f:
            line = line.strip().split(" ")
            ph_counts[line[0]] = int(line[1])

    
    # normalize counts to get frequencies
    total = sum(list(ph_counts.values()))
    ph_freq = {}
    for ph in ph_counts:
        ph_freq[ph] = ph_counts[ph]/total
    
    x = []
    y = []
    mean = []
    stdev = []
    for ph in path_lengths:
        # log transform the count and repeat the value
        # x += [np.log10(ph_freq[ph])] * len(path_lengths[ph]) 
        x.append(np.log10(ph_freq[ph]))
        y = path_lengths[ph]
        mean.append(np.mean(y))
        stdev.append(np.std(y))
    order = np.argsort(x)
    x = np.array(x)[order]
    mean = np.array(mean)[order]
    stdev = np.array(stdev)[order]
    ax.errorbar(x, mean, stdev, **kwargs, linestyle="None")
    # ax.fill_between(x, mean - 1.96 * stdev, mean + 1.96 * stdev, alpha=0.5)

    


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--lengths", help="Path to neutral path lengths "
                        "files", required=True)
    parser.add_argument("-f", "--frequencies", help="Path to file containing "
                        "phenotype frequencies", required=True)
    parser.add_argument("-r", "--ref_lengths", help="Reference neutral path "
                        "lengths ", required=True)
    parser.add_argument("-k", "--ref_frequencies", help="Path to file "
                        "containing reference ph. frequencies", required=True)
    parser.add_argument("-b", "--bp_rule", help="bp rule number ", required=True)
    parser.add_argument("-m", "--ref_bp_rule", help="ref bp rule number ", required=True)
    parser.add_argument("-n", "--ranking", help="Ranking number ", required=True)
    # parser.add_argument("-d", "--distributions", help="Path to file "
    #                     "containing path length distributions", required=True,
    #                     nargs="+")
    # parser.add_argument("-l", "--ref_distribution", help="Path to file "
    #                     "containing reference path length distributions", 
    #                     required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .png)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    # define a color iterator
    color = iter(cm.rainbow(np.linspace(0, 1, len(args.lengths))))
    linestyles=["solid", "dashed", "dotted"]
    
    fig, ax1 = plt.subplots()
    # fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, 
    #                                sharey=True, 
    #                                width_ratios=[3, 1])

    ax1.set_axisbelow(True)  # so grid appears below points

    ax1.set_xlabel("Phenotype frequency (log10)")
    ax1.set_ylabel("Neutral path length")
    ax1.set_xlim(-7, 0)
    ax1.set_ylim(0, 35)
    
    all_path_lengths = []
    plot_length_over_counts(args.lengths, 
                            args.frequencies, 
                            ax1,
                            all_path_lengths=all_path_lengths,
                            alpha=.5,
                            marker=".",
                            elinewidth=1,
                            errorevery=1,
                            zorder=10,
                            label=f"Base-pairing {args.bp_rule}, Random Ranking {args.ranking}", color="tab:blue") #, markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])
        
    all_ref_path_lengths = []
    plot_length_over_counts(args.ref_lengths, 
                            args.ref_frequencies, 
                            ax1,
                            all_path_lengths=all_ref_path_lengths,
                            alpha=.5,
                            marker=".",
                            errorevery=1,
                            elinewidth=1,
                            zorder=0,
                            label=f"Base-pairing {args.ref_bp_rule}, Random Ranking {args.ranking}", color="black") # markevery=(i*3, 10), markersize=5) #, linestyle=linestyles[i])
        
    # # manually create handles for the legend
    # rank_handles = []
    # for i in range(len(args.lengths)):
    #     rank_handles.append(Line2D([0], [0], 
    #                         label=f'Random ranking {i+1}', 
    #                         marker='', 
    #                         markersize=5, 
    #                         color="black",
    #                         linestyle=linestyles[i]))
    
    # bp_rule_handles = []
    # bp_rule_handles.append(Line2D([0], [0], 
    #                     label=f"Base-pairing {args.bp_rule}", 
    #                     marker='', 
    #                     markersize=5, 
    #                     color="tab:blue",
    #                     linestyle="-", linewidth=10))

    # bp_rule_handles.append(Line2D([0], [0], 
    #                     label=f"Base-pairing {args.ref_bp_rule}", 
    #                     marker='', 
    #                     markersize=5, 
    #                     color="black",
    #                     linestyle="-", linewidth=10))
    
    # ax1.legend(handles=rank_handles+bp_rule_handles)
    
    # plot count histogram for random rankings
    # for i, file in enumerate(args.distributions):
    #     data = np.loadtxt(file, delimiter=" ", dtype=int).T
    #     ax2.barh(data[0], data[1], color="orange", alpha=0.5)

    # # # plot count histogram for reference
    # data = np.loadtxt(args.ref_distribution, delimiter=" ", dtype=int).T
    # ax2.barh(data[0], data[1], color="black", alpha=0.5)


    # ax2.hist(all_path_lengths, orientation="horizontal", color="red", alpha=0.5)
    # ax2.hist(all_ref_path_lengths, orientation="horizontal", color="black", alpha=0.5)
    # ax2.set_xlabel("Counts")
    # ax2.set_yticks([])
    # ax2.set_yticks([], minor=True)  # for minor ticks

    plt.legend()
    plt.subplots_adjust(wspace=0.03)
    plt.tight_layout()
    plt.grid()
    # plt.savefig(args.output, format="pdf", dpi=8)
    plt.savefig(args.output, format="pdf", dpi=30)
