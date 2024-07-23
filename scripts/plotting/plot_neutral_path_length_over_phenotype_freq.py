#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import cm
from matplotlib.lines import Line2D


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

    for ph in path_lengths:
        # log transform the count and repeat the value
        x = [np.log10(ph_counts[ph])] * len(path_lengths[ph])  
        y = path_lengths[ph]
        ax.scatter(x, y, **kwargs)


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--lengths", help="Path to neutral path lengths "
                        "files", required=True, nargs='+')
    parser.add_argument("-f", "--frequencies", help="Path to file containing "
                        "phenotype frequencies", required=True, nargs='+')
    parser.add_argument("-r", "--ref_lengths", help="Reference neutral path "
                        "lengths ", required=True)
    parser.add_argument("-k", "--ref_frequencies", help="Path to file "
                        "containing reference ph. frequencies", required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    # define a color iterator
    color = iter(cm.rainbow(np.linspace(0, 1, len(args.lengths))))
    
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, 
                                   sharey=True, 
                                   width_ratios=[3, 1])

    ax1.set_xlabel("Phenotype count (log10)")
    ax1.set_ylabel("Neutral path lengths")
    
    all_path_lengths = []
    for i, (path_length, ph_counts) in enumerate(zip(args.lengths, args.frequencies)):
        c = next(color)  # pick one color for each file
        
        plot_length_over_counts(path_length, 
                                ph_counts, 
                                ax1,
                                all_path_lengths=all_path_lengths,
                                color=c,
                                alpha=0.2)

    all_ref_path_lengths = []
    plot_length_over_counts(args.ref_lengths, 
                            args.ref_frequencies,
                            ax1,
                            all_path_lengths=all_ref_path_lengths,
                            color="black",
                            alpha=0.2)
    
    # manually create handles for the legend
    rank_handle = Line2D([0], [0], 
                         label='Random ranking', 
                         marker='.', 
                         markersize=15, 
                         color='r',
                         linestyle='')
    
    vienna_handle = Line2D([0], [0], 
                         label='ViennaRNA', 
                         marker='.', 
                         markersize=15, 
                         color='black',
                         linestyle='')
    
    ax1.legend(handles=[vienna_handle, rank_handle])

    ax2.hist(all_path_lengths, orientation="horizontal", color="red", alpha=0.5)
    ax2.hist(all_ref_path_lengths, orientation="horizontal", color="black", alpha=0.5)
    ax2.set_xlabel("Counts")
    ax2.set_yticks([])
    ax2.set_yticks([], minor=True)  # for minor ticks

    plt.subplots_adjust(wspace=0.03)
    # plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
