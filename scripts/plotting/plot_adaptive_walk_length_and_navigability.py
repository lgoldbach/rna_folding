#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.ticker import MaxNLocator
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--low", help="Input files with adaptive "
                        "walk length for low selection pressure of following "
                        "format: (((...))) 5 2 1 4 2 4 ", 
                        required=True, nargs='+')
    parser.add_argument("-m", "--medium", help="Input files with adaptive "
                        "walk length for medium selection pressure of following "
                        "format: (((...))) 5 2 1 4 2 4 ", 
                        required=True, nargs='+')
    parser.add_argument("-g", "--high", help="Input files with adaptive "
                        "walk length for high selection pressure of following "
                        "format: (((...))) 5 2 1 4 2 4 ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

    def get_walk_lengths_and_navigability(input):
        walk_lengths=[]
        aborted = []
        navigability = []
        for i, file_path in enumerate(input, start=2):
            # if i == 10:
            #     continue
            walk_lengths.append([])
            aborted.append(0)
            with open(file_path, "r") as file:
                for line in file:
                    lengths = [int(length) for length in line.strip().split(" ")[1:]]
                    for l in lengths:
                        if l > -1:
                            walk_lengths[-1].append(l)
                        else:
                            aborted[-1] += 1

            navigability.append(int(np.round(len(walk_lengths[-1])/(aborted[-1]+len(walk_lengths[-1]))*100, 0)))
        return walk_lengths, navigability

    labels = ["low", "medium", "high"]
    cmap = colormaps['Oranges']
    colors = [cmap(.5), cmap(.75), cmap(.99)]
    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]


    for i, inp in enumerate([args.low, args.medium, args.high]):
        walk_lengths, navigability = get_walk_lengths_and_navigability(inp)
      
        walk_length_means = [np.mean(lengths) for lengths in walk_lengths]
        walk_length_medians = [np.median(lengths) for lengths in walk_lengths]

        ax1.scatter(range(0, 10), [navigability[i] for i in new_order_idx], color=colors[i], label=labels[i])
        ax1.xaxis.set_major_locator(MaxNLocator(integer=True))

        ax2.scatter(range(0, 10), [walk_length_means[i] for i in new_order_idx], color=colors[i], label=labels[i])
        # ax2.scatter(range(2, 12), walk_length_medians, color=colors[i], label=labels[i])
        ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax1.set_ylabel("Navigability (% of successful walks)")
    ax1.set_xlabel("Base-pairing rule")
    ax2.set_ylabel("Adaptive walk length (# steps)")
    ax2.set_xlabel("Base-pairing rule")

    ax1.set_xticks(list(range(0, 10)))
    ax2.set_xticks(list(range(0, 10)))
    ax1.set_xticklabels([str(i) for i in new_order])
    ax2.set_xticklabels([str(i) for i in new_order])

    l = ax1.legend(title="Selection pressure", prop={'size': 9})
    plt.setp(l.get_title(),fontsize=9)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
