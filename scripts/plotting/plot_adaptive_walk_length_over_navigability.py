#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from rna_folding.utils import ranked_ph_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input files with adaptive "
                        "walk lengths of following format: "
                        "(((...))) 5 2 1 4 2 4 ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    walk_lengths=[]
    aborted = []
    success_perc = []
    for i, file_path in enumerate(args.input, start=2):
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

        success_perc.append(int(np.round(len(walk_lengths[-1])/(aborted[-1]+len(walk_lengths[-1]))*100, 0)))
        
        # t = ax.text(i, 100, f"{}%", fontsize=8, weight="bold", backgroundcolor="white", ha="center", bbox=dict(facecolor='white', alpha=0, edgecolor="none"))
        # ax.bar(i, np.round(len(walk_lengths[-1])/(aborted[-1]+len(walk_lengths[-1])), 2), zorder=10)
        # print(aborted[-1]/(aborted[-1]+len(walk_lengths[-1])), flush=True)
        # print(aborted[-1]+len(walk_lengths[-1]), flush=True)
    # ax.set_ylim(0, 1)     
    walk_length_means = [np.mean(lengths) for lengths in walk_lengths]


    ax.scatter(success_perc, walk_length_means)

    for x, y, bp in zip(success_perc, walk_length_means, list(range(2, 12))):
        ax.annotate(str(bp), xy=(x,y))

    ax.set_xlabel("Navigability")
    ax.set_ylabel("Adaptive walk lengths")
    # ax.set_ylabel("Navigability\n Perc. of adaptive walks reaching target", fontsize=10)

    # plt.xticks(range(2, 2+len(args.input)))

    # ax.grid(axis='y', zorder=1)
    plt.tight_layout()
    # plt.xticks([0, 1, 2, 3, 4, 5, 6, 7, 8], labels=)    
    plt.savefig(args.output, format="pdf", dpi=30)
    # plt.savefig(args.output, format="pdf", dpi=30)
