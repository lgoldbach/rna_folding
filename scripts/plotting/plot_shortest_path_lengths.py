#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from rna_folding.utils import ranked_ph_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input files with shortest path "
                        "lengths. SImple one line space-separated file, e.g. "
                        " 5 2 1 4 2 4 ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    path_lengths=[]
    means, medians = [], []
    for i, file_path in enumerate(args.input, start=2):
        with open(file_path, "r") as file:
            for line in file:
                # path_lengths.append([int(l) for l in line.strip().split(" ")])
                mean, median = line.strip().split(" ")
                means.append(float(mean))
                medians.append(float(median))

    x = [i for i in range(2, 2+len(args.input))]
    print(x, means, medians)
    ax.scatter(x, means, marker="x", label="mean")
    ax.scatter(x, medians, marker="o", label="median")
        
    # # ax.boxplot(walk_lengths, positions=[2, 3, 4, 5, 6, 7, 8, 9, 10, 11], showfliers=False)
    # ax.violinplot(path_lengths)
    plt.xticks(x)
    ax.set_xlabel("Base-pairing rule")
    ax.set_ylabel("Shortes path lengths between neutral components")
    # ax.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
