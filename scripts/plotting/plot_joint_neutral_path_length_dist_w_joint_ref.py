#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--lengths", help="Path to neutral path length "
                        "distribution files", required=True, nargs='+')
    parser.add_argument("-r", "--ref_lengths", help="Reference neutral path "
                        "lengths distribution files ", required=True, nargs="+")
    parser.add_argument("-b", "--bp_rule", help="bp rule number ", required=True)
    parser.add_argument("-k", "--ref_bp_rule", help="ref bp rule number ", required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    linestyles=["solid", "dashed", "dotted"]


    for i, file in enumerate(args.lengths):
        data = np.loadtxt(file, delimiter=" ", dtype=int).T
        path_lenghts = data[0]
        counts = data[1]
        count_sum = np.sum(counts)
        
        # fill up the missing entries for path lengths that weren't sampled 
        # with zeros
        path_lenghts_all = np.arange(path_lenghts[-1]+2)
        counts_all = np.zeros_like(path_lenghts_all)  # create all zero array
        for c, l in zip(counts, path_lenghts):
            counts_all[l] = c  # add the counts we do have

        # get cdf of path counts
        cdf = np.empty_like(counts_all, dtype=np.float32)
        for length in range(len(counts_all)):
            cdf[length] = np.sum(counts_all[length:]) / count_sum
        ax.plot(path_lenghts_all, cdf, marker="", label=f"Base-pairing {args.bp_rule}, Ranking {i + 1}", color="tab:blue", markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])

    for i, file in enumerate(args.ref_lengths):
        data = np.loadtxt(file, delimiter=" ", dtype=int).T
        path_lenghts = data[0]
        counts = data[1]
        count_sum = np.sum(counts)
        
        # fill up the missing entries for path lengths that weren't sampled 
        # with zeros
        path_lenghts_all = np.arange(path_lenghts[-1]+2)
        counts_all = np.zeros_like(path_lenghts_all)  # create all zero array
        for c, l in zip(counts, path_lenghts):
            counts_all[l] = c  # add the counts we do have

        # get cdf of path counts
        cdf = np.empty_like(counts_all, dtype=np.float32)
        for length in range(len(counts_all)):
            cdf[length] = np.sum(counts_all[length:]) / count_sum
        ax.plot(path_lenghts_all, cdf, marker="", label=f"Base-pairing {args.ref_bp_rule}, Ranking {i + 1}", color="black", markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])

    ax.set_xlabel("Path Length")
    ax.set_ylabel("Probability of Finding a Neutral Path")
    
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
