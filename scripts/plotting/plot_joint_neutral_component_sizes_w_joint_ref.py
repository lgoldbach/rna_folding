#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nc", help="Path to neutral component files", required=True, nargs='+')
    parser.add_argument("-r", "--ref_nc", help="Reference neutral component "
                        "files ", required=True, nargs="+")
    parser.add_argument("-c", "--rank_cutoff", help="How many ranks to plot max.",
                        type=int, required=False)
    parser.add_argument("-b", "--bp_rule", help="bp rule number ", required=True)
    parser.add_argument("-k", "--ref_bp_rule", help="ref bp rule number ", required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    parser.add_argument("-l", "--log", action="store_true")
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    linestyles=["solid", "dashed", "dotted"]
   
    nc_sizes_all = []
    for i, filename in enumerate(args.nc):
        nc_sizes_all.append([])  # add new list for new file
        with open(filename, "r") as file:
            nc_sizes = []
            for line in file:
                for size in line.strip().split(" "):
                    nc_sizes_all[-1].append(int(size))

    nc_sizes_all_sort = [np.sort(l)[::-1] for l in nc_sizes_all]

    for i, nc_sizes in enumerate(nc_sizes_all_sort):
        if args.rank_cutoff:
            try:
                nc_sizes = nc_sizes[:args.rank_cutoff]
            except IndexError:
                pass

        x = range(nc_sizes.shape[0])
        if args.log:
            ax.plot(x, np.log10(nc_sizes), marker="", label=f"Base-pairing {args.bp_rule}, Ranking {i + 1}", color="tab:blue", markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])
        else:
            ax.plot(x, nc_sizes, marker="", label=f"Base-pairing {args.bp_rule}, Ranking {i + 1}", color="tab:blue", markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])

    nc_sizes_all = []
    for i, filename in enumerate(args.ref_nc):
        nc_sizes_all.append([])  # add new list for new file
        with open(filename, "r") as file:
            nc_sizes = []
            for line in file:
                for size in line.strip().split(" "):
                    nc_sizes_all[-1].append(int(size))

    nc_sizes_all_sort = [np.sort(l)[::-1] for l in nc_sizes_all]

    for i, nc_sizes in enumerate(nc_sizes_all_sort):
        if args.rank_cutoff:
            try:
                nc_sizes = nc_sizes[:args.rank_cutoff]
            except IndexError:
                pass

        x = range(nc_sizes.shape[0])
        if args.log:
            ax.plot(x, np.log10(nc_sizes), marker="", label=f"Base-pairing {args.ref_bp_rule}, Ranking {i + 1}", color="black", markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])
        else:
            ax.plot(x, nc_sizes, marker="", label=f"Base-pairing {args.ref_bp_rule}, Ranking {i + 1}", color="black", markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])

    ax.set_xlabel("Rank")
    if args.log:
        ax.set_ylabel("Neutral component size (log10)")
    else:
        ax.set_ylabel("Neutral component size")
        
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
