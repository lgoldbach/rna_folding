#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nc", help="Path to neutral component files", required=True, nargs='+')
    parser.add_argument("-r", "--ref_nc", help="Reference neutral component "
                        "file ", required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    parser.add_argument("-l", "--log", action="store_true")
    
    args = parser.parse_args()
    fig, ax = plt.subplots()
    nc_sizes_all = []

    for i, filename in enumerate(args.nc):
        nc_sizes_all.append([])  # add new list for new file
        with open(filename, "r") as file:
            nc_sizes = []
            for line in file:
                for size in line.strip().split(" "):
                    nc_sizes_all[-1].append(int(size))

    nc_sizes_all_sort = [np.sort(l)[::-1] for l in nc_sizes_all]

    for nc_sizes in nc_sizes_all_sort:
        x = range(nc_sizes.shape[0])
        if args.log:
            ax.scatter(x, np.log10(nc_sizes), color="blue")
        else:
            ax.scatter(x, nc_sizes, color="blue")

    nc_ref_sizes = []
    with open(args.ref_nc, "r") as file:
        for line in file:
            for size in line.strip().split(" "):
                nc_ref_sizes.append(int(size))
    nc_ref_sizes_sort = np.sort(nc_ref_sizes)[::-1]

    x = range(nc_ref_sizes_sort.shape[0])
    if args.log:
        ax.plot(x, np.log10(nc_ref_sizes_sort), color="black", marker=".")
    else:
        ax.plot(x, nc_ref_sizes_sort, color="black", marker=".")

    ax.set_xlabel("Rank")
    if args.log:
        ax.set_ylabel("Neutral component size (log10)")
    else:
        ax.set_ylabel("Neutral component size")
        
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
