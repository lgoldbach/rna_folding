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

    for i, file in enumerate(args.nc):
        nc = pickle.load(open(file, "rb"))
        nc_sizes_all.append([])  # add new list for new file
        for ph in nc:
            for c in ph:
                nc_sizes_all[-1].append(len(c))
    nc_sizes_all_sort = [np.sort(l)[::-1] for l in nc_sizes_all]

    for nc_sizes in nc_sizes_all_sort:
        x = range(nc_sizes.shape[0])
        if args.log:
            ax.scatter(x, np.log10(nc_sizes), label=i, color="blue")
        else:
            ax.scatter(x, nc_sizes, label=i, color="blue")


    nc_ref = pickle.load(open(args.ref_nc, "rb"))
    nc_ref_sizes = []
    for ph in nc_ref:
        for c in ph:
            nc_ref_sizes.append(len(c))
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
