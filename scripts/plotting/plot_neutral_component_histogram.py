#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nc", help="Path to neutral component files", required=True, type=str)
    parser.add_argument("-b", "--bp_rule", help="bp rule number ", required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    parser.add_argument("-l", "--log", action="store_true")

    args = parser.parse_args()
    fig, ax = plt.subplots()
   
    # store all neutral component sizes for all phenotypes in one list
    nc_sizes = []
    with open(args.nc, "r") as file:
        for line in file:
            for size in line.strip().split(" "):
                if size:
                    nc_sizes.append(int(size))

    if args.log:
        ax.hist(np.log10(nc_sizes), label=f"Base-pairing {args.bp_rule}")
    else:
        ax.hist(nc_sizes, label=f"Base-pairing {args.bp_rule}")

    print(sum(nc_sizes), 4**12)

    ax.set_ylabel(f"count, sum: {len(nc_sizes)}")
    ax.set_ylim(0, 1000)

    if args.log:
        ax.set_xlabel("Neutral component size (log10)")
    else:
        ax.set_xlabel("Neutral component size")
        
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
