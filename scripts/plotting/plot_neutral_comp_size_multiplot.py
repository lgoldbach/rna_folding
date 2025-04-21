#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np


def read_nc_file(file):
        nc = []
        with open(file, "r") as f:
            for line in f:
                nc += [int(s) for s in line.strip().split(" ")]

        nc = sorted(nc)[::-1]
        nc = nc[1:]  # ignore largest one for unfolded

        return nc

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--nc", help="Paths to neutral component size "
                        "files", required=True, nargs='+')
    parser.add_argument("-r", "--ref", help="Reference bp rule",
                        required=True, type=int)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    args.ref = 4
    
    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(20, 20), sharey=True, sharex=False)
    axes = axes.flatten()
            #   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  
    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]
    nc_files = np.array(args.nc)[new_order_idx]

    adjust_i = 0  # absolutely disgusting hack to skip reference bp rule but still grep correct ax
    for i, file in enumerate(nc_files):  # 2 because rules start at 2
        if i == args.ref-1:
            adjust_i = -1
            continue

        hacky_i = i+adjust_i

        ax = axes[hacky_i]
        if hacky_i in [3]:
            ax.set_ylabel("Neutral component size (log10)", fontsize=40)
        else:
            for tick in ax.yaxis.get_major_ticks():
                tick.tick1line.set_visible(False)
                tick.tick2line.set_visible(False)
                # tick.label1.set_visible(False)
                # tick.label2.set_visible(False)     
        if hacky_i == 7:
            ax.set_xlabel("Neutral components", fontsize=40)

        # ax.set_ylim([-7, -1])
        if hacky_i != 2:
            ax.set_xlim([-10, 250])
        else:
            ax.set_xlim([-10, 1100])
        ax.grid(axis='y')
        ax.grid(axis='x')

        nc = read_nc_file(file)
    
        x = range(1, len(nc)+1)
        ax.plot(x, np.log10(nc), label=f"Base-pairing {i+1}", color="black", linewidth=10, zorder=0)

    ref_file = nc_files[args.ref-1]
    
    ref_nc = read_nc_file(ref_file)
    # only take every step'th data point
    step = 20

    ref_nc_coarse = [ref_nc[i] for i in range(0, len(ref_nc), step)]  
    ref_x = range(1, len(ref_nc)+1, step)
    for ax in axes:
        ax.plot(ref_x, np.log10(ref_nc_coarse), color="0.4", linewidth=10, label="Natural base-pairing", zorder=1, linestyle="dotted")
        ax.legend(loc="upper right", prop={'size': 25}, frameon=False)
    
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=8)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
