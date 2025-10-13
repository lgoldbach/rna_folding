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
    
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(21, 7))
    for ax in axes:
        ax.set_box_aspect(1)
    
            #   [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  
    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]
    nc_files = np.array(args.nc)[new_order_idx]

    linewidth = 2
    scatter_s = 20
    for i, file in enumerate(nc_files, start=1):

        nc = read_nc_file(file)
        x = range(1, len(nc)+1)
        if i in [1, 2, 6, 7]:
            axes[0].plot(x, np.log10(nc), label=f"{i}", color=f"C{i-1}", linewidth=linewidth, zorder=-1, linestyle="-")
            axes[0].scatter(x[0], np.log10(nc)[0], color=f"C{i-1}", s=scatter_s, zorder=-i)
        elif i in [5, 8, 9, 10]:
            axes[1].plot(x, np.log10(nc), label=f"{i}", color=f"C{i-1}", linewidth=linewidth, zorder=-i)
            axes[1].scatter(x[0], np.log10(nc)[0], color=f"C{i-1}", s=scatter_s, zorder=-i)
        elif i == 3:
            axes[2].plot(x, np.log10(nc), label=f"{i}", color=f"C{i-1}", linewidth=linewidth, zorder=10)
            axes[2].scatter(x[0], np.log10(nc)[0], color=f"C{i-1}", s=scatter_s, zorder=-i)
        elif i == 4:
            axes[0].plot(x, np.log10(nc), label=f"canon.", color=f"C{i-1}", linewidth=linewidth+1, zorder=-12, linestyle=":")
            axes[1].plot(x, np.log10(nc), label=f"canon.", color=f"C{i-1}", linewidth=linewidth+1, zorder=-12, linestyle=":")
            axes[2].plot(x, np.log10(nc), label=f"canon.", color=f"C{i-1}", linewidth=linewidth+1, zorder=-12, linestyle=":")
            axes[0].scatter(x[0], np.log10(nc)[0], color=f"C{i-1}", s=scatter_s, zorder=-12)
            axes[1].scatter(x[0], np.log10(nc)[0], color=f"C{i-1}", s=scatter_s, zorder=-12)
            axes[2].scatter(x[0], np.log10(nc)[0], color=f"C{i-1}", s=scatter_s, zorder=-12)

    
    # make sure that canonical legend handle comes first
    for ax in axes:
        new_labels = []
        new_handels = []    
        new_handels_num = []
        handles, labels = ax.get_legend_handles_labels()
        for h, l in zip(handles, labels):
            if l == "canon.":
                new_handels.insert(0, h)
                new_labels.insert(0, l)
            else:
                new_labels.append(l)
                new_handels.append(h)

        ax.legend(new_handels, new_labels, loc="upper right", prop={'size': 18}, frameon=False, title="g-p map", title_fontsize=18)
   

    axes[0].set_xlim(-5, 275)
    axes[1].set_xlim(-5, 275)
    axes[2].set_xlim(-20, 1050)

    for ax in axes:
        ax.set_ylim(1.9, 5.7)
        ax.set_ylabel("Neutral component size (log10)", fontsize=22)
        ax.set_xlabel("Rank", fontsize=22)
        ax.tick_params(axis='both', which='major', labelsize=22)
        ax.tick_params(axis='both', which='minor', labelsize=22)

    # plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
