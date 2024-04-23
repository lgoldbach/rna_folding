#!/usr/bin/env python

import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--files", nargs="+")
    parser.add_argument("-o", "--output", help="File output for neutral components",
                        required=True)
    parser.add_argument("-n", "--num_of_ranks", help="Number of rankings per bp rule",
                        required=True)
    
    args = parser.parse_args()

    fig, axes = plt.subplots(nrows=3, ncols=int(len(args.files)/int(args.num_of_ranks)), 
                             figsize=(2*int(args.num_of_ranks), 10), sharex="col", sharey="row")

    rn = int(args.num_of_ranks)
    for bp_rule, i in enumerate(range(0, len(args.files), rn)):
        for rank, file in enumerate(args.files[i:i+rn]): 
    # for i in range(0, int(len(args.files)/int(args.num_of_ranks))):
    #     files_per_bprule = args.files[i*int(args.num_of_ranks):i*int(args.num_of_ranks)+int(args.num_of_ranks)]
    #     for rank, file in enumerate(files_per_bprule):  # loop over rankings
            # AX1
            ncs_per_ph = pickle.load(open(file, "rb"))
            ncs_count_per_ph = [len(ncs) for ncs in ncs_per_ph]
            for count in ncs_count_per_ph:
                axes[0][bp_rule].scatter(rank, count, alpha=0.5)
            
            if rank == 0:
                axes[0][bp_rule].set_title(f"bp rule {bp_rule+1}")
            if bp_rule == 0:
                axes[0][bp_rule].set_ylabel("Number of connected components\nfor different phenotypes")
            
            # AX2
            ncs_max_per_ph = []
            ncs_sizes_all = []
            for ncs in ncs_per_ph:
                ncs_sizes = []
                for nc in ncs:
                    ncs_sizes.append(len(nc))
                    ncs_sizes_all.append(len(nc))
                ncs_max_per_ph.append(max(ncs_sizes))

            if bp_rule == 0:
                axes[1][bp_rule].set_ylabel("Size of largest of connected components\nfor different phenotypes")

            for max_nc in ncs_max_per_ph:
                axes[1][bp_rule].scatter(rank, max_nc, alpha=0.5)

            # AX3
            axes[2][bp_rule].scatter(rank, np.mean(ncs_sizes_all), marker="o", label="mean" if rank==0 else "", color="black")
            axes[2][bp_rule].scatter(rank, np.median(ncs_sizes_all), marker="x", label="median" if rank==0 else "", color="red")
            axes[2][bp_rule].legend()
            if bp_rule == 0:
                axes[2][bp_rule].set_ylabel("Mean/median neutral component size")
    
    # set ticklabels from 1 to number of rankings
    labels = [str(j) for j in range(1, int(args.num_of_ranks)+1)]
    for ax in axes[-1]:
        ax.set_xlabel("Rankings")
        ax.set_xticks(list(range(0, int(args.num_of_ranks))), labels=labels)

        

    plt.tight_layout()
    plt.savefig(args.output, format="pdf")


        

    
    