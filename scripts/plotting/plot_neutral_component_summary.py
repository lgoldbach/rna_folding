#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--files", nargs="+")
    parser.add_argument("-o", "--output", help="File output for neutral components",
                        required=True)
    

    args = parser.parse_args()

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10,5))

    for i, file in enumerate(args.files):  # loop over rankings
        ncs_per_ph = pickle.load(open(file, "rb"))
        ncs_count_per_ph = [len(ncs) for ncs in ncs_per_ph]
        for count in ncs_count_per_ph:
            ax1.scatter(i, count, alpha=0.5)

        ncs_max_per_ph = []
        for ncs in ncs_per_ph:
            ncs_sizes = []
            for nc in ncs:
                ncs_sizes.append(len(nc))
            ncs_max_per_ph.append(max(ncs_sizes))

        for max_nc in ncs_max_per_ph:
            ax2.scatter(i, max_nc, alpha=0.5)

    ax1.set_title("Number of neutral components per phenotype")
    ax1.set_xlabel("Rankings")
    ax1.set_ylabel("Number of connected components\nfor different phenotypes")
    ax2.set_title("Largest neutral components per phenotype")
    ax2.set_xlabel("Rankings")
    ax2.set_ylabel("Size of largest of connected components\nfor different phenotypes")
    plt.savefig(args.output, format="pdf")


        

    
    