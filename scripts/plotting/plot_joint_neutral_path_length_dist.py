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
                        "lengths distribution file ", required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    for i, file in enumerate(args.lengths):
        data = np.loadtxt(file, delimiter=" ", dtype=int).T
        ax.plot(data[0], data[1], label=f"{i+2}") #, color="blue")

    # plot reference data
    data = np.loadtxt(args.ref_lengths, delimiter=" ", dtype=int).T
    ax.plot(data[0], data[1], color="black", marker=".")
    

    ax.set_xlabel("Count")
    ax.set_ylabel("Neutral path length")
    
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
