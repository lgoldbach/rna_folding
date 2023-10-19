#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--phenotype_dist", help="Path to phenotype "
                        "distribution file", required=True, type=str)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()

    data = []
    # load data and get second column (fist only contains phenotype IDs)
    distr = np.loadtxt(args.phenotype_dist)[:,1]
    print(distr)
    distr = np.sort(np.log10(distr))[::-1]
    print(distr)
    x = range(distr.shape[0])

    fig, ax = plt.subplots()
    ax.bar(x, distr)

    ax.set_xlabel("Rank")
    ax.set_ylabel("Phenotype count (log10)")
    ax.set_title("Phenotype distribution")

    plt.savefig(args.output, format="pdf", dpi=30)
