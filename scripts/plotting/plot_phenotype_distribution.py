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
    file_data = np.loadtxt(args.phenotype_dist, dtype=str)
    phenotypes = file_data[:,0]
    distr = file_data[:,1].astype(int)
    distr = np.log10(distr)
    order = np.argsort(distr)[::-1]
    distr = distr[order]
    phenotypes = phenotypes[order]


    x = range(distr.shape[0])

    fig, ax = plt.subplots()
    ax.bar(x, distr)
    ax.set_xticks(x, phenotypes, rotation=45, ha='right')
    ax.set_xlabel("Phenotypes")
    ax.set_ylabel("Phenotype count (log10)")
    ax.set_title("Phenotype distribution")
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
