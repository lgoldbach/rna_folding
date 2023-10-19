#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--robustness", help="Path to phenotype "
                        "robustness file", required=True, type=str)
    parser.add_argument("-d", "--phenotype_dist", help="Path to phenotype "
                        "distribution file", required=True, type=str)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()

    data = []
    robus = np.loadtxt(args.robustness)
    distr = np.loadtxt(args.phenotype_dist)

    # sort based on first column (pheno. IDs) so both arrays are in same order
    robus = robus[robus[:, 0].argsort()]
    distr = distr[distr[:, 0].argsort()]

    # get relevant column with robustness/ph. distr values respectively
    robus = robus[:, 1]
    distr = distr[:, 1]

    # turn phenotype counts into log frequency
    distr /= np.sum(distr)
    distr = np.log10(distr)

    fig, ax = plt.subplots()

    ax.scatter(distr, robus, s=5, alpha=0.5)

    ax.set_xlabel("Phenotype frequency (log10)")
    ax.set_ylabel("Phenotype robustness")
    ax.set_title("Phenotype robustness over freqency plot")

    plt.savefig(args.output, format="pdf", dpi=30)