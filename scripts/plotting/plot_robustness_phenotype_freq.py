#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

def list_of_strings(arg):
    return arg.split(',')

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--robustness", help="Path to phenotype "
                        "robustness file", required=True, type=list_of_strings)
    parser.add_argument("-d", "--phenotype_dist", help="Path to phenotype "
                        "distribution file", required=True, type=list_of_strings)
    parser.add_argument("-l", "--label", help="Label for data, e.g. graph ids "
                        , required=True, type=list_of_strings)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()

    fig, ax = plt.subplots()

    for robust_in, dist_in, l in zip(args.robustness, args.phenotype_dist, args.label):
        robus = {}
        with open(robust_in, "r") as file:
            for line_ in file:
                line = line_.strip().split(" ")
                robus[line[0]] = float(line[1])
        
        distr = {}
        with open(dist_in, "r") as file:
            for line_ in file:
                line = line_.split(" ")
                distr[line[0]] = int(line[1])

        phenotype_count_sum = sum(distr.values())
        x = []
        y = []
        for ph in robus:
            y.append(robus[ph])
            # turn phenotype counts into log frequency
            d = distr[ph]/phenotype_count_sum
            d = np.log10(d)
            x.append(d)
        ax.scatter(x, y, s=5, alpha=0.5, label=l)
        
        
    ax.set_xlabel("Phenotype frequency (log10)")
    ax.set_ylabel("Phenotype robustness")
    ax.set_title("Phenotype robustness over freqency plot")
    ax.legend()

    plt.savefig(args.output, format="pdf", dpi=30)