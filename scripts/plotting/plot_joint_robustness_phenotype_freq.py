#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

def list_of_strings(arg):
    return arg.split(',')

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--robust_distr", help="Path to robustness x phenotype data where column 1 is rob. and column 2 is phenotype freq ",
                         required=True, type=list_of_strings)
    parser.add_argument("-l", "--label", help="Label for data, e.g. graph ids "
                        , required=True, type=list_of_strings)
    parser.add_argument("-e", "--plot_null_expectation", help="Plot null "
                        "expectation of robustness ", required=False, 
                        type=bool, default=False)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True, type=str)
    
    args = parser.parse_args()

    fig, ax = plt.subplots()

    for file, l in zip(args.robust_distr, args.label):
        x = []
        y = []
        with open(file, "r") as f:
            for line in f:
                rob, dist = line.strip().split(" ")
                y.append(float(rob))
                x.append(np.log10(float(dist)))
        
        coeff = np.polyfit(x, y, 2)
        p = np.poly1d(coeff)
        xp = np.linspace(np.min(x)*0.95, np.max(x)*1.05, 100)
        ax.plot(xp, p(xp), label=l)
        # ax.scatter(x, y, s=5, alpha=0.5, label=l)
    ax.set_ylim(bottom=-0.05, top=1)
    ax.set_xlim(left=-8, right=0)
    if args.plot_null_expectation:
        expec = []
        x_expec = []
        exp_of_xlim_max = 10**ax.get_xlim()[1]
        exp_of_xlim_min = 10**ax.get_xlim()[0]
        for freq in np.arange(exp_of_xlim_min, 1, 0.01):
            x_expec.append(np.log10(freq))
            expec.append(freq)
        ax.plot(x_expec, expec, color="grey", ls="--", lw="1")
        
    ax.set_xlabel("Phenotype frequency (log10)")
    ax.set_ylabel("Phenotype robustness")
    ax.set_title("Phenotype robustness over freqency plot")
    ax.legend()

    plt.savefig(args.output, format="pdf", dpi=30)