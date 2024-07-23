#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

def list_of_strings(arg):
    return arg.split(',')

if __name__ ==  "__main__":
    # assumes the first element is the reference data
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--robustness", help="Path to phenotype "
                        "robustness file", required=True, nargs="+")
    parser.add_argument("-d", "--phenotype_dist", help="Path to phenotype "
                        "distribution file", required=True, nargs="+")
    parser.add_argument("-l", "--label", help="Label for data, e.g. graph ids "
                        , required=True, nargs="+")
    parser.add_argument("-e", "--plot_null_expectation", help="Plot null "
                        "expectation of robustness ", required=False, 
                        type=bool, default=False)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True, type=str)
    
    args = parser.parse_args()

    fig, ax = plt.subplots()

    for i, (robust_in, dist_in, l) in enumerate(zip(args.robustness, args.phenotype_dist, args.label)):
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
            # turn phenotype counts into log frequency
            try: 
                a = distr[ph]
            except KeyError:
                print(ph, "\n", len(distr.keys()), len(robus.keys()))
                continue
            y.append(robus[ph])
            d = distr[ph]/phenotype_count_sum
            d = np.log10(d)
            x.append(d)
        
        # assume the first element is the reference and plot that differently
        if i == 0:  
            ax.scatter(x, y, s=20, label=l, marker="x", color="0.4")
        else:
            ax.scatter(x, y, s=5, alpha=1, label=f"Random ranking {l}")

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
        ax.plot(x_expec, expec, color="grey", ls="--", lw="1", label="Null expectation")
        
    ax.set_xlabel("Phenotype frequency (log10)")
    ax.set_ylabel("Phenotype robustness")
    # ax.set_title("Phenotype robustness over freqency plot")
    ax.legend(loc="upper left")
    
    ax.set_axisbelow(True)  # needed to set grid below points
    plt.grid()
    
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)