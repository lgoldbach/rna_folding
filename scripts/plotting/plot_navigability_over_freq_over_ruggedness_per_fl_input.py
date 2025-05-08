#!/usr/bin/env python

import argparse
import pickle
import RNA
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.parsing import load_phenotype_and_metric_from_file, read_ruggedness_per_ph_file
from rna_folding.parsing import read_navigability_per_ph_per_fl_file


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--navigability", help="Navigability files", nargs="+", required=True)
    parser.add_argument("-r", "--ruggedness", help="Ruggedness", nargs="+", required=True)
    parser.add_argument("-f", "--ph_dist", help="phenotype distribution", nargs="+", required=True)
    parser.add_argument("-k", "--ruggedness_sample_size", help="Sample size for ruggedness", type=int, required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), sharey=False)
    
    
    for i, (ph_dist, navigability, ruggedness) in enumerate(zip(args.ph_dist, args.navigability, args.ruggedness)):
        phenotypes, counts = load_phenotype_and_metric_from_file(ph_dist)
        navig = read_navigability_per_ph_per_fl_file(navigability)

        # local peak sizes
        peak_sizes = read_ruggedness_per_ph_file(ruggedness, n=args.ruggedness_sample_size)

        rugged_av = {}  # compute average local peak size = ruggedness
        for ph in peak_sizes:
            # sum peak sizes and take average over sums
            peak_size_sums = [sum(ps) for ps in peak_sizes[ph]]
            rugged_av[ph] = np.mean(peak_size_sums)  # average size of local peaks
        
        x = []
        y = []
        cs = []
        color = []
        for ph, c in zip(phenotypes, counts):
            x.append(c/(rugged_av[ph]+c))
            y.append(np.mean(navig[ph]))  # mean navigability   

        im = ax1.scatter(x, y, label=f"GP map {i+2}", s=10, alpha=.7, linewidths=0)
                    
        x_m = np.mean(x)
        x_std = np.std(x)
        y_m = np.mean(y)
        y_std = np.std(y)
        ax2.errorbar(x_m, y_m, xerr=x_std, yerr=y_std, elinewidth=1, marker="s", markersize=5, label=f"GP map {i+2}")
        

    ax1.set_xlabel("Target NN size / <Sum of all peak sizes>")
    ax1.set_ylabel("Average Navigability")
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)


    ax2.set_xlabel("<(Target NN size / <Sum of all peak sizes>)> (Averaged over all phenotypes)")
    ax2.set_ylabel("Average Navigability averaged over phenotypes")
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)    
    
    plt.legend(fontsize=6, loc="lower right")
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
