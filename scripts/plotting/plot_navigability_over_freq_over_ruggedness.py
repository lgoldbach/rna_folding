#!/usr/bin/env python

import argparse
import pickle
import RNA
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.stats import pearsonr

from rna_folding.parsing import load_phenotype_and_metric_from_file, read_ruggedness_per_ph_file



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--navigability", help="Navigability files", nargs="+", required=True)
    parser.add_argument("-r", "--ruggedness", help="Ruggedness", nargs="+", required=True)
    parser.add_argument("-f", "--ph_dist", help="phenotype distribution", nargs="+", required=True)
    parser.add_argument("-k", "--ruggedness_sample_size", help="Sample size for ruggedness", type=int, required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5), sharey=False)
    
    rug_x_m_all = []
    y_m_all = []
    x_m_all = []
    for i, (ph_dist, navigability, ruggedness) in enumerate(zip(args.ph_dist, args.navigability, args.ruggedness)):
        phenotypes, counts = load_phenotype_and_metric_from_file(ph_dist)

        # read navigability file that has one ph and multiple navigability
        # values per line
        ph_to_navig = {}
        with open(navigability, "r") as file:
            for line_ in file:
                line = line_.strip().split(" ")
                ph = line[0]
                ph_to_navig[ph] = [float(n) for n in line[1:]]

        # local peak sizes
        peak_sizes = read_ruggedness_per_ph_file(ruggedness, n=args.ruggedness_sample_size)

        rugged_av = {}  # compute average local peak size = ruggedness
        for ph in peak_sizes:
            # sum peak sizes and take average over sums
            peak_size_sums = [sum(ps) for ps in peak_sizes[ph]]
            rugged_av[ph] = np.mean(peak_size_sums)  # average size of local peaks
        
        rug_x = []
        x = []
        y = []
        cs = []
        color = []
        for ph, c in zip(phenotypes, counts):
            if c > 0 and ph != "............":
                x.append(c/(rugged_av[ph]+c))
                rug_x.append(rugged_av[ph]+c)
                y.append(np.mean(ph_to_navig[ph]))  # mean navigability   

        im = ax1.scatter(x, y, s=15, alpha=.5, linewidths=0)

        rug_x_m = np.mean(rug_x)
        rug_x_std = np.std(rug_x)

        x_m = np.mean(x)
        x_std = np.std(x)
    
        y_m = np.mean(y)
        y_std = np.std(y)

        rxq1 = [rug_x_m-np.percentile(rug_x, q=25)]
        rxq2 = [max(0, np.percentile(rug_x, q=75)-rug_x_m)]

        xq1 = [x_m-np.percentile(x, q=25)]
        xq2 = [np.percentile(x, q=75)-x_m]
        
        yq1 = [y_m-np.percentile(y, q=25)]
        yq2 = [np.percentile(y, q=75)-y_m]

        ax2.errorbar(rug_x_m, y_m, xerr=(rxq1, rxq2), yerr=(yq1, yq2), elinewidth=1, marker="s", markersize=5, label=f"{i+2}")
    
        ax1.errorbar(x_m, y_m, xerr=(xq1, xq2), yerr=(yq1, yq2), elinewidth=1.5, marker="s", markersize=5, label=f"{i+2}")

        rug_x_m_all.append(rug_x_m)
        x_m_all.append(x_m)
        y_m_all.append(y_m)

    r, p = pearsonr(rug_x_m_all, y_m_all)
    p_str = "%.3g" % p
    ax2.text(.025, .9, f'r = {np.round(r, 2)}\np = {p_str}', transform=ax2.transAxes, horizontalalignment='left', size=10)

    r, p = pearsonr(x_m_all, y_m_all)
    p_str = "%.3g" % p
    ax1.text(.025, .9, f'r = {np.round(r, 2)}\np = {p_str}', transform=ax1.transAxes, horizontalalignment='left', size=10)
        
    ax1.set_xlabel("Target phenotype size / Ruggedness", size=15)
    ax1.set_ylabel("Phenotype accessibility", size=15)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 1)

    ax2.set_xlabel("Ruggedness", size=15)
    ax2.set_ylabel("Phenotype accessibility", size=15)

    ax3.set_xlabel("Target phenotype size / Ruggedness", size=15)
    ax3.set_ylabel("Phenotype accessibility", size=15)
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)    
    
    ax1.legend(fontsize=10, title="g-p map", loc="lower right", frameon=False)
    ax2.legend(fontsize=8, title="g-p map", loc="lower left", frameon=False)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
