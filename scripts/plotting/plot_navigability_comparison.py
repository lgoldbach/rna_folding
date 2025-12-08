#!/usr/bin/env python

import argparse
import pickle
import RNA
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.stats import pearsonr, linregress

from rna_folding.parsing import load_phenotype_and_metric_from_file, read_ruggedness_file



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--navigability", help="Navigability with normal adaptive walks", required=True)
    parser.add_argument("-n", "--nc_navigability", help="Navigability with neutral component adaptive walks", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))

    navs = []
    for file in [args.navigability, args.nc_navigability]:
        navs.append([])
        with open(file, "r") as f:
            for line in f:
                nav = float(line.strip().split(" ")[1])  # get navigability value
                navs[-1].append(nav)

    ax.scatter(navs[0], navs[1], s=2)

    ax.set_xlabel("Navigability")
    ax.set_ylabel("Navigability (NC walks)")

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
