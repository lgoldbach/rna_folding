#!/usr/bin/env python

import argparse
import pickle
import RNA
import matplotlib.pyplot as plt
import numpy as np


from rna_folding.parsing import load_phenotype_and_metric_from_file



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--phenotype_distribution", help="phenotype distribution", required=True)
    parser.add_argument("-r", "--peak_sizes", help="File with peak sizes", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()

    peak_sizes = {}
    for i in range(1, 100):
        path = f"{args.peak_sizes}/peak_sizes_{str(i)}.txt"
        try:
            with open(path, "r") as f:
                lines = list(f)
                ph = lines[0].strip()
                peak_sizes[ph] = []
                for line in lines[1:]:
                    try:
                        peak_sizes[ph].append([int(s) for s in line.strip().split(" ")])
                    except ValueError:
                        pass
        except FileNotFoundError:
            break
    
    ph, count = load_phenotype_and_metric_from_file(args.phenotype_distribution)

    ph_freq=dict(zip(ph, count))

    fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10,5))

    x = []
    y = []
    for p in ph:
        if p == "............":
            continue
        f = ph_freq[p]
        if f > 0:
            for peaks in peak_sizes[p]:
                y.append(sum(peaks))
                x.append(f)
    ax1.scatter(x, y)
    ax1.set_ylabel("Combined size of local peaks")
    ax1.set_xlabel("Global peak size")
    ax1.set_ylim(0, 1000000)
    x = []
    y = []
    for p in ph:
        if p == "............":
            continue
        f = ph_freq[p]
        if f > 0:
            for peaks in peak_sizes[p]:
                y.append(len(peaks))
                x.append(f)
    ax2.scatter(x, y)
    
    ax2.set_ylabel("Local peak count")
    ax2.set_xlabel("Global peak size")

    plt.savefig(args.output, format="pdf", dpi=30)
