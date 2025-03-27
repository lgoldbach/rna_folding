#!/usr/bin/env python

import argparse
import pickle
import RNA
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.parsing import load_phenotype_and_metric_from_file



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--walk_lengths", help="Adaptive walk length ", required=True)
    parser.add_argument("-r", "--ruggedness", help="Ruggedness", required=True)
    parser.add_argument("-f", "--ph_dist", help="phenotype distribution", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5), sharey=True)

    def read_walk_file(filename):
        d = {}
        with open(filename, "r") as f:
            for line_ in f:
                line = line_.strip().split(" ")
                p = line[0]
                d[p] = 0  # start a count
                # loop over sets of walk lengths, each set coming from one
                # random fitness landscape instance
                for walk_length in line[1:]:
                    if int(walk_length) != -1:  # if not -1 which stands for unsuccessful walk
                        d[p] += 1
                if d[p] > 0:
                    # compute fraction of successful walks (navig.)
                    d[p] /= len(line)-1
                    d[p] *= 100
        return d
    
    def read_rugged_file_sum(filename, n=10):
        r = {}
        with open(filename, "r") as f:
            lines = list(f)
            for i in range(0, len(lines), n+1):
                ph = lines[i].strip().split(" ")[0]
                r[ph] = []
                for j in range(i+1, i+1+n):  # loop over next n lines
                    # print(lines[j-1], lines[j], lines[j+1])
                    peaks_sizes = lines[j].strip().split(" ")
                    try:
                        r[ph].append(sum([int(k) for k in peaks_sizes]))  # sum of peak sizes
                    except ValueError:
                        print(peaks_sizes)
        return r

    def read_rugged_file_len(filename, n=10):
        r = {}
        with open(filename, "r") as f:
            lines = list(f)
            for i in range(0, len(lines), n+1):
                ph = lines[i].strip().split(" ")[0]
                r[ph] = []
                for j in range(i+1, i+1+n):  # loop over next n lines
                    # print(lines[j-1], lines[j], lines[j+1])
                    peaks_sizes = lines[j].strip().split(" ")
                    try:
                        r[ph].append(len([int(k) for k in peaks_sizes if int(k) > 0]))  # sum of peak sizes
                    except ValueError:
                        print(peaks_sizes)
        return r
    
    phenotypes, counts = load_phenotype_and_metric_from_file(args.ph_dist)
    walk_success = read_walk_file(args.walk_lengths)

    # peak size
    rugged = read_rugged_file_sum(args.ruggedness)
    rugged_av = {}
    for ph in rugged:
        rugged_av[ph] = np.mean(rugged[ph])
    
    x = []
    y = []
    for ph, c in zip(phenotypes, counts):
        if c > 0 and ph != "............":
            x.append(np.log10(c/rugged_av[ph]))
            y.append(walk_success[ph])

    ax1.scatter(x, y)

    # peak count
    rugged = read_rugged_file_len(args.ruggedness)
    rugged_av = {}
    for ph in rugged:
        rugged_av[ph] = np.mean(rugged[ph])
    
    x = []
    y = []
    for ph, c in zip(phenotypes, counts):
        if c > 0 and ph != "............":
            x.append(np.log10(c/rugged_av[ph]))
            y.append(walk_success[ph])

    ax2.scatter(x, y)

    ax1.set_xlabel("Target NN size / <Sum of local peak sizes> (log10)")
    ax1.set_ylabel("Average Navigability")

    ax2.set_xlabel("Target NN size / <No. local peaks> (log10)")
    # ax.set_ylim(0, 100)
    # ax.set_xlim(-5, 3.2)

    plt.savefig(args.output, format="pdf", dpi=30)
