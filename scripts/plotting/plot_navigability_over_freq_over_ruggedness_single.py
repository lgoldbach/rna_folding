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
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5), sharey=False)

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
    
    # compute for each target phenoype

    x = []
    y = []
    x3 = []
    y3 = [] 
    c3 = []
    ax3.plot([-7, -1], [-7, -1], color="grey")
    for ph, c in zip(phenotypes, counts):
        if c > 0 and ph != "............":
            x.append(c/(rugged_av[ph]+c))
            y.append(walk_success[ph])

            x3.append(c/sum(counts))
            y3.append(rugged_av[ph]/sum(counts))
            c3.append(walk_success[ph])


    ax1.scatter(x, y)
    ax1.set_ylim(0, 100)
    ax1.set_xlim(0, 1)

    # target size over ruggedness color by navig
    cm = plt.cm.get_cmap('YlGnBu')
    im=ax3.scatter(np.log10(x3), np.log10(y3), c=c3, cmap=cm)
    fig.colorbar(im, ax=ax3, label="Navigablity")
    ax3.set_xlabel("Target frequency")
    ax3.set_ylabel("Local peak frequency")

    ax3.set_xlim(-7, -1)
    ax3.set_ylim(-7, -1)
    # peak count
    rugged_len = read_rugged_file_len(args.ruggedness)
    rugged_len_av = {}
    for ph in rugged_len:
        rugged_len_av[ph] = np.mean(rugged_len[ph])
    
    x = []
    y = []
    for ph, c in zip(phenotypes, counts):
        if c > 0 and ph != "............":
            x.append(np.log10(c/rugged_len_av[ph]))
            y.append(walk_success[ph])

    ax2.scatter(x, y)
    ax2.set_ylim(0, 100)
    



    ax1.set_xlabel("Target NN size / <Sum of local peak sizes> (log10)")
    ax1.set_ylabel("Average Navigability")

    ax2.set_xlabel("Target NN size / <No. local peaks> (log10)")

    # ax.set_ylim(0, 100)
    # ax.set_xlim(-5, 3.2)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
