#!/usr/bin/env python

import argparse
import pickle
import RNA
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.parsing import load_phenotype_and_metric_from_file



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--walk_lengths", help="Adaptive walk length ", nargs="+", required=True)
    parser.add_argument("-r", "--ruggedness", help="Ruggedness", nargs="+", required=True)
    parser.add_argument("-p", "--peak_heights", help="Peak heights", nargs="+", required=True)
    parser.add_argument("-f", "--ph_dist", help="phenotype distribution", nargs="+", required=True)
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

    def read_peak_height_file_avg(filename, n=10):
        h = {}
        with open(filename, "r") as f:
            lines = list(f)
            for i in range(0, len(lines), n+1):
                ph = lines[i].strip().split(" ")[0]
                for j in range(i+1, i+1+n):  # loop over next n lines
                    # print(lines[j-1], lines[j], lines[j+1])
                    peak_heights = lines[j].strip().split(" ")                    
                    h[ph] = np.mean([float(k) for k in peak_heights])  # sum of peak sizes

        return h

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
    
    for i, (ph_dist, walk_lengths, ruggedness, heights) in enumerate(zip(args.ph_dist, args.walk_lengths, args.ruggedness, args.peak_heights)):
        phenotypes, counts = load_phenotype_and_metric_from_file(ph_dist)
        walk_success = read_walk_file(walk_lengths)

        # peak size
        rugged = read_rugged_file_sum(ruggedness)
        rugged_av = {}
        for ph in rugged:
            rugged_av[ph] = np.mean(rugged[ph])
        
        peak_h = read_peak_height_file_avg(heights)

        x = []
        y = []
        cs = []
        color = []
        for ph, c in zip(phenotypes, counts):
            if c > 0 and ph != "............":
                cs.append(c)
                x.append(c/(rugged_av[ph]+c))
                y.append(walk_success[ph])
                color.append(peak_h[ph])
   
        cm = plt.cm.get_cmap('YlGnBu') 
        im = ax1.scatter(x, y, label=f"GP map {i+2}", s=10, alpha=.7, linewidths=0, c=color)
    
        # average over targets
        # rugged = read_rugged_file_len(ruggedness)
        # rugged_av = {}
        # for ph in rugged:
        #     rugged_av[ph] = np.mean(rugged[ph])
        
        # x = []
        # y = []
        # for ph, c in zip(phenotypes, counts):
        #     if c > 0 and ph != "............":
        #         x.append(c/(rugged_av[ph]+c))
        #         y.append(walk_success[ph])
        
        count_sum = 0
        for ph, c in zip(phenotypes, counts):
            if ph != "............":
                count_sum += c
        count_m = count_sum/(len(phenotypes)-1)
        r_m = np.mean([rugged_av[ph] for ph in rugged_av if ph != "............"])
        
        x_ = count_m/(r_m+count_m)
        y_ = np.mean(y)

        ax3.scatter([x_], [y_], label=f"GP map {i+2}")


        
        x_m = np.mean(x)
        x_std = np.std(x)
        y_m = np.mean(y)
        y_std = np.std(y)
        ax2.errorbar(x_m, y_m, xerr=x_std, yerr=y_std, elinewidth=1, marker="s", markersize=5, label=f"GP map {i+2}")
        
        ### for median NN size
        # x = np.median(cs)
        # ax2.errorbar(x, y_m, xerr=x_std, yerr=y_std, elinewidth=1, marker="s", markersize=5, label=f"GP map {i+2}")
        # ax2.err(x, y, label=f"GP map {i+1}", s=10, alpha=.7, linewidths=0)

    fig.colorbar(im, ax=ax1, label="Average fitness of a local peak genotype")

    ax1.set_xlabel("Target NN size / <Sum of all peak sizes>")
    ax1.set_ylabel("Average Navigability")

    ax3.set_xlabel("<NN size> / <<Sum of all peak sizes>>")
    ax3.set_ylabel("Average Navigability averaged over phenotypes")
    ax3.set_xlim(0, 1)

    ax2.set_xlabel("<(Target NN size / <Sum of all peak sizes>)> (Averaged over all phenotypes)")
    ax2.set_ylabel("Average Navigability averaged over phenotypes")
    # ax2.set_xlim(0, 1)
    
    ax1.set_ylim(0, 100)
    ax2.set_ylim(0, 100)
    ax3.set_ylim(0, 100)
    
    plt.legend(fontsize=6, loc="lower right")
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
