#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.parsing import load_phenotype_and_metric_from_file

def list_of_strings(arg):
    return arg.split(',')

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reffreq", help="phenotype frequencies of reference", 
                        required=True, 
                        type=str)
    parser.add_argument("--freq", help="phenotype frequencies of query", 
                        required=True, 
                        type=str)
    parser.add_argument("--walks", help="Walk success of query", 
                        required=True, 
                        type=str)
    parser.add_argument("--refwalks", help="Walk success of query", 
                        required=True, 
                        type=str)
    parser.add_argument("-o", "--output", help="pdf file", required=True, type=str)
    
    args = parser.parse_args()

    ref_ph, ref_freq = load_phenotype_and_metric_from_file(args.reffreq, dtype=int)
    ph, freq = load_phenotype_and_metric_from_file(args.freq, dtype=int)

    unf_idx = np.where(ph=="............")
    ph = np.delete(ph, unf_idx)
    freq = np.delete(freq, unf_idx)


    sum_f = np.sum(ref_freq)
    ref_freq = [fre/sum_f for fre in ref_freq]
    sum_f = np.sum(freq)
    freq = [fre/sum_f for fre in freq]

    ref_freq_sort, ref_ph_sort = zip(*sorted(zip(ref_freq, ref_ph)))  # sort both lists by frequency
    freq_sort, ph_sort = zip(*sorted(zip(freq, ph)))  # sort both lists by frequency

    ph_sort = np.array(ph_sort)
    ref_ph_sort = np.array(ref_ph_sort)

    ref_d = dict(zip(ref_ph_sort, ref_freq_sort))  # make dict by ph
    query_d = dict(zip(ph, freq))  # make dict by ph

    def read_walk_file(filename):
        d = {}
        with open(filename, "r") as f:
            for line_ in f:
                line = line_.strip().split(" ")
                p = line[0]
                d[p] = 0
                for walk_length in line[1:]:
                    if int(walk_length) != -1:  # if not -1 which stands for unsuccessful walk
                        d[p] += 1
                if d[p] > 0:
                    d[p] /= len(line[1:])

        return d
    
    walk_success = read_walk_file(args.walks)
    ref_walk_success = read_walk_file(args.refwalks)

    
    y_ref = [ref_walk_success[p] for p in ref_ph_sort]
    
    x_query = []
    y_query = []
    for p in ref_ph_sort:
        if p in walk_success:
            y_query.append(walk_success[p])  # get walk success
            x_query.append(query_d[p])  # get freq


    fig, ax = plt.subplots()

    ax.scatter(np.log10(ref_freq_sort), y_ref, label="ViennaRNA", marker="x")
    ax.scatter(np.log10(x_query), y_query, label="Ranking approach", marker="x")

    ax.legend()

    ax.set_xlabel("Phenotype frequency (log10)")
    ax.set_ylabel("Navigability")

    plt.savefig(args.output, format="pdf", dpi=30)



    x = np.log10(ref_freq_sort)
    y = np.log10(freq_sort)
    
    # x = np.arange(len(ref_ph_sort))
    y = []
    y_f = []

    ax.set_ylim([-6.5, -0.8])
    ax.set_xlim([-6.5, -0.8])
    for j, p in enumerate(ref_ph_sort[::-1]):
        if p in ph_sort:        
            i = np.where(ph_sort[::-1]==p)[0][0]
            print(p, j, i)
            y.append(i)
            y_f.append(freq_sort[i])
        else:
            y.append(0)
            y_f.append(0)
    ax.scatter(x, np.log10(y_f))
    ax.plot([-1,-6], [-1,-6], linestyle="--", color="darkgrey")
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.gca().set_aspect('equal')
    ax.set_xlabel("Phenotype frequency ViennaRNA")
    ax.set_ylabel("Phenotype frequency ranking approach")
    plt.savefig("frequency_correlation.pdf", format="pdf", dpi=30)

    fig, ax = plt.subplots()
    ax.set_xticks(np.arange(2, 12))
    ax.set_ylim(0, 100)
    ax.set_ylabel("Navigability\n% of successful adaptive walks")
    ax.set_xlabel("Nucleotide alphabets")
    ax.bar(np.arange(2, 12), [64, 58, 89, 60, 51, 45, 37, 81, 66, 37])
    plt.savefig("navigability.pdf", format="pdf", dpi=30)