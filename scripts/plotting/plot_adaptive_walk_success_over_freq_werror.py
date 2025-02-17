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
    parser.add_argument("--sample_size", help="How big is one set of samples per phenotype, i.e. how many samples are takend for a given fitness landscape instance of a phenotype",
                        required=True, type=int)
    parser.add_argument("-o", "--output", help="pdf file", required=True, type=str)
    
    args = parser.parse_args()

    ref_ph, ref_freq = load_phenotype_and_metric_from_file(args.reffreq, dtype=int)
    ph, freq = load_phenotype_and_metric_from_file(args.freq, dtype=int)

    unf_idx = np.where(ph=="............")
    ph = np.delete(ph, unf_idx)
    freq = np.delete(freq, unf_idx)

    unf_idx = np.where(ref_ph=="............")
    ref_ph = np.delete(ref_ph, unf_idx)
    ref_freq = np.delete(ref_freq, unf_idx)


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
                d[p] = []  # a list for all navigabilities for a phenotypes
                # loop over sets of walk lengths, each set coming from one
                # random fitness landscape instance
                for sample_start in range(0, len(line[1:]), args.sample_size):
                    d[p].append(0)  # init a counter for this sample
                    for walk_length in line[1:][sample_start:sample_start+args.sample_size]:
                        if int(walk_length) != -1:  # if not -1 which stands for unsuccessful walk
                            d[p][-1] += 1
                    if d[p][-1] > 0:
                        # compute fraction of successful walks (navig.)
                        d[p][-1] /= args.sample_size
        return d
    
    walk_success = read_walk_file(args.walks)
    ref_walk_success = read_walk_file(args.refwalks)

    y_ref = []
    x_ref = []
    y_err_ref = []
    for i, p in enumerate(ref_ph_sort):
        if p in ref_walk_success:
            y_ref.append(np.mean(ref_walk_success[p]))
            x_ref.append(ref_freq_sort[i])
            y_err_ref.append(np.std(ref_walk_success[p]))
    # y_ref = [ref_walk_success[p] for p in ref_ph_sort if p in ref_walk_success]
    
    x_query = []
    y_query = []
    y_err = []
    for p in ph_sort:
        if p in walk_success and p in query_d:
            y_query.append(np.mean(walk_success[p]))  # get walk success
            x_query.append(query_d[p])  # get freq
            y_err.append(np.std(ref_walk_success[p]))


    fig, ax = plt.subplots()

    # ax.scatter(np.log10(x_ref), y_ref, label="Ref", marker="x")
    # ax.scatter(np.log10(x_query), y_query, label="Query", marker="x")

    ax.errorbar(np.log10(x_ref), y_ref, yerr=y_err_ref, label="Ref", linestyle='', marker='x', elinewidth=.2)
    ax.errorbar(np.log10(x_query), y_query, yerr=y_err, label="Query", linestyle='', marker='x', elinewidth=.2)

    # ax.errorbar(x_ref, y_ref, yerr=y_err_ref, label="Ref", linestyle='', marker='x', elinewidth=.2)
    # ax.errorbar(x_query, y_query, yerr=y_err, label="Query", linestyle='', marker='x', elinewidth=.2)

    order = 1
    def log_fit(x, y):
        p = np.polyfit(np.log10(x), y, order)

        func = lambda x: [p[0] * np.log10(i) + p[1] for i in x]
        
        x = np.logspace(min(np.log10(x)), max(np.log10(x)), 250)

        return func, x

    func, x = log_fit(x_ref, y_ref)
    ax.plot(np.log10(x), func(x), color="blue")
    func, x = log_fit(x_query, y_query)
    ax.plot(np.log10(x), func(x), color="orange")

    # p = np.poly1d(np.polyfit(x_query, y_query, order))
    # t = np.logspace(min(np.log10(x_query)), max(np.log10(x_query)), 250)
    # ax.plot(np.log10(t), p(t))

    ax.legend()

    ax.set_xlabel("Phenotype frequency (log10)")
    ax.set_ylabel("Navigability")

    plt.savefig(args.output, format="pdf", dpi=30)



    # x = np.log10(ref_freq_sort)
    # y = np.log10(freq_sort)
    
    # # x = np.arange(len(ref_ph_sort))
    # y = []
    # y_f = []

    # ax.set_ylim([-6.5, -0.8])
    # ax.set_xlim([-6.5, -0.8])
    # for j, p in enumerate(ref_ph_sort[::-1]):
    #     if p in ph_sort:        
    #         i = np.where(ph_sort[::-1]==p)[0][0]
    #         y.append(i)
    #         y_f.append(freq_sort[i])
    #     else:
    #         y.append(0)
    #         y_f.append(0)
    # ax.scatter(x, np.log10(y_f))
    # ax.plot([-1,-6], [-1,-6], linestyle="--", color="darkgrey")
    # plt.gca().invert_xaxis()
    # plt.gca().invert_yaxis()
    # plt.gca().set_aspect('equal')
    # ax.set_xlabel("Phenotype frequency reference")
    # ax.set_ylabel("Phenotype frequency query")
    # plt.savefig("frequency_correlation.pdf", format="pdf", dpi=30)

    # fig, ax = plt.subplots()
    # ax.set_xticks(np.arange(2, 12))
    # ax.set_ylim(0, 100)
    # ax.set_ylabel("Navigability\n% of successful adaptive walks")
    # ax.set_xlabel("Nucleotide alphabets")
    # ax.bar(np.arange(2, 12), [64, 58, 89, 60, 51, 45, 37, 81, 66, 37])
    # plt.savefig("navigability.pdf", format="pdf", dpi=30)