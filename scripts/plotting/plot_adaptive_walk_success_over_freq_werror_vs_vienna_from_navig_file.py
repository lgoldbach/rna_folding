#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from rna_folding.parsing import load_phenotype_and_metric_from_file, read_navigability_per_ph_per_fl_file


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--reffreq", help="phenotype frequencies of reference", 
                        required=True, 
                        type=str)
    parser.add_argument("--freq", help="phenotype frequencies of query", 
                        required=True, 
                        type=str)
    parser.add_argument("--navig", help="Navigability file",
                        required=True, 
                        type=str)
    parser.add_argument("--refwalks", help="Walk success of query", 
                        required=True, 
                        type=str)
    parser.add_argument("--sample_size", help="How big is one set of samples per phenotype, i.e. how many samples are takend for a given fitness landscape instance of a phenotype",
                        required=True, type=int)
    parser.add_argument("-o", "--output", help="pdf file", required=True, type=str)
    parser.add_argument("-l", "--label_ref_bp", help="Legend labels for base pairing rule", 
                        type=str, 
                        required=True)
    
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
    # sum_f = 4**12
    print(sum_f)
    ref_freq = [fre/sum_f for fre in ref_freq]
    print(ref_freq)
    sum_f = np.sum(freq)
    print(sum_f)
    # sum_f = 4**12
    freq = [fre/sum_f for fre in freq]
    print(freq)

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
                        d[p][-1] *= 100
        return d
    
    ref_walk_success = read_walk_file(args.refwalks)

    walk_success = read_navigability_per_ph_per_fl_file(args.navig)
    with open(args.navig, "r") as file:
        for line_ in file:
            line = line_.strip().split(" ")
            ph = line[0]
            walk_success[ph] = [float(n)*100 for n in line[1:]]
    
    y_ref = []
    x_ref = []
    y_err_ref = []
    y_err_2d_ref = [[],[]]
    for i, p in enumerate(ref_ph_sort):
        if p in ref_walk_success:
            s = ref_walk_success[p]
            m = np.mean(s)
            p_low = np.abs(m-np.percentile(s, q=25))
            p_high = np.abs(m-np.percentile(s, q=75))

            y_ref.append(m)
            x_ref.append(ref_freq_sort[i])
            y_err_ref.append(np.std(s))

            y_err_2d_ref[0].append(p_low)
            y_err_2d_ref[1].append(p_high)
    # y_ref = [ref_walk_success[p] for p in ref_ph_sort if p in ref_walk_success]
    
    x_query = []
    y_query = []
    y_err = []
    y_err_2d = [[],[]]
    for p in ph_sort:
        if p in walk_success and p in query_d:
            s = walk_success[p]
            m = np.mean(s)
            
            p_low = np.abs(m-np.percentile(s, q=25))
            p_high = np.abs(m-np.percentile(s, q=75))

            y_query.append(m)
            x_query.append(query_d[p])
            y_err.append(np.std(s))

            y_err_2d[0].append(p_low)
            y_err_2d[1].append(p_high)
            

    fig, ax = plt.subplots()

    # ax.scatter(np.log10(x_ref), y_ref, label="Ref", marker="x")
    # ax.scatter(np.log10(x_query), y_query, label="Query", marker="x")

    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]
    
    ref_l = "ViennaRNA"
    query_l = "Global ranking"

    ax.errorbar(np.log10(x_ref), y_ref, yerr=y_err_2d_ref, label=ref_l, linestyle='', marker='s', elinewidth=0.2, color="black", alpha=0.8, markeredgewidth=0)
    ax.errorbar(np.log10(x_query), y_query, yerr=y_err_2d, label=query_l, linestyle='', marker='s', elinewidth=0.2, color="orange", alpha=0.8, markeredgewidth=0)
    


    order = 1
    def log_fit(x, y):
        p = np.polyfit(np.log10(x), y, order)

        func = lambda x: [p[0] * i + p[1] for i in x]
        
        x = np.logspace(min(np.log10(x)), max(np.log10(x)), 250)

        return func, x
    
    def logarithm(x, a, b, c):
        return a*np.log10(x+b)+c
    
    def sigmoid(x, L ,x0, k, b):
        y = L / (1 + np.exp(-k*(x-x0))) + b
        return (y)
    
    def fit(x, y):
        p = np.polyfit(x, y, 5)
        func = np.poly1d(p)
        # try: 
        #     p0 = [max(y), np.median(x),1,min(y)] # this is an mandatory initial guess
        #     popt, popcov = curve_fit(sigmoid, x, y, p0, method='dogbox')
        #     func = lambda x_: sigmoid(x_, *popt)
            
        # except RuntimeError:  # sigmoid fit failed
        #     p0 = [2, np.median(x), min(y)] # this is an mandatory initial guess
        #     popt, popcov = curve_fit(logarithm, x, y, p0, method='dogbox')
        #     func = lambda x_: logarithm(x_, *popt)

            # func, x_trash = log_fit(x, y)
            
        x__ = np.logspace(min(x), max(x), 250)

        return func, x__

    # func, x = fit(np.log10(x_query), y_query)
    # ax.plot(np.log10(x), func(np.log10(x)), color="orange")

    # func, x = fit(np.log10(x_ref), y_ref)
    # ax.plot(np.log10(x), func(np.log10(x)), color="black")
    

    # func, x = log_fit(np.log10(np.array(x_ref)), y_ref)
    # ax.plot(np.log10(x), func(x), color="black")
    # func, x = log_fit(x_query, y_query)
    # ax.plot(np.log10(x), func(x), color="orange")

    # p = np.poly1d(np.polyfit(x_query, y_query, order))
    # t = np.logspace(min(np.log10(x_query)), max(np.log10(x_query)), 250)
    # ax.plot(np.log10(t), p(t))

    ax.set_xlabel("Phenotype frequency (log10)", fontsize=15)
    ax.set_ylabel("Navigability (%)", fontsize=15)

    plt.tight_layout()

    plt.yticks([0, 20, 40, 60, 80, 100])

    ax.grid(axis="y", zorder=-1)

    ax.set_ylim(-3, 103)
    ax.legend(loc="upper left", frameon=False, fancybox=False, prop={'size': 15})

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