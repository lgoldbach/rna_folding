#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from rna_folding.parsing import load_phenotype_and_metric_from_file
from rna_folding.utils import ranked_ph_distribution


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--phenotype_dist", help="Path to phenotype "
                        "distribution files", required=True)
    parser.add_argument("-r", "--ref", help="Reference phenotype distribution  "
                        "files ", required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    parser.add_argument("-n", "--nc", help="Path to neutral component file", required=True)
    parser.add_argument("-m", "--ref_nc", help="Reference neutral component "
                        "file ", required=True)
    parser.add_argument("-c", "--rank_cutoff", help="How many ranks to plot max.",
                        type=int, required=False)
    parser.add_argument("--walks", help="Walk success of query", 
                        required=True, 
                        type=str)
    parser.add_argument("--refwalks", help="Walk success of query", 
                        required=True, 
                        type=str)
    parser.add_argument("--sample_size", help="How big is one set of samples per phenotype, i.e. how many samples are takend for a given fitness landscape instance of a phenotype",
                        required=True, type=int)


    args = parser.parse_args()

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

    ### phenotype bias part
    phenotypes, distr = ranked_ph_distribution(ph_distr_file=args.ref,
                                                log=True)

    distr = distr[1:]  # remove unfolded phenotype frequency
    x = range(1, distr.shape[0]+1)
    
    sc = ax1.plot(x, distr, marker="", label=f"ViennaRNA", color="black", linewidth=7)

    phenotypes, distr = ranked_ph_distribution(ph_distr_file=args.phenotype_dist,
                                                log=True)
    distr = distr[1:]  # remove unfolded phenotype frequency
    x = range(1, distr.shape[0]+ 1)
    
    sc = ax1.plot(x, distr, marker="", label=f"Global ranking", color="0.4", linewidth=7, linestyle="dotted", zorder=10)

    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.tick_params(axis='both', which='minor', labelsize=8)

    ax1.set_xlabel("Phenotypes", fontsize=15)
    ax1.set_ylabel("Phenotype frequency (log10)", fontsize=15)
    ax1.grid()
    ax1.legend(loc="lower left", prop={'size': 15}, frameon=False)

    ### neutral component part
    nc_sizes = []
    with open(args.ref_nc, "r") as file:
        for line in file:
            for size in line.strip().split(" "):
                if size:
                    nc_sizes.append(int(size))

    nc_sizes_sort = np.sort(nc_sizes)[::-1]

    if args.rank_cutoff:
        try:
            nc_sizes = nc_sizes_sort[:args.rank_cutoff]
        except IndexError:
            pass

    x = range(nc_sizes.shape[0])
    ax2.plot(x[1:], np.log10(nc_sizes)[1:], marker="", label=f"ViennaRNA", color="black", linewidth=7, markersize=5)

    nc_sizes = []
    with open(args.nc, "r") as file:
        nc_sizes = []
        for line in file:
            for size in line.strip().split(" "):
                if size:
                    nc_sizes.append(int(size))

    nc_sizes_sort = np.sort(nc_sizes)[::-1]

    if args.rank_cutoff:
        try:
            nc_sizes = nc_sizes_sort[:args.rank_cutoff]
        except IndexError:
            pass


    step = 20
    nc_sizes = nc_sizes[1:]

    step_nc_sizes = [nc_sizes[i] for i in range(0, len(nc_sizes), step)]  
    
    step_x = range(0, len(nc_sizes), step)
    ax2.plot(step_x, np.log10(step_nc_sizes), marker="", label=f"Global ranking", color="0.4", linewidth=7, zorder=10, linestyle="dotted")

    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.tick_params(axis='both', which='minor', labelsize=8)

    ax2.set_xlabel("Neutral components", fontsize=15)
    ax2.grid()
    ax2.legend(loc="upper right", prop={'size': 15}, frameon=False)


    ax2.set_ylabel("Neutral component size (log10)", fontsize=15)
    
    ### navigability over frequency part
    ref_ph, ref_freq = load_phenotype_and_metric_from_file(args.ref, dtype=int)
    ph, freq = load_phenotype_and_metric_from_file(args.phenotype_dist, dtype=int)

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
                        d[p][-1] *= 100
        return d
    
    walk_success = read_walk_file(args.walks)
    ref_walk_success = read_walk_file(args.refwalks)

    y_ref = []
    x_ref = []
    y_err_ref = []
    y_err_2d_ref = [[],[]]
    for i, p in enumerate(ref_ph_sort):
        if p in ref_walk_success:
            s = ref_walk_success[p]
            m = np.median(s)
            p_low = np.abs(m-np.percentile(s, q=25))
            p_high = np.abs(m-np.percentile(s, q=75))

            y_ref.append(m)
            x_ref.append(ref_freq_sort[i])
            y_err_ref.append(np.std(s))

            y_err_2d_ref[0].append(p_low)
            y_err_2d_ref[1].append(p_high)
    
    x_query = []
    y_query = []
    y_err = []
    y_err_2d = [[],[]]
    for p in ph_sort:
        if p in walk_success and p in query_d:
            s = walk_success[p]
            m = np.median(s)
            p_low = np.abs(m-np.percentile(s, q=25))
            p_high = np.abs(m-np.percentile(s, q=75))

            y_query.append(m)
            x_query.append(query_d[p])
            y_err.append(np.std(s))

            y_err_2d[0].append(p_low)
            y_err_2d[1].append(p_high)
            

    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]
    
    ref_l = "ViennaRNA"
    query_l = "Global ranking"
    
    ax3.errorbar(np.log10(x_ref), y_ref, yerr=y_err_2d_ref, label=ref_l, linestyle='', marker='s', elinewidth=0.2, color="black", alpha=1, markeredgewidth=0, markersize=7)
    ax3.errorbar(np.log10(x_query), y_query, yerr=y_err_2d, label=query_l, linestyle='', marker='o', elinewidth=0.2, color="0.6", alpha=0.8, markeredgewidth=0, markersize=7)
  

    ax3.set_xlabel("Phenotype frequency (log10)", fontsize=15)
    ax3.set_ylabel("Navigability (%)", fontsize=15)

    plt.tight_layout()

    plt.yticks([0, 20, 40, 60, 80, 100])
    ax3.tick_params(axis='both', which='major', labelsize=12)
    ax3.tick_params(axis='both', which='minor', labelsize=8)

    ax3.grid(zorder=-1)
    # ax3.grid(axis="y", zorder=-1)

    # ax3.set_ylim(-3, 103)
    ax3.legend(loc="upper left", frameon=False, fancybox=False, prop={'size': 15})
        
    plt.savefig(args.output, format="pdf", dpi=30)
