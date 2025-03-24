#!/usr/bin/env python

import argparse
import pickle
import RNA
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.adaptive_walks import kimura_fixation
from rna_folding.analysis import get_peaks



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--walk_lengths", help="Adaptive walk length ", nargs="+", required=True)
    parser.add_argument("-n", "--nc_graph", help="neutral component graph in "
                        ".pickle format ", required=True)
    parser.add_argument("-f", "--fl", help="fitness landscapes, i.e. phenotype "
                        "to fitness map", nargs="+", required=True)
    parser.add_argument("-p", "--target_phenotype", help="Phenotype that was "
                        "targeted by adaptive walks", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()
    
    # load neutral component nx.Graph object  
    nc_graph = pickle.load(open(args.nc_graph, "rb"))

    navig = []
    peak_count = []
    peaks_size = []
    for (fl_f, wl_f) in zip(args.fl, args.walk_lengths):
        # read fitness landscape
        ph_to_f = {}
        with open(fl_f, "r") as f:
            for line in f:
                data = line.strip().split(" ")
                ph_to_f[data[0]] = float(data[1])

        # compute peaks
        peaks_nc, peaks_f = get_peaks(nc_graph, ph_to_f)

        peak_count.append(len(peaks_nc))

        peaks_size.append(sum([nc_graph.nodes[peak]["size"] for peak in peaks_nc if nc_graph.nodes[peak]["phenotype"] != args.target_phenotype]))

        # read adaptive walk lengths
        wl = np.loadtxt(wl_f, dtype=int)
        fail = len(np.where(wl == -1)[0])

        navig.append(1-(fail/len(wl)))

    # count nc for target phenotype
    ph_nc_count = len([node for node,attr in nc_graph.nodes(data=True) if attr['phenotype']==args.target_phenotype])
    nn_size = sum([attr["size"] for node,attr in nc_graph.nodes(data=True) if attr['phenotype']==args.target_phenotype])
    all_genotype_count = sum([attr["size"] for node,attr in nc_graph.nodes(data=True)])
    nc_count = len(nc_graph.nodes)

    ph_freq = np.round(np.log10(nn_size/all_genotype_count), 1)

    fig, (ax1,ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10,5), sharey=True)

    
    x = peaks_size
    ax1.scatter(np.log10(x), navig)

    ax1.set_xlim(2, 5.8)
    ax1.set_ylim(0, 1)
    ax1.set_title(f"Target phenotype: {args.target_phenotype}\n"
                 f"Phenotype freq.: 10^{ph_freq}\n"
                 f"Number of target phenotype NC: {ph_nc_count}\n"
                 f"NC count total: {nc_count}")
    
    x = [p/nc_count for p in peak_count]
    ax2.scatter(x, navig)
    ax2.set_xlabel("Fraction of NC that are peaks")
    ax2.set_ylabel("Navigability")
    ax2.set_xlim(0, .3)
    plt.subplots_adjust(top=0.8)
    plt.savefig(args.output, format="pdf", dpi=30)
