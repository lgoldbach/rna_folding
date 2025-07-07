#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import matplotlib.pyplot as plt
import datetime
import time
import numpy as np
import os

from rna_folding.analysis import get_peaks


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nc_graph", help="pickled networkx Graph of "
                        "neutral components", required=True)
    parser.add_argument("-s", "--sample_size", help="sample size", type=int, required=True)
    parser.add_argument("-l", "--low_f", help="Lower fitness limit of an open"
                        "fitness interval, for all f: min_f < f < max_f", 
                        type=float, required=True)
    parser.add_argument("-u", "--upp_f", help="Upper fitness limit of an open"
                        "fitness interval, for all f: min_f < f < max_f", 
                        type=float, required=True)
    parser.add_argument("-d", "--lethal_ph", help="Lethal phenotype to which "
                        "low_f fitness limit will be assigned, i.e. it will "
                        "have the lowest fitness of all phenotypes exlusively", 
                        required=False)                  
    parser.add_argument("-o", "--output", help="Output file name", 
                        type=str, required=True)

    args = parser.parse_args()
    path = os.path.dirname(os.path.abspath(args.output))

    # load neutral component nx.Graph object  
    nc_graph = pickle.load(open(args.nc_graph, "rb"))
    
    phenotypes = set([nc_graph.nodes[node]["phenotype"] for node in nc_graph.nodes])

    rng = np.random.default_rng()
    
    peaks_size = {}
    count = 0
    for target_ph in phenotypes:
        if target_ph == args.lethal_ph:
            continue
        count += 1
        peak_count = []
        peaks_size[target_ph] = []
        for i in range(args.sample_size):
            ph_to_f = {}
            # randomly assign fitness from the open interval (low_f, upp_f), i.e. for
            # all f: low_f < f < upp_f
            for ph in phenotypes:
                # make sure f is never =low_f because numpy uses half open intervals [,)
                # f=low_f is reserved for lethal phenotypes which is applied below
                f = args.low_f
                while f == args.low_f:
                    f = rng.uniform(args.low_f, args.upp_f)
                ph_to_f[ph] = f

            # assign lowest fitness to respective phenotype if applicable
            if args.lethal_ph:
                ph_to_f[args.lethal_ph] = args.low_f

            ph_to_f[target_ph] = args.upp_f
            
            peaks_nc, peaks_f = get_peaks(nc_graph, ph_to_f)
            peak_count.append(len(peaks_nc)-1)  # -1 to substract global
            peaks_size[target_ph].append([str(nc_graph.nodes[peak]["size"]) for peak in peaks_nc if nc_graph.nodes[peak]["phenotype"] != target_ph])
            
    
    with open(args.output, "w") as f:
        for ph in peaks_size:
            f.write(ph + "\n")
            for peaks in peaks_size[ph]:
                if peaks:  
                    f.write(" ".join(peaks) + "\n") 
                else:  # no peaks
                    f.write("0\n")
                    

            
