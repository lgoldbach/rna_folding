#!/usr/bin/env python

import argparse
import pickle
import RNA
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.parsing import load_phenotype_and_metric_from_file, read_ruggedness_file



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--navigability", help="Navigability per peak ", nargs=True, required=True)
    parser.add_argument("-n", "--nc_graph", help="NC graph", nargs=True, required=True)
    parser.add_argument("-r", "--ruggedness", help="Peak sizes from sample", nargs=True, required=True)
    parser.add_argument("-p", "--ph_dist", help="phenotype distribution", nargs=True, required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    for i, nc_graph, navig, rugg, ph_dist in enumerate(zip(args.nc_graph, args.navigability, args.ruggedness, args.ph_dist)): 
        if i == 1:
            break

        nc_graph = pickle.load(open(nc_graph, "rb"))

        ph_to_nc = {nc_graph.nodes[nc]["phenotype"]: nc for nc in nc_graph.nodes}        

        ph, counts = load_phenotype_and_metric_from_file(ph_dist)
        ph_to_count = dict(ph, counts)

        peak_sizes = read_ruggedness_file(rugg)

        rugg_mean = np.mean([sum(s) for s in peak_sizes])  # average of peak size 

        y_nav = []
        x = []
        with open(navig, "r") as f:
            for line_ in f:
                line = line_.strip().split(" ")
                global_peak_ph = line[0]
                nav = int(line[1])
                y_nav.append(nav)

                ph_count = ph_to_count[global_peak_ph]
                x.append(ph_to_count/rugg_mean)

        ax.scatter(x, y_nav)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
