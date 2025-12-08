#!/usr/bin/env python

import argparse
import numpy as np
import pickle
from rna_folding.parsing import load_phenotype_and_metric_from_file
from rna_folding.analysis import get_peaks



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--nc_paths", help="Files that contain " \
    "adaptive walk lengths for each phenotype. One file per fitness landscape", 
    nargs="+", required=True)       
    parser.add_argument("-f", "--fitness_landscapes", help="Fitness landscapes", nargs="+", required=True)
    parser.add_argument("-n", "--nc_graph", help="Neutral component graphs", required=True)
    parser.add_argument("-o", "--output", help="Output file name",
                        type=str, required=True)

    args = parser.parse_args()

    nc_graph = pickle.load(open(args.nc_graph, "rb"))

    fl_success_count = {}
    g_peak_ph = {}
    ruggedness = {}
    for i, (paths, fl) in enumerate(zip(args.nc_paths, args.fitness_landscapes)):  # i enumerates fitn. landsc.
        ph, f = load_phenotype_and_metric_from_file(fl)
        ph_to_f = dict(zip(ph, f))
        g_peak_ph[i] = ph[np.argmax(f)]

        ruggedness[i] = 0
        peaks_nc, peaks_f = get_peaks(nc_graph, ph_to_f)
        for nc in peaks_nc:
            ruggedness[i] += nc_graph.nodes[nc]["size"]
        
        path_count = 0
        fl_success_count[i] = 0

        with open(paths, "r") as file:
            for j, line in enumerate(file):
                path = line.strip().split()
                end_nc = int(path[-1])
                # only count path if it ended at either local or global peak
                # ignore paths that are still at a plateau
                if end_nc in peaks_nc:
                    path_count += 1
                # check if global peak by checking if 
                if nc_graph.nodes[end_nc]["phenotype"] == g_peak_ph[i]:
                    fl_success_count[i] += 1  # count as success
        fl_success_count[i] /= path_count  # get fraction of successful paths

    with open(args.output, "w") as f:
        for i in fl_success_count:
            f.write(f"{g_peak_ph[i]} {str(fl_success_count[i])} {ruggedness[i]} \n")  # new line after every fl block
    
        

                        
                        



    
                    

            
