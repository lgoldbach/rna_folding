#!/usr/bin/env python

import argparse
import numpy as np
import pickle
from rna_folding.parsing import load_phenotype_and_metric_from_file
from rna_folding.analysis import get_peaks



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("-f", "--fitness_landscapes", help="Fitness landscapes", nargs="+", required=True)
    parser.add_argument("-n", "--nc_graph", help="Neutral component graphs", required=True)
    parser.add_argument("-o", "--output", help="Output file name",
                        type=str, required=True)

    args = parser.parse_args()

    nc_graph = pickle.load(open(args.nc_graph, "rb"))

    global_peak_ratio = {}
    for i, fl in enumerate(args.fitness_landscapes):  # i enumerates fitn. landsc.
        ph, f = load_phenotype_and_metric_from_file(fl)
        ph_to_f = dict(zip(ph, f))
        max_f = max(f)

        peaks, peaks_f = get_peaks(nc_graph, ph_to_f)

        global_peaks = []
        for j, peak in enumerate(peaks):
            if peaks_f[j] == max_f:
                global_peaks.append(peak)
        
        global_peak_ratio[i] = len(global_peaks)/len(peaks)
        

    with open(args.output, "w") as f:
        for i in global_peak_ratio:
            f.write(f"{global_peak_ratio[i]}\n")  # new line after every fl block
    
        

                        
                        



    
                    

            
