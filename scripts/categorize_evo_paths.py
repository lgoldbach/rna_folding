#!/usr/bin/env python

import argparse
import pickle
import RNA
from rna_folding.adaptive_walks import load_fl_file_to_dict, genotype_path_to_fitness_path, contains_downhill_steps
from rna_folding.analysis import get_peaks
from rna_folding.parsing import read_adaptive_walks_w_ph_headers_to_dict

import matplotlib.pyplot as plt


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--paths", help="Adaptive walk paths ", nargs="+", required=True)
    parser.add_argument("-f", "--fls", help="fitness landscape, i.e. phenotype "
                        "to fitness map", nargs="+", required=True)
    parser.add_argument("-n", "--nc_graph", help="neutral component graph in "
                        ".pickle format ", type=str, required=True)
    parser.add_argument("-g", "--gp_map", help="GP map in .pickle format ", 
                        type=str, required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)

    args = parser.parse_args()

    nc_graph = pickle.load(open(args.nc_graph, "rb"))
    gpmap = pickle.load(open(args.gp_map, "rb"))
    
    # init dictionary with categories
    categ = {"successful_mono": 0, 
            "unsuccessful_mono": 0, 
            "successful_nonmono": 0,
            "unsuccessful_nonmono": 0}

    # read in path and fitness landscapes
    for path_file, fl_file in zip(args.paths, args.fls):
        ph_to_f_general = load_fl_file_to_dict(fl_file)
        phenotypes = list(ph_to_f_general.keys())
        ph_to_paths = read_adaptive_walks_w_ph_headers_to_dict(path_file, phenotypes=phenotypes)
        
        for ph in ph_to_paths:
            # We have to set the fitness of the target phenotype to 1
            ph_to_f = ph_to_f_general.copy() 
            ph_to_f[ph] = 1
            f_paths = genotype_path_to_fitness_path(ph_to_paths[ph], gpmap, ph_to_f, ignore_neutral=True)

            for path in ph_to_paths[ph]:
                if gpmap.map(path[-1]) == ph:  # reached target -> successful
                    if contains_downhill_steps(path, gpmap, ph_to_f):  # non-monotonic
                        categ["successful_nonmono"] += 1
                    else:  # monotonic
                        categ["successful_mono"] += 1   
                else:
                    if contains_downhill_steps(path, gpmap, ph_to_f):  # non-monotonic
                        categ["unsuccessful_nonmono"] += 1
                    else:  # monotonic
                        categ["unsuccessful_mono"] += 1 

    print(categ)
    pickle.dump(categ, open(args.output, "wb"))
 