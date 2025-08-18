#!/usr/bin/env python

import argparse
import pickle
from rna_folding.adaptive_walks import load_fl_file_to_dict, read_genotype_paths_from_file

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--gp_map", help="gp map pickle", required=True)
    parser.add_argument("-p", "--paths", help="adaptive walk paths", nargs="+", required=True)
    parser.add_argument("-f", "--fitness_landscapes", help="fitness_landscape", nargs="+", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()

    gp_map = pickle.load(open(args.gp_map, "rb"))

    with open(args.output, "w") as out:
        for paths, fl in zip(args.paths, args.fitness_landscapes):
            ph_to_f = load_fl_file_to_dict(fl)
            ph_global_peak = max(ph_to_f, key=ph_to_f.get)

            paths = read_genotype_paths_from_file(file=paths, delimiter=" ", gt_type=str, map_to=gp_map._map)
            nav_count = 0
            for path in paths:
                if path[-1] == ph_global_peak:
                    nav_count += 1

            nav = nav_count/len(paths)

            out.write(f"{ph_global_peak} {str(nav)} \n")

    