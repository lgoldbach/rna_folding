#!/usr/bin/env python

import argparse
import json

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="File with nc endpoints", nargs="+", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()

    with open(args.output, "w") as out:
        for inp in args.input:  # for results from each fitness landscape
            nc_path_count = {}
            with open(inp, "r") as file:
                for line_ in file:
                    line = line_.strip()
                    end_nc = line[0]
                    # count how many paths end in a given nc
                    if end_nc not in nc_path_count:
                        nc_path_count[end_nc] = 1
                    else:
                        nc_path_count[end_nc] += 1

            all_paths = sum(list(nc_path_count.values()))  # number of paths
            for nc in nc_path_count:
                nc_nav = nc_path_count[nc]/all_paths  # get fraction of paths
                out.write(f"{nc} {str(nc_nav)}\n")  # write to out





