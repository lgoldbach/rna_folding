#!/usr/bin/env python

import argparse
import pickle
import networkx as nx


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-o", "--output", help="File output for neutral components",
                        required=True)
    

    args = parser.parse_args()
    
    gpm = pickle.load(open(args.file, "rb"))
    # if edges are not added, add them and save updated gpm in pickle file
    if not gpm.edges:
        gpm.add_hamming_edges()
        pickle.dump(gpm, open(args.file, "wb")) 

    nc = gpm.neutral_components(return_ids=True)

    nc_l = [list(i) for i in nc]
    pickle.dump(nc_l, open(args.output, "wb"))
    