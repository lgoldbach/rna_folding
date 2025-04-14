#!/usr/bin/env python

import argparse
import pickle
from datetime import datetime


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-i", "--ignore", help="Phenotype to ignore, e.g "
                        "unfolded", type=str, required=False)
    parser.add_argument("-o", "--output", help="File output for neutral components",
                        required=True)
    

    args = parser.parse_args()
    
    print("start loading", datetime.now())
    gpm = pickle.load(open(args.file, "rb"))
    print("done", datetime.now())

    if args.ignore:
        phenotypes = [ph for ph in gpm.phenotype_set if ph != args.ignore]
    else:
        phenotypes = None

    nc_counts = gpm.neutral_components(phenotypes=phenotypes)
    
    # pickle.dump(gpm, open(args.file, "wb"))
    with open(args.output, "w") as file:
        for counts in nc_counts:
            file.write(" ".join([str(c) for c in counts]) + "\n")
    