#!/usr/bin/env python

import argparse
import pickle
from datetime import datetime


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-o", "--output", help="File output for neutral components",
                        required=True)
    

    args = parser.parse_args()
    
    gpm = pickle.load(open(args.file, "rb"))

    nc_counts = gpm.neutral_components()

    with open(args.output, "w") as file:
        for counts in nc_counts:
            file.write(" ".join([str(c) for c in counts]) + "\n")
    