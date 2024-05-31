#!/usr/bin/env python

import argparse
import pickle
from datetime import datetime
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", 
                        help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-n", "--num_of_genotypes", 
                        help="From how many genotypes to sample neutral paths",
                        required=True, type=int)
    parser.add_argument("-s", "--num_of_paths", 
                        help="How many neutral paths to sample per genotype",
                        required=True, type=int)
    parser.add_argument("-o", "--output", 
                        help="File output for neutral components",
                        required=True)
    
    args = parser.parse_args()
    
    gpm = pickle.load(open(args.file, "rb"))
    # if edges are not added, add them and save updated gpm in pickle file

    neutral_paths = []
    for i in range(args.num_of_genotypes):
        genotype = np.random.choice(gpm.genotypes)
        neutral_paths.append(gpm.neutral_paths(genotype, n=args.num_of_paths))

    with open(args.output, "w") as file:
        for path_set in neutral_paths:
            for path in path_set:
                file.write(",".join(path) + "\n")
    