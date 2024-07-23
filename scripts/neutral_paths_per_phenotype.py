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
    parser.add_argument("-p", "--paths", 
                        help="File output for neutral paths",
                        required=True)
    
    args = parser.parse_args()

    print("A", flush=True)
    gpm = pickle.load(open(args.file, "rb"))
    print("B", flush=True)
    with open(args.paths, "w") as file:
        for ph in gpm.phenotype_set:
            print(ph, flush=True)
            file.write(ph)  # add first line to csv
            genotypes = gpm.genotypes_of_phenotype(ph)
            # if larger sample requested than genotypes use len(genotype)
            sample_size = min(len(genotypes), args.num_of_genotypes)
            genotype_sample = np.random.choice(genotypes, 
                                               size=sample_size, 
                                               replace=False)
            neutral_paths = []
            for g in genotype_sample:
                print(g, flush=True)
                neutral_paths.append(gpm.neutral_paths(g, n=args.num_of_paths))
                print(neutral_paths[-1], flush=True)
            
            for path_set in neutral_paths:
                for path in path_set:
                    # add path length to output file
                    file.write(" " + str(len(path)))
            file.write("\n")  # new line for next phenotypes
        file.close()
