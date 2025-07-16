#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.stats import pearsonr
import copy
import networkx as nx

from rna_folding.adaptive_walks import read_genotype_paths_from_file, write_paths_to_file
from rna_folding.parsing import many_to_one_map_from_file_to_dict


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--paths", help="File with adaptive walks", required=True)
    parser.add_argument("-g", "--nc_to_gt", help="Map from neutral components to genotypes", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)", required=True)
    
    args = parser.parse_args()

    # Dictionary to map genotypes to neutral component id
    gt_to_nc = many_to_one_map_from_file_to_dict(file=args.nc_to_gt, 
                                                 source_type=str, 
                                                 target_type=str, 
                                                 delimiter=" ")

    paths = read_genotype_paths_from_file(file=args.paths,
                                  delimiter=" ",
                                  gt_type="str",
                                  map=gt_to_nc)
    
    write_paths_to_file(paths, args.output, delimiter=" ")
    