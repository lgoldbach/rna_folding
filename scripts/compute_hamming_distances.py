#!/usr/bin/env python

import argparse
import numpy as np
from scipy.spatial.distance import hamming


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, type=str, 
                    help="genotype file, one genotype per line")
parser.add_argument("-o", "--out", required=False, type=str,
                    help="Name of output file, ideally something like "
                    "phenotypes.txt. Contains all phenotypes in order, one "
                    "per line")

args = parser.parse_args()

# read data and create empty array for hamming distances
seq_arr_num = np.load(input)
num_of_seq = seq_arr_num.shape[0]
seq_len = seq_arr_num.shape[1]
hamming_dist = np.zeros((num_of_seq, num_of_seq), dtype=int)

# compute pairwise hamming distances
for i in range(num_of_seq):
    for j in range(i):
        # multiply by sequence length to get hamming distance not as fraction
        # but as mutational step. Default of scipy hamming is fraction [0, 1].
        hamming_dist[i, j] = int(hamming(seq_arr_num[i], seq_arr_num[j]) * seq_len)

np.save(output, hamming_dist)
