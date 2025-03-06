#!/usr/bin/env python

import argparse
import numpy as np
from typing import Iterable

from rna_folding.utils import sequence_to_integers

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, type=str, 
                    help="genotype file, one genotype per line")
parser.add_argument("-l", "--letters", required=True, type=Iterable, 
                    help="Letters used in the input genotypes")
parser.add_argument("-o", "--out", required=False, type=str,
                    help=".npy file holding array of integer genotypes")

args = parser.parse_args()

seq_arr = np.loadtxt(args.input, dtype=str)
# turn genotypes strings to list of letters
seq_arr = np.array([list(seq) for seq in seq_arr])

seq_arr_num = sequence_to_integers(seq_arr, letters=args.letters)  # encode letters as integers

np.save(args.output, seq_arr_num)
