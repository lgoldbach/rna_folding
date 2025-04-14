#!/usr/bin/env python
import argparse
import numpy as np

from rna_folding.parsing import gpmap_to_lists

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, type=str, 
                    help="GP map file in p-g format")
parser.add_argument("-o", "--out", required=False, type=str,
                    help="Name of output file, ideally something like "
                    "phenotypes.txt. Contains all phenotypes in order, one "
                    "per line")

args = parser.parse_args()

gt, ph = gpmap_to_lists(args.input)

with open(args.out, "r") as f:
    for p in ph:
        f.write(p + "\n")
