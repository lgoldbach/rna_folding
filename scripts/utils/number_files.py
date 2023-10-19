#!/usr/bin/env python

"""Rename files to prefix + ascending number

e.g.:
input files: a b ...

number_files.py --prefix gt_ --suffix .txt --start_one --files a b ...

output files: gt_1.txt g1_2.txt ...

"""

import sys
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--prefix", required=True, type=str, 
                        help="file name prefix before the suffix {i}.txt, e.g."
                        "prefix 'RNA12_gt' will become RNA12_gt1.txt")

parser.add_argument("-s", "--suffix", required=False, type=str, default="",
                        help="file name suffix")

parser.add_argument("-f", "--files", required=True, nargs='+', type=str, 
                        help="List of files to be renamed")

parser.add_argument("-o", "--start_one", required=False,
                    action='store_true',
                    default=False, help="Start list at one")

args = parser.parse_args()
for i, j in enumerate(args.files, args.start_one):
    os.rename(j, args.prefix + str(i) + ".txt")
