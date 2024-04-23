#!/usr/bin/env python

"""Count the number of genotypes per phenotype and vice versa and save results
as pickled dictionary.
"""

import argparse
import pickle

from rna_folding.analysis import count_gt_per_ph_and_ph_per_gt


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--gp_map_file", help="Path to gp_map file.", 
                        required=True)
    parser.add_argument("-s", "--separator", help="Column separator", 
                        required=False, default=" ")
    parser.add_argument("--output_gt_per_ph", help="Path to output file", 
                        required=True)
    parser.add_argument("--output_ph_per_gt", help="Path to output file", 
                        required=True)
  
    args = parser.parse_args()

    gt_per_ph, ph_per_gt = count_gt_per_ph_and_ph_per_gt(gp_map_file=args.gp_map_file,
                                  sep=args.separator)
    
    pickle.dump(gt_per_ph, open(args.output_gt_per_ph, "wb"))
    pickle.dump(ph_per_gt, open(args.output_ph_per_gt, "wb"))
    