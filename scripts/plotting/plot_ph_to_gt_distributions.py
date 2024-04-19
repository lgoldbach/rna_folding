#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_gt_per_ph", help="pickle file with dict containing genotype per phenotype counts",
                        required=True)
    parser.add_argument("--input_ph_per_gt", help="pickle file with dict containing phenotype per genotype counts",
                        required=True)
    parser.add_argument("-s", "--separator", help="Column separator", 
                        required=False, default=" ")
    parser.add_argument("-r", "--bp_rule", help="Number of base pair rule", 
                        required=True)
    parser.add_argument("-o", "--output", help="Path to output file", 
                        required=True)
  
    
    args = parser.parse_args()

    gt_per_ph = pickle.load(open(args.input_gt_per_ph, "rb"))
    ph_per_gt = pickle.load(open(args.input_ph_per_gt, "rb"))

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)

    ax1.hist(gt_per_ph.values())
    ax2.hist(ph_per_gt.values())

    ax1.set_xlabel("Genotypes per phenotype")
    ax1.set_ylabel("Count")
    ax2.set_xlabel("Phenotypes per genotypes")
    ax2.set_ylabel("Count")
    ax2.title.set_text(f"bp_rule{args.bp_rule}")

    plt.tight_layout()
    plt.savefig(args.output, format='pdf')

    