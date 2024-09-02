#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from rna_folding.utils import ranked_ph_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input files with navigability "
                        "value in last column (space-serparated)", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    navig=[]
    for i, file in enumerate(args.input, start=2):
        with open(file, "r") as f:            
            navig.append([float(line.strip().split(" ")[-1]) for line in f])  # get value from last columns
    sns.boxplot(navig)
    
    ax.set_xlabel("Base-pairing rule")
    ax.set_ylabel("Navigability")
    
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
