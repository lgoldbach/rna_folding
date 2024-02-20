#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.utils import load_phenotype_and_metric_from_file

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", help="Files containing phenotype distribution and robustness "
                        "file", required=True, nargs="+")
    parser.add_argument("-o", "--output", help="File output for plot",
                        required=True)

    args = parser.parse_args()

    fig, (ax1, ax2, ax3)  = plt.subplots(nrows=3, figsize=(4, 5))

    for file_name in args.files:
        data = np.loadtxt(file_name, delimiter=' ')
        ax1.hist(data[:, 0], alpha=0.5)
        ax2.hist(np.log10(data[:, 1]), alpha=0.5)
        ratio = data[:, 0]/data[:, 1]
        ratio = np.log10(ratio[ratio!=0])
        ax3.hist(ratio, alpha=0.5)

    ax1.set_xlabel("Phenotype robustness")
    ax2.set_xlabel("Phenotype frequency (log10)")
    ax3.set_xlabel("Phenotype robustness/Phenotype frequency (log10)")
    for ax in (ax1, ax2, ax3):
        ax.set_ylabel("Count")

    plt.tight_layout()
    plt.savefig(args.output, format="pdf")






            
