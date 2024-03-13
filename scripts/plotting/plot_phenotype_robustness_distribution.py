#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.utils import load_phenotype_and_metric_from_file

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="File with phenotype distribution and robustness "
                        "file", required=True)
    parser.add_argument("-r", "--ref_file", help="File with reference phenotype distribution and robustness "
                        "file", required=True)
    parser.add_argument("-b", "--bp_rule", help="Name of the bp_rule",
                        required=True)
    parser.add_argument("-i", "--ref_bp_rule", help="Name of the ref bp_rule",
                        required=True)
    parser.add_argument("-o", "--output", help="File output for plot",
                        required=True)

    args = parser.parse_args()

    fig, (ax1, ax2, ax3)  = plt.subplots(nrows=3, figsize=(4, 5))
    phenotype_counts = []
    for bp_rule, file_name in zip([args.bp_rule, args.ref_bp_rule], [args.file, args.ref_file]):
        data = np.loadtxt(file_name, delimiter=' ')
        n, bins, patches = ax1.hist(data[:, 0], alpha=0.5, label=bp_rule)
        
        phenotype_counts.append(data[:, 0].shape[0])

        ax1.vlines(np.mean(data[:, 0]), ymin=0, ymax=10, color=patches[0].get_facecolor())
        ax1.vlines(np.median(data[:, 0]), ymin=0, ymax=10, color=patches[0].get_facecolor(), linestyles="dotted")
        ax2.hist(np.log10(data[:, 1]), alpha=0.5)
        ratio = data[:, 0]/data[:, 1]
        ratio = np.log10(ratio[ratio!=0])
        ax3.hist(ratio, alpha=0.5)
    
    ax1.set_title(f"Number of distinct phenotypes:\n{args.bp_rule}: {phenotype_counts[0]}\n{args.ref_bp_rule}: {phenotype_counts[1]}", 
              fontsize=7)
        
    ax1.set_xlabel("Phenotype robustness")
    ax2.set_xlabel("Phenotype frequency (log10)")
    ax3.set_xlabel("Phenotype robustness/Phenotype frequency (log10)")
    for ax in (ax1, ax2, ax3):
        ax.set_ylabel("Phenotype count")

    ax1.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")






            
