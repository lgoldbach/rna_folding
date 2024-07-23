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

    fig, ax = plt.subplots()

    phenotype_counts = []
    for bp_rule, file_name in zip([args.bp_rule], [args.file]):
        data = np.loadtxt(file_name, delimiter=' ')

        freq = np.log10(data[:, 1])
        ax.scatter(freq, data[:, 0], s=10, alpha=0.5, label=f"Base-pairing {args.bp_rule}", color="tab:blue")
    
    for bp_rule, file_name in zip([args.ref_bp_rule], [args.ref_file]):
        data = np.loadtxt(file_name, delimiter=' ')

        freq = np.log10(data[:, 1])
        ax.scatter(freq, data[:, 0], s=10, alpha=0.5, label=f"Base-pairing {args.ref_bp_rule}", color="black")
    
    ax.set_xlabel("Phenotype frequency (log10)")
    ax.set_ylabel("Phenotype robustness")

    ax.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf")






            
