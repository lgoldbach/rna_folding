#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

from rna_folding.parsing import load_phenotype_and_metric_from_file


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Phenotypes and scores in two "
                        "space separated columns", 
                        required=True)
    parser.add_argument("-d", "--ignore", help="Phenotype to ignore", type=str,
                        required=False)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()

    phenotypes, scores = load_phenotype_and_metric_from_file(args.input)
    phenotypes = list(phenotypes)
    scores = list(scores)
    if args.ignore:
        for i, (p, s) in enumerate(zip(phenotypes, scores)):
            if p == args.ignore:
                p = phenotypes.pop(i)
                scores.pop(i)
    fig, ax = plt.subplots()

    ma =  max(scores)
    scores = [s/ma for s in scores]
    # ax.bar(range(len(scores)), np.log10(scores))
    ax.bar(range(len(scores)), scores)

    ax.set_xticks(range(len(scores)))
    ax.set_xticklabels(phenotypes)
    # ax.tick_params(axis='x', labelrotation=60, labelsize=8)
    plt.xticks(fontsize=6, rotation=55, ha="right")

    ax.set_xlabel("Phenotypes (dot-bracket notation)")
    ax.set_ylabel("Bradley-Terry scores (a.u.)")
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
