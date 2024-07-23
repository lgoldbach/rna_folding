#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt

from rna_folding.utils import ranked_ph_distribution


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--phenotype_dist", help="Path to phenotype "
                        "distribution file", required=True, type=str)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    parser.add_argument("-l", "--log", action="store_true")
    
    args = parser.parse_args()
    phenotypes, distr = ranked_ph_distribution(ph_distr_file=
                                                   args.phenotype_dist,
                                                   log=args.log)


    x = range(distr.shape[0])


    print(phenotypes)
    fig, ax = plt.subplots()
    ax.bar(x, distr)
    ax.set_xticks(x, phenotypes, rotation=45, ha='right')
    ax.set_xlabel("Phenotypes")
    if args.log:
        ax.set_ylabel("Phenotype frequency (log10)")
    else:
        ax.set_ylabel("Phenotype frequency")
        
    ax.set_title("Phenotype distribution")
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
