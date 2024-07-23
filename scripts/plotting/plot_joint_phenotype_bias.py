#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt

from rna_folding.utils import ranked_ph_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--phenotype_dist", help="Path to phenotype "
                        "distribution files", required=True, nargs='+')
    parser.add_argument("-r", "--ref", help="Reference phenotype distribution  "
                        "files ", required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    parser.add_argument("-l", "--log", action="store_true")
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    phenotypes, distr = ranked_ph_distribution(ph_distr_file=args.ref,
                                                    log=args.log)
    x = range(distr.shape[0])

    ax.plot(x, distr, color="black", marker="x", markersize=4, label="ViennaRNA")

    for i, file in enumerate(args.phenotype_dist):
        phenotypes, distr = ranked_ph_distribution(ph_distr_file=file,
                                                    log=args.log)


        x = range(distr.shape[0])
        
        sc = ax.plot(x, distr, marker=".", label=f"Random ranking {i + 1}")
        print(sc[0].get_color())

    
    ax.set_xlabel("Rank")
    if args.log:
        ax.set_ylabel("Phenotype frequency (log10)")
    else:
        ax.set_ylabel("Phenotype frequency")
        
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
