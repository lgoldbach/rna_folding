#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt

from rna_folding.utils import ranked_ph_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--phenotype_dist", help="Path to phenotype "
                        "distribution files", required=True, nargs='+')
    parser.add_argument("-r", "--ref", help="Reference phenotype distribution  "
                        "files ", required=True, nargs="+")
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    parser.add_argument("-b", "--bp_rule", help="bp rule number ", required=True)
    parser.add_argument("-k", "--ref_bp_rule", help="ref bp rule number ", required=True)
    parser.add_argument("-l", "--log", action="store_true")
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    marker_types=[".", "s", "D"]
    linestyles=["solid", "dashed", "dotted"]

    for i, file in enumerate(args.phenotype_dist):
        phenotypes, distr = ranked_ph_distribution(ph_distr_file=file,
                                                    log=args.log)


        x = range(distr.shape[0])
        
        sc = ax.plot(x, distr, marker="", label=f"Base-pairing {args.bp_rule}, Ranking {i + 1}", color="tab:blue", markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])

    
    for i, file in enumerate(args.ref):
        phenotypes, distr = ranked_ph_distribution(ph_distr_file=file,
                                                    log=args.log)


        x = range(distr.shape[0])
        
        sc = ax.plot(x, distr, marker="", label=f"Base-pairing {args.ref_bp_rule}, Ranking {i + 1}", color="black", markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])

    
    ax.set_xlabel("Rank")
    if args.log:
        ax.set_ylabel("Phenotype frequency (log10)")
    else:
        ax.set_ylabel("Phenotype frequency")
        
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
