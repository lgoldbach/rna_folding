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
    parser.add_argument("-b", "--bp_rule", help="bp rule number ", required=False)
    parser.add_argument("-k", "--ref_bp_rule", help="ref bp rule number ", required=False)
    parser.add_argument("-l", "--log", action="store_true")
    
    args = parser.parse_args()
    fig, ax = plt.subplots(figsize=(5, 4))

    marker_types=[".", "s", "D"]
    linestyles=["solid", "dashed", "dotted"]

    for i, file in enumerate(args.phenotype_dist):
        phenotypes, distr = ranked_ph_distribution(ph_distr_file=file,
                                                    log=args.log)
        distr = distr[1:]  # remove unfolded phenotype frequency
        x = range(1, distr.shape[0]+ 1)
        
        sc = ax.plot(x, distr, marker="", label=f"Base-pairing 4", color="tab:blue", linewidth=2, linestyle=linestyles[i], zorder=10)

    
    for i, file in enumerate(args.ref):
        phenotypes, distr = ranked_ph_distribution(ph_distr_file=file,
                                                    log=args.log)

        distr = distr[1:]  # remove unfolded phenotype frequency
        x = range(1, distr.shape[0]+1)
        
        sc = ax.plot(x, distr, marker="", label=f"ViennaRNA", color="black", linewidth=2, linestyle=linestyles[i])

    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.tick_params(axis='both', which='minor', labelsize=8)
    ax.set_ylim([-7.35, -1])
    ax.set_xlim([1, 75])
    ax.set_xlabel("Rank")
    ax.grid()
    ax.legend(loc="upper right", prop={'size': 10}, frameon=True)

    if args.log:
        ax.set_ylabel("Phenotype frequency (log10)")
    else:
        ax.set_ylabel("Phenotype frequency")
        
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
