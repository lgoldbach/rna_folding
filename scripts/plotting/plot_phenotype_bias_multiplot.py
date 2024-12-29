#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt

from rna_folding.utils import ranked_ph_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ph_dist", help="Paths to phenotype "
                        "distribution files", required=True, nargs='+')
    parser.add_argument("-r", "--ref", help="Reference phenotype distribution  "
                        "files ", required=True, type=int)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, axes = plt.subplots(nrows=1, ncols=len(args.ph_dist)-1)

    ref_file = args.ph_dist[args.ref-2]  # -2 because rules start counting at 2
    phenotypes, ref_distr = ranked_ph_distribution(ph_distr_file=ref_file,
                                                log=args.log)
    ref_x = range(ref_distr.shape[0])

    for i, file in enumerate(args.ph_dist, start=2):  # 2 because rules start at 2
        if i == args.ref:
            continue

        phenotypes, distr = ranked_ph_distribution(ph_distr_file=file,
                                                    log=args.log)

        x = range(distr.shape[0])
        
        sc = axes[i].plot(x, distr, marker="", label=f"Base-pairing {args.bp_rule}, Ranking {i + 1}", color="tab:blue", markevery=(i*3, 10), markersize=5, linestyle=linestyles[i])

    
    ax.set_xlabel("Rank")
    if args.log:
        ax.set_ylabel("Phenotype frequency (log10)")
    else:
        ax.set_ylabel("Phenotype frequency")
        
    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
