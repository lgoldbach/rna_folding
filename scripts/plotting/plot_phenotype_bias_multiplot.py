#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.utils import ranked_ph_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ph_dist", help="Paths to phenotype "
                        "distribution files", required=True, nargs='+')
    parser.add_argument("-r", "--ref", help="Reference bp rule", required=True,
                        type=int)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    

    args = parser.parse_args()
    args.ref = 4

    fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(20, 20), sharey=True, sharex=True)
    axes = axes.flatten()
    
    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]
    ph_dist_files = np.array(args.ph_dist)[new_order_idx]

    adjust_i = 0  # absolutely disgusting hack to skip reference bp rule but still grep correct ax
    for i, file in enumerate(ph_dist_files):  # 2 because rules start at 2
        if i == args.ref-1:
            adjust_i = -1
            continue

        hacky_i = i+adjust_i

        ax = axes[hacky_i]
        if hacky_i in [3]:
            ax.set_ylabel("Phenotype frequency (log10)", fontsize=40)
        else:
            for tick in ax.yaxis.get_major_ticks():
                tick.tick1line.set_visible(False)
                tick.tick2line.set_visible(False)
                # tick.label1.set_visible(False)
                # tick.label2.set_visible(False)     
        if hacky_i == 7:
            ax.set_xlabel("Phenotypes", fontsize=40)

        ax.set_ylim([-7, -1])
        ax.set_xlim([1, 37])
        ax.grid(axis='y')
        ax.grid(axis='x')

        phenotypes, distr = ranked_ph_distribution(ph_distr_file=file,
                                                    log=True)
        distr = distr[1:]
        x = range(1, distr.shape[0]+1)
     
        ax.plot(x, distr, label=f"Base-pairing {i+1}", color="black", linewidth=7, zorder=0)

    ref_file = ph_dist_files[args.ref-1]
    phenotypes, ref_distr = ranked_ph_distribution(ph_distr_file=ref_file,
                                                   log=True)
    ref_distr = ref_distr[1:]   # ignore unfolded
    ref_x = range(1, ref_distr.shape[0]+1)
    for ax in axes:
        ax.plot(ref_x, ref_distr, color="0.4", linewidth=7, label="Natural base-pairing", linestyle="dotted", zorder=1)
        ax.legend(loc="lower left", prop={'size': 25}, frameon=False)
    
        ax.tick_params(axis='both', which='major', labelsize=20)
        ax.tick_params(axis='both', which='minor', labelsize=8)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
