#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.ticker import MaxNLocator
import numpy as np
import seaborn as sns

from rna_folding.parsing import read_navigability_per_ph_per_fl_file


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--navigability", help="Input files with " \
                        "navigability values for each phenotype for each fitness landscape"
                        "format: (((...))) 5 2 1 4 2 4 ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    
    cmap = colormaps['Oranges']
    colors = [cmap(.5), cmap(.75), cmap(.99)]
    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]
    
    navigs = [read_navigability_per_ph_per_fl_file(file) for file in args.navigability]
    
    # here we just assemble all navigabilities for each phenotype on each fl
    nav = []
    nav_per_bp = []
    x = []
    # for i, id in enumerate(new_order_idx, start=1):
    #     nav_per_bp.append([])
    #     for ph in navigs[id]:
    #         for n in navigs[id][ph]:
    #             nav.append(n)  # append value for every phenotype for every fl
    #             x.append(i)

    #             nav_per_bp[-1].append(n)  # make list of lists 


    # Take the mean of navigabilities for each phenotype over every fl
    nav_per_bp_per_fl = []
    nav_per_bp = []
    for i, id in enumerate(new_order_idx, start=1):
        nav_per_bp.append([])
        for navigs_over_fls in navigs[id].values():  # per ph navigs
            navig_per_ph_mean_over_fl = np.mean(navigs_over_fls)  # mean per ph over fls
            nav_per_bp_per_fl.append(navig_per_ph_mean_over_fl)
            x.append(i)

            nav_per_bp[-1].append(navig_per_ph_mean_over_fl)
        

    palette = ["0.9" if i != 4 else "0.6" for i in range(1, len(x))]
    ax = sns.boxplot(x=x,
                y=nav_per_bp_per_fl, 
                ax=ax,
                color="0.9",
                linewidth=0.2,
                linecolor="black",
                legend=False,
                showfliers=False,
                zorder=3,
                palette=palette,
                whis=[0, 95])    
    
    # ax = sns.violinplot(x=x,
    #                y=nav,
    #                palette=palette,
    #                ax=ax,
    #                legend=False,
    #                density_norm="area",
    #                width=0.95,
    #                common_norm=True,
    #                cut=0,
    #                inner="quart",
    #                linewidth=0.2,
    #                linecolor="black",
    #                zorder=3,
    #                inner_kws={"zorder": 4})
    
    for l in ax.lines:
        l.set_linestyle('-')
        l.set_color('black')
        l.set_linewidth(0.5)


    for i, id in enumerate(new_order_idx):
        d = nav_per_bp[i]
        q1 = np.percentile(d, q=25)
        q3 = np.percentile(d, q=75)
        mean = np.mean(d)
        median = np.median(d)

        # ax1.vlines(i, q1, q3, color="black", linewidth=1, zorder=5)
        if i == 1:  # only add label once for legend
            ax.scatter(i, mean, color="black", marker="s", s=4, zorder=10, label="mean") 
        else:
            ax.scatter(i, mean, color="black", marker="s", s=4, zorder=10) 



    ax.set_ylabel("Phenotype navigability\naveraged over fitness lanscapes (%)", fontsize=15)
    ax.set_xlabel("Base-pairing rule", fontsize=15)

    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.tick_params(axis='both', which='minor', labelsize=8)

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xticks(list(range(0, 10)))
    ax.set_xticklabels(list(range(1, 11)))

    ax.grid(axis="y", zorder=-1)

    # l = ax2.legend(title="Selection pressure", prop={'size': 12})
    # plt.setp(l.get_title(),fontsize=12)

    # ax.legend(loc="lower center", frameon=False, fancybox=False)

    ax.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
