#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.ticker import MaxNLocator
import numpy as np
import seaborn as sns


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--walks", help="Input files with adaptive "
                        "walk length for low selection pressure of following "
                        "format: (((...))) 5 2 1 4 2 4 ", 
                        required=True, nargs='+')
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

    def get_walk_lengths_and_navigability(input):
        walk_lengths_all = []
        navigability = []
        for i, file_path in enumerate(input, start=2):
            navigability.append([])  # add a list for each bp rule
            walk_lengths_all.append([])
            # if i == 10:
            #     continue
            with open(file_path, "r") as file:
                for line in file:
                    walk_lengths_per_ph = []
                    aborted = 0
                    navigability
                    lengths = [int(length) for length in line.strip().split(" ")[1:]]
                    for l in lengths:
                        if l > -1:
                            walk_lengths_per_ph.append(l)
                            walk_lengths_all[-1].append(l)
                        else:
                            aborted += 1

                    navigability[-1].append(int(np.round(len(walk_lengths_per_ph)/(aborted+len(walk_lengths_per_ph))*100, 0)))  # for each phenotype add a success percentage to the bp rules sublist
        return walk_lengths_all, navigability
    
    cmap = colormaps['Oranges']
    colors = [cmap(.5), cmap(.75), cmap(.99)]
    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]
    
    walk_lengths, navigability = get_walk_lengths_and_navigability(args.walks)

    walk_length_means = [np.mean(lengths) for lengths in walk_lengths]
    walk_length_medians = [np.median(lengths) for lengths in walk_lengths]
    
    # for i, id in enumerate(new_order_idx, start=1):
    #     ax1.scatter([i]*len(navigability[id]), navigability[id])

    # ax1.violinplot([navigability[i] for i in new_order_idx], positions=range(1, 11), showmeans=True, showextrema=False, showmedians=False)
    nav = []
    x = []
    for i, id in enumerate(new_order_idx, start=1):
        for j in navigability[id]:
            nav.append(j)
            x.append(i)

    ax = sns.boxplot(x=x,
                y=nav, 
                ax=ax1,
                color="0.9",
                linewidth=0.2,
                linecolor="black",
                legend=False,
                showfliers=False,
                zorder=3,
                whis=[0, 95])
    
    # ax = sns.violinplot(x=x,
    #                y=nav, 
    #                ax=ax1,
    #                legend=False,
    #                density_norm="area",
    #                width=0.95,
    #                common_norm=True,
    #                cut=0,
    #                inner="quart",
    #                linewidth=0.2,
    #                linecolor="black",
    #                color="0.9",
    #                zorder=3,
    #                inner_kws={"zorder": 4})
    
    for l in ax.lines:
        l.set_linestyle('-')
        l.set_color('black')
        l.set_linewidth(0.5)


    for i, id in enumerate(new_order_idx):
        d = navigability[id]
        q1 = np.percentile(d, q=25)
        q3 = np.percentile(d, q=75)
        mean = np.mean(d)
        median = np.median(d)

        # ax1.vlines(i, q1, q3, color="black", linewidth=1, zorder=5)
        if i == 1:  # only add label once for legend
            ax1.scatter(i, mean, color="black", marker="s", s=4, zorder=10, label="mean") 
        else:
            ax1.scatter(i, mean, color="black", marker="s", s=4, zorder=10) 


    # ax2.scatter(range(1, 11), [walk_length_means[i] for i in new_order_idx], alpha=0.8)
    # ax2.scatter(range(2, 12), walk_length_medians, color=colors[i], label=labels[i])
    # ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

    walks = []
    x = []
    for i, id in enumerate(new_order_idx, start=1):
        for j in walk_lengths[id]:
            walks.append(j)
            x.append(i)

    ax = sns.boxplot(x=x,
                y=walks, 
                ax=ax2,
                color="0.9",
                linewidth=0.2,
                linecolor="black",
                legend=False,
                showfliers=False,
                zorder=3,
                whis=[0, 95])
    
    for i, id in enumerate(new_order_idx):
        d = walk_lengths[id]
        mean = np.mean(d)
        median = np.median(d)
        ax2.scatter(i, mean, color="black", marker="s", s=4, zorder=10) 


    ax1.set_ylabel("Navigability", fontsize=15)
    ax1.set_xlabel("Base-pairing rule", fontsize=15)
    ax2.set_ylabel("Adaptive walk length of successful walks", fontsize=15)
    ax2.set_xlabel("Base-pairing rule", fontsize=15)

    ax1.tick_params(axis='both', which='major', labelsize=13)
    ax1.tick_params(axis='both', which='minor', labelsize=8)
    ax2.tick_params(axis='both', which='major', labelsize=13)
    ax2.tick_params(axis='both', which='minor', labelsize=8)

    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_xticks(list(range(0, 10)))
    ax1.set_xticklabels(list(range(1, 11)))

    ax1.grid(axis="y", zorder=-1)
    ax2.grid(axis="y", zorder=-1)

    # l = ax2.legend(title="Selection pressure", prop={'size': 12})
    # plt.setp(l.get_title(),fontsize=12)

    # ax1.legend()

    # ax1.set_ylim(0, 100)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
