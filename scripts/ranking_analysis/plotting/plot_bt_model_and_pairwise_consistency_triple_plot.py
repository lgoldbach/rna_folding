#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import pickle
import seaborn as sns

from rna_folding.parsing import load_phenotype_and_metric_from_file


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--bt_scores", help="Phenotypes and scores in two "
                        "space-separated columns", 
                        required=True)
    parser.add_argument("-d", "--ignore", help="Phenotype to ignore", type=str,
                        required=False)
    parser.add_argument("-c", "--consensus_matrix", help="2D consensus ranking matrix in "
                        "pickle format (numpy 2D array) ", required=True, 
                        type=str)
    parser.add_argument("-p", "--phenotypes", help="List of phenotypes", required=True)    
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()


    
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=1, ncols=3, figsize=(21, 7), width_ratios=[1, 1, 1])
    # plt.subplots(layout="constrained")
    
    plt.subplots_adjust(top=.9, bottom=.1, left=.1, right=.9)  # make space for labels

    A = pickle.load(open(args.consensus_matrix, "rb"))

    # consensus matrix plot
    phenotypes = np.loadtxt(args.phenotypes, dtype=str)


    if args.ignore:
        # find index of phenotype to ignore
        ignore_idx = np.where(phenotypes==args.ignore)[0][0]
        A = np.delete(A, (ignore_idx), axis=0)  # delete respective row
        A = np.delete(A, (ignore_idx), axis=1)  # delete respective column
        # remove from phenotypes
        phenotypes = np.array([ph for ph in phenotypes if ph != args.ignore])

    A_ratio = np.zeros_like(A)
    for i in range(A.shape[0]):
        for j in range(i+1, A.shape[1]):
            sum_ = A[i, j] + A[j, i]
            if A[i, j]:
                A_ratio[i, j] = A[i, j]/sum_
            if A[j, i]:
                A_ratio[j, i] = A[j, i]/sum_
    
    order = np.argsort(np.sum(A_ratio, axis=1))[::-1]
    
    ph_ordered = phenotypes[order]

    # bring both axes in order
    A_ratio = A_ratio[order, :]
    A_ratio = A_ratio[:, order[::-1]]

    A_ratio = np.flip(A_ratio, axis=1)  # mirror horizontally


    sns.heatmap(A_ratio, ax=ax1, square=True, cmap="YlGnBu", cbar=True)
    
    # make outlines visible
    for _, spine in ax1.spines.items():
        spine.set_visible(True)

    ax1.tick_params(axis='x', which='major', labelrotation=30, labelsize=15, labelbottom = False, bottom=False, top = True, labeltop=True)
    ax1.tick_params(axis='y', which='major', labelsize=5, labelrotation=0)

    ph_to_show = [0, 4, 8, 14, 22, 29, 38, 46]
    ax1.set_xticks([i + .5 for i in ph_to_show])
    ax1.set_yticks([i + .5 for i in ph_to_show])

    ph_labels = [ph_ordered[i] for i in ph_to_show]
    ax1.set_xticklabels(ph_labels, horizontalalignment="left", verticalalignment="bottom", fontsize=12)
    ax1.set_yticklabels(ph_labels, fontsize=12)

    # ax1.set_xticks(fontsize=4, rotation=45, ha="right")
    # ax1.set_yticks(fontsize=6)

    # pairwise consistency

    # keep track of qualitative properties
    balanced = 0  # perfectly balanced pairs of phenotypes
    unbalanced = 0  # unbalanced pairs of phenotypes
    unmatched = 0  # phenotypes that never appear in the same suboptimal set.
    # count number of matchups for each of the phenotype pairs
    balanced_list_count = []
    unbalanced_list = []
    unbalanced_list_count = []
    # loop over upper triangle of matrix
    for i in range(A.shape[0]):
        for j in range(i+1, A.shape[1]):
            # if the two phenotypes never matched up both entries will be 0
            if A[i,j] == 0 and A[j,i] == 0:
                unmatched += 1
            # if either of them is 0 that means that they are balanced
            # note that this would require and XOR but since we use an elif
            # and checked for both being 0 before it behaves like an XOR
            elif A[i,j] == 0 or A[j,i] == 0:
                balanced += 1
                # count how many times this balanced matchup happened
                balanced_list_count.append(sum([A[j,i], A[i, j]]))
            # both entries are non-zero, now we keep track of the ratio.
            else:
                v = max([A[i,j],A[j, i]])  # get max value
                frac = v/sum([A[i,j],A[j, i]]) # get fraction
                c = (frac-0.5)/0.5  # normalize to [0,1]
                # compute the consistency as 1 minus the ratio
                unbalanced_list.append(c)
                unbalanced_list_count.append(sum([A[i,j],A[j, i]]))
                unbalanced += 1

    # combine all except the unmatches pairs
    total_balance = unbalanced_list + ([1]*balanced)  
    count_above_50 = 0
    for c in total_balance:
        if c > 0.75:
            count_above_50 += 1
    print(count_above_50, len(total_balance), count_above_50/len(total_balance))


    # s = sum(total_balance)
    # total_balance_norm = [i/s for i in total_balance]
    ax2.hist(total_balance, color=".4", edgecolor='white')
    ax2.set_xlabel("Pairwise consistency", fontsize=18)
    ax2.set_ylabel("Count", fontsize=18)

    ax2.tick_params(axis='both', which='major', labelsize=12)

    # BT scores
    # ignore these phenotypes because we already load phenotypes above
    # and remove the ignored phenotype if applicable

    phenotypes_, scores = load_phenotype_and_metric_from_file(args.bt_scores)

    phenotypes_ = phenotypes_[1:]
    scores = scores[1:]

    # scores = list(scores)
    # if args.ignore:
    #     for i, (p, s) in enumerate(zip(phenotypes, scores)):
    #         print(p, s)
    #         if p == args.ignore:
    #             p = phenotypes.pop(i)
    #             scores.pop(i)
    #             print("I")

    ma =  max(scores)
    scores = [s/ma for s in scores]
    # ax.bar(range(len(scores)), np.log10(scores))
    ax3.bar(range(len(scores)), scores, color=".4")

    left, bottom, width, height = [0.81, 0.3, 0.08, 0.23]
    ax4 = fig.add_axes([left, bottom, width, height])
    ax4.bar(range(30, len(scores)), scores[30:], color=".4")
    ax4.tick_params(axis='both', which='major', labelsize=20)
    # ax3.set_xticks(range(len(scores)))
    # ax3.set_xticklabels(phenotypes)
    # ax.tick_params(axis='x', labelrotation=60, labelsize=8)
    # ax3.xticks(fontsize=10, rotation=55, ha="right")

    ax3.tick_params(axis='both', which='major', labelsize=20)


    ax3.set_xlabel("Rank", fontsize=28)
    ax3.set_ylabel("mfe-score (a.u.)", fontsize=28)

    plt.tight_layout(rect=[0.1, 0.1, 0.9, 0.9])
    plt.savefig(args.output, format="pdf", dpi=30)

    plt.close()


     # plot colorbar separately
    fig_c, ax_c = plt.subplots(figsize=(7, 7))
    sns.heatmap(A_ratio, ax=ax_c, square=True, cmap="YlGnBu", cbar_kws={'shrink': 1, 'label': 'Fraction of genotypes where phenotype i ranks above phenotype j'}, linecolor="white", linewidth=.05)
    ax_c.set_xlabel("Phenotype i    Phenotype j")
    plt.savefig("colorbar.pdf", format="pdf", dpi=30)