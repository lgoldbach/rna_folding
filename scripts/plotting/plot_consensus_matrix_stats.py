#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt
import pickle

def list_of_strings(arg):
    return arg.split(',')

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="2D consensus ranking matrix in "
                        "pickle format (numpy 2D array) ", required=True, 
                        type=str)
    parser.add_argument("-o", "--output", help="Output file name for plot"
                        "(should end in .pdf)", required=True, type=str)
    
    args = parser.parse_args()

    A = pickle.load(open(args.input, "rb"))

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
                v = sorted([A[i,j],A[j, i]])  # sort so consistency is in [0,1]
                # compute the consistency as 1 minus the ratio
                unbalanced_list.append(1-(v[0]/v[1]))
                unbalanced_list_count.append(sum(v))
                unbalanced += 1

    # combine all except the unmatches pairs
    total_balance = unbalanced_list + ([1]*balanced)  


    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))
    axes[0].hist(total_balance)
    axes[0].set_xlabel("Consistency of pairwise ranking")
    axes[0].set_ylabel("Count")
    axes[0].title.set_text('Histogram of consistency scores')

    axes[1].hist(total_balance, cumulative=1)
    axes[1].set_xlabel("Consistency")
    axes[1].set_ylabel("Cumulative count")
    axes[1].title.set_text('Cumulative count of consistency scores')

    axes[2].scatter(unbalanced_list_count + balanced_list_count, unbalanced_list + ([1]*balanced))
    axes[2].set_xlabel("How often the two phenotypes appear\nin same suboptimal set")
    axes[2].set_ylabel("Consistency of pairwise ranking")
    axes[2].title.set_text('Consistency over frequency of match-ups')

    plt.tight_layout()

    plt.savefig(args.output, format="pdf", dpi=30)
