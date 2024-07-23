#!/usr/bin/env python

import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle

def list_of_strings(arg):
    return arg.split(',')

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="2D consensus ranking matrix in "
                        "pickle format (numpy 2D array) ", required=True, 
                        type=str)
    parser.add_argument("-p", "--phenotypes", help="List of phenotypes", required=True)
    parser.add_argument("-o", "--output", help="Output file name for plot"
                        "(should end in .pdf)", required=True, type=str)
    
    args = parser.parse_args()

    A = pickle.load(open(args.input, "rb"))


    fig, ax = plt.subplots()
    
    phenotypes = np.loadtxt(args.phenotypes, dtype=str)
    

    A_ratio = np.zeros_like(A)
    for i in range(A.shape[0]):
        for j in range(i+1, A.shape[1]):
            sum_ = A[i, j] + A[j, i]
            if A[i, j]:
                A_ratio[i, j] = A[i, j]/sum_
            if A[j, i]:
                A_ratio[j, i] = A[j, i]/sum_
    
    order = np.argsort(np.sum(A_ratio, axis=1))[::-1]
    
    A_ratio = A_ratio[order, :]
    A_ratio = A_ratio[:, order[::-1]]
    
    sns.heatmap(A_ratio, ax=ax, square=True,  cmap="YlGnBu")
    ax.set_xticks(np.arange(len(order))+.5)
    ax.set_yticks(np.arange(len(order))+.5)
    # print(len(phenotypes))
    # print(phenotypes[order[::-1]])
    ax.set_xticklabels(phenotypes[order[::-1]])
    ax.set_yticklabels(phenotypes[order])
    plt.xticks(fontsize=4, rotation=45, ha="right")
    plt.yticks(fontsize=6)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
