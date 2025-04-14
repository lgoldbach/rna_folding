#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input files with phenotype in "
                        "first column followed by mfe values in following "
                        "columns. One row for each phenotypes", 
                        required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    mfes = {}
    with open(args.input, "r") as f:
        for line_ in f:
            line = line_.strip().split(" ")
            mfes[line[0]] = [float(i) for i in line[1:]]  #line[0] is phenotype

    phenotypes = []
    means = []
    for ph in mfes:
        phenotypes.append(ph)
        means.append(np.mean(mfes[ph]))

    phenotypes = np.array(phenotypes)
    means = np.array(means)
    order = np.argsort(means)

    means_sort = means[order]
    phenotypes_sort = phenotypes[order]

    data = []
    for ph in phenotypes_sort:
        data.append(mfes[ph])
        
    sns.violinplot(data)  #, showmeans=True, showmedians=False)

    print(phenotypes_sort)

    ax.set_xlabel("Phenotypes")
    ax.set_ylabel("Minimum free energy")

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
