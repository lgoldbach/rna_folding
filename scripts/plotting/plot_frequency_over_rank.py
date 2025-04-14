#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.utils import ranked_ph_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--phenotype_dist", help="Path to phenotype "
                        "distribution files", required=True)
    parser.add_argument("-r", "--ranking", help="Phenotype ranking", 
                        required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    ph, freq = ranked_ph_distribution(ph_distr_file=args.phenotype_dist)
    ph_to_freq = dict(zip(ph, freq))

    ranking = np.loadtxt(args.ranking, dtype=str)

    freq_by_rank = [ph_to_freq[ph] for ph in ranking[:-1]]

    x = range(len(freq_by_rank))

    ax.bar(x, freq_by_rank)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
