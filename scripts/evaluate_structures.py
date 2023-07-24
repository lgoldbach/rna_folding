import argparse
import numpy as np
import matplotlib.pyplot as plt
import pickle

from rna_folding.evaluate import f1_score
from rna_folding.parsing import gpmap_to_dict


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genotypes", required=True, type=str, help="List of genotypes")
    parser.add_argument("-p", "--phenotypes", required=True, type=str, help="List of phenotypes and genotype IDs that map to them")
    parser.add_argument("-r", "--reference", required=True, type=str, help="Reference RNA secondary structure in dot-bracket notation")
    parser.add_argument("-a", "--abstract", required=False, type=int, help="Level of abstraction on which to compare structures")

    args = parser.parse_args()
    
    gp_map = gpmap_to_dict(args.phenotypes, args.genotypes)
    
    with open(args.reference, "r") as ref:
        ref_dict = dict([line.strip().split() for line in ref])
    
    f1_scores = []
    for query_seq in gp_map:
        ref_ss = ref_dict[query_seq]
        # loop over suboptimal strucs and compare to ref
        f1 = [f1_score(ref_ss, query_ss) 
                for query_ss in gp_map[query_seq]]  
        f1_scores.append(f1)
            
    f1_score_avg = [np.mean(scores) for scores in f1_scores]

    f1_score_avg_no_zero = [i for i in f1_score_avg if i > 0]
    print(len(f1_score_avg) - len(f1_score_avg_no_zero))
    plt.hist(f1_score_avg_no_zero)
    plt.savefig('hist_no_zeros.png')
        