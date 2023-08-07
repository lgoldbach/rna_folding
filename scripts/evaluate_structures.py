#!/usr/bin/env python

import argparse
import numpy as np
import matplotlib.pyplot as plt

from rna_folding.evaluate import f1_score
from rna_folding.parsing import gpmap_to_dict
from rna_folding.utils import count_bp

def folding_confusion_matrix(gp_map, f1_scores, ref, L):
    """Compute confusion matrix for the two classes folded and unfolded

    Args:
        gp_map (dict): gp_map dict
        f1_scores (dict): dict of f1_scores
        ref (dict): reference dict
        L (int): len of sequence

    
    Returns:
        2x2 np.array
        [[q_uf_r_uf, q_f_r_uf],
         [q_uf_r_f, q_f_r_f]]
        legend: q = query, r = reference, uf = unfolded, f = folded
    
    """
    unfolded = "." * L

    # q = query, r = reference, uf = unfolded, f = folded
    q_uf_r_uf = 0
    q_f_r_uf = 0
    q_f_r_f = 0
    q_uf_r_f = 0

    f1_scores_by_folding = [[[], []],[[], []]]

    for seq in f1_scores:
        if gp_map[seq][0] == unfolded:
            if ref[seq] == unfolded:
                f1_scores_by_folding[0][0].append(f1_scores[seq])
                q_uf_r_uf += 1
            else:
                f1_scores_by_folding[1][0].append(f1_scores[seq])
                q_uf_r_f += 1
        else:
            if ref[seq] == unfolded:
                f1_scores_by_folding[0][1].append(f1_scores[seq])
                q_f_r_uf += 1
            else:
                f1_scores_by_folding[1][1].append(f1_scores[seq])
                q_f_r_f += 1

    folding_matches = np.array([[q_uf_r_uf, q_f_r_uf],[q_uf_r_f, q_f_r_f]])
    return folding_matches, f1_scores_by_folding


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genotypes", required=True, type=str, 
                        help="List of genotypes")
    parser.add_argument("-p", "--phenotypes", required=True, type=str, 
                        help="List of phenotypes and genotype IDs that map to "
                        "them")
    parser.add_argument("-r", "--reference", required=True, type=str, 
                        help="Reference RNA secondary structure in dot-bracket"
                        " notation")
    parser.add_argument("-a", "--abstract", required=False, type=int, 
                        help="Level of abstraction on which to compare "
                        "structures")
    
    print("Parse files ...")
    args = parser.parse_args()
    gp_map = gpmap_to_dict(args.phenotypes, args.genotypes)
    
    L = len(list(gp_map.keys())[0])

    with open(args.reference, "r") as ref:
        ref_dict = dict([line.strip().split() for line in ref])
    
    gp_map = gpmap_to_dict(args.phenotypes, args.genotypes)
    
    with open(args.reference, "r") as ref:
        ref_dict = dict([line.strip().split() for line in ref])

    print("Done!")
    print("Compute F1 scores ...\n(can take a while for suboptimal folding or " 
          "large maps)")

    # compute f1 scores and f1 score averages
    f1_scores = {}
    f1_scores_avg = {}
    for query_seq in gp_map:
        ref_ss = ref_dict[query_seq]
        # loop over suboptimal strucs and compare to ref
        f1 = [f1_score(ref_ss, query_ss) 
                for query_ss in gp_map[query_seq]]
        f1_scores[query_seq] = f1
        f1_scores_avg[query_seq] = np.mean(f1)

    # get f1 scores without zeros entries
    f1_score_avg_no_zero = [i for i in f1_scores_avg.values() if i > 0]
    print("Done!")
    
    print("Compute confusion matrix for folded/unfolded ...")
    folding_cm, f1_scores_by_folding = folding_confusion_matrix(gp_map, 
                                                                f1_scores_avg, 
                                                                ref_dict, L)
    
    print(folding_cm)
    print("Done!")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    ax1.hist([np.mean(i) for i in f1_scores_by_folding[0][0]])
    ax2.hist([np.mean(i) for i in f1_scores_by_folding[0][1]])
    ax3.hist([np.mean(i) for i in f1_scores_by_folding[1][0]])
    ax4.hist([np.mean(i) for i in f1_scores_by_folding[1][1]])
    ax1.set_title("both unfolded")
    ax2.set_title("ref. unfolded, query folded")
    ax3.set_title("query unfolded, ref folded")
    ax4.set_title("both folded")
    fig.suptitle("Mean F1 score histograms for confusion matrix.")
    fig.savefig("mean_f1_score_conf_matrix.pdf")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    ax1.hist([np.max(i) for i in f1_scores_by_folding[0][0]])
    ax2.hist([np.max(i) for i in f1_scores_by_folding[0][1]])
    ax3.hist([np.max(i) for i in f1_scores_by_folding[1][0]])
    ax4.hist([np.max(i) for i in f1_scores_by_folding[1][1]])
    ax1.set_title("both unfolded")
    ax2.set_title("ref. unfolded, query folded")
    ax3.set_title("query unfolded, ref folded")
    ax4.set_title("both folded")
    fig.suptitle("Max F1 score histograms for confusion matrix.")
    fig.savefig("max_f1_score_conf_matrix.pdf")

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    ax1.hist([np.min(i) for i in f1_scores_by_folding[0][0]])
    ax2.hist([np.min(i) for i in f1_scores_by_folding[0][1]])
    ax3.hist([np.min(i) for i in f1_scores_by_folding[1][0]])
    ax4.hist([np.min(i) for i in f1_scores_by_folding[1][1]])
    ax1.set_title("both unfolded")
    ax2.set_title("ref. unfolded, query folded")
    ax3.set_title("query unfolded, ref folded")
    ax4.set_title("both folded")
    fig.suptitle("Min F1 score histograms for confusion matrix.")
    fig.savefig("min_f1_score_conf_matrix.pdf")


    # Create heat map for number of bp in ref and query structure
    ref_bp_counts_sorted = []
    pred_bp_counts_sorted = []
    for seq in gp_map:
        ref_bp_counts_sorted.append(count_bp(ref_dict[seq]))
        pred_bp_counts_sorted.append(count_bp(gp_map[seq][0]))

    # maximum number of base-pairs in any structure
    max_bp = max(max(ref_bp_counts_sorted, pred_bp_counts_sorted))

    bp_count_pairs = zip(ref_bp_counts_sorted, pred_bp_counts_sorted)
    counts_matrix_ = np.zeros((max_bp,max_bp))
    for count_pair in bp_count_pairs:
            counts_matrix_[count_pair] += 1
    counts_matrix = np.log10(counts_matrix_)
    
    fig, ax = plt.subplots()
    ax.imshow(counts_matrix, cmap='Oranges', interpolation='nearest', origin='lower')
    for i in range(counts_matrix.shape[0]):
        for j in range(counts_matrix.shape[1]):
            if counts_matrix[i, j] < 0:
                ax.text(j, i, "0")
            else:
                ax.text(j, i, "{:.1f}".format(counts_matrix[i, j]), va='center', ha='center')
    ax.set_xlabel("#base-pairs (predicted)")
    ax.set_ylabel("#base-pairs (reference)")
    ax.set_title("Base-pair counts (log10)")
    fig.savefig("basepair_counts_heatmap.pdf")
