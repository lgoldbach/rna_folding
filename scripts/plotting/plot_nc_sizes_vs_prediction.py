#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import pickle
from scipy.stats.stats import pearsonr

from rna_folding.utils import count_bp

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--nc_graphs", help="Neutral component graph pickle "
                        "objects", required=True, nargs="+")    
    parser.add_argument("-o", "--output", help="pdf file", required=True, type=str)
    
    args = parser.parse_args()

mut_graph_n = {2: 2,
               3: 2,
               4: 2,
               5: 4,
               6: 1,
               7: 2,
               8: 1,
               9: 2,
               10: 1,
               11: 1}

mut_graph_s = {2: 1,
               3: 2,
               4: 3,
               5: 1,
               6: 6,
               7: 3,
               8: 8,
               9: 4,
               10: 10,
               11: 12}

max_bp = {2: 2,
            3: 2,
            4: 2,
            5: 2,
            6: 2,
            7: 2,
            8: 3,
            9: 3,
            10: 4,
            11: 4}

# max_bp = {2: 4,
#             3: 4,
#             4: 4,
#             5: 4,
#             6: 4,
#             7: 4,
#             8: 4,
#             9: 4,
#             10: 4,
#             11: 4}

num_red_mut = {2: 2,
                3: 2,
                4: 6,
                5: 0,
                6: 0,
                7: 0,
                8: 0,
                9: 4,
                10: 2,
                11: 0}


fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

ax1.set_box_aspect(1)
ax2.set_box_aspect(1)

x = []
y = []

x2 = []
y2 = []

# predict 
new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
new_order_idx = [i - 2 for i in new_order]
for i, l in enumerate(new_order, start=1):
    nc_graph_f = args.nc_graphs[l-2]
    nc_graph = pickle.load(open(nc_graph_f, "rb"))
    ph_all = list(nx.get_node_attributes(nc_graph, name="phenotype").values())
    phenotypes = np.unique(ph_all)
    nc_n_est_all = []
    for bp_n in range(5):  # assume all phenotypes have either 2, 3 or 4 bp
        nc_n_est = 0
        for j in range(len(phenotypes)):  # loop over number of phenotypes
            bp = np.random.choice([2, 3, 4])
            nc_n_est += mut_graph_n[l]**bp
        
        nc_n_est_all.append(nc_n_est)
    nc_n_est_avg = np.mean(nc_n_est_all)
    nc_n_truth = len(nc_graph.nodes)

    folded_gt = sum(list(nx.get_node_attributes(nc_graph, name="size").values()))

    gt_total = 4**12

    nc_ests = []
    for nc_n_est in nc_n_est_all:
        nc_est = []
        frac = folded_gt/gt_total
        # nc_size = 1000000000
        for j in range(nc_n_est):
            stack = mut_graph_s[l] ** max_bp[l]  # 4: number of pairs
            loop = 4**(12 - (max_bp[l]*2)) # 6: loop size, 4: alphabet size
            nc_size = stack * loop * frac
            if nc_size < 100 or j > nc_n_est:
                break
            nc_est.append(nc_size)
            frac = frac - (nc_est[-1] / gt_total)

        nc_ests.append(nc_est)
    nc_n_std = np.std([len(nc_est) for nc_est in nc_ests])
    nc_n_mean = np.mean([len(nc_est) for nc_est in nc_ests])
    if i == 4:
        label = "canon."
    else:
        label = f"{i}"
    sc = ax2.scatter(nc_n_truth, nc_n_mean, label=label)
    col = sc.get_facecolors()[0].tolist()  # get color
    # ax.scatter(nc_n_truth, nc_n_est_all[1], color=col)
    # ax.scatter(nc_n_truth, nc_n_est_all[2], color=col)
    x.append(nc_n_truth)
    y.append(nc_n_mean)
    ax2.set_xlim(0, 1050)
    ax2.set_ylim(0, 1050)
    ax2.plot([0, 1050], [0, 1050], zorder=-5, linestyle="--", color="black", linewidth=0.5)

    ### Estimate neutral component sizes
    nc_est = []
    frac = folded_gt/gt_total
    for j in range(nc_n_est):
        stack = mut_graph_s[l] ** max_bp[l]  # 4: number of pairs
        loop = 4**(12 - (max_bp[l]*2)) # 6: loop size, 4: alphabet size
        nc_size = stack * loop * frac
        if nc_size < 100 or j > nc_n_est:
            break
        nc_est.append(nc_size)
        frac = frac - (nc_est[-1] / gt_total)

    print(i)
    for p in range(2, 5):
        frac = folded_gt/gt_total
        stack = mut_graph_s[l] ** p  # 4: number of pairs
        loop = 4**(12 - (p*2)) # 6: loop size, 4: alphabet size
        nc_size = stack * loop * frac
        print(p, nc_size)

    largest_nc_truth = max(list(nx.get_node_attributes(nc_graph, name="size").values()))
    print(largest_nc_truth)
    x2.append(np.log10(largest_nc_truth))
    y2.append(np.log10(nc_est[0]))
    ax1.scatter(np.log10(largest_nc_truth), np.log10(nc_est[0]), label=label, color=col)
    ax1.set_xlim(3, 6)
    ax1.set_ylim(3, 6)
    ax1.plot([3, 6], [3, 6], zorder=-5, linestyle="--", color="black", linewidth=0.5)
    # if i != 3:
    #     ax1.scatter(np.log10(largest_nc_truth), nc_n_truth, label=label)

fontsize = 14
r, p = pearsonr(x, y)
p_str = "%.3g" % p
ax2.text(30, 950, f'r = {np.round(r, 2)}\np = {p_str}')
ax2.set_xlabel("No. neutral networks", fontsize=fontsize)
ax2.set_ylabel("Predicted no. neutral networks", fontsize=fontsize)
# ax2.legend(loc="lower right", prop={'size': 8}, frameon=False, title="Base-pairing")

ax1.tick_params(axis='both', which='major', labelsize=fontsize)
ax1.tick_params(axis='both', which='minor', labelsize=fontsize)
ax2.tick_params(axis='both', which='major', labelsize=fontsize)
ax2.tick_params(axis='both', which='minor', labelsize=fontsize)

r, p = pearsonr(x2, y2)
p_str = "%.3g" % p
ax1.text(3.08, 5.7, f'r = {np.round(r, 2)}\np = {p_str}')
ax1.set_xlabel("Size of largest\nneutral network (log10)", fontsize=fontsize)
ax1.set_ylabel("Predicted size of\nlargest neutral network (log10)", fontsize=fontsize)
ax1.legend(loc="lower right", prop={'size': 10}, frameon=False, title="GP map", title_fontsize=fontsize)

plt.tight_layout()
plt.savefig(args.output, format="pdf", dpi=30)
