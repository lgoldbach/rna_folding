#!/usr/bin/env python

import argparse
import numpy as np
from scipy.stats import f_oneway, tukey_hsd
import matplotlib.pyplot as plt

from rna_folding.gp_map import GenotypePhenotypeGraph

import networkx as nx

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file ", required=True)
    parser.add_argument("-n", "--sample_size", help="Number of rankings "
                        "considered per bp rule ", 
                        required=True, type=int)
    parser.add_argument("-o", "--output", help="Output file with plot (pdf)",
                        required=True)
    

    args = parser.parse_args()
    
    data_full = np.loadtxt(args.input)
    data_split = np.split(data_full, data_full.size/args.sample_size)
    data_resh = np.reshape(data_full, (data_full.size//args.sample_size, args.sample_size)).T
    means = np.mean(data_resh, axis=0)
    means_argsort = np.argsort(means)


    fig, ax = plt.subplots()
    for x, i in enumerate(means_argsort):
        c = plt.cm.tab10(x)
        ax.scatter(np.repeat(i, args.sample_size), data_resh[:, i], color=c, s=10, label=f"{i+2}")
        ax.scatter(i, means[i], marker="_", color="black", zorder=10, s=120)
    ax.set_ylim(0, 1)

    handles, labels = ax.get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: int(t[0])))
    legend = ax.legend(handles, labels, title="Base-pair rules", fontsize="small", ncol=2, loc=4)
    legend._legend_box.align = "left"
    plt.setp(legend.get_title(),fontsize='small')
    ax.set_xticks([])
    ax.set_ylabel("Genotype robustness")
    ax.set_title("Average genotype robustness and scatter of samples")
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=15)

    F, p = f_oneway(*data_split)
    tukey_res = tukey_hsd(*data_split)
    x, y = np.where(tukey_res.pvalue>.05)
    i_j = list(zip(x,y))
    x2, y2 = np.where(tukey_res.pvalue<=.05)
    i_j2 = list(zip(x2,y2))

    A = tukey_res.pvalue.copy()
    # C = A.copy()
    # C[tukey_res.pvalue < 0.05] = 0
    # C[tukey_res.pvalue >= 0.05] = 1
    np.fill_diagonal(A, 0)

    fig2, ax2 = plt.subplots()
    G = nx.from_numpy_array(A)
    G.edges(data=True)
    edges = G.edges()
    weights = [G[u][v]['weight']*3 for u,v in edges]
    colors = []
    for u, v in edges:
        if A[u,v] > 0.5:
            colors.append("red")
        else:
            colors.append("black")

    pos = nx.spring_layout(G)
    nx.draw_networkx_nodes(G, pos=pos)
    nx.draw_networkx_labels(G, pos=pos, labels=dict(zip(G.nodes, range(2, 12))))
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=colors, width=weights, ax=ax2)
    fig2.savefig("analysis/test.pdf", format="pdf")