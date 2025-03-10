#!/usr/bin/env python

import argparse
import pickle
import networkx as nx

import matplotlib.pyplot as plt
import datetime
import time
import numpy as np


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-i", "--ignore", help="Phenotype to ignore, e.g "
                        "unfolded", type=str, required=False)
    parser.add_argument("-o", "--output", help="File output for phenotype "
                        "Should end in .pickle",
                        required=True)
    

    args = parser.parse_args()
    
    # print("Start loading gpm", flush=True)
    a = datetime.datetime.now()
    gpm = pickle.load(open(args.file, "rb"))
    b = datetime.datetime.now()
    c = b-a
    # print(f"Done in ${np.round(c.seconds/60, 2), c.seconds}", flush=True)

    if args.ignore:
        phenotypes = [ph for ph in gpm.phenotype_set if ph != args.ignore]
    else:
        phenotypes = None
    
    # print("Compute neutral components", flush=True)
    a = datetime.datetime.now()

    ncs, boundaries = gpm.get_neutral_components(phenotypes=phenotypes,
                                                return_boundaries=True,
                                                add_labels=True)

    # b = datetime.datetime.now()
    # c = b-a
    # print(f"Done in ${np.round(c.seconds/60, 2), c.seconds}", flush=True)

    nc_graph = nx.Graph()
    # add a node for each neutral component and its phenotype value
    for ph in ncs:
        for nc in ncs[ph]:
            nc_graph.add_node(nc, phenotype=ph, size=len(ncs[ph][nc]))
    
    # a = datetime.datetime.now()

    for (i, j) in boundaries:
        nc_i = gpm.nodes[i]["neutral_component"]  # get nc id
        nc_j = gpm.nodes[j]["neutral_component"]  # get nc id

        # edge weights are added in .5 increments because all nc boundaries are
        # counted twice, once from each direction. We do not care about 
        # direction and thus divide it by 2. This can lead to problems in 
        # incomplete gp maps but that's a story for another day
        if not nc_graph.has_edge(nc_i, nc_j):
            nc_graph.add_edge(nc_i, nc_j, weight=.5)  # add edge to nc graph
        else:
            # increase weight if the edge already exists. In undirected graph 
            # the order of nodes does not matter when referencing the edge,
            # (ni, nj) or (nj, ni) are the same edge
            nc_graph.edges[(nc_i, nc_j)]["weight"] += .5 
        



    pickle.dump(nc_graph, open(args.output, "wb"))
    
    # arr = nx.to_numpy_array(ph_graph)
    # attr = nx.get_node_attributes(ph_graph, "phenotype")
    # pickle.dump((arr, attr), open(args.output, "wb"))

    labels = {node: str(node) + " " + nc_graph.nodes[node]["phenotype"] for node in nc_graph}
    pos = nx.spring_layout(nc_graph)
    
    weights = [nc_graph[u][v]['weight']*3 for u,v in nc_graph.edges()]

    color_by_ph = {ph: np.random.choice(range(256), size=3)/256 for ph in ncs}
    colors = [color_by_ph[nc_graph.nodes[node]["phenotype"]] for node in nc_graph.nodes]

    nx.draw(nc_graph, pos=pos, labels=labels, node_color=colors, node_size=[s*1000 for s in nx.get_node_attributes(nc_graph, "size").values()], width=weights)

    edge_labels = {(u, v): int(nc_graph[u][v]['weight']) for u,v in nc_graph.edges()}
    nx.draw_networkx_edge_labels(nc_graph, pos=pos, edge_labels=edge_labels)
    plt.savefig("phenotype_graph.pdf", dpi=30)
    
    plt.clf()

    gpm.add_hamming_edges()
    labels = nx.get_node_attributes(gpm, 'neutral_component')
    nx.draw(gpm, labels=labels)
    plt.savefig("neutral_components.pdf", dpi=30)
    
    plt.clf()

    labels = {node: str(gpm.nodes[node]["neutral_component"]) + " " + gpm.nodes[node]["phenotype"] for node in gpm}
    nx.draw(gpm, labels=labels)
    plt.savefig("gpm_ph.pdf", dpi=30)
    