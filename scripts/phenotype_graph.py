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
    
    print("Start loading gpm", flush=True)
    a = datetime.datetime.now()
    gpm = pickle.load(open(args.file, "rb"))
    b = datetime.datetime.now()
    c = b-a
    print(f"Done in ${np.round(c.seconds/60, 2), c.seconds}", flush=True)

    if args.ignore:
        phenotypes = [ph for ph in gpm.phenotype_set if ph != args.ignore]
    else:
        phenotypes = None
    
    print("Compute neutral components", flush=True)
    a = datetime.datetime.now()

    # get list of generators of set, one generator for each phenotype
    nc_sizes, boundaries = gpm.neutral_component_sizes(phenotypes=phenotypes,
                                                return_boundaries=True,
                                                add_labels=True)

    # pickle.dump(gpm, open(args.file, "wb"))

    b = datetime.datetime.now()
    c = b-a
    print(f"Done in ${np.round(c.seconds/60, 2), c.seconds}", flush=True)
    
    ph_graph = nx.Graph()
        
    a = datetime.datetime.now()

    for (i, j) in boundaries:
        nc_i = gpm.nodes[i]["neutral_component"]  # get nc id
        nc_j = gpm.nodes[j]["neutral_component"]  # get nc id
       
        phenotype_dict = {nc_i: gpm.nodes[i]["phenotype"], nc_j: gpm.nodes[j]["phenotype"]}
        ph_graph.add_nodes_from([nc_i, nc_j])
        nx.set_node_attributes(ph_graph, phenotype_dict, "phenotype")
        nc_dict = {nc_i: gpm.nodes[i]["neutral_component"], nc_j: gpm.nodes[j]["neutral_component"]}
        nx.set_node_attributes(ph_graph, nc_dict, "neutral_component")
        ph_graph.add_edge(nc_i, nc_j)  # add to phenotype graph


    print(ph_graph.nodes, len(ph_graph.edges), nx.get_node_attributes(ph_graph, "phenotype"))
    pickle.dump(ph_graph, open(args.output, "wb"))
    
    # arr = nx.to_numpy_array(ph_graph)
    # attr = nx.get_node_attributes(ph_graph, "phenotype")
    # pickle.dump((arr, attr), open(args.output, "wb"))
    

    # labels = nx.get_node_attributes(ph_graph, 'neutral_component')
    # nx.draw(ph_graph, labels=labels)
    # plt.savefig("phenotype_graph.pdf", dpi=30)
    
    # # plt.clf()

    # # gpm.add_hamming_edges()
    # # labels = nx.get_node_attributes(gpm, 'neutral_component')
    # # nx.draw(gpm, labels=labels)
    # # plt.savefig("neutral_components.pdf", dpi=30)
    
    # # plt.clf()

    # # labels = nx.get_node_attributes(gpm, "phenotype")
    # # nx.draw(gpm, labels=labels)
    # # plt.savefig("gpm_ph.pdf", dpi=30)
    