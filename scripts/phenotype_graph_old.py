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
        phenotypes = [ph for ph in gpm.phenotypes if ph != args.ignore]
    else:
        phenotypes = None

    
    print("Compute neutral components", flush=True)
    a = datetime.datetime.now()

    # get list of generators of set, one generator for each phenotype
    nn_per_ph_gen = gpm.neutral_components(phenotypes=phenotypes)
    b = datetime.datetime.now()
    c = b-a
    print(f"Done in ${np.round(c.seconds/60, 2), c.seconds}", flush=True)
    
    ph_graph = nx.Graph()

        
    print("Start labeling nodes with their neutral components", flush=True)
    a = datetime.datetime.now()
    # label nodes with their neutral component:
    nc_counter = 0
    neutral_components = []  # fill a list of lists of genotypes from the iterator
    for nn in nn_per_ph_gen:
        for nc in nn:
            nc_counter += 1
            neutral_components.append([gt for gt in nc])
            nc_dict = {gt: nc_counter for gt in neutral_components[-1]}  # map each genotype to its nc num.
            nx.set_node_attributes(gpm, nc_dict, "neutral_component")

            # add a node to the phenotype fraph for each neutral component (nc)
            ph = gpm.nodes[neutral_components[-1][0]]["phenotype"]
            ph_graph.add_node(nc_counter, phenotype=ph)  # add node for neutral component
    print(nc_counter)
    b = datetime.datetime.now()
    c = b-a
    print(f"Done in ${np.round(c.seconds/60, 2), c.seconds}", flush=True)
    print("Number of neutral components", len(neutral_components))

        
    print("Start building neutral component graph", flush=True)
    a = datetime.datetime.now()
    for i, nc in enumerate(neutral_components, start=1):  # i is nc id
        for genotype in nc:
            for neighbor in gpm._neighbors(genotype):
                try:
                    j = gpm.nodes[neighbor]["neutral_component"]  # get nc id
                except KeyError:
                    # node without "neutral_component" attribute is part of 
                    # "--ignored" phenotype
                    continue  

                if j != i:  # different nc
                    ph_graph.add_edge(i, j)  # add to phenotype graph
    b = datetime.datetime.now()
    c = b-a
    print(f"Done in ${np.round(c.seconds/60, 2), c.seconds}", flush=True)
    
    print("Number of nodes in ph. graph:", len(ph_graph.nodes))
    
        
    print("Pickle dump ph_graph", flush=True)
    a = datetime.datetime.now()
    pickle.dump(ph_graph, open(args.output, "wb"))
    b = datetime.datetime.now()
    c = b-a
    print(f"Done in ${np.round(c.seconds/60, 2), c.seconds}", flush=True)

    labels = nx.get_node_attributes(ph_graph, 'phenotype')
    nx.draw(ph_graph, labels=labels)
    plt.savefig("phenotype_graph.pdf", dpi=30)
    
    plt.clf()

    labels = nx.get_node_attributes(gpm, 'phenotype')
    nx.draw(gpm, labels=labels)
    plt.savefig("gpm.pdf", dpi=30)
    