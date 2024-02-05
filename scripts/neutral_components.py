#!/usr/bin/env python

import argparse
import pickle
import networkx as nx


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-o", "--output", help="File output for robustness",
                        required=True)
    

    args = parser.parse_args()
    
    gpm = pickle.load(open(args.file, "rb"))

    for ph in gpm.phenotype_set:
        neutral_set = gpm.nodes_with_phenotype(ph)
        G_neutral_set = nx.induced_subgraph(gpm, neutral_set)
        print(list(nx.connected_components(G_neutral_set)))

    # with open(args.output, "w") as outfile:
    #     for ph, r in zip(gpm.phenotype_set, robustnesses):
    #         outfile.write(ph + " " + str(r) + "\n")
    