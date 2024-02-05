#!/usr/bin/env python

import argparse
import pickle

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map, "
                        "i.e. file containing a pickle GenoypePenotypeGraph "
                        "object", required=True)
    parser.add_argument("-o", "--output", help="File output for robustness",
                        required=True)
    

    args = parser.parse_args()
    
    gpm = pickle.load(open(args.file, "rb"))

    ncs = []
    for ph in gpm.phenotype_set:
        ncs.append(gpm.nodes_with_phenotype(ph))

    robustnesses = []
    for i, nc in enumerate(ncs):
        robustnesses.append(gpm.phenotype_robustness(nc))

    with open(args.output, "w") as outfile:
        for ph, r in zip(gpm.phenotype_set, robustnesses):
            outfile.write(ph + " " + str(r) + "\n")
    