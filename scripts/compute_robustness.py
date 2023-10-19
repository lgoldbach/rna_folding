#!/usr/bin/env python

import argparse
from itertools import chain
from datetime import datetime

from rna_folding.gp_map import GenotypePhenotypeGraph

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-a", "--alphabet", help="Define the alphabet used in "
                        "the genotype map, e.g. -a AUGC", required=True)
    parser.add_argument("-o", "--output", help="File output for robustness",
                        required=True)
    

    args = parser.parse_args()
    
    # currentDateAndTime = datetime.now()
    # print(f"Build Genotype-Phenotype Map ({currentDateAndTime.strftime('%H:%M:%S')})")
    gpm = GenotypePhenotypeGraph.read_from_file(path=args.file, 
                                              alphabet=args.alphabet)

    # currentDateAndTime = datetime.now()
    # print(f"Compute neutral components ({currentDateAndTime.strftime('%H:%M:%S')})")
    ncs = []
    for ph in gpm.phenotype_set:
        ncs.append(gpm.nodes_with_phenotype(ph))

    robustnesses = []
    for i, nc in enumerate(ncs):
        # currentDateAndTime = datetime.now()
        # print(f"Compute robustness of NC {i} ({currentDateAndTime.strftime('%H:%M:%S')})")
        robustnesses.append(gpm.phenotype_robustness(nc))

    with open(args.output, "w") as outfile:
        for ph, r in zip(gpm.phenotype_set, robustnesses):
            outfile.write(ph + " " + str(r) + "\n")
    