#!/usr/bin/env python

import argparse
from itertools import chain
from datetime import datetime

from rna_folding.gp_map import GenotypePhenotypeGraph

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-g", "--genotypes", help="Reference file that "
                        "contains the ith genotype on the ith line. Used to "
                        "map the numeric genotypes back to the sequences ", 
                        required=True)
    parser.add_argument("-a", "--alphabet", help="Define the alphabet used in "
                        "the genotype map, e.g. -a AUGC", required=True)
    parser.add_argument("-o", "--output", help="File output for robustness",
                        required=True)
    

    args = parser.parse_args()
    
    gpm = GenotypePhenotypeGraph.read_from_ph_to_gt_file(path=args.file, 
                                              genotype_ref_path=args.genotypes,
                                              alphabet=args.alphabet)

    ncs = []
    for ph in gpm.phenotype_set:
        ncs.append(gpm.nodes_with_phenotype(ph))

    robustnesses = []
    for i, nc in enumerate(ncs):
        robustnesses.append(gpm.phenotype_robustness(nc))

    with open(args.output, "w") as outfile:
        for ph, r in zip(gpm.phenotype_set, robustnesses):
            outfile.write(ph + " " + str(r) + "\n")
    