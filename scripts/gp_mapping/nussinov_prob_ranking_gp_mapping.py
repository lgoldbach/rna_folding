#!/usr/bin/env python

import argparse
import numpy as np

from rna_folding.base_pairing import BasePairing
from rna_folding.mapping_functions import gp_mapper, nussinov_with_probabilistic_scoring
from rna_folding.parsing import load_phenotype_and_metric_from_file


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file with genotypes")
    parser.add_argument("-o", "--output", help="File output for phenotypes")
    parser.add_argument("-p", "--phenotype_scores", required=True, type=str, help="Allowed phenotypes and their scores. Phenotypes that don't appear are assumed to have a score of 0")
    parser.add_argument("-b", "--base_pairing", required=False, type=int, default=-1, help="Which base-pairing to choose. I.e. from the base-pairing simple graphs, which one to pick "
                        "e.g. for 4 bases there are 11 possible base-pairings, so possible input is any number between 1 and 11, If given -1 then it uses canonical base-pairing and AUGC bases")
    parser.add_argument("-a", "--alphabet", required=False, type=str, default="AUGC", help="Which bases do the genotypes contain, e.g. 'AUGC' for canonical RNA")
    parser.add_argument("-g", "--graph_path", required=True, type=str,
                        help="Path to folder containing the base-pairing "
                        "graphs files, e.g. graph4.adj. Check base_pairing.py "
                        "for info on where these graphs come from.")

    args = parser.parse_args()

    rng = np.random.default_rng(858292)
    
    pairing = BasePairing(bases=args.alphabet,
                          graph_path=args.graph_path,
                          id=args.base_pairing)
    
    phenotypes, scores = load_phenotype_and_metric_from_file(args.phenotype_scores)
    ph_scores = dict(zip(phenotypes, scores))

    mapping = lambda seq: nussinov_with_probabilistic_scoring(seq,
                                                              scores=ph_scores,
                                                              base_pairing=pairing, 
                                                              rng=rng)

    # generate g-p map and save to output file
    gp_mapper(input=args.input, output=args.output, 
              mapping_function=mapping)
