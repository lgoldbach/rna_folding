#!/usr/bin/env python

"""Compute the MFE for a gp_map with alternative alphabets. This includes 
the following steps:
For each sequence:
    For each phenotype:
1) Turn the genotype into an AU or GC genotype. This is done by assigning
the positions of all opening and closing brackets a GC or AU base-pair and the
unpaired position get either a random or the same base assigned.
Then the MFE is computed for that sequence-structure pair.
2) Take the minimum MFE

e.g.:

input: 
LLKKMMMOO (((...))) ()().... 
LLLKKKMMM (.......) ()()()()
...

step 1 (for line 1):
ph        gt         mfe 
(((...))) GGGGGGCCC  -10
()()..... GCGCGGGGG  -5

output:
LLKKMMMOO (((...)))

"""
import argparse

from rna_folding.base_pairing import BasePairing
from rna_folding.mapping_functions import gp_mapper, nussinov_canonical_fe


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input file with genotypes")
    parser.add_argument("-o", "--output", help="File output for phenotypes")
    parser.add_argument("-m", "--min_loop_size", required=True, type=int, default=1, help="Minimum size for loop")
    parser.add_argument("-s", "--suboptimal", type=int, required=True,
                        help="Create all suboptimal structures with number of base-pairs in the range"
                             "of max - s, where s is an integer. Without the flag, only one structure "
                             "is computed.")
    parser.add_argument("-z", "--structures_max", required=False, type=int, help="Limit on how many suboptimal structures to generate")
    parser.add_argument("-a", "--alphabet", required=False, type=str, default="AUGC", help="Which bases do the genotypes contain, e.g. 'AUGC' for canonical RNA")

    args = parser.parse_args()
    
    # use canonical base-pairing. we use -1 to indicate that
    pairing = BasePairing(bases=args.alphabet,
                          graph_path=None, 
                          id=-1)
    
    mapping = lambda seq: nussinov_canonical_fe(seq, 
                                   base_pairing=pairing, 
                                   min_loop_size=args.min_loop_size, 
                                   suboptimal=args.suboptimal, 
                                   structures_max=args.structures_max)

    # generate g-p map and save to output file
    gp_mapper(input=args.input, output=args.output, 
              mapping_function=mapping)

