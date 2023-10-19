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

import time
import numpy as np
import argparse
import RNA

from rna_folding.utils import dotbracket_to_genotype
from rna_folding.parsing import gpmap_to_dict


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genotypes", required=True, type=str, 
                        help="List of genotypes")
    parser.add_argument("-m", "--gp_map", required=True, type=str, 
                        help="List of phenotypes and genotype IDs that map to "
                        "them")
    parser.add_argument("-r", "--phenotype_ref", required=True, type=str, 
                        help="Reference that defines the phenotype numbering, "
                        "Simple file where each row contains one phenotype")
    parser.add_argument("-b", "--basepair", type=str,
                        help="Which base-pair to convert to, max. 1 pair "
                        "(2 letters), e.g. 'GC'", default="GC")
    parser.add_argument("-d", "--deterministic",
                        action='store_true',
                        help="Do not assign bases randomly.")
    parser.add_argument("-o", "--output", required=True, type=str, 
                        help="Output file")
    
    args = parser.parse_args()
    gp_map = gpmap_to_dict(args.gp_map, args.genotypes)
    L = len(list(gp_map.keys())[0])

    # read in phenotype ref. Simple list of phenotypes
    with open(args.phenotype_ref, "r") as file:
        phenos = [line.strip() for line in file]
        ph_ref = dict(zip(phenos, range(1, len(phenos)+1)))

    g_fe_map = {}
    g_mfe_map = {}

    start_time = time.time()
    for seq in gp_map:
        g_fe_map[seq] = []
        for struc in gp_map[seq]:
            # turn into a canonical alphabet
            seq_canon = dotbracket_to_genotype(dotbracket=struc,
                                               base_pair=args.basepair,
                                               random=True,
                                               seed=1996)
            g_fe_map[seq].append(RNA.eval_structure_simple(seq_canon, struc))
        
        mfe_struc_id = np.argmin(g_fe_map[seq])  # get index of mfe structure
        mfe_struc = gp_map[seq][mfe_struc_id]  # get mfestructure
        
        g_mfe_map[seq] = ph_ref[mfe_struc]  # get ph id of structure
    
    with open(args.output, "w") as file:
        for seq in g_mfe_map:
            file.write(seq + " " + str(g_mfe_map[seq]) + "\n")

