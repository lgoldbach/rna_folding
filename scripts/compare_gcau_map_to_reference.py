#!/usr/bin/env python

import argparse
from rna_folding.mapping_functions import nussinov, debug_nussinov_mfe
from rna_folding.base_pairing import BasePairing
from rna_folding.utils import combinatorically_complete_genotypes
import RNA
import numpy as np
import time


seed = np.random.randint(0, 1000)
seed = 23
print(f"random seed: {seed}")
np.random.seed(seed)
L=12
bases = "GACU"  # order matters to match base-pairing!
bases_alt = "JKLM"  # order matters to match base-pairing!
bases_trans = {"J":"G", 'K': 'A', 'L': 'C', 'M': 'U'}
graph_id = 7  # for seq of len 4, base-pairing graph is the canonical one
graph_path = "/home/lgold/phd/research/projects/connectivity/rna_folding/data/graphs/"
min_loop_size=3
base_pairing_alt = BasePairing(bases=bases_alt,
                     graph_path=graph_path,
                     id=graph_id)

sample_size = 100
genotype_tuples = combinatorically_complete_genotypes(l=L, a=bases_alt)
genotypes = ["".join(gt) for gt in genotype_tuples]

# keep track of phenotypes in order have a consistent numbering
ph_ids = {}
phenotypes = []
ph_counter = -1

# store phentype ids and mfe for each sequence in list of lists
gp_map_ph_2 = []
gp_map_mfe_2 = []

subset_counter = 0

canon_in_s2 = 0
mfe_match_s2 = 0
folded_2_counter = 0

start = time.time()

counter = 0
while counter < sample_size:
    seqID = counter
    seq = np.random.choice(genotypes, size=1)[0]
    seq_canon = ''.join([bases_trans[b] for b in seq])
    canon_ph, canon_mfe = RNA.fold(seq_canon)
    if canon_ph == '.' * 12:
        continue
    # canon_bp_count = canon_ph.count("(")
    print(seq, seq_canon)
    print(canon_ph, canon_mfe)
    counter += 1
    gp_map_ph_2.append([])  # start new list for current seq
    gp_map_mfe_2.append([])  # start new list for current seq


    s2, phenos, mfes = debug_nussinov_mfe(genotype=seq, base_pairing=base_pairing_alt,
                    min_loop_size=min_loop_size, suboptimal=2, 
                    structures_max=1000, seed=np.random.randint(0, 10000),
                      deterministic=True)
    print(s2, phenos)
    idx = phenos.index(canon_ph)
    print('MFE of canon in nussinov:', mfes[idx], 'nussinov min mfe:', np.min(mfes))
    print('Size of suboptimal:', len(phenos))

    min_bp=1000
    max_bp=0
    correct_pred = False
    for ph in s2:
        if ph in ph_ids:
            gp_map_ph_2[-1].append(ph_ids[ph])  # append phenotype
        else:
            phenotypes.append(ph)
            ph_counter += 1
            ph_ids[ph] = ph_counter
            gp_map_ph_2[-1].append(ph_counter)
        
        e = RNA.eval_structure_simple(seq, ph)
        gp_map_mfe_2[-1].append(e)
        if canon_ph == ph:
            correct_pred = True
            canon_in_s2 += 1
            print(f"# Num of correct ph in suboptimal 2: {canon_in_s2} out of {counter}\n")
            if np.round(canon_mfe,2) == np.round(e, 2):
                mfe_match_s2 += 1
                print(f"# Num of correct mfe: {mfe_match_s2} out of {counter}\n\n")
                if canon_ph == "."*L:
                    folded_2_counter += 1
        print('\n\n')
    # if correct_pred == False:
    #     if canon_ph == 12*"." and np.min(gp_map_mfe_2[-1]) >=0:
    #         canon_in_s2 +=1
    #     else:
    #         print("X", seq, canon_ph, gp_map_ph_2[-1])
    #     print(seq, min_bp, max_bp, canon_bp_count, canon_ph, s2, np.round(gp_map_mfe_2[-1],2))
    # else:
    #     print("AAAA")

    # intersect = len(list(set(gp_map_ph_0[-1]) & set(gp_map_ph_2[-1])))
    # if intersect == len(gp_map_ph_0[-1]):
    #     subset_counter += 1
    # else:
    #     print(seq, canon_ph)
    #     print(s0)
    #     print(list(set(gp_map_ph_0[-1]) & set(gp_map_ph_2[-1])))
    #     print(intersect, len(gp_map_ph_0[-1]))
    #     print(gp_map_ph_0[-1], gp_map_ph_2[-1])
    #     break

end = time.time()
#Subtract Start Time from The End Time
total_time = end - start
print("\n"+ str(total_time))

with open("phenotypes.txt", "w") as ph_out:
    for ph in phenotypes:
        ph_out.write(ph + "\n")

with open("results.txt", "w") as out:
    out.write(f"In {subset_counter} out of {counter} cases is the nussinov suboptimal 0 result a subset of suboptimal 2\n\n")
    # out.write(f"Correct structure in suboptimal 0: {canon_in_s0}, correct mfe also: {mfe_match_s0}, unfolded: {folded_0_counter}\n")
    out.write(f"Correct structure in suboptimal 2: {canon_in_s2}, correct mfe also: {mfe_match_s2}, unfolded: {folded_2_counter}\n") 

# with open("gp_map_ph_suboptimal0.txt", "w") as out:
#     for line in gp_map_ph_0:
#         out.write(" ".join([str(l) for l in line])+"\n")
with open("gp_map_ph_suboptimal2.txt", "w") as out:
    for line in gp_map_ph_2:
        out.write(" ".join([str(l) for l in line])+"\n")
# with open("gp_map_mfe_suboptimal0.txt", "w") as out:
#     for line in gp_map_mfe_0:
#         out.write(" ".join([str(l) for l in line])+"\n")
with open("gp_map_mfe_suboptimal2.txt", "w") as out:
    for line in gp_map_mfe_2:
        out.write(" ".join([str(l) for l in line])+"\n")

end = time.time()
#Subtract Start Time from The End Time
total_time = end - start
print("\n"+ str(total_time))