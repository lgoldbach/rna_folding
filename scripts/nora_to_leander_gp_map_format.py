#!/usr/bin/env python

import pickle
import numpy as np

from rna_folding.gp_map import GenotypePhenotypeGraph
from rna_folding.utils import combinatorically_complete_genotypes
from rna_folding.parsing import dict_to_gpmap


gp = pickle.load(open("gp_map_Nora.pkl", "rb"))
alphabet = "10"

print("Load done", flush=True)

gt_gen = combinatorically_complete_genotypes(20, alphabet)

phenotypes = []
genotypes = []
for g in gt_gen:
    gt_str = "".join(g)
    genotypes.append(gt_str)
    
    g_tup = tuple(int(s) for s in g) # turn str genotype into tuple
    phenotypes.append(str(gp[g_tup]))

print("gp map gen start", flush=True)
gp_map = GenotypePhenotypeGraph(genotypes=genotypes, phenotypes=phenotypes, alphabet=alphabet)
print("gp map gen Done", flush=True)

pickle.dump(gp_map, open("gp_map.pickle", "bw"))
print("gp map dump Done", flush=True)

with open("genotypes.txt", "w") as f:
    for gt in genotypes:
         f.write(gt + "\n")

ph_uniq = np.unique(phenotypes)
with open("phenotypes.txt", "w") as f:
    for ph in ph_uniq:
         f.write(ph + "\n")

gt_id = range(len(genotypes))
pg_map = {}
for (gt, ph) in zip(gt_id, phenotypes):
    if ph in pg_map:
        pg_map[ph].append(gt)
    else:
        pg_map[ph] = [gt]

with open("gp_map.txt", "w") as file_out:
    for p in pg_map:
        line = p + " " + " ".join(map(str, pg_map[p])) + "\n"
        file_out.write(line)

print("gp map txt Done", flush=True)

print(genotypes[:10], flush=True)
print(phenotypes[:10], flush=True)
print(ph_uniq, type(ph_uniq[5]), flush=True)
print(gp_map.genotypes[:10], flush=True)
print(gp_map.phenotypes[:10], flush=True)
print("type", type(gp_map.phenotypes[10]), flush=True)
print(gp_map.map(genotypes[10000]), flush=True)
print(gp_map.phenotype_set, flush=True)