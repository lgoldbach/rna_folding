#!/usr/bin/env python

import pickle
import numpy as np

from rna_folding.gp_map import GenotypePhenotypeGraph
from rna_folding.utils import combinatorically_complete_genotypes


gp = pickle.load(open("gp_map_Nora.pkl", "rb"))
alphabet = "10"

print("Load done", flush=True)
g_t = (0,0,0,0,0,1,1,1,0,1,0,1,1,1,0,1,0,1,1,1,1,0,1,1,1)

gt_gen = combinatorically_complete_genotypes(25, alphabet)

phenotypes = []
genotypes = []
for g in gt_gen:
    gt_str = "".join(g)
    genotypes.append(gt_str)
    
    g_tup = tuple(int(s) for s in g) # turn str genotype into tuple
    phenotypes.append(gp[g_tup])

print("gp map gen start", flush=True)
gp_map = GenotypePhenotypeGraph(genotypes=genotypes, phenotypes=phenotypes, alphabet=alphabet)
print("gp map gen Done", flush=True)

pickle.dump(gp_map, open("gp_map.pickle", "bw"))
print("gp map dump Done", flush=True)

print(genotypes[:10], flush=True)
print(phenotypes[:10], flush=True)
print(gp_map.genotypes[:10], flush=True)
print(gp_map.phenotypes[:10], flush=True)
print(gp_map.map(genotypes[10000]), flush=True)
print(gp_map.phenotype_set, flush=True)