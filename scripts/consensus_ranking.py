import numpy as np
import random
import pickle
import time
import matplotlib.pyplot as plt
from rna_folding.parsing import gpmap_to_dict
from rna_folding.utils import dotbracket_to_bp

### Initial data loading process
# gp_map_nussi7 = gpmap_to_dict(gpmap_file="gp_map_nussi7.txt", genotype_file="genotypes_nussi7.txt")
# # gp_map_nussi7_1000 = dict(random.sample(gp_map_nussi7.items(), 1000))


# trans_dict = {"L": "A", "J": "U", "M": "G", "K": "C"}
# gp_map_nussi7_augc = {}
# for g in gp_map_nussi7:
#     new_g = "".join([trans_dict[i] for i in g])
#     gp_map_nussi7_augc[new_g] = gp_map_nussi7[g]

# vienna_gp_map = gpmap_to_dict("gp_map_vienna.txt", "genotypes_vienna.txt")

# pickle.dump(vienna_gp_map, open("vienna_gp_map_dict.pickle", "wb"))
# pickle.dump(gp_map_nussi7_augc, open("gp_map_nussi7_augc_dict.pickle", "wb"))




# D = {}
# for gt in gp_map_nussi7_augc:
#     for ph in gp_map_nussi7_augc[gt]:
#         if ph in D:
#             D[ph].append(gt)
#         else:
#             D[ph] = [gt]

# pickle.dump(D, open("pg_map_nussi7_augc.pickle", "wb"))


# unfolded = "."*12
# gp_map_nussi7_augc_folded = {}
# for gt in gp_map_nussi7_augc:
#     if vienna_gp_map[gt] != unfolded:
#         gp_map_nussi7_augc_folded[gt] = gp_map_nussi7_augc[gt]

# pickle.dump(gp_map_nussi7_augc_folded, open("gp_map_nussi7_augc_folded_dict.pickle", "wb"))


# pg_map_nussi7_folded_augc = {}
# for gt in gp_map_nussi7_folded_augc:
#     for ph in gp_map_nussi7_folded_augc[gt]:
#         if ph in pg_map_nussi7_folded_augc:
#             pg_map_nussi7_folded_augc[ph].append(gt)
#         else:
#             pg_map_nussi7_folded_augc[ph] = [gt]

# pickle.dump(pg_map_nussi7_folded_augc, open("pg_map_nussi7_augc_folded_dict.pickle", "wb"))


# phenos = set([i[0] for i in vienna_gp_map.values()])
# print(phenos, len(phenos))
# for i in list(pg_map_nussi7_folded_augc.keys()):
#     if i not in phenos:
#         pg_map_nussi7_folded_augc.pop(i)

start_time = time.time()

vienna_gp_map = pickle.load(open("vienna_gp_map_dict.pickle", "rb"))
gp_map_nussi7_folded_augc = pickle.load(open("gp_map_nussi7_augc_folded_dict.pickle", "rb"))
pg_map_nussi7_folded_augc = pickle.load(open("pg_map_nussi7_augc_folded_dict.pickle", "rb"))

print("--- %s seconds ---" % (time.time() - start_time))


phenos_v = list(set([i[0] for i in vienna_gp_map.values()]))
phenos_v.remove("."*12)
# phenos_n = list(pg_map_nussi7_folded_augc.keys())
# inters = set(phenos_v).intersection(set(phenos_n))

def pairwise_consensus_matrix(phenotypes, pg_map, ref_gp_map):
    A = np.zeros(shape=(len(phenotypes), len(phenotypes)))
    # Loop over upper triangle of matrix
    for i in range(len(phenotypes)):
        for j in range(i+1,len(phenotypes)):
            gt_intersect = set(pg_map[phenotypes[i]]).intersection(pg_map[phenotypes[j]])  # get genotypes that map to both phenos
            for gt in gt_intersect:
                if ref_gp_map[gt][0] == phenotypes[i]:
                    A[i, j] += 1
                elif ref_gp_map[gt][0] == phenotypes[j]:
                    A[j, i] += 1
    return A
             

A = pairwise_consensus_matrix(phenotypes=phenos_v, 
                          pg_map=pg_map_nussi7_folded_augc,
                          ref_gp_map=vienna_gp_map)

# A = pickle.load(open("pairwise_rank_matrix.pickle", "rb"))
pickle.dump(A, open("pairwise_rank_matrix.pickle", "wb"))

balanced = 0
unbalanced = 0
unmatched = 0
unbalanced_list = []
unbalanced_list_count = []
for i in range(len(phenos_v)):
    for j in range(i+1,len(phenos_v)):
        if A[i,j] == 0 and A[j,i] == 0:
            unmatched += 1
        elif A[i,j] == 0 or A[j,i] == 0:
            balanced += 1
        else:
            v = sorted([A[i,j],A[j, i]])
            unbalanced_list.append(1-(v[0]/v[1]))
            unbalanced_list_count.append(sum(v))
            # unbalanced_list.append(abs()/(A[i,j] + A[j, i]))
            unbalanced += 1

print(A.sum(), (len(phenos_v)*(len(phenos_v)-1))/2)
# print(balanced, unbalanced, unmatched)
# print(unbalanced_list)
# plt.hist(unbalanced_list)
# plt.savefig("unbalanced_phenotypes_fraction_diff.pdf", format="pdf", dpi=15)
plt.scatter(unbalanced_list_count, unbalanced_list)
plt.xlabel("Number of match-ups")
plt.ylabel("Consistency of match-ups")
plt.savefig("unbalanced_phenotypes_fraction_diff_over_count.pdf", format="pdf", dpi=15)

for i in range(len(phenos_v)):
    for j in range(i+1,len(phenos_v)):
        if A[i,j] > A[j,i]:
            A[i, j] = 1
            A[j, i] = 0
        elif A[i,j] < A[j,i]:
            A[j, i] = 1
            A[i, j] = 0
        else:
            A[i, j] = 0
            A[j, i] = 0


def position(element, rank, ranking, pair_ranking, inconsistent):
    shift = 0
    if rank+shift == len(ranking):
        ranking.append([element])

    if pair_ranking[element, ranking[rank]].sum() > 0:  # if ranking above at least one of the current rank
        shift = -1

    if pair_ranking[ranking[rank], element].sum() > 0:
        if shift == -1:
            inconsistent.append((element, ranking[rank]))
            return
            # print( ValueError(f"Inconsistent ranking between element {element} and elements {ranking[rank]}")
        else:
            shift = 1

    if shift == 0:
        ranking[rank].append(element)
    
    elif shift == -1:       # insert as new rank above current rank
        ranking.insert(rank, [element])
    
    
    elif shift == 1:   # continue search at next position
        if rank+shift == len(ranking):  # reached end, so append
            ranking.append([element])
        else:
            position(element, rank=rank+1, ranking=ranking, pair_ranking=pair_ranking, inconsistent=inconsistent)  # search next pos
    
ranking = [[0]]
inconsistent = []
for ph in range(1, len(phenos_v)):
    position(ph, rank=0, ranking=ranking, pair_ranking=A, inconsistent=inconsistent)

# print(ranking, "\n", inconsistent)
pickle.dump(ranking, open("ph_ranking_listlist.pickle", "wb"))
pickle.dump(inconsistent, open("ph_inconsistent_listtuple.pickle", "wb"))

ranking = pickle.load(open("ph_ranking_listlist.pickle", "rb"))
inconsistent =pickle.load(open("ph_inconsistent_listtuple.pickle", "rb"))

ranking_ph = []
for rank in ranking:
    ranking_ph.append([phenos_v[ph] for ph in rank])

inconsistent_ph = []    
for tup in inconsistent:
    inconsistent_ph.append([])
    inconsistent_ph[-1].append(phenos_v[tup[0]])
    inconsistent_ph[-1].append([phenos_v[i] for i in tup[1]])
    inconsistent[-1] = tuple(inconsistent[-1])

pickle.dump(ranking_ph, open("ph_ranking_ph_listlist.pickle", "wb"))
pickle.dump(inconsistent_ph, open("ph_inconsistent_ph_listtuple.pickle", "wb"))
print(ranking_ph)
print(inconsistent_ph)
