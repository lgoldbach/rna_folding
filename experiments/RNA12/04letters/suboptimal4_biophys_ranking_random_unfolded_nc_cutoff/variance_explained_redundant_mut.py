import pickle
import numpy as np
import matplotlib.pyplot as plt
from rna_folding.parsing import many_to_one_map_from_file_to_dict

# load gp maps

# for each gp map:
# find 10(00) foldable genotypes
# generate all possible point mutants for site 1, 3 and 5
# categorize them into redundant and non redundant mutations
# count fraction of mutations that produce new phenotype
# compare overall fraction and without redundant mut
# show that without redundant mut, they are similar in evolvability

fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))
bp_ids = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
# bp_ids = [4]

allowed_mut = {1: {'J': ['K', 'L', 'M'], 'K': ['J', 'M'], 'L': ['J', 'M'], 'M': ['J', 'L', 'K']},
               2: {'J': ['L', 'M'], 'K': ['L', 'M'], 'L': ['J', 'K', 'M'], 'M': ['J', 'L', 'K']},
               3: {'J': ['K', 'L', 'M'], 'K': ['J', 'L', 'M'], 'L': ['J', 'K', 'M'], 'M': ['J', 'L', 'K']},
               4: {'J': ['K', 'L', 'M'], 'K': ['J', 'L', 'M'], 'L': ['J', 'K', 'M'], 'M': ['J', 'L', 'K']},
               5: {'J': ['K', 'L', 'M'], 'K': ['J', 'L', 'M'], 'L': ['J', 'K', 'M'], 'M': ['J', 'L', 'K']},
               6: {'J': ['K', 'L', 'M'], 'K': ['J'], 'L': ['J'], 'M': ['J']},
               7: {'J': ['L', 'M'], 'K': ['L', 'M'], 'L': ['J', 'K'], 'M': ['J', 'K']},
               8: {'J': ['K', 'L', 'M'], 'K': ['J', 'L', 'M'], 'L': ['J', 'K', 'M'], 'M': ['J', 'L', 'K']},
               9: {'J': ['L', 'M'], 'K': ['L', 'M'], 'L': ['J', 'K', 'M'], 'M': ['J', 'L', 'K']},
               10: {'J': ['K', 'L', 'M'], 'K': ['J', 'L', 'M'], 'L': ['J', 'K', 'M'], 'M': ['J', 'L', 'K']}}

# allowed_mut = {'J': [], 'K': [], 'L': [], 'M': []}
def create_all_neighbors_allowed_mut(genotype, allowed_mut: dict):
        neighbors = []
        for site, state in enumerate(genotype):
            for l in allowed_mut[state]:  # get allowed mutations from this nucleotide (state)
                neigh = gt[:site] + l + gt[site + 1:]
                neighbors.append(neigh)
        return neighbors

def get_neighbors(gp_map, genotype):
        neighbors = []
        for site, state in enumerate(genotype):
            for l in gp_map.alphabet:  # get allowed mutations from this nucleotide (state)
                if l != state:
                    neigh = gt[:site] + l + gt[site + 1:]
                    neighbors.append(neigh)
        return neighbors

for i, bpid in enumerate(bp_ids, start=1):
    gp_map_path = f"bp_graph{bpid}/ranking1/gp_map.pickle"
    
    gp_map = pickle.load(open(gp_map_path, "rb"))

    nc_to_gt_path = f"bp_graph{bpid}/ranking1/nc_to_gt.txt"
    # gt_to_nc = many_to_one_map_from_file_to_dict(file=nc_to_gt_path, 
    #                                              source_type=str,
    #                                              target_type=int,
    #                                              delimiter=" ", skip_first=True)

    
    nc_to_gt = {}
    with open(nc_to_gt_path, "r") as f:
        for line_ in f:
            line = line_.strip().split(" ")
            nc = line[0]
            gts = line[2:]
            nc_to_gt[nc] = gts
    
    # compute nc evo factor per nc which is the factor by which we downsize the nc according to the fract
    nc_evolvs_all = []
    nc_evolvs_allowed = []
    nc_evolv_all = {nc: {} for nc in list(nc_to_gt.keys())}
    nc_evolv_allowed = {nc: {} for nc in list(nc_to_gt.keys())}
    for nc, gts in nc_to_gt.items():
        if len(gts) < 50:
            continue
        neigh_count_all = 0
        neigh_count_allowed = 0
    gt_sample = np.random.choice(gts, size=50, replace=False)
    nc_ph = gp_map.map(gts[0])  # get phenotype of the nc
    for gt in gt_sample:
        neighbors_all = get_neighbors(gp_map, gt)
        neighbors_allowed = create_all_neighbors_allowed_mut(gt, allowed_mut = allowed_mut[i])
        for n in neighbors_all:
            neigh_count_all += 1
            try:
                p = gp_map.map(n)
            except KeyError:  # unfolded neighbor, does not matter
                continue
            if p not in nc_evolv_all[nc]:  # add phenotype if it is not in there yet
                nc_evolv_all[nc][p] = 1
        for n in neighbors_allowed:
            neigh_count_allowed += 1
            try:
                p = gp_map.map(n)
            except KeyError:  # unfolded neighbor, does not matter
                continue
            if p not in nc_evolv_allowed[nc]:  # add phenotype if it is not in there yet
                nc_evolv_allowed[nc][p] = 1
            

        nc_evolvs_all.append(len(nc_evolv_all[nc].keys())/neigh_count_all)  # count unqiue phenotypes
        nc_evolvs_allowed.append(len(nc_evolv_allowed[nc].keys())/neigh_count_allowed)  # count unqiue phenotypes

    ax1.scatter(i, np.mean(nc_evolvs_all))
    ax2.scatter(i, np.mean(nc_evolvs_allowed))


    
    # def create_all_neighbors(gt, alphabet, sites):
    #     neighbors = []
    #     for site in sites:
    #         for l in alphabet:
    #             if gt[site] != l:
    #                 neigh = gt[:site] + l + gt[site + 1:]
    #                 neighbors.append(neigh)
    #     return neighbors

    

    # nc_size = {nc: 0 for nc in list(gt_to_nc.values())}
    # nc_neighs = {nc: 0 for nc in list(gt_to_nc.values())}
    # nc_evolv = {nc: {} for nc in list(gt_to_nc.values())}
    # fraction_of_new_per_gt_all = []
    # fraction_of_new_nc_avg = {nc: 0 for nc in list(gt_to_nc.values())}
    # for j, gt in enumerate(gt_to_nc): 
    #     nc = gt_to_nc[gt]
    #     p_ref = gp_map.map(gt)

    #     diff = 0
    #     neighbors = gp_map._neighbors(gt)
    #     for n in neighbors:
    #         p = gp_map.map(n)
    #         if p_ref != p:
    #             diff += 1
                
    #         if p not in nc_evolv[nc]:  # add phenotype if it is not in there yet
    #             nc_evolv[nc][p] = 1
    #         nc_neighs[nc] += 1
    #     fraction_of_new_per_gt = diff/len(neighbors)

    #     fraction_of_new_per_gt_all.append(fraction_of_new_per_gt)

    #     fraction_of_new_nc_avg[nc] += fraction_of_new_per_gt
    #     nc_size[nc] += 1

    # for nc in fraction_of_new_nc_avg:
    #     fraction_of_new_nc_avg[nc] = fraction_of_new_nc_avg[nc]/nc_size[nc]
    

    # fract_of_new_nc_avg_avg = np.mean(list(fraction_of_new_nc_avg.values()))
    # ax1.scatter(i, fract_of_new_nc_avg_avg)
    # ax2.scatter(i, np.mean(fraction_of_new_per_gt_all))

plt.savefig("mut_impact.pdf", format="pdf")
        
for i in range(10):
    print("AAASDASDAS" * 1)



# evolvs = []
# bp_ids = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
# for i, bpid in enumerate(bp_ids):
#     gp_map_path = f"bp_graph{bpid}/ranking1/gp_map.pickle"

#     gp_map = pickle.load(open(gp_map_path, "rb"))

#     gts = []

#     gts_ = np.random.choice(gp_map.genotypes, size=10000, replace=False)

#     def create_all_neighbors(gt, alphabet, sites):
#         neighbors = []
#         for site in sites:
#             for l in alphabet:
#                 if gt[site] != l:
#                     neigh = gt[:site] + l + gt[site + 1:]
#                     neighbors.append(neigh)
#         return neighbors

#     diff = 0
#     folded = 0
#     total = 0
#     same = 0
#     for gt in gts_:
#         p_ref = gp_map.map(gt)
#         # print("\n\nref:", p_ref)
#         neighbors = create_all_neighbors(gt, gp_map.alphabet, sites=range(12))
#         total += len(neighbors)
#         folded += len(neighbors)
#         for n in neighbors:
#             try:
#                 # print(gp_map.map(n))
#                 p = gp_map.map(n)
#                 if p_ref != p:
#                     diff += 1
#                 elif p_ref == p:
#                     same += 1
#             except KeyError:
#                 folded -= 1

#     try:
#         evol = diff/same
#     except ZeroDivisionError:
#         evol = 0
#         print("AAA")
#     evolvs.append(evol)
#     print(i, folded/total)

# fig, ax = plt.subplots()
# ax.scatter(np.arange(1, 11), evolvs)

# plt.savefig("mut_impact.pdf", format="pdf")
        



