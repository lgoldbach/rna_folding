

import numpy as np
import pickle
import networkx as nx
import matplotlib.pyplot as plt
from rna_folding.utils import count_bp

nc_graph = {}
for i in range(2, 12):
    file = f"bp_graph{i}/ranking1/nc_graph.pickle"
    nc_graph[i] = pickle.load(open(file, "rb"))

mut_graph_n = {2: 2,
               3: 2,
               4: 2,
               5: 4,
               6: 1,
               7: 2,
               8: 1,
               9: 2,
               10: 1,
               11: 1}

mut_graph_s = {2: 1,
               3: 2,
               4: 3,
               5: 1,
               6: 6,
               7: 3,
               8: 8,
               9: 4,
               10: 10,
               11: 12}

max_bp = {2: 2,
            3: 2,
            4: 2,
            5: 2,
            6: 2,
            7: 2,
            8: 3,
            9: 3,
            10: 4,
            11: 4}

ph_n = 32
# nc_n_est_2p = {}
# for i in range(2, 12):
#     nc_n_per_ph = mut_graph_n[i]**2
#     nc_n_est_2p[i] = nc_n_per_ph * ph_n

# nc_n_est_3p = {}
# for i in range(2, 12):
#     nc_n_per_ph = mut_graph_n[i]**3
#     nc_n_est_3p[i] = nc_n_per_ph * ph_n

# nc_n_est_4p = {}
# for i in range(2, 12):
#     nc_n_per_ph = mut_graph_n[i]**4
#     nc_n_est_4p[i] = nc_n_per_ph * ph_n

nc_n_est = {}
for i in range(2, 12):
    ph_all = list(nx.get_node_attributes(nc_graph[i], name="phenotype").values())
    phenotypes = np.unique(ph_all)
    nc_n_est[i] = 0
    for ph in phenotypes:
        bp_n = count_bp(ph)
        nc_n_est[i] += mut_graph_n[i]**bp_n

# power_law_est = {}
# for i in range(2, 12):
#     ph_all = list(nx.get_node_attributes(nc_graph[i], name="phenotype").values())
#     phenotypes = np.unique(ph_all)
#     nc_n = 0
#     for ph in phenotypes:
#         bp_n = count_bp(ph)
#         nc_n += mut_graph_n[i]**bp_n
#     A = 1/np.log(nc_n)
#     power_law_est[i] = []
#     for r in range(1, len(nc_graph[i].nodes)+1):
#         power_law_est[i].append(A/r * 1000000)

largest_nc_est = {}
for i in range(2, 12):
    stack = mut_graph_s[i] ** max_bp[i]  # 4: number of pairs
    loop = 4**4 # 4: alphabet size, 4: loop size
    nc_size = stack * loop
    largest_nc_est[i] = nc_size * 0.1



nc_est = {}
for i in range(2, 12):
    nc_est[i] = []
    frac = 0.1

    nc_size = 1000000000
    for j in range(nc_n_est[i]+100):
        stack = mut_graph_s[i] ** max_bp[i]  # 4: number of pairs
        loop = 4**(12 - (max_bp[i]*2)) # 6: loop size, 4: alphabet size
        nc_size = stack * loop * frac
        if nc_size < 50 or j > nc_n_est[i]:
            print(i)
            break
        nc_est[i].append(nc_size)
        frac = frac - (nc_est[i][-1] / 4**12)

fig, axes = plt.subplots(nrows=2, ncols=10, figsize=(50, 10), sharey="row")

cutoff = 50

for ax, bp in zip(axes[0], nc_graph):
    # print([nc_graph[bp].nodes[node]["phenotype"] for node in nc_graph[bp] if nc_graph[bp].nodes[node] == "............"])
    count = 0
    count_all = 0
    for nc in nc_graph[bp].nodes:
        sss = nc_graph[bp].nodes[nc]["size"]
        if sss > cutoff:
            count += nc_graph[bp].nodes[nc]["size"]
        count_all += sss

    y = list(reversed(sorted([nc_graph[bp].nodes[nc]["size"] for nc in nc_graph[bp].nodes if nc_graph[bp].nodes[nc]["size"] > cutoff])))

    ax.scatter(range(len(y)), np.log10(y))
    # ax.plot([0, nc_n_est_2p[bp]], np.log10([y[0], 0]), label="2")
    # ax.plot([0, nc_n_est_3p[bp]], [y[0], 0], label="3")
    # ax.plot([0, nc_n_est_4p[bp]], [y[0], 0], label="4")
    # ax.plot([0, nc_n_est[bp]], [np.log10(y[0]), np.log10(cutoff)], label="exact", color="green")
    # ax.plot([0, len(y)],[np.log10(largest_nc_est[bp]), np.log10(largest_nc_est[bp])], color="red")
    ax.scatter(range(len(nc_est[bp])), np.log10(nc_est[bp]), color="red", s = 3, marker="x")

    ax.set_title(f"bp rule {bp}", size=20)

bins = range(0, 100000, 100)
for ax, bp in zip(axes[1], nc_graph):
    sizes = []
    for nc in nc_graph[bp].nodes:
        sss = nc_graph[bp].nodes[nc]["size"]
        if sss > cutoff:
            sizes.append(nc_graph[bp].nodes[nc]["size"])
    
    logbins = np.geomspace(np.min(sizes), np.max(sizes), 20)

    ax.hist(sizes, bins=logbins)

    ax.set_xscale('log')
    ax.set_yscale('log')
plt.legend()

plt.savefig("explore_nc.pdf", format="pdf", dpi=30)
