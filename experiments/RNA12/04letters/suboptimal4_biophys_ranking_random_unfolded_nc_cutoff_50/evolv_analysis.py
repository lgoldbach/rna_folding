
import pickle
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from scipy.stats import gaussian_kde

from rna_folding.parsing import load_phenotype_and_metric_from_file
from rna_folding.utils import count_bp

fig, (axes) = plt.subplots(ncols=10, nrows=4, figsize=(50, 20)) #, sharey="row", sharex="row")

# x4 = []
# y4 = []
# for i in range(2, 12):
#     f = f"bp_graph{i}/ranking1/nc_graph.pickle"
#     ph_freq_f = f"bp_graph{i}/ranking1/phenotype_distribution.txt"
#     # navigability = f"bp_graph{i}/ranking1/productive_walk_lengths_1000000_seed7339311_seldif1.csv"
#     rug_d_f = f"bp_graph{i}/ranking1/peak_sizes.txt"

#     phenotypes, ph_count = load_phenotype_and_metric_from_file(ph_freq_f)
#     ph_freqs = [c/sum(ph_count) for c in ph_count if c > 0]
#     ph_to_freqs = dict(zip(phenotypes, ph_freqs))
#     ph_to_count = dict(zip(phenotypes, ph_count))
#     nc_graph = pickle.load(open(f, "br"))
#     ph_freqs_sort = np.array(sorted(ph_freqs))[::-1]
    
#     ph_freqs_sort_no_unf = ph_freqs_sort[1:]

#     nc_sizes = []
#     out_deg = []
#     evolvab = []
#     ph_frs = []

#     for node in nc_graph.nodes:
#         nc_sizes.append(nc_graph.nodes[node]["size"])  # size
#         out_deg.append(sum([nc_graph.edges[e]["weight"] for e in nc_graph.edges([node])]))  # edge weight sum
#         # if out_deg[-1] > 10**7:
#         #     print(out_deg[-1], nc_graph.nodes[node]["phenotype"], nc_sizes[-1])
#         evolvab.append(len(np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(node)]))) # unique neighbors


#     x = nc_sizes
#     y = evolvab
#     xy = np.vstack([x,y])
#     z = gaussian_kde(xy)(xy)

#     ax = axes[0, i-2]
#     ax.scatter(np.log10(x), y, c=z, s=100)
#     ax.set_xlabel("Size of NC (log10)")
#     ax.set_ylabel("Evolvability of NC")
#     ax.set_title(f"bp rule {i}", size=20)
 
nc_graph = {}
for i in range(2, 12):
    file = f"bp_graph{i}/ranking1/nc_graph.pickle"
    nc_graph[i] = pickle.load(open(file, "rb"))

ph_count = {}
for i in range(2, 12):
    file = f"bp_graph{i}/ranking1/phenotype_distribution.txt"
    ph, counts = load_phenotype_and_metric_from_file(file)

    ph_count[i] = dict([(ph, c) for ph, c in zip(ph, counts) if c > 0])

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

num_red_mut = {2: 2,
                3: 2,
                4: 6,
                5: 0,
                6: 0,
                7: 0,
                8: 0,
                9: 4,
                10: 2,
                11: 0}

def estimate_evol(freqs, phenos, n, bp):
    ch = np.random.choice(phenos, p=freqs, replace=True, size=n)  # freqs
    uni = np.unique(ch)
    evol = len(uni) # * ((12-num_red_mut[bp])/12)
    return evol

nc_n_est = {}
for i in range(2, 12):
    ph_all = list(nx.get_node_attributes(nc_graph[i], name="phenotype").values())
    phenotypes = np.unique(ph_all)
    nc_n_est[i] = 0
    for ph in phenotypes:
        bp_n = count_bp(ph)
        nc_n_est[i] += mut_graph_n[i]**bp_n

nc_est = {}
for i in range(2, 12):
    nc_est[i] = []
    frac = 0.1
    nc_size = 1000000000
    for j in range(nc_n_est[i]):
        stack = mut_graph_s[i] ** max_bp[i]  # 4: number of pairs
        loop = 4**(12 - (max_bp[i]*2)) # 6: loop size, 4: alphabet size
        nc_size = stack * loop * frac
        if nc_size < 50 or j > nc_n_est[i]:
            print(i)
            break
        nc_est[i].append(nc_size)
        frac = frac - (nc_est[i][-1] / 4**12)


navigs_mean = {2: 73,
                3: 76,
                4: 52,
                5: 79,
                6: 84,
                7: 82,
                8: 73,
                9: 53,
                10: 60,
                11: 66}

for i in range(2, 12):
    nc_g = nc_graph[i]
    nc_sizes = []
    evolvab = []
    outdeg = []
    for node in nc_g.nodes:
        nc_sizes.append(nc_g.nodes[node]["size"])  # size
        evolvab.append(len(np.unique([nc_g.nodes[neigh]["phenotype"] for neigh in nc_g.neighbors(node)]))) # unique neighbors
        outdeg.append(sum([nc_g.edges[e]["weight"] for e in nc_g.edges([node])]))  # edge weight sum

    nc_sizes_sort, evo_sort = zip(*[(nc, e) for nc, e in sorted(zip(nc_sizes, evolvab), reverse=True)])

    x = []
    y = []
    rug = 0
    rug_es = 0
    for c, s in enumerate(nc_est[i]):
        rug_es += s/((np.log10(s)+1)*10)
        try:
            y.append(evo_sort[c])
            x.append(s)
            rug += s/(evo_sort[c]+1)
            rug_es += s/((np.log10(s)+1)*10)

        except IndexError:
            rug_es += s/((np.log10(s)+1)*10)  # always estimate this
    axes[0][0].scatter([rug], [navigs_mean[i]], label=str(i))
    axes[0][1].scatter([rug_es], [navigs_mean[i]], label=str(i))
    
    # xy = np.vstack([x,y])
    # z = gaussian_kde(xy)(xy)
    # axes[1][i-2].scatter(np.log10(x), y, c=z, s=100)

    # print(outdeg, nc_sizes)
    axes[2][i-2].scatter(np.log10(outdeg), np.log10(nc_sizes), label=str(i))
    axes[2][i-2].legend()
    # axes[2][i-2].set_xlabel("Evolvability")
    # axes[2][i-2].set_xlabel("Out degree (log10)")

    # x = outdeg
    # y = evolvab
    # xy = np.vstack([x,y])
    # z = gaussian_kde(xy)(xy)
    # cm = plt.cm.get_cmap('RdYlBu')
    sc = axes[1][i-2].scatter(np.log10(nc_sizes), evolvab)
    # plt.colorbar(sc)
    freqs = [c/4**12 for c in ph_count[i].values()]
    phenos = list(ph_count[i].keys())
    evol_es = []

    rug_es = 0
    # for nc_size in nc_est[i]:
    for nc_size in nc_sizes:
        # freqs_node = [f for f in freqs if f != freq_ph]
        evol_es.append(estimate_evol(freqs, phenos, n=int(nc_size), bp=i))
        rug_es += nc_size/evol_es[-1]

    axes[0][3].scatter([rug_es], [navigs_mean[i]], label=str(i))
    axes[0][3].set_xlabel("Ruggedness")
    axes[0][3].set_ylabel("Navigability")
    
    axes[3][i-2].scatter(evolvab, evol_es, label=str(i))
    axes[3][i-2].set_aspect('equal', adjustable='box')
    axes[3][i-2].plot([3, 35], [3, 35], color="grey")
    axes[3][i-2].set_title(f"Bp rule {i}")

axes[0][0].legend()
axes[0][3].legend()
axes[2][0].legend()
axes[2][1].legend()


plt.tight_layout()
plt.savefig("ana_evolv.png", format="png", dpi=100)



