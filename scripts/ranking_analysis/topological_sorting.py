#!/usr/bin/env python

import argparse
import pickle
import networkx as nx
import numpy as np

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input consensus matrix, 2D "
                        "numpy array in pickle format", required=True)
    parser.add_argument("-o", "--output", help="The topologically sorted "
                        "list ", required=True)
    parser.add_argument("-p", "--phenotypes", help="phenotypes ", required=True)

    args = parser.parse_args()
    
    with open(args.phenotypes, "r") as file:
        phenotypes = np.array([line.strip() for line in file])

    A = pickle.load(open(args.input, "rb"))

    # def position(element, rank, ranking, pair_ranking, inconsistent):
    #     shift = 0
    #     if rank+shift == len(ranking):
    #         ranking.append([element])

    #     if pair_ranking[element, ranking[rank]].sum() > 0:  # if ranking above at least one of the current rank
    #         shift = -1

    #     if pair_ranking[ranking[rank], element].sum() > 0:
    #         if shift == -1:
    #             inconsistent.append((element, ranking[rank]))
    #             return
    #             # print( ValueError(f"Inconsistent ranking between element {element} and elements {ranking[rank]}")
    #             # also check if the inconsistency is with two others. In that case this element can be used to sort the other two
    #         else:
    #             shift = 1

    #     if shift == 0:  # here I have to see if it is outranked by lower ones, so actually continue the search
    #         ranking[rank].append(element)
        
    #     elif shift == -1:       # insert as new rank above current rank
    #         ranking.insert(rank, [element])
        
        
    #     elif shift == 1:   # continue search at next position
    #         if rank+shift == len(ranking):  # reached end, so append
    #             ranking.append([element])
    #         else:
    #             position(element, rank=rank+1, ranking=ranking, pair_ranking=pair_ranking, inconsistent=inconsistent)  # search next pos
    
    
 
    # from random import sample
    # phs = list(range(1, A.shape[0]))
    # print(phs, "\n\n")
    # for i in range(1):
    #     phs_r = sample(phs, len(phs))
    #     ranking = [[phs_r[0]]]
    #     inconsistent = []
    #     print(f"Phenotype order: {phs_r}\n")
    #     for ph in phs_r:
    #         position(ph, rank=0, ranking=ranking, pair_ranking=A, inconsistent=inconsistent)

    #     print(f"Ranking: {ranking}\nInconsistent: {inconsistent}\n\n")

    # beats = A.sum(axis=1)
    # beaten = A.sum(axis=0)
    # win_ratio = beats/beaten
    # game_count = beats+beaten
    # rank = np.argsort(win_ratio)
    # # np.set_printoptions(formatter={'float': lambda x: "{0:0.3f}".format(x)})
    # print(game_count[rank], "\n\n")
    # print(beats[rank], "\n", beaten[rank], "\n\n")
    # print(rank, "\n", np.sort(win_ratio), "\n\n")
    # with open("win_ratio_ranking.txt", "w") as file:
    #     for i, p in enumerate(phenotypes[rank][::-1]):
    #         file.write(p + "\n")

    # # B = np.tril(A)
    # B = A
    # idx = np.where(B == 1)
    # with open("consensus_edges.csv", "w") as f:
    #     f.write("SOURCE,TARGET\n")
    #     for e in list(zip(*idx)):
    #         f.write(phenotypes[e[0]] + "," + phenotypes[e[1]] + "\n")
    
    # print(A[9], A[9].sum())
    
    

     ### DEBUG CODE
    # import numpy as np
    # A = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]])

    # G = nx.from_numpy_array(A, create_using=nx.DiGraph)


    # import matplotlib.pyplot as plt
    # nx.draw(G, pos=nx.spring_layout(G))
    # plt.draw()

    # plt.savefig("digraph.pdf", format="pdf", dpi=20)
    
    ## DEBUG CODE

    G = nx.from_numpy_array(A, create_using=nx.DiGraph)

    cycle = True
    cycles = []
    while cycle:
        try:
            c = nx.find_cycle(G)
            cycles.append(c)
            # for edge in c:
                # print(phenotypes[edge[0]], phenotypes[edge[1]])
            # print(c)
            for edge in c:
                G.remove_edge(edge[0], edge[1])
        except nx.exception.NetworkXNoCycle:
            cycle = False
    
    nodes_in_cycle = set([node for cycle in cycles for edge in cycle for node in edge])
    print(nodes_in_cycle, len(nodes_in_cycle))
    
    print(A.shape)
    initial_sort = list(nx.topological_sort(G))
    
    for i, node in enumerate(initial_sort):
        if node not in nodes_in_cycle:
            G_sub = nx.induced_subgraph(G, initial_sort[0:i+1])
            all_sorts = list(nx.all_topological_sorts(G_sub))
        
            print(all_sorts)  
        else:
            break  

    # top17 = [36, 46, 40, 29, 0, 19, 44, 9, 1, 7, 6, 20, 26, 24, 41, 42, 30]
    # G_17 = nx.induced_subgraph(G, top17)
    # s_gen = nx.all_topological_sorts(G_17)
    # s = 1
    # count = 0
    # while s:
    #     try:
    #         s = next(s_gen)
    #         print(s)
    #         count += 1
    #     except StopIteration:
    #         print(count)
    #         break
    
    # s_gen = nx.all_topological_sorts(G_17)
    # print(s)
    # edges = []
    # for i, j in enumerate(s[:-1]):
    #     edges.append((j, s[i+1]))
    # edges = [(36, 46), (46, 40), (40, 29), (29, 0), (0, 19), (19, 44), (44, 9), (9, 7), (9, 6), (9, 1), (7, 20), (6, 20), (1, 20), (20, 26), (26, 24), (24, 41), (41, 42), (42, 30)]

    # print(cycles)
    # cycles_merged = []
    # for c in cycles:
    #     for e in c:
    #         cycles_merged.append(e)
    # edges_cycles = cycles_merged + edges
    # with open("consensus_edges_w_cycles.csv", "w") as f:
    #     f.write("SOURCE,TARGET\n")
    #     for e in edges_cycles:
    #         f.write(phenotypes[e[0]] + "," + phenotypes[e[1]] + "\n")

    # with open("consensus_edges.csv", "w") as f:
    #     f.write("SOURCE,TARGET\n")
    #     for e in edges:
    #         f.write(phenotypes[e[0]] + "," + phenotypes[e[1]] + "\n")


    # # with open(args.output, "w") as file:
    # #     for s in S:
    # #         file.write(" ".join(s))
        