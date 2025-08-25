

import numpy as np
from rna_folding.base_pairing import BasePairing

def dfs(pair, stack, visited, pair_neighbors): 
    visited[pair] = True
    for neigh in pair_neighbors[pair]:
        if not visited[neigh]:
            stack.append(neigh)
    return stack
    
def connected(pair1, pair2):
    return (pair1[0] == pair2[0]) != (pair1[1] == pair2[1])  #XOR

def find_mut_group(init_pair, possible_pairs, pair_neighbors):
    visited = {pair: False for pair in possible_pairs}
    stack = [init_pair]

    mut_group = []
    while stack:
        pair = stack.pop()
        if not visited[pair]:
            mut_group.append(pair)
            dfs(pair, stack, visited, pair_neighbors)

    return mut_group


graph_path = "../../../../data/graphs/"



def compute_mut_groups(rule_id, bases):
    print("BP: ", rule_id, " bases: ", bases)
    BP = BasePairing(bases=bases, id=rule_id, graph_path=graph_path)
    possible_pairs = []
    for b1 in bases:
        for b2 in bases:
            if BP.pairs(b1, b2):
                possible_pairs.append((b1, b2))

    pair_neighbors = {}
    for pair1 in possible_pairs:
        pair_neighbors[pair1] = []
        for pair2 in possible_pairs:
            if connected(pair1, pair2):
                pair_neighbors[pair1].append(pair2)

    mut_groups = []
    for pair in possible_pairs:
        mut_group = sorted(find_mut_group(pair, possible_pairs, pair_neighbors))
        if mut_group not in mut_groups:
            mut_groups.append(mut_group)

    return mut_groups

bases = "JKLM"

for rule_id in range(1, 12):
    mut_groups = compute_mut_groups(rule_id, bases)
    for mg in mut_groups:
        print(mg, len(mg), "\n")


    

