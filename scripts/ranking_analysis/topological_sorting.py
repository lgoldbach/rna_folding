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

    # create a networkx directed graph using adjacecy matrix A
    G = nx.from_numpy_array(A, create_using=nx.DiGraph)
    # Some nodes may not be connected, i.e. some phenotypes have not appeared
    # in a suboptimal set together. This will lead to non-unique topological
    # sortings later on. To prevent that we add two edges connecting the nodes
    # in opposite directions to artificially create a cycle between them, so
    # that they will be lumped together as unsortable later in the algorithm
    for n_i in G.nodes:
        for n_j in G.nodes:
            if (n_i, n_j) not in G.edges and (n_j, n_i) not in G.edges and n_i != n_j:
                print(n_i, n_j)
                G.add_edges_from([(n_i, n_j), (n_j, n_i)])

    # Find all the cycles and remove edges of these cycles (not nodes!) until
    # there are no more cycles and we can do topological sort
    cycle = True
    cycles = []
    while cycle:
        try:
            c = nx.find_cycle(G)
            cycles.append(c)
            for edge in c: # remove all the edges of the cycle
                G.remove_edge(edge[0], edge[1])  
        except nx.exception.NetworkXNoCycle:
            cycle = False
    
    # get all the nodes that are part of a cycle
    nodes_in_cycle = set([node for cycle in cycles for edge in cycle for node in edge])
    
    # Get a topolical sorting on G. This is probably not unique all the way 
    # through but we do not care since we only care about the unique parts
    initial_sort = list(nx.topological_sort(G))

    # extract the parts of the topoligcal sorting that are unique
    # The logic is: Loop over the sorting from top to bottom. Continue for as
    # long as nodes have not been part of a cycle. These nodes have to 
    # have a topological sorting. Once we hit a node that was part of a cycle
    # we loop until we hit a nodes wasn't part of a cycle. All the nodes we
    # found until then were part of a cycle and thus cannot be sorted. 
    # All these nodes are just lumped into one cluster.

    sort = []  # the overall sorting
    all_sub_sorts = [[]]  # track the sorts for subsets of nodes
    no_sort = []  # track the nodes that are not sortable
    k = 0  # track the start of each new interval
    for i, node in enumerate(initial_sort, start=1):
        if node not in nodes_in_cycle:  # this starts an interval of sortable nodes
            if len(no_sort) > 0:  # if we had an unsortable interval before, first append that intervall to the main sort
                sort.append(no_sort)
                no_sort = []  # reset no_sort
            G_sub = nx.induced_subgraph(G, initial_sort[k:i])  # induce a subgraph on the subinterval. 
     
            all_sub_sorts = list(nx.all_topological_sorts(G_sub))
            # we have more than one way to sort this, which should not happen 
            # in this case as all the pairs node are connected in the beginning 
            # and thus cycle free sets of nodes should have a unique sorting.
            # alternatively one can say that if the nodes form a hamiltonian 
            # path, there is a unique sorting
            if len(all_sub_sorts) > 1:
                print(f"Non-unique topoligical sorting detected. Should not "
                      f"happen if there is a hamiltonian path including all the "
                      f"nodes, which is the case in a fully connected graph.\n "
                      f"Only plausible case is if the two structures never met "
                      f"in a suboptimal set\n"
                      f"All sortings: {all_sub_sorts}.\nList of unsorted "
                      f"previous interval if applicable {no_sort}\n "
                      f"Inital sort for this interval: {initial_sort[k:i+1]}\n "
                      f"Whole initial sort: {initial_sort}")
            
                
                            

                print(G_sub.edges, 32 in nodes_in_cycle, 1 in nodes_in_cycle, 41 in nodes_in_cycle)
                break
        else:
            k = i
            if len(all_sub_sorts[0]) > 0:
                for n in all_sub_sorts[0]:
                    sort.append([n])
            no_sort.append(node)
            all_sub_sorts = [[]]  # reset
    if len(no_sort) > 0:
        sort.append(no_sort)

    with open(args.output, "w") as file:
        for cluster in sort:
            file.write(" ".join([str(node) for node in cluster]) + "\n")
