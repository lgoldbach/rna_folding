import numpy as np
from rna_folding.gp_map import GenotypePhenotypeGraph


def kimura_fixation(s: float, N: int):
    """Formula to compute kimuara's fixation probability

    Args:
        s (float): Selection coefficient
        N (int): Population size

    Returns:
        float: fixation probability (float in [0, 1]).
    
    """
    if s == 0:
        s = -10**-10

    p = (1-np.exp(-2*s))/(1-np.exp(-2*N*s))
    return p


def adaptive_walk(gpmap: GenotypePhenotypeGraph, 
                  starting_genotype,
                  fitness_function,
                  max_steps,
                  population_size,
                  fixation_function,
                  rng) -> list:
    path = [starting_genotype]
    if fitness_function[gpmap.nodes[path[-1]]["phenotype"]] == 1:
        return path
    
    while len(path) < max_steps:
        candidate = rng.choice(gpmap._neighbors(path[-1]))
        f1 = fitness_function[gpmap.nodes[path[-1]]["phenotype"]]
        f2 = fitness_function[gpmap.nodes[candidate]["phenotype"]]
        if f2 == 1:  # found target phenotype
            path.append(candidate)  # append and break
            break
        s = f2-f1
        p = fixation_function(s, N=population_size)
        # this ignores waiting times for mutations, in a sense that we only
        # track the occurences of mutatons and if the mutation is rejected
        # we instead append the same genotype again. From this sequence we
        # can infer waiting times from a realistic mutation probability
        # with high population size it will be impossible to traverse neutral
        # nets with this approach
        if rng.uniform() < p:
            path.append(candidate)
        else:
            path.append(path[-1])
    return path

def productive_adaptive_walk(gpmap: GenotypePhenotypeGraph, 
                  starting_genotype,
                  fitness_function,
                  max_steps,
                  population_size,
                  fixation_function,
                  rng) -> list:
    path = [starting_genotype]
    if fitness_function[gpmap.nodes[path[-1]]["phenotype"]] == 1:
        return path
    
    while len(path) < max_steps:
        probs = []
        f1 = fitness_function[gpmap.nodes[path[-1]]["phenotype"]]
        neighbors = gpmap._neighbors(path[-1])
        for neigh in neighbors:
            f2 = fitness_function[gpmap.nodes[neigh]["phenotype"]]
            s = f2-f1
            p = fixation_function(s, N=population_size)
            probs.append(p)
        
        if sum(probs) == 0:  # no way to go
            break
        normed_probs = np.array(probs) / sum(probs)
        candidate = rng.choice(neighbors, p=normed_probs)

        path.append(candidate)
        if fitness_function[gpmap.nodes[candidate]["phenotype"]] == 1:  # found target phenotype
            break
    return path


def greedy_adaptive_walk(gpmap: GenotypePhenotypeGraph, 
                  starting_genotype,
                  fitness_function,
                  max_steps,
                  rng) -> list:
    path = [starting_genotype]
    if fitness_function[gpmap.nodes[path[-1]]["phenotype"]] == 1:
        return path
    
    while len(path) < max_steps:
        s_coeffs = []
        f1 = fitness_function[gpmap.nodes[path[-1]]["phenotype"]]
        neighbors = gpmap._neighbors(path[-1])
        for neigh in neighbors:
            f2 = fitness_function[gpmap.nodes[neigh]["phenotype"]]
            s = f2-f1
            s_coeffs.append(s)
        
        s_coeffs = np.array(s_coeffs)
        m = np.max(s_coeffs)  # find max
        if m >= 0:
            maxima = np.where(s_coeffs==m)[0]  # find all maxima
            next_gt = rng.choice(maxima)  # pick one at random
            path.append(neighbors[next_gt])
            if fitness_function[gpmap.nodes[path[-1]]["phenotype"]] == 1:  # found target phenotype
                break
        else:
            break  # no higher or equal fitness found, end path
        
    return path


def greedy_adaptive_walk_no_neutral(gpmap: GenotypePhenotypeGraph, 
                  starting_genotype,
                  fitness_function,
                  max_steps,
                  rng) -> list:
    path = [starting_genotype]
    if fitness_function[gpmap.nodes[path[-1]]["phenotype"]] == 1:
        return path
    
    while len(path) < max_steps:
        s_coeffs = []
        f1 = fitness_function[gpmap.nodes[path[-1]]["phenotype"]]
        neighbors = gpmap._neighbors(path[-1])
        for neigh in neighbors:
            f2 = fitness_function[gpmap.nodes[neigh]["phenotype"]]
            s = f2-f1
            s_coeffs.append(s)
        
        s_coeffs = np.array(s_coeffs)
        m = np.max(s_coeffs)  # find max
        if m > 0:
            maxima = np.where(s_coeffs==m)[0]  # find all maxima
            next_gt = rng.choice(maxima)  # pick one at random
            path.append(neighbors[next_gt])
            if fitness_function[gpmap.nodes[path[-1]]["phenotype"]] == 1:  # found target phenotype
                break
        else:
            break  # no higher fitness found, end path
        
    return path

