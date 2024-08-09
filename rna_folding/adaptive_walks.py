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
        # np.random.seed(12343124*len(path)**4)
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

