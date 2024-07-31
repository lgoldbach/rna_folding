import numpy as np
from itertools import product
import networkx as nx


def remove_nonadaptive_edges(gp_graph: nx.graph) -> nx.graph:
    """Remove edges that connect two nodes i and j where node j has higher 
    fitness than node i.

    Args:
        gp_graph (nx.graph):    A networkx graph where nodes have the attribute 
                                "fitness"
    
    Returns:
        gp_graph (nx.graph):    The input graph object without non-adaptive 
                                edges

    """
    for (u, v) in gp_graph.edges:
        if gp_graph.nodes[u]["fitness"] > gp_graph.nodes[v]["fitness"]:
            gp_graph.remove_edge(u, v)

    return gp_graph


def random_fitness_landscape_from_nx_graph(gp_graph: nx.graph, 
                                           peak_phenotype: str) -> nx.graph:
    """Create a random fitness landscape. Assign fitnesses to phenotypes
    from a uniform distribution between 0 and 1. peak_phenotypes always gets
    maximum fitness (1) assigned.

    Args:
        gp_graph (nx.graph):    A networkx graph where nodes have the attribute 
                                "phenotype"
        peak_phenotype (str):   The phenotype to which maximum fitness is 
                                assigned
    
    Returns:
        gp_graph (nx.graph):    The input graph object with "fitness: attribute
                                added to each node.

    """
    phenotypes = set(nx.get_node_attributes(gp_graph, "phenotype").values())

    for ph in phenotypes:
        nodes = [g for g, attr in gp_graph.nodes(data=True) 
                 if attr['phenotype']==ph]  # get all nodes for given phenotype
        F = np.random.uniform(low=0, high=1, size=1)  # [0, 1) interval (1 ex.)
        # assign fitness value to all nodes
        nx.set_node_attributes(gp_graph, dict(zip(nodes, [F]*len(nodes))), "fitness")

    # assign 1 to all peak_phenotype nodes
    nodes = [g for g, attr in gp_graph.nodes(data=True) 
                 if attr['phenotype']==peak_phenotype]
    nx.set_node_attributes(gp_graph, dict(zip(nodes, [1]*len(nodes))), "fitness")

    return gp_graph


def shuffle_array(A, k):
    """Shuffle an array k times. In each step two random elements are swaped
    See Janzer et al., 2023, DOI. 10.1137/22M1530677.

    Args:
        A (list-like): An array or list like mutable object.
    
    Returns:
        A (list-like): Shuffled version of A
    """
    for shuffle in range(k):
        i = np.random.choice(len(A), size=1)
        j = np.random.choice(len(A), size=1)
        A[i], A[j] = A[j], A[i]
    
    return A


def canonical_base_pairs(A, B):
    """Define all possible base-pairs. Use dictionary for quick look-up with hashing

    Args:
        A (str): First base
        B (str): Second base

    Returns:
        bool: True if A and B can pair, False otherwise
    """
    pairs = {("A", "A"): 0, ("A", "U"): 1, ("A", "G"): 0, ("A", "C"): 0,
         ("U", "A"): 1, ("U", "U"): 0, ("U", "G"): 1, ("U", "C"): 0,
         ("G", "A"): 0, ("G", "U"): 1, ("G", "G"): 0, ("G", "C"): 1,
         ("C", "A"): 0, ("C", "U"): 0, ("C", "G"): 1, ("C", "C"): 0}

    return pairs[A, B]

def canonical_adjacency_matrix():
    """Just a hard coded adjacency matrix for canonical base-pairing
    A-U, U-A, G-C, C-G, G-U, U-G,
    The matrix entries are as folloews: (A, U, G, C)x(A, U, G, C),
    so M[0,0] = 0 because A and A can't pair,
    M[2, 3] = 1 because G and C can pair, etc.
    """
    M = np.array([[0, 1, 0, 0],
                  [1, 0, 1, 0],
                  [0, 1, 0, 1],
                  [0, 0, 1, 0]])
    return M

def bp_to_dotbracket(bp: list, l: int) -> str:
    """Turn list of base-pairs into dot-bracket notation

    Args:
        bp (list): List of tuples defining base-pairs
        l (int): Length of sequence

    Returns:
        db (str): Secondary structure in dot-bracket notation

    """
    db = ["."] * l
    for pair in bp:
        db[pair[0] - 1] = "("
        db[pair[1] - 1] = ")"
    db = ''.join(db)
    return db


def dotbracket_to_bp(db: str) -> set:
    """Convert RNA secondary structure from dot-bracket format to list of 
    base-pairs format

    Args:
        db (str): dot-bracket string

    Returns:
        bp (str): list of base-pair tuples (numeric).

    """
    opening_stack = []
    bp = []
    for i, site in enumerate(db):
        if site == "(":
            opening_stack.append(i)
        elif site == ")":
            bp.append((opening_stack.pop(), i))
    return set(bp)


def count_bp(seq):
    """Count number of base-pairs in sequence

    Args:
        seq (str): dot-bracket notation string, e.g. "(..)(..)"

    Returns:
        int: Number of base-pairs in string
    """
    bp_count = 0
    for i in seq:
        if i == "(":
            bp_count += 1
    return bp_count


def combinatorically_complete_genotypes(l, a):
    """

    Args:
        l (int): Sequence length.
        a (str): Alphabet as a continuous string, e.g. "AUGC".

    Returns:
        g (list): List of genotypes (str).

    """
    g = product(list(a), repeat=l)
    return g


def dotbracket_to_genotype(dotbracket: str, 
                           base_pair: str = "GC",
                           random: bool = True,
                           seed: int = None) -> str:
    """Take a dot-bracket RNA phenotype and translate it into a sequence of
    that is consistent with the given base-pairing.
    E.g. (((...))) -> GGG...CCC

    Args:
        dotbracket (str): A RNA phenotype in dot-bracket notation, e.g. ((..))
        base_pair (str): What base-pair to translate to. Limited to one pair.
                         Default = "GC"
        random (bool): If true, bases are picked at random. If False, bases are
                       picked so that the first base given in <base_pair> will
                       be assigned to all opening brackets and the second one
                       to all closing brackets. The unpaired regions will be 
                       assigned the first base.
        seed (int): Provide random seed for assigning bases to positons.

    Returns:
        str: Genotype that is consistent with the base-pairing.
    
    """
    if seed:
        np.random.seed(seed)
    stack = []
    genotype = []
    for site, db in enumerate(dotbracket):
        if db == "(":
            stack.append(site)
            if random:
                base = base_pair[np.random.choice([0, 1])]
            else:
                base = base_pair[0]
        elif db == ")":
            if random:
                paired_site = stack.pop()
                paired_base = genotype[paired_site]
                base = base_pair[paired_base == base_pair[0]]  # get matching base
            else:
                base = base_pair[1]
        elif db == ".":
            if random:
                base = base_pair[np.random.choice([0, 1])]
            else:
                base = base_pair[0]
        genotype.append(base)
    
    return ''.join(genotype)

def dotbracket_to_genotype_random(dotbracket: str,
                           base_pairs: str = ["AU", "GC"],
                           seed: int = None) -> str:
    """Take a dot-bracket RNA phenotype and translate it into a sequence of
    that is consistent with the given base-pairing.
    E.g. (((...))) -> GGG...CCC

    Args:
        dotbracket (str): A RNA phenotype in dot-bracket notation, e.g. ((..))
        base_pairs (list): List of base-pairs, default: ["AU", "GC", "GU"].
        seed (int): Provide random seed for assigning bases to positons.

    Returns:
        str: Genotype that is consistent with the base-pairing.

    """
    if seed:
        np.random.seed(seed)

    bases = list("".join(base_pairs))
    complement_stack = []

    genotype = ""
    for site, db in enumerate(dotbracket):
        if db == "(":
            base_pair = np.random.choice(base_pairs)  # pick pair
            b1 = np.random.choice([0, 1])  # pick base from pair
            b2 = b1 == 0  # get complementary
            base = base_pair[b1]
            complement_stack.append(base_pair[b2])  # store complementary on stack
        elif db == ")":
            base = complement_stack.pop()
        elif db == ".":
            base = np.random.choice(bases)
        else:
            raise ValueError(f"Unknown character in dotbracket: {db}")
        genotype += base
    return genotype


def dict_to_gpmap(ph_to_gt: dict, file: str) -> None:
    """Take a dict that maps phenotype to list of genotypes and save it
    as a space-separated "c"sv file, where each line looks like this:
    "{ph} {gt_id} {gt_id} {gt_id}"

    Args:
        ph_to_gt (dict): _description_
        file (str): _description_
    """
    # Write to output file (
    with open(file, "w") as file_out:
        for p in ph_to_gt:
            line = p + " " + " ".join(map(str, ph_to_gt[p])) + "\n"
            file_out.write(line)
    file_out.close()


def ranked_ph_distribution(ph_distr_file, log=False) -> tuple:
    """Rank phenotypes by their count and return log10 frequency

    Args:
        ph_distr_file (str): file containing phenotypes and their count, e.g.:
                                ph1 19229
                                ph2 123123123
                                ph3 212
                                ...
    Returns:
        (phenotypes, distr): phenotypes: Phenotypes in descending order of 
                                         their frequency
                             distr: Frequencies in descending order. (matches
                                    order of phenotype

    """
    # load data and get second column (fist only contains phenotype IDs)
    
    phenotypes, distr = load_phenotype_and_metric_from_file(ph_distr_file)
    distr = distr / np.sum(distr)
    if log:
        distr = np.log10(distr)

    order = np.argsort(distr)[::-1]
    distr = distr[order]
    phenotypes = phenotypes[order]
    return phenotypes, distr


def load_phenotype_and_metric_from_file(file: str, dtype=float):
    """Take a file in the common phenotype (col1) metric (col2) data-type 
    I am using and reat it as two array.
    Example file:
    ((...)) 0.8
    (.....) 0.7
    ...

    Args:
        file (str): Path to the file

    Retruns:
        phentypes, data
    """
    file_data = np.loadtxt(file, dtype=str)
    if file_data.ndim == 1:  # in case there is only one phenotype
        file_data = np.expand_dims(file_data, axis=0)
    phenotypes = file_data[:,0]
    distr = file_data[:,1].astype(dtype)

    return phenotypes, distr
