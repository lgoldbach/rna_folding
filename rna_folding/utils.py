import numpy as np
from itertools import product


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
