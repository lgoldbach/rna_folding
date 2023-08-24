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

