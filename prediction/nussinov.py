"""
Implementation of Nussinov's algorithm for RNA structure prediction
Author: Leander Goldbach
GitHub: lgolbach

"""

import sys
import numpy as np


# Define all possible base-pairs. Use dictionary for quick look-up with hashing
pairs = {("A", "A"): 0,
         ("A", "U"): 1,
         ("A", "G"): 0,
         ("A", "C"): 0,
         ("U", "A"): 1,
         ("U", "U"): 0,
         ("U", "G"): 1,
         ("U", "C"): 0,
         ("G", "A"): 0,
         ("G", "U"): 1,
         ("G", "G"): 0,
         ("G", "C"): 1,
         ("C", "A"): 0,
         ("C", "U"): 0,
         ("C", "G"): 1,
         ("C", "C"): 0}


def initialize_matrix(L: int):
    """Initialize a LxL+1 matrix with zeros on diagonal and the diagonal below.

    Args:
        L (int) : Sequence length.

    Returns:
        M (np.ndarray) : 2d array of integers

    """
    P = np.zeros((L+1, L+1), dtype=int)

    # x, y = np.diag_indices(L)
    # P[x, y] = 0    # set diagonal to zero
    # P[x[1:], y[1:]-1] = 0    # set lower one-off diagonal to zero

    return P


def fill_matrix(P: np.ndarray, S: str, q: int = 1):
    """Main step of Nussinov's algorithm, i.e. filling the P matrix to find maximum base-pairing

    Args:
        P (np.ndarray): 2D array that should be empty except for zeros along diagonal and lower off-diagonal
        S (string): The RNA sequence comprised of the letters A, U, G or C.

    Returns:
        P (np.ndarray): 2D array where entry i, j defines the maximum number of base-pairs for segment [i, j].

    """
    L = len(S)

    for k in range(1, L):  # loop over segment sizes
        for i in range(1, L - k + 1):  # loop over starting index of segment
            j = i + k
            j_unpaired = P[i, j-1]
            l_j_paired = [P[i, l-1] + P[l+1, j-1] + 1 for l in range(i, j-q) if pairs[S[l-1], S[j-1]]]

            if l_j_paired:
                P[i, j] = max(j_unpaired, *l_j_paired)
            else:
                P[i, j] = j_unpaired
    return P


def traceback(P: np.ndarray, S: str):
    """Find an base-pairs of an optimal secondary structure, i.e. a structure with maximum number of base-pairs

    Args:
        P (np.ndarray): 2D array that should be empty except for zeros along diagonal and lower off-diagonal
        S (string): The RNA sequence comprised of the letters A, U, G or C.

    Returns:
        B (list): List of base pairs in tuple format

    Remarks:
        Only returns one structure. There may be more structures with the same number of base-pairs.

    """
    L = len(S)
    sigma = [(1, L)]
    B = []

    while sigma:
        i, j = sigma.pop()
        if i >= j:  # ignore segments too small for a base-pair (produced when two neighboring sites pair or first)
            continue
        if P[i, j] == P[i, j-1]:
            sigma.append((i, j-1))
        else:
            for l in range(i, j):
                if pairs[S[l-1], S[j-1]]:
                    if P[i, j] == P[i, l-1] + P[l+1, j-1] + 1:
                        B.append((l, j))
                        sigma.extend([(i, l-1), (l+1, j-1)])
                        break
    return B


def traceback_subopt(P: np.ndarray, S: str, d: int = 0):
    """Find all suboptimal structures within a certain number of base-pairs from the maximum.

    Args:
        P (np.ndarray): 2D array that should be empty except for zeros along diagonal and lower off-diagonal
        S (string): The RNA sequence comprised of the letters A, U, G or C.
        d (int): allowed difference in number of base-pairs between optimal and suboptimal structures. Default d = 0
                     generates all possible optimal structures.

    Returns:
        B (list): List of base pairs in tuple format

    """
    # L = len(S)
    # sigma = [(1, L)]
    # B = []
    #
    # while sigma:
    #     i, j = sigma.pop()
    #     if i >= j:  # ignore segments too small for a base-pair (produced when two neighboring sites pair or first)
    #         continue
    #     if P[i, j] == P[i, j-1]:
    #         sigma.append((i, j-1))
    #     else:
    #         for l in range(i, j):
    #             if pairs[S[l-1], S[j-1]]:
    #                 if P[i, j] == P[i, l-1] + P[l+1, j-1] + 1:
    #                     B.append((l, j))
    #                     sigma.extend([(i, l-1), (l+1, j-1)])
    #                     break
    return B


def bp_to_dotbracket(bp: list, L: int) -> str:
    """Turn list of base-pairs into dot-bracket notation

    Args:
        bp (list): List of tuples defining base-pairs
        L (int): Length of sequence

    Returns:
        db (str): Secondary structure in dot-bracket notation

    """
    db = ["."] * L
    for pair in bp:
        db[pair[0] - 1] = "("
        db[pair[1] - 1] = ")"
    db = ''.join(db)
    return db

if __name__ ==  "__main__":
    S = sys.argv[1]
    P = initialize_matrix(len(S))
    P = fill_matrix(P, S, q=0)
    print(P)
    B = traceback(P, S)
    print(B)
    db = bp_to_dotbracket(B, L=len(S))
    print(db)
