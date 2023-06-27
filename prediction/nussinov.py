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


def fill_matrix(P: np.ndarray, S: str):
    """Main step of Nussinov's algorithm, i.e. filling the P matrix to find maximum base-pairing

    Args:
        P (np.ndarray): 2D array that should be empty except for zeros along diagonal and lower off-diagonal
        S (string): The RNA sequence comprised of the letters A, U, G or C.

    Returns:
        P (np.ndarray):

    """
    L = len(S)

    for k in range(1, L):  # loop over segment sizes
        for i in range(1, L - k + 1):  # loop over starting index of segment
            j = i + k
            j_unpaired = P[i, j-1]
            l_j_paired = [P[i, l-1] + P[l+1, j-1] + 1 for l in range(i, j) if pairs[S[l-1], S[j-1]]]

            # for l in range(i, j):
            #     if pairs[S[l - 1], S[j - 1]]:
            #         t = P[i, l - 1] + P[l + 1, j - 1] + 1
            #         print(i, j, l, t, P[i, l - 1], P[l + 1, j - 1])
            #         print("REFS:", (i, l-1), (l+1, j-1))

            P[i, j] = max(j_unpaired, *l_j_paired)

    return P


if __name__ ==  "__main__":
    S = sys.argv[1]
    print(S)
    P = initialize_matrix(len(S))
    print(P)
    print("START")
    P = fill_matrix(P, S)
    print("END")
    print(P)
