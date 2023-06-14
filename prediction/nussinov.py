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
    M = np.zeros((L, L), dtype=int)

    x, y = np.diag_indices(L)
    M[x, y] = 0    # set diagonal to zero
    M[x[1:], y[1:]-1] = 0    # set lower one-off diagonal to zero

    return M


def fill_matrix(M: np.ndarray, S: str):
    """Main step of Nussinov's algorithm, i.e. filling the M matrix to find maximum base-pairing

    Args:
        M (np.ndarray): 2D array that should be empty except for zeros along diagonal and lower off-diagonal
        S (string): The RNA sequence comprised of the letters A, U, G or C.

    Returns:
        M (np.ndarray):

    """
    L = len(S)

    for j in range(1, L):
        for i in range(j-1, -1, -1):
            T = np.zeros((L, L))
            T[i, j] = 1
            if pairs[(S[i], S[j])]:
                M[i, j] = M[i + 1, j - 1] + 1
                T[i + 1, j - 1] = 2
                print(T)

            else:
                M[i, j] = max(M[i, j - 1], M[i + 1, j - 1], M[i + 1, j])
                T[i, j - 1] = 2
                T[i + 1, j - 1] = 3
                T[i + 1, j] = 4
                print(T)



    return M


if __name__ ==  "__main__":
    S = sys.argv[1]
    print(S)
    M = initialize_matrix(len(S))
    print(M)
    print("START")
    M = fill_matrix(M, S)
    print("END")
    print(M)
