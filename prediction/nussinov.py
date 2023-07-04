"""
Implementation of Nussinov's algorithm for RNA structure prediction
Author: Leander Goldbach
GitHub: lgolbach

"""

import sys
import numpy as np


# Define all possible base-pairs. Use dictionary for quick look-up with hashing
pairs = {("A", "A"): 0, ("A", "U"): 1, ("A", "G"): 0, ("A", "C"): 0,
         ("U", "A"): 1, ("U", "U"): 0, ("U", "G"): 1, ("U", "C"): 0,
         ("G", "A"): 0, ("G", "U"): 1, ("G", "G"): 0, ("G", "C"): 1,
         ("C", "A"): 0, ("C", "U"): 0, ("C", "G"): 1, ("C", "C"): 0}


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
    """Find a base-pairs of an optimal secondary structure, i.e. a structure with maximum number of base-pairs

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
    L = len(S)
    s_init = SecondaryStructure(sigma=[(1, L)], B=[])  # initiate first suboptimal structure
    R = [s_init]  # initiate stack of suboptimal structures
    final_structures = []  # where we collect suboptimal structures
    p_max = P[1, -1]  # maximum possible number of base-pairs

    while R:
        added_to_R = False  # track whether something has been put on R stack since popping s
        s = R.pop()
        if s.is_folded():
            final_structures.append(s)
            continue
        while not s.is_folded():
            i, j = s.pop()
            s_ = SecondaryStructure(sigma=[(i, j-1), *s.sigma], B=s.B)
            if s_.maximum_bp(P) >= p_max - d:
                R.append(s_)
                added_to_R = True
            for l in range(i, j):
                if pairs[S[l - 1], S[j - 1]]:
                    s_ = SecondaryStructure(sigma=[(i, l-1), (l+1, j-1), *s.sigma], B=[*s.B, (l, j)])
                    if s_.maximum_bp(P) >= p_max - d:
                        R.append(s_)
                        added_to_R = True
        if not added_to_R:  # nothing has been put on stack since popping s
            R.append(s)  # continue with s next iteration (no infinite loop because each iteration we pop from s.sigma)
    return final_structures


class SecondaryStructure:
    """Implements data structure for backtracking an RNA secondary structure. Holds identified base-pairs and sequence
    segments that have yet to be explored by backtracking.

    """
    def __init__(self, sigma: list, B: list):
        """Initialize a stack of segments (sigma) and a list of base-pairs

        Args:
            sigma (list): Stack of segments (tuples of size 2)
            B (list): Stack of base-pairs (tuples of size 2)
        """
        self._sigma = []
        self.__add_segments(sigma)
        self._B = []
        self.__add_base_pairs(B)

    def maximum_bp(self, P: np.ndarray) -> int:
        """Determine the maximum possible number of base pairs this structure can have

        Args:
            P (np.ndarray): Matrix containing the maximum number of base-pairs possible for each segment [i, j] in P_i,j

        Returns:
            max_bps (int): Maximum number of base-pairs possible for this structure.

        """
        determined_bps = len(self._B)
        potential_bps = sum([P[seg] for seg in self._sigma])
        max_bps = determined_bps + potential_bps
        return max_bps

    @property
    def B(self):
        return self._B

    def __add_base_pairs(self, bp: list):
        self._B.extend(bp)

    @property
    def sigma(self):
        return self._sigma

    def __add_segments(self, segments: list):
        for s in segments:
            if s[0] < s[1]:  # ignore segments too small for a base-pair (produced when two neighboring sites pair or first)
                self._sigma.append(s)

    def pop(self):
        """Return segment last added to the stack

        Returns:
            (tuple): Interval (i, j), with i, j in [0, L] (L=seq. length) and i<j, which defines a sequence segment.

        """
        if not self.is_folded():
            return self._sigma.pop()

    def is_folded(self):
        if self._sigma:
            return False
        else:
            return True

    def update(self, segments: list, base_pairs: list):
        """Updates segment stack and base pairs simultaneously to avoid unsynched access to either of the two

        Args:
            segments (list): list of tuples defining segments to be added to the stack
            base_pairs (list): list of tuples defining new base pairs to be added

        Returns:
            None

        """
        self.__add_segments(segments)
        self.__add_base_pairs(base_pairs)


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
    d = int(sys.argv[2])
    P = initialize_matrix(len(S))
    P = fill_matrix(P, S, q=0)
    subopt_strucs = traceback_subopt(P=P, S=S, d=d)
    for s in subopt_strucs:
        print(bp_to_dotbracket(s.B, L=len(S)))