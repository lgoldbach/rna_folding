"""
Nussinov-Jacobson python algorithm implementation
	Predicts the secondary RNA structure from an RNA sequence.
	The minimal loop length is set to a default of 0
"""

import numpy as np
import pandas as pd  # you can remove this... it's just for nice output (matrix formatting)


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

def fill(nm, rna):
    """
    Fill the matrix as per the Nussinov algorithm
    """
    minimal_loop_length = 0

    for k in range(1, len(rna)):
        for i in range(len(rna) - k):
            j = i + k
            print(i, j)
            if j - i >= minimal_loop_length:
                down = nm[i + 1][j]  # 1st rule
                left = nm[i][j - 1]  # 2nd rule
                diag = nm[i + 1][j - 1] + pairs[(rna[i], rna[j])]  # 3rd rule

                try:
                    rc = max([nm[i][t] + nm[t + 1][j] for t in range(i+1, j)])  # 4th rule
                except ValueError:
                    rc = 0

                nm[i][j] = max(down, left, diag, rc)  # max of all
            else:
                nm[i][j] = 0

    return nm


def traceback(nm, rna, fold, i, L):
    """
    Traceback through complete Nussinov matrix to find optimial RNA secondary structure solution through max base-pairs
    """
    j = L

    if i < j:
        if nm[i][j] == nm[i + 1][j]:  # 1st rule
            traceback(nm, rna, fold, i + 1, j)
        elif nm[i][j] == nm[i][j - 1]:  # 2nd rule
            traceback(nm, rna, fold, i, j - 1)
        elif nm[i][j] == nm[i + 1][j - 1] + pairs[(rna[i], rna[j])]:  # 3rd rule
            fold.append((i, j))
            traceback(nm, rna, fold, i + 1, j - 1)
        else:
            for k in range(i + 1, j - 1):
                if nm[i][j] == nm[i, k] + nm[k + 1][j]:  # 4th rule
                    traceback(nm, rna, fold, i, k)
                    traceback(nm, rna, fold, k + 1, j)
                    break

    return fold


def dot_write(rna, fold):
    dot = ["." for i in range(len(rna))]

    for s in fold:
        # print(min(s), max(s))
        dot[min(s)] = "("
        dot[max(s)] = ")"

    return "".join(dot)


def init_matrix(rna):
    M = len(rna)

    # init matrix
    nm = np.empty([M, M])
    nm[:] = np.NAN

    # init diaganols to 0
    # few ways to do this: np.fill_diaganol(), np.diag(), nested loop, ...
    nm[range(M), range(M)] = 0
    nm[range(1, len(rna)), range(len(rna) - 1)] = 0

    return nm


if __name__ == "__main__":
    rna = "AUGC"

    nm = init_matrix(rna)
    nm = fill(nm, rna)

    fold = []
    sec = traceback(nm, rna, fold, 0, len(rna) - 1)

    res = dot_write(rna, fold)

    names = [_ for _ in rna]
    df = pd.DataFrame(nm, index=names, columns=names)
    print(df, "\n", rna, res)