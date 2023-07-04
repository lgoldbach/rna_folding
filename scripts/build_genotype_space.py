import sys
import argparse
import numpy as np
from itertools import product


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


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", help="Output file", type=str)
    parser.add_argument("-l", "--seq_len", help="Sequence length", type=int)
    parser.add_argument("-a", "--alphabet", help="Define the alphabet as an unseparated string, e.g: AUGC", type=str)

    args = parser.parse_args()
    genotypes = combinatorically_complete_genotypes(args.seq_len, args.alphabet)

    with open(args.output, "w") as f:
        for g in genotypes:
            f.write("".join(g + ("\n",)))
    f.close()
