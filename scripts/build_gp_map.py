#!/usr/bin/env python

import argparse

from rna_folding.nussinov import BasePairMatrixNussinov
from rna_folding.utils import bp_to_dotbracket
from rna_folding.base_pairing import BasePairing


"""
This is a dirty hack but the graph data is too big to be hosted on github, so
you need to hard-code the full path to the graph data here. If this ever
becomes a real package I will find a solution to this but for now it doesn't 
matter
"""
graph_path = "/home/lgold/phd/research/projects/connectivity/rna_folding/data/graphs/"


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input file with genotypes")
    parser.add_argument("-o", "--output", help="File output for phenotypes")
    parser.add_argument("-m", "--min_loop_size", required=True, type=int, default=1, help="Minimum size for loop")
    parser.add_argument("-s", "--suboptimal", type=int, required=True,
                        help="Create all suboptimal structures with number of base-pairs in the range"
                             "of max - s, where s is an integer. Without the flag, only one structure "
                             "is computed.")
    parser.add_argument("-z", "--structures_max", required=False, type=int, help="Limit on how many suboptimal structures to generate")
    parser.add_argument("-p", "--base_pairing", required=False, type=int, default=-1, help="Which base-pairing to choose. I.e. from the base-pairing simple graphs, which one to pick "
                        "e.g. for 4 bases there are 11 possible base-pairings, so possible input is any number between 1 and 11, If given -1 then it uses canonical base-pairing and AUGC bases")
    parser.add_argument("-b", "--bases", required=False, type=str, default="AUGC", help="Which bases do the genotypes contain, e.g. 'AUGC' for canonical RNA")
    parser.add_argument("-g", "--graph_path", required=True, type=str, 
                        default=graph_path, 
                        help="Path to folder containing the base-pairing "
                        "graphs files, e.g. graph4.adj. Check base_pairing.py "
                        "for info on where these graphs come from.")

    args = parser.parse_args()

    pairing = BasePairing(bases=args.bases,
                          graph_path=graph_path, 
                          id=args.base_pairing)
    phenotypes = {}
    with open(args.file, "r") as file_in:
        for i, sequence in enumerate(file_in):
            seq = sequence.strip()
            P = BasePairMatrixNussinov(n=len(seq), base_pairing=pairing)
            P.fill_matrix(seq=seq, min_loop_size=args.min_loop_size)
            strucs = P.traceback_subopt(seq=seq, d=args.suboptimal, structures_max=args.structures_max)
            for s in strucs:
                db = bp_to_dotbracket(s.B, l=len(seq))
                try:
                    phenotypes[db].append(i)
                except KeyError:
                    phenotypes[db] = [i]

    with open(args.output, "w") as file_out:
        for p in phenotypes:
            line = p + " " + " ".join(map(str, phenotypes[p])) + "\n"
            file_out.write(line)
    file_out.close()