import sys
import argparse

from rna_folding.nussinov import BasePairMatrixNussinov
from rna_folding.utils import bp_to_dotbracket


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--sequence", required=True, type=str, help="RNA sequence")
    parser.add_argument("-m", "--min_loop_size", required=True, type=int, default=1, help="Minimum size for loop")
    parser.add_argument("-s", "--suboptimal", type=int, required=True,
                        help="Create all suboptimal structures with number of base-pairs in the range"
                             "of max - s, where s is an integer. Without the flag, only one structure "
                             "is computed.")

    args = parser.parse_args()

    P = BasePairMatrixNussinov(n=len(args.sequence))
    P.fill_matrix(seq=args.sequence, min_loop_size=args.min_loop_size)
    strucs = P.traceback_subopt(seq=args.sequence, d=args.suboptimal)
    for s in strucs:
        db = bp_to_dotbracket(s.B, l=len(args.sequence))
        print(db)