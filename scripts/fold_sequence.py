#!/usr/bin/env python

import argparse

from rna_folding.nussinov import BasePairMatrixNussinov
from rna_folding.utils import bp_to_dotbracket
from rna_folding.base_pairing import BasePairing


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--sequence", required=True, type=str, help="RNA sequence")
    parser.add_argument("-m", "--min_loop_size", required=True, type=int, default=1, help="Minimum size for loop")
    parser.add_argument("-s", "--suboptimal", type=int, required=True,
                        help="Create all suboptimal structures with number of base-pairs in the range"
                             "of max - s, where s is an integer. Without the flag, only one structure "
                             "is computed.")
    parser.add_argument("-z", "--structures_max", required=False, type=int, help="Limit on how many suboptimal structures to generate")

    args = parser.parse_args()

    pairing = BasePairing("AUGC", -1)
    P = BasePairMatrixNussinov(n=len(args.sequence), base_pairing=pairing)
    P.fill_matrix(seq=args.sequence, min_loop_size=args.min_loop_size)
    strucs = P.traceback_subopt(seq=args.sequence, d=args.suboptimal, structures_max=args.structures_max)
    for s in strucs:
        db = bp_to_dotbracket(s.B, l=len(args.sequence))
        print(db)

    #### DEBUG CODE START #####
   # S = "UUG"
   # mls = 1
   # sub = 0
   # P = BasePairMatrixNussinov(n=len(S))
   # P.fill_matrix(seq=S, min_loop_size=mls)
   #  strucs = P.traceback_subopt(seq=S, d=sub)
   # 
   # for s in strucs:
   #     db = bp_to_dotbracket(s.B, l=len(S))
   #     print(db)
    #### DEBUG CODE END #####

