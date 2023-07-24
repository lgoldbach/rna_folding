import argparse

from rna_folding.nussinov import BasePairMatrixNussinov
from rna_folding.utils import bp_to_dotbracket


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="File input for genotypes")
    parser.add_argument("-o", "--output", help="File output for phenotypes")
    parser.add_argument("-m", "--min_loop_size", required=True, type=int, default=1, help="Minimum size for loop")
    parser.add_argument("-s", "--suboptimal", type=int, required=True,
                        help="Create all suboptimal structures with number of base-pairs in the range"
                             "of max - s, where s is an integer. Without the flag, only one structure "
                             "is computed.")
    parser.add_argument("-z", "--structures_max", required=False, type=int, help="Limit on how many suboptimal structures to generate")

    args = parser.parse_args()

    phenotypes = {}
    with open(args.file, "r") as file_in:
        for i, sequence in enumerate(file_in):
            seq = sequence.strip()
            P = BasePairMatrixNussinov(n=len(seq))
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