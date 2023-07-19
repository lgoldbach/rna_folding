import argparse
from rna_folding.evaluate import f1_score


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-s1", "--structure1", required=True, type=str, help="Reference RNA secondary structrue in dot-bracket notation")
    parser.add_argument("-s2", "--structure2", required=True, type=str, help="Query RNA secondary structrue in dot-bracket notation")
    parser.add_argument("-a", "--abstract", required=False, type=int, help="Level of abstraction on which to compare structures")

    args = parser.parse_args()

    f1 = f1_score(A=args.structure1, B=args.structure2, 
                  abstract_level=args.abstract)
    
    print(f1)
    