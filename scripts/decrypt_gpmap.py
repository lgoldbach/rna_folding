"""Decrypt the sgreenbury format for storing g-p map info into a single file
of the format (see: https://github.com/sgreenbury/gp-maps-nav):
<seq> <secondary structure>

e.g.:
AUGC ()()
CCCC ....
...

"""


import argparse
import numpy as np

def rna_numeric_to_letter(seq_num):
    code = {0: 'A', 1: 'U', 2: 'C', 3: 'G'}

    seq = []
    for i in seq_num:
        seq.append(code[int(i)])
        
    return ''.join(seq)


def base_convert(n, b):
    result = np.base_repr(n, b, padding=12)[::-1][:12]
    return result
    

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gpmap", required=True, type=str, help="List of genotypes")
    parser.add_argument("-p", "--phenotypes", required=True, type=str, help="List of phenotypes and genotype IDs that map to them")
    parser.add_argument("-o", "--output", required=True, type=str, help="List of phenotypes and genotype IDs that map to them")


    args = parser.parse_args()
    ph = np.loadtxt(args.phenotypes, dtype=str)

    seq_ss_map = []
    with open(args.gpmap, "r") as gp:
        for i, line in enumerate(gp, start=1):
            b4 = base_convert(i - 1, 4)

            seq = rna_numeric_to_letter(b4)
            
            seq_ss_map.append(seq + " " + ph[int(line)] + "\n")

    with open(args.output, "w") as outfile:
         for line in seq_ss_map:
              outfile.write(line)
              