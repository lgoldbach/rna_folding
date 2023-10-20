#!/usr/bin/env python

"""Merge the output files of a genotype-phenotype map that was generated in
separate jobs. 
Assumes that all files are in same directory and follow the following format:

<ph> <gt_id> <gt_id> ...
<ph> <gt_id> <gt_id> ...

e.g.
(..) 1 3 5
()() 4 2

Assumes that each file starts countin genotypes from 1 which will be corrected
by the script by adding the number of already processed genotypes to each 
genotype ID.

"""

import argparse


def count_lines(file: str) -> int:
    """Count lines of file

    Args:
        file (str): Filename

    Returns:
        int: Number of lines in file
    """
    with open(file , "r") as f:
        line_count = len(f.readlines())
    return line_count


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--filename", required=True, type=str, 
                        help="Filename of all g-p maps. Has to be identical" 
                        "except for number")
    parser.add_argument("-r", "--range", required=True, type=str, 
                        help="Number range (closed interval) for filename"
                        "suffix. E.g: 1-50. No leading zeros")
    parser.add_argument("-g", "--genotype_reference", required=True, type=str, 
                        help="Name of corresponding genotypes files (needed to"
                        "get correct genotype ids)")
    parser.add_argument("-o", "--out", required=False, type=str, 
                        help="Name of output file")

    args = parser.parse_args()

    # turn range argument string into tuple
    suffix_range = list((int(i) for i in args.range.split("-")))

    gp_map = {}  # dictionary that will hold merged data

    # Each file starts counting genotypes from one. So we need to keep
    # track of how many genotypes have been processed already and then
    # add that number to the genotype ID.
    genotype_count = 0

    for i in range(suffix_range[0], suffix_range[1]+1):  # +1->closed interval
        filename = args.filename + str(i)  # make filename
        with open(filename , "r") as file:
            for line in file:
                line_split = line.strip().split(" ")
                phenotype = line_split[0]
                genotypes = line_split[1:]
                
                for gt in genotypes:
                    gt_id = int(gt) + genotype_count  # correct id by count
                    if phenotype in gp_map:
                        gp_map[phenotype].append(gt_id)
                    else:
                        gp_map[phenotype] = [gt_id]

        # Counts how many genotypes were processed in this round
        genotype_count += count_lines(args.genotype_reference + str(i))

    with open(args.out, "w") as file_out:
        for p in gp_map:
            line = p + " " + " ".join(map(str, gp_map[p])) + "\n"
            file_out.write(line)
