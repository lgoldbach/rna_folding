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
# import sys
# from pympler.asizeof import asizeof
# import os, psutil

def count_lines(file: str) -> int:
    """Count lines of file

    Args:
        file (str): Filename

    Returns:
        int: Number of lines in file
    """
    with open(file , "r") as f:
        line_count = len(f.readlines())
    f.close()
    return line_count


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputs", required=True, type=str, 
                        help="Filenames of all g-p maps", nargs="+")
    parser.add_argument("-g", "--genotype_files", required=True, type=str,
                        help="List of genotype files in order", nargs="+")
    parser.add_argument("-o", "--out", required=False, type=str,
                        help="Name of output file")

    args = parser.parse_args()

    gp_map = {}  # dictionary that will hold merged data

    # Each file starts counting genotypes at 0, so every file will have
    # genotype 0, 1, 2, ... . So we need to keep
    # track of how many genotypes have been processed already and then
    # add that number to the genotype ID.
    genotype_count = 0
    line_counts = [count_lines(fi) for fi in args.genotype_files]

    for i, filename in enumerate(args.inputs): 
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
            # print("A", asizeof(gp_map))
            file.close()
        # Counts how many genotypes were processed in this round
        genotype_count += line_counts[i]
        # process = psutil.Process()
        # print("B", process.memory_info().rss)  # in bytes

    with open(args.out, "w") as file_out:
        for p in gp_map:
            line = p + " " + " ".join(map(str, gp_map[p])) + "\n"
            file_out.write(line)
