import numpy as np
from itertools import product


def count_gt_per_ph_and_ph_per_gt(gp_map_file, sep=" "):
    """Open a gp_map file and count phenotypes per genotype and vice versa

    Args:
        gp_map_file (str): gp_map file in format: 
                            <ph1> <gtX> <gtY> ...
                            <ph2> <gtY> <gtX> <gtZ>...
        sep (str): Column seperator. Defaults to " ".

    Returns:
        gt_per_ph (dict): Maps genotypes to number of phenotypes they map to
        ph_per_gt (dict): Maps phenotypes to number of genotypes that map to them

    """
    gt_per_ph = {}
    ph_per_gt = {}

    with open(gp_map_file, "r") as file:
        for line_ in file:  # at start of each loop the garbage collector should delete previous line from memory
            line = line_.strip().split(sep)
            gt_per_ph[line[0]] = len(line)-1
            for gt in line[1:]:
                if gt not in ph_per_gt:
                    ph_per_gt[gt] = 1
                else:
                    ph_per_gt[gt] += 1

    return gt_per_ph, ph_per_gt