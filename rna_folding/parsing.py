import argparse
import numpy as np


def gpmap_to_dict(gpmap_file: str, genotype_file: str = None) -> dict:
    """Takes a file that stores genotype-phenotype mapping and a list of 
    genotypes and parses it into dictionary 
    {<genotype> (str): [phenotypes (str)]}. Intended for mappings to multiple 
    phenotypes.    

    Args:
        gpmap_file (str): Path to file in following format:
                            <phenotype> <genotypeID_x> <genotypeID_y>
                            <phenotype> <genotypeID_y>
                            ...
                            Example for RNA secondary structure (dot-bracket):
                            ((..)) 2 1
                            (....) 3
                            ().... 4 1
                            ...

        genotype_file (str): file path to a file that contains simple list of 
                                genotypes (one per line)

    Returns:
        gpmap (dict): Dictionary that maps genotype (str) to list of one or 
                        more phenotypes (str).

    """
    # read in genotypes as list
    if genotype_file:
        with open(genotype_file, "r") as g_file:    
            genotype_list = [line.strip() for line in g_file]
         
    gp_map = {}
    with open(gpmap_file, "r") as gp_file:
        for line in gp_file:
            l = line.split()
            db = l[0]
            for gt in l[1:]:
                if genotype_file:
                    gt = genotype_list[int(gt)]  # get genotype using id (i)
                if gt in gp_map:
                    gp_map[gt].append(db)
                else:
                    gp_map[gt] = [db]

    return gp_map


def viennarna_to_gp_map_file(viennarna_output: str) -> dict:
    """Takes an output file and parses it into dictionary 
    {<genotype> (str): [phenotypes (str)]}. Intended for mappings to multiple 
    phenotypes.    

    Args:
        viennarna_output (str): path to file with format:
                                <seq1>
                                <dot-bracket> <free energy>
                                <seq2>
                                <dot-bracket> <free energy>
                                ...
                                e.g:
                                GGGAAACCC
                                (((...))) (-1.20)
                                AAAAAAAAA
                                ......... (0.00)
    Returns:
        gpmap (dict): Dictionary that maps genotype (str) to list of one or 
                        more phenotypes (str).

    """
    gp_map = {}
    with open(viennarna_output, "r") as file:
        for line in file:
            gt = line.strip()
            ph = next(file).split(" ")[0]
            gp_map[gt] = ph

    return gp_map

        