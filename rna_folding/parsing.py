import argparse
import numpy as np

def gpmap_pgdict(gpmap_file: str, genotype_file: str = None) -> dict:
    """Takes a file that stores genotype-phenotype mapping and a list of 
    genotypes and parses it into dictionary {<phe> (str): [genotypes (str)]}. 
    Intended for many-to-many mappings

    Args:
        gpmap_file (str):       Path to file in following format:
                                <phenotype> <genotypeID_x> <genotypeID_y>
                                <phenotype> <genotypeID_y>
                                ...
                                Example for RNA secondary structure:
                                ((..)) 2 1
                                (....) 3
                                ().... 4 1
                                ...

        genotype_file (str):    file path to a file that contains simple list of 
                                genotypes (one per line)

    Returns:
        pgmap (dict):           Dictionary that maps phenotype (str) to list of
                                one or more genotypes (str).

    """
    if genotype_file:
        with open(genotype_file, "r") as g_file:    
            genotype_list = [line.strip() for line in g_file]
         
    pg_map = {}
    with open(gpmap_file, "r") as gp_file:
        for line in gp_file:
            l = line.strip().split(" ")
            pg_map[l[0]] = [genotype_list[int(gt)] for gt in l[1:]]
    
    return pg_map


def gpmap_to_lists(gpmap_file: str) -> tuple:
    """Takes a gp map text file and returns a tuple which contains a list of
    genotypes and a list of phenotypes where the ith genotype of the first list
    maps to the iths phenotype of the second list

    Args:
        gpmap_file (str):       Path to file in following format:
                                <phenotype> <genotypeID_x> <genotypeID_y>
                                <phenotype> <genotypeID_y>
                                ...
                                Example for RNA secondary structure:
                                ((..)) 2 1
                                (....) 3
                                ().... 4 1
                                ...

    Returns:
        tuple (list, list):     One list with genotypes and another phenotypes
                                Both same length

    """
    genotypes = []
    phenotypes = []
    with open(gpmap_file, "r") as gp_file:
        for line in gp_file:
            l = line.split()
            ph = l[0]
            for gt in l[1:]:
                genotypes.append(gt)
                phenotypes.append(ph)
    
    return genotypes, phenotypes


def lists_to_gp_map(genotypes, phenotypes, output_filename) -> None:
    """Take a list of genotypes and phenotypes and saves them in standard
    gp map format.

    Args:
        genotypes (list): list of genotype strings
        phenotypes (list): list of phenotypes strings
        output_filename (str):    Path to file in following format:
                                    <phenotype> <genotypeID_x> <genotypeID_y>
                                    <phenotype> <genotypeID_y>
                                    ...
                                    Example for RNA secondary structure:
                                    ((..)) 2 1
                                    (....) 3
                                    ().... 4 1
                                    ...

    """
    




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


def dict_to_gpmap(ph_to_gt: dict, file: str) -> None:
    """Take a dict that maps phenotype to list of genotypes and save it
    as a space-separated "c"sv file, where each line looks like this:
    "{ph} {gt_id} {gt_id} {gt_id}"

    Args:
        ph_to_gt (dict): _description_
        file (str): _description_
    """
    # Write to output file (
    with open(file, "w") as file_out:
        for p in ph_to_gt:
            line = p + " " + " ".join(map(str, ph_to_gt[p])) + "\n"
            file_out.write(line)
    file_out.close()


def load_phenotype_and_metric_from_file(file: str, dtype=float):
    """Take a file in the common phenotype (col1) metric (col2) data-type 
    I am using and reat it as two array.
    Example file:
    ((...)) 0.8
    (.....) 0.7
    ...

    Args:
        file (str): Path to the file

    Retruns:
        phentypes, data
    """
    file_data = np.loadtxt(file, dtype=str)
    if file_data.ndim == 1:  # in case there is only one phenotype
        file_data = np.expand_dims(file_data, axis=0)
    phenotypes = file_data[:,0]
    distr = file_data[:,1].astype(dtype)

    return phenotypes, distr
