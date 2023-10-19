"""Contains all functions fo different genotype to phenotype mappings
All function should take a single sequence as input and output a list
of phenotypes (can be of size one as well)

"""
import numpy as np
from typing import Callable
from abc import ABC, abstractmethod

from rna_folding.base_pairing import BasePairing
from rna_folding.nussinov import BasePairMatrixNussinov
from rna_folding.utils import bp_to_dotbracket, dotbracket_to_genotype
import RNA



def gp_mapper(input: str, output: str, mapping_function: Callable):
    """Takes file with genotypes, maps them to phenotypes and saves them in
    output file

    Args:
        input (str): Path to input file.
        output (str): Path to output file.
        mapping_function (function): A function takes a genotype (str) as 
        single positional argument and returns a list of phenotypes (str).

    Returns:
        None
    
    """
    phenotypes = {}

    # Read genotypes and map to phenotypes
    with open(input, "r") as file_in:
        for i, sequence in enumerate(file_in):
            seq = sequence.strip()
            phenotypes_ = mapping_function(seq)
            # add sequence ID to the phenotype that they map to
            for ph in phenotypes_:
                try:
                    phenotypes[ph].append(i)
                except KeyError:
                    phenotypes[ph] = [i]

    # Write to output file (line example: "{ph} {gt_id} {gt_id} {gt_id}\n"
    with open(output, "w") as file_out:
        for p in phenotypes:
            line = p + " " + " ".join(map(str, phenotypes[p])) + "\n"
            file_out.write(line)
    file_out.close()


def nussinov(genotype: str, 
             base_pairing: BasePairing, 
             min_loop_size: int, 
             suboptimal: int, 
             structures_max: int) -> list:
    """Nussinov genotype-phenotype mapping wrapper

    Args:
        genotype (str): genotype to be mapped
        base_pairing (BasePairing): An BasePairing object defining pairing ules
        min_loop_size (int): minimum size that RNA loops must have
        suboptimal (int): How many base-pairs off from optimum are allowed
        structures_max (int): How many structures to generate at most

    Returns:
        list: List of phenotypes that the genotypes maps to

    """
    P = BasePairMatrixNussinov(n=len(genotype), base_pairing=base_pairing)
    P.fill_matrix(seq=genotype, min_loop_size=min_loop_size)
    strucs = P.traceback_subopt(seq=genotype, d=suboptimal,
                                structures_max=structures_max)
    
    phenotypes = [bp_to_dotbracket(s.B, l=len(genotype)) for s in strucs]

    return phenotypes


def nussinov_mfe(genotype: str, 
                 base_pairing: BasePairing, 
                 min_loop_size: int, 
                 suboptimal: int, 
                 structures_max: int,
                 seed: int,
                 base_pair: str = "GC",
                 deterministic: bool = False) -> list:
    """Nussinov + mfe ranking genotype-phenotype mapping wrapper.
    Candidate phenotypes are generated using Nussinov's algorithm which are
    then mapped to a canonical genotype and scored using viennaRNA package

    Args:
        genotype (str): genotype to be mapped
        base_pairing (BasePairing): An BasePairing object defining pairing ules
        min_loop_size (int): minimum size that RNA loops must have
        suboptimal (int): How many base-pairs off from optimum are allowed
        structures_max (int): How many structures to generate at most
        base_pair (str): Which base-pair is used for MFE calc., either "GC" or "AU".
        seed (int): Random seed to use for generation of canonical genotype
        deterministic (bool): If True, always use G to for unpaired sites, else
        randomly pick G or C, if "GC" is given as base-pair

    Returns:
        list: List of phenotypes that the genotypes maps to
        
    """
    phenotypes = nussinov(genotype=genotype, base_pairing=base_pairing, 
                            min_loop_size=min_loop_size, suboptimal=suboptimal, 
                            structures_max=structures_max)
    
    g_fe_map = []

    for ph in phenotypes:
        # turn into a canonical alphabet
        seq_canon = dotbracket_to_genotype(dotbracket=ph,
                                            base_pair=base_pair,
                                            random=deterministic,
                                            seed=seed)
        g_fe_map.append(RNA.eval_structure_simple(seq_canon, ph))
    
    mfe_ph_id = np.argmin(g_fe_map)  # get index of mfe phenotype
    mfe_ph = phenotypes[mfe_ph_id]  # get mfe phenotype
    
    return [mfe_ph]


def viennaRNA_mfe(genotype: str) -> list:
    """Predict RNA secondary structure using default ViennaRNA mfe function.
    Args:
        genotype (str): Input genotype

    Returns:
        list: List of phenotype the genotype maps to.
        
    """
    mfe_ph, mfe = RNA.fold(genotype)

    return [mfe_ph]
