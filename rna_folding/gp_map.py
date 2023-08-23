import numpy as np
import networkx as nx


class GenotypePhenotypeMap:
    """Storing genotype-phenotype map data as a graph and wrap 
    networkx functionalities

    """
    def __init__(self, genotypes, phenotypes) -> None:
        # build networkx
        pass
        
    @classmethod
    def read_from_file(cls, path: str):
        genotypes = None
        phenotypes = None
        gpm = cls(genotypes, phenotypes)
        return gpm

    @classmethod
    def read_from_dict(cls, dict: dict):
        genotypes = None
        phenotypes = None
        gpm = cls(genotypes, phenotypes)
        
        