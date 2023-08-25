import numpy as np
import networkx as nx
from rna_folding.utils import combinatorically_complete_genotypes


class GenotypePhenotypeMap(nx.Graph):
    """Storing genotype-phenotype map data as a graph and wrap 
    networkx functionalities

    """
    def __init__(self, genotypes: list = None, phenotypes: list = None, alphabet: list = None) -> None:
        """Read genotype-phenotype data (optional) and generate hamming graph

        Args:
            genotypes (list, optional): List of genotypes. Defaults to None.
            phenotypes (list, optional): List of phenotypes
                        Assumes that the ith genotype maps to the ith 
                        phenotype. Defaults to None.
            alphabet (list, optional): List of the alphabet used for the 
                        genotypes. Needed to built hamming graph. 
                        Defaults to None.

        """
        super().__init__(self)
        self.genotypes = genotypes
        self.phenotypes = phenotypes
        self.alphabet = alphabet

        if genotypes:
            self.phenotype_set = set(self.phenotypes)
            for g, p in zip(self.genotypes, self.phenotypes):
                self.add_node(g, phenotype=p)
            self.add_hamming_edges()

    @classmethod
    def read_from_file(cls, path: str, alphabet: list):
        """Read genotype-phenotype data from file

        Args:
            path (str): Path to g-p map file, assumes a csv file with one 
                        comma-separated genotype-phenotype mapping per line, 
                        e.g.:   AUGC ()()
                                CCCC ....
                                GUUC (..)
            alphabet (list): list of all letters used in the genotype space.
                             Order does not matter.

        Returns:
            GenotypePhenotypeMap: Class instance

        """
        genotypes = None
        phenotypes = None
        gpm = cls(genotypes, phenotypes, alphabet)
        return gpm

    @classmethod
    def read_from_dict(cls, dict: dict, alphabet: list):
        """Read genotype-phenotype data from file

        Args:
            dict (dict): Dictionary of genotype-phenotype mapping, 
                         e.g.: { "AUGC": "()()",
                                 "CCCC": "....",
                                 "GUUC": "(..)" }

            alphabet (list): list of all letters used in the genotype space.
                             Order does not matter.

        Returns:
            GenotypePhenotypeMap: Class instance
            
        """
        genotypes = list(dict.keys())
        phenotypes = list(dict.values())
        gpm = cls(genotypes, phenotypes, alphabet)
        return gpm

    def add_hamming_edges(self):
        """Compute all hamming edges for each genotype, i.e. add an edge 
        between a genotype and all genotypes that differ by one letter.
        
        Note:
            Current implementation assumes combinatorically complete g-p map
        """
        for g in self.genotypes:
            for site, l in enumerate(g):
                for m in self.alphabet:
                    if m != l:
                        neighbor = g[:site] + m + g[site + 1:]
                        self.add_edge(g, neighbor)
    
    def neutral_components(self, phenotypes: list = []) -> list:
        """Compute all neutral components for given phenotypes. A neutral 
        component is defined as a connected set of nodes that all map to the
        same phenotype. A phenotype can have between one and #(phenotype) 
        neutral components. 

        Args:
            phenotypes (list, optional): List of phenotypes for which neutral
            components will be returned. If none are given, neutral components 
            for all phenotypes will be returned. Defaults to [].

        Returns:
            list: list of sets. The ith element in the list contains the 
            neutral components (sets) of the ith phenotype

        """
        if not phenotypes:
            phenotypes = self.phenotype_set
        
        neutral_components = []
        for ph in phenotypes:
            nodes = [node for node, attr in self.nodes(data=True) if 
                     attr['phenotype'] == ph]
            G_sub = self.subgraph(nodes)
            cc = nx.connected_components(G_sub)
            neutral_components.append(cc)
        return neutral_components

    def phenotype_robustness(self, nodes: list) -> float:
        """Compute robustness of given set of nodes. Defined as fraction of
        neighboring nodes with a different phenotype averaged over the whole
        set of nodes provided.
        Note that this method does not check whether all provided nodes 
        actually have the same phenotype or whether they are in the same
        connected component but simply checks whether nodes have the same 
        phenotypes as their neighbor.

        Args:
            nodes (list): List of nodes (str) to consides, e.g. ["AA", "AU"]

        Returns:
            float: Phenotype robustness.

        """
        fractions_of_identical_neighbors = []
        for node in nodes:
            total_neighb = 0
            same_ph = 0
            ref_ph = self.nodes[node]["phenotype"]
            for neighbor in self.neighbors(node):
                total_neighb += 1
                if self.nodes[neighbor]["phenotype"] == ref_ph:
                    same_ph += 1
            fractions_of_identical_neighbors.append(same_ph / total_neighb)

        robustness = np.mean(fractions_of_identical_neighbors)
        return robustness
        