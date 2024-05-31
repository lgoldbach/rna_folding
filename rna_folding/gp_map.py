import numpy as np
import networkx as nx
from rna_folding.utils import combinatorically_complete_genotypes


class GenotypePhenotypeGraph(nx.Graph):
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
        self.genotypes = np.array(genotypes)
        self.phenotypes = np.array(phenotypes)
        self.alphabet = alphabet

        if genotypes:
            self.phenotype_set = list(set(self.phenotypes))
            for i, (g, p) in enumerate(zip(self.genotypes, self.phenotypes)):
                self.add_node(g, phenotype=p, id=i)

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
        genotypes = []
        phenotypes = []

        with open(path, "r") as file:
            for line in file:
                genotype, phenotype = line.strip().split(' ')
                genotypes.append(genotype)
                phenotypes.append(phenotype)

        gpm = cls(genotypes, phenotypes, alphabet)
        return gpm

    @classmethod
    def read_from_ph_to_gt_file(cls, path: str, genotype_ref_path: str, 
                                alphabet: list):
        """Read genotype-phenotype data from file

        Args:
            path (str): Path to g-p map file, assumes a csv file where each 
                        phenotype is followed by all the genotypes that map to
                        it, separated by spaces. Assumes one-to-one mapping.
                        e.g.:   ()() 10 4 2 1
                                .... 3 5 6
                                (..) 7 8 9
            genotype_ref_path (str): Path to genotype reference where the genotype 
                        on the ith line corresponds to the numbering of the 
                        gp_map file (path).
            alphabet (list): list of all letters used in the genotype space.
                             Order does not matter.
        Returns:
            GenotypePhenotypeMap: Class instance

        """
        genotypes = []
        phenotypes = []

        with open(genotype_ref_path, "r") as gt_file:
            genotype_ref = [line.strip() for line in gt_file]
            

        with open(path, "r") as file:
            for line in file:
                line_ = line.strip().split(' ')
                phenotype = line_[0]
                for genotype in line_[1:]:
                    genotypes.append(genotype_ref[int(genotype)])  # num to seq
                    phenotypes.append(phenotype)

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

    def nodes_with_phenotype(self, phenotype: str) -> list:
        """Get all nodes with a given phenotype

        Args:
            phenotype (str): Phenotype string

        Returns:
            list: list of node names
        
        """
        nodes = self.genotypes[np.where(self.phenotypes == phenotype)[0]]
        return nodes

    def add_hamming_edges(self):
        """Compute all hamming edges for each genotype, i.e. add an edge 
        between a genotype and all genotypes that differ by one letter.
        
        Note:
            Current implementation assumes combinatorically complete g-p map
        """
        for g in self.genotypes:
            for neighbor in self._neighbors(g):
                self.add_edge(g, neighbor)
    
    def _neighbors(self, node):
        """Get all neighbors of a node. Either looks them up from existing
        edges or generates them one the fly. Neighbors are defined as nodes
        whose genotype differ in one positon.

        Args:
            node (str): Node descriptor (genotype), e.g. "AGCA"

        Returns:
            list: list of neighboring nodes (str)

        """
        neighbors = []
        for site, wt_l in enumerate(node):
            for l in self.alphabet:
                if l != wt_l:
                    neighbors.append(node[:site] + l + node[site + 1:])
        return neighbors

    def neutral_components(self, phenotypes: list = [], return_ids=False) -> list:
        """Compute all neutral components for given phenotypes. A neutral 
        component is defined as a connected set of nodes that all map to the
        same phenotype. A phenotype can have between one and #(phenotype) 
        neutral components. 

        Args:
            phenotypes (list, optional): List of phenotypes for which neutral
            components will be returned. If none are given, neutral components 
            for all phenotypes will be returned. Defaults to [].
            return_ids (bool): If true, genotypes will be returnes as numerical
                ids and not as full sequence. This saves memory for large maps.

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
            # translate the components from full sequences to numeric id
            # for memory efficiency
            if return_ids:
                final_cc = []
                for c in cc:
                    final_cc.append({self.nodes[node]["id"] for node in c})
            else:
                final_cc = cc
                
            neutral_components.append(final_cc)
        return neutral_components
    
    def neutral_paths(self, source_genotype: str, n: int) -> list:
        """Start n random walks from <genotype> along genotypes with the same
        phenotype

        Args:
            genotype (str): Starting genotype for random walks
            n (int): Number of random walks to take

        Returns:
            list: List of lists of random walks.

        """
        phenotype = self.nodes[source_genotype]["phenotype"]  # get ref phenotype
        neutral_paths = []
        
        for i in range(n):            
            # create mutation generator for each iteration
            mutation_gen = self.mutation_space(source_genotype)  
            neutral_path = [source_genotype]  # init list
            
            genotype = source_genotype
            for mutation in mutation_gen:  # loop through possible mutations
                new_genotype = self.mutate(genotype, mutation)
                # only accept neutral mutants
                if self.nodes[new_genotype]["phenotype"] == phenotype:
                    neutral_path.append(new_genotype)
                    genotype = new_genotype
                else:
                    continue
            neutral_paths.append(neutral_path)

        return neutral_paths
    
    def mutation_space(self, genotype: str):
        """Generate new mutations for genotype. Each mutation can only appear
        once.

        Args:
            genotype (str): A genotype.

        Yields:
            mutation (tuple): (<site, int>, <letter, str>), The first element
                            defines the site of the mutation and the second
                            the letter.

        """
        # create a space of all possible mutations for the genotype
        # one list per site. Contains all mutations except the wildtype.
        mutation_space = []
        for site in genotype:
            mutation_space.append([])
            for letter in self.alphabet:
                if letter != site:
                    mutation_space[-1].append(letter)
        
        # as long as there are possible mutations left, yield a random one
        while any(mutation_space):
            site = np.random.choice(len(genotype))
            # this could be done more sophisticated. I am brute-forcing
            # random.choice until I find a site that still has available mut.
            try:  
                mut_i = np.random.choice(len(mutation_space[site]))
            except ValueError:
                continue
            mut = mutation_space[site].pop(mut_i)
            yield((site, mut))
    
    @staticmethod
    def mutate(genotype: str, m: tuple) -> str:
        """Take a genotype and mutate is

        Args:
            genotype (str): A genotype, e.g. ACGUA
            m (tuple): A tuple (<site, int>, <letter, str>) defining
                                the mutation.

        Returns:
            mutant (str): The mutant genotype

        """
        new_genotype = genotype[:m[0]] + m[1] + genotype[m[0]+1:]
        return new_genotype

    def phenotype_robustness(self, nodes: list) -> float:
        """Compute robustness of given set of nodes. Defined as fraction of
        neighboring nodes with a different phenotype averaged over the whole
        set of nodes provided.
        Note that this method does not check whether all provided nodes 
        actually have the same phenotype or whether they are in the same
        connected component but simply checks whether nodes have the same 
        phenotypes as their neighbor.

        Args:
            nodes (list): List of nodes (str) to consider, e.g. ["AA", "AU"]

        Returns:
            float: Phenotype robustness.

        """
        fractions_of_identical_neighbors = []
        for node in nodes:
            total_neighb = 0
            same_ph = 0
            ref_ph = self.nodes[node]["phenotype"]
            for neighbor in self._neighbors(node):
                total_neighb += 1
                if self.nodes[neighbor]["phenotype"] == ref_ph:
                    same_ph += 1
            fractions_of_identical_neighbors.append(same_ph / total_neighb)

        robustness = np.mean(fractions_of_identical_neighbors)
        return robustness
        