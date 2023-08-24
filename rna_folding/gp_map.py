import numpy as np
import networkx as nx
from rna_folding.utils import combinatorically_complete_genotypes


class GenotypePhenotypeMap(nx.Graph):
    """Storing genotype-phenotype map data as a graph and wrap 
    networkx functionalities

    """


    def __init__(self, genotypes: list = None, phenotypes: list = None, alphabet: list = None) -> None:
        # build networkx
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
        genotypes = None
        phenotypes = None
        gpm = cls(genotypes, phenotypes, alphabet)
        return gpm

    @classmethod
    def read_from_dict(cls, dict: dict, alphabet: list):
        genotypes = list(dict.keys())
        phenotypes = list(dict.values())
        gpm = cls(genotypes, phenotypes, alphabet)
        return gpm

    def add_hamming_edges(self):
        for g in self.genotypes:
            for site, l in enumerate(g):
                for m in self.alphabet:
                    if m != l:
                        neighbor = g[:site] + m + g[site + 1:]
                        self.add_edge(g, neighbor)
    
    def neutral_components(self, phenotypes: list = []) -> list:
        if not phenotypes:
            phenotypes = self.phenotype_set
        
        neutral_components = []
        for ph in phenotypes:
            print(ph)
            nodes = [node for node, attr in self.nodes(data=True) if 
                     attr['phenotype'] == ph]
            print(nodes)
            G_sub = self.subgraph(nodes)
            cc = nx.connected_components(G_sub)
            neutral_components.append(cc)
        return neutral_components


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    l = 3  # sequence length
    alphabet = "AU"
    np.random.seed(12)
    colors = ["blue", "red", "green", "yellow", "pink", "cyan", "purple"]

    genotypes_tuple = combinatorically_complete_genotypes(l=l, a=alphabet)
    genotypes = ["".join(g) for g in genotypes_tuple]
    phenotypes = []

    # make up some nonsense phenotypes
    num_of_ph = 3
    phenotype_set = []
    ph_to_color = {}    
    while len(phenotype_set) < num_of_ph:
        ph = "".join(np.random.choice(["(", ".", ")"], l))
        if ph not in phenotype_set:
            ph_to_color[ph] = colors[len(phenotype_set)]
            phenotype_set.append(ph)

    p_distr = []
    p = 1

    for i in range(len(phenotype_set)):
        p = p/2
        p_distr.append(p)



    p_sum = sum(p_distr)    
    p_distr_norm = [p/p_sum for p in p_distr]
    

    phenotypes = np.random.choice(phenotype_set, len(genotypes),
                                  p=p_distr_norm)
    print(list(zip(genotypes, phenotypes)))
    gp_dict = dict(zip(genotypes, phenotypes))
    GP = GenotypePhenotypeMap.read_from_dict(gp_dict, alphabet)

    c_map = []
    for node in GP.nodes:
        c_map.append(ph_to_color[GP.nodes[node]["phenotype"]])

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    labels = nx.get_node_attributes(GP, "phenotype")
    print("X", labels)
    nx.draw(GP, with_labels=True, ax=ax1, pos=nx.circular_layout(GP), node_color=c_map)
    nx.draw(GP, labels=labels, ax=ax2, pos=nx.circular_layout(GP), node_color=c_map)
    plt.savefig("test.png")

    cc = GP.neutral_components()
    print([list(c) for c in cc])