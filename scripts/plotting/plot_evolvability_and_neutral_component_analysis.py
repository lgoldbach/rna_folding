#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

from rna_folding.parsing import many_to_one_map_from_file_to_dict


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--nc_graph", help="Neutral component graph", required=True)
    parser.add_argument("-m", "--nc_to_gt", help="Map from neutral components to genotypes", required=True)
    parser.add_argument("-g", "--gp_map", help="whole gp map", required=True)
    parser.add_argument("-u", "--unviable", help="Phenotype to ignore", type=str, required=False)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()

    # Dictionary to map genotypes to neutral component id
    gt_to_nc = many_to_one_map_from_file_to_dict(file=args.nc_to_gt, 
                                                 source_type=str,
                                                 target_type=int,
                                                 delimiter=" ")
    
    nc_graph = pickle.load(open(args.nc_graph, "rb"))
    gp_map = pickle.load(open(args.gp_map, "rb"))
    
    viable_non_neutral_fraction = []
    viable_fraction = []
    non_neutral_fraction_of_viable = []
    for gt in gt_to_nc:  # loop over gt_to_nc because it should only contain genotypes that map to viable phenotypes
        ph_of_gt = nc_graph.nodes[gt_to_nc[gt]]["phenotype"]  # get phenotype of gt
        neighbors = gp_map._neighbors(gt)  # get its neighbor's genotypes

        count_viable = 0
        count_viable_non_neutral = 0
        for neighbor in neighbors:
            ph = gp_map.map(neighbor)  # get each neighbors phenotype
            # if neighbor phenotype is viable and not the same as that of gt
            if ph != args.unviable:
                count_viable += 1
                if ph != ph_of_gt:
                    count_viable_non_neutral += 1  # count the neighbor that is a different and viable phenotype

        
        viable_non_neutral_fraction.append(count_viable_non_neutral)
        viable_fraction.append(count_viable)
        try:
            non_neutral_fraction_of_viable.append(count_viable_non_neutral/count_viable)
        except ZeroDivisionError:
            pass


    # compute the outdegree in this gp map
    neighbor_count = (len(gp_map.alphabet)-1) * len(gp_map.genotypes[0])

    # Normalize count of different and viable phenotype neighbors by outdegree
    viable_non_neutral_fraction = np.array(viable_non_neutral_fraction)/neighbor_count
    viable_fraction = np.array(viable_fraction)/neighbor_count

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5))

    n, bins, p = axes[0].hist(viable_fraction)
    axes[0].set_xlabel("Fraction of neighbors with viable phenotype")
    axes[0].set_ylabel(f"Genotype count\nSum: {len(viable_fraction)}")
    axes[0].vlines([np.mean(viable_fraction)], [0], [max(n)], label="mean", color="black")

    n, bins, p = axes[1].hist(viable_non_neutral_fraction)
    axes[1].set_xlabel("Fraction of neighbors with viable, different phenotype")
    axes[1].set_ylabel(f"Genotype count\nSum: {len(viable_non_neutral_fraction)}")
    axes[1].vlines([np.mean(viable_non_neutral_fraction)], [0], [max(n)], label="mean", color="black")

    n, bins, p = axes[2].hist(non_neutral_fraction_of_viable)
    axes[2].set_xlabel("Fraction of viable neighbors with different phenotype")
    axes[2].set_ylabel(f"Genotype count\nSum: {len(non_neutral_fraction_of_viable)}")
    axes[2].vlines([np.mean(non_neutral_fraction_of_viable)], [0], [max(n)], label="mean", color="black")

    axes[0].set_xlim(0, 1)
    axes[1].set_xlim(0, 1)
    axes[2].set_xlim(0, 1)

    plt.legend()
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)

    