#!/usr/bin/env python

import argparse
import pickle
from datetime import datetime


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="Input genotype-phenotype map "
                        "file", required=True)
    parser.add_argument("-o", "--output", help="File output for neutral components",
                        required=True)
    
    """
    #######
    ### Implement own neutral componenent algo to prevent memory issues
    ### write it parallizable. (possible starting from the stage of phenotype
    ### induces subgraphs, then can run BFS for each phenotype in parallel
    ### or do I even have to? 
    Could just get induced subgraph one by one and then use networkx connected_component
    because that is implemented as a generator.
    """
    
    