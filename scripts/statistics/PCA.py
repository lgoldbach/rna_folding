#!/usr/bin/env python

import argparse
import numpy as np
from sklearn.decomposition import PCA

def list_of_strings(arg):
    return arg.split(',')

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", type=str, 
                        required=True, help="output file name")
    parser.add_argument("-i", "--inputs", type=list_of_strings,
                        required=True, help="List of files each containing a "
                        "single 1D data set")

    args = parser.parse_args()
    print(args.inputs)
    X = np.array([np.loadtxt(file) for file in args.inputs])
    pca = PCA()
    pca.fit(X)

    print(pca.components_)
    print(pca.explained_variance_)
    
