#!/usr/bin/env python

import argparse
import pickle

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input consensus matrix, 2D "
                        "numpy array in pickle format", required=True)
    parser.add_argument("-o", "--output", help="Binary consensus matrix, "
                        "2D numpy array", required=True)

    args = parser.parse_args()
        
    A = pickle.load(open(args.input, "rb"))

    for i in range(A.shape[0]):
        for j in range(i+1, A.shape[1]):
            if A[i,j] > A[j,i]:
                A[i, j] = 1
                A[j, i] = 0
            elif A[i,j] < A[j,i]:
                A[j, i] = 1
                A[i, j] = 0
            else:
                A[i, j] = 0
                A[j, i] = 0
    
    pickle.dump(A, open(args.output, "wb"))
