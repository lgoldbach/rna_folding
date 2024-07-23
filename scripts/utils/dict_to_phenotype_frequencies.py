#!/usr/bin/env python

import argparse
import pickle


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, type=str, 
                        help="Pickle file containing a dict that maps ph to "
                        "its count")
    parser.add_argument("-o", "--output", required=False, type=str,
                        help="File to save phenotype counts")

    args = parser.parse_args()

    
    gt_per_ph = pickle.load(open(args.input, "rb"))
    
    with open(args.output, "w") as file:
        for ph in gt_per_ph:
            file.write(ph + " " + str(gt_per_ph[ph]) + "\n")
