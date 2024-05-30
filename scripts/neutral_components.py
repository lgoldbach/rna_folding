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
    

    args = parser.parse_args()
    
    print("start loading", datetime.now())
    gpm = pickle.load(open(args.file, "rb"))
    print("done", datetime.now())
    # if edges are not added, add them and save updated gpm in pickle file
    if not gpm.edges:
        print("start edge addition", datetime.now())
        gpm.add_hamming_edges()
        pickle.dump(gpm, open(args.file, "wb")) 
        print("Done", datetime.now())

    print("start neutral components", datetime.now())
    nc = gpm.neutral_components(return_ids=True)
    print("Done", datetime.now())

    nc_l = [list(i) for i in nc]
    pickle.dump(nc_l, open(args.output, "wb"))
    