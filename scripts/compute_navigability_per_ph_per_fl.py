#!/usr/bin/env python

import argparse
import numpy as np



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--walk_lengths", help="Files that contain " \
    "adaptive walk lengths for each phenotype. One file per fitness landscape", 
    nargs="+", required=True)       
    parser.add_argument("-o", "--output", help="Output file name", 
                        type=str, required=True)
    parser.add_argument("-s", "--sample_size", help="How many walks were taken" \
    "per phenotype per fitness landscape", type=int, required=True)

    args = parser.parse_args()

    ph_success_count = {}
    print(args.sample_size)
    for i, path in enumerate(args.walk_lengths):  # i enumerates fitn. landsc.
        with open(path, "r") as file:
            for j, line_ in enumerate(file):
                line = line_.strip().split()
                print(j, line, path)
                # The data for one phenotype will contain args.sample_size many
                # lines +1 line for the header
                if j % (args.sample_size + 1) == 0:  # phenotype header
                    ph = line[0]
                    if not ph in ph_success_count:
                        ph_success_count[ph] = {}
                    ph_success_count[ph][i] = 0  # init entry for this fl for this ph
          
                else:
                    # successful walk
                    walk_length = int(line[0])
                    if walk_length > 0:  # -1 would indicate unsuccessful walk
                        ph_success_count[ph][i] += 1  # count as success

    with open(args.output, "w") as f:
        for ph in ph_success_count:
            f.write(ph + " ")
            for fl in ph_success_count[ph]:
                # divide count by sample size to get value between 0 and 1
                navig = ph_success_count[ph][fl] / args.sample_size 
                f.write(str(navig) + " ")  # space after ever fl success count 
            f.write("\n")  # new line after every phenotype block
    
        

                        
                        



    
                    

            
