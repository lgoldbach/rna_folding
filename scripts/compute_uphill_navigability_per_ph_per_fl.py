#!/usr/bin/env python

import argparse
import pickle
from rna_folding.adaptive_walks import contains_downhill_steps, load_fl_file_to_dict



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--paths", help="Files that contain " \
    "adaptive walk lengths for each phenotype. One file per fitness landscape", 
    nargs="+", required=True)       
    parser.add_argument("-g", "--gp_map", help="gp map pickle object", required=True)
    parser.add_argument("-f", "--fl", help="Fitness landscape files", required=True, nargs="+")
    parser.add_argument("-s", "--sample_size", help="How many walks were taken" \
    "per phenotype per fitness landscape", type=int, required=True)
    parser.add_argument("-o", "--output", help="Output file name", type=str, required=True)

    args = parser.parse_args()

    gpmap = pickle.load(open(args.gp_map, "rb"))  # load gp map
    
    ph_success_count = {}
    for i, (path, fl) in enumerate(zip(args.paths, args.fl)):  # i enumerates fitn. landsc.
        with open(path, "r") as file:
            for j, line_ in enumerate(file):
                line = line_.strip().split()
                # The data for one phenotype will contain args.sample_size many
                # lines +1 line for the header
                if j % (args.sample_size + 1) == 0:  # phenotype header
                    ph = line[0]
                    
                    if not ph in ph_success_count:
                        ph_success_count[ph] = {}
                    ph_success_count[ph][i] = 0  # init entry for this fl for this ph
                    ph_to_f = load_fl_file_to_dict(fl)  # load fitness landscape
                    ph_to_f[ph] = 1
                    continue  # enter lines that contain paths

                # we only reach this code after hitting ph header
                path = line
                if gpmap.map(path[-1]) != ph:  # not a successful path because it didnt reach target
                    continue
                else:
                    # has to contain only uphill or neutral steps
                    if not contains_downhill_steps(path, gpmap, ph_to_f):
                        ph_success_count[ph][i] += 1  # count as success

    with open(args.output, "w") as f:
        for ph in ph_success_count:
            f.write(ph + " ")
            for fl in ph_success_count[ph]:
                # divide count by sample size to get value between 0 and 1
                navig = ph_success_count[ph][fl] / args.sample_size 
                f.write(str(navig) + " ")  # space after ever fl success count 
            f.write("\n")  # new line after every phenotype block
    
        

                        
                        



    
                    

            
