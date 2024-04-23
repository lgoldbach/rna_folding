#!/usr/bin/env python

import argparse


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--distribution", help="Phenotype distribution file "
                        "file", required=True)
    parser.add_argument("-r", "--robustness", help="Phenotype robustness file", 
                        required=True)
    parser.add_argument("-o", "--output", help="Joined output", 
                        required=True)
    
    args = parser.parse_args()
    
    def file_to_dict(file):
        with open(file, "r") as f:
            dic = {}
            for line in f:
                ph, v = line.strip().split(" ")
                dic[ph] = float(v)
        return dic
    

    rob_dic = file_to_dict(args.robustness)
    distr_dic = file_to_dict(args.distribution)

    ph_sum = sum(list(distr_dic.values()))
    with open(args.output, "w") as f:
        for ph in rob_dic:
            ph_freq = distr_dic[ph]/ph_sum
            f.write(str(rob_dic[ph]) + " " + str(ph_freq) + "\n")
