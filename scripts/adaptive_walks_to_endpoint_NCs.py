#!/usr/bin/env python

import argparse
import json

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--nc_paths", help="File with adaptive walks along neutral components", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()

end_ncs = []
with open(args.nc_paths, "r") as file:
    for line_ in file:
        line = line_.strip().split(" ")
        end_nc = line[-1]
        end_ncs.append(end_nc)

with open(args.output, "w") as out:
    for nc in end_ncs:
        out.write(nc + "\n")
