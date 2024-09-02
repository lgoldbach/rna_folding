#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from rna_folding.utils import ranked_ph_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input files with adaptive walk "
                        "stats triplets: <genotype>,<nc_id>,<fitness> separated "
                        " by spaces. One line represents one walk and  ", 
                        required=True)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots()

    xs = []
    ys = []
    count = 0
    count_w = 0
    with open(args.input, "r") as file:
        for c, line_ in enumerate(file):
            # if c == 100:
            #     break
            line = line_.strip().split(" ")
            if line[0] == "ph":  # it's just a phenotype header
                continue
            count_w += 1
            x = []
            y = []  # init with impossible fitness value to match
            gts = []
            for i, triplet_str in enumerate(line):
                triplet = triplet_str.split(",")
                gt = triplet[0]
                nc = int(triplet[1])
                f = float(triplet[2])
               
                y.append(f)
                x.append(i)
               
                gts.append(gt)
            
                new_x = [x[0]]
                new_y = [y[0]]
                for p, e in enumerate(x[1:-1], start=1):  # start iteration at second element
                    if gts[p] != gts[p-1] or gts[p] != gts[p+1]:
                        # if y[p] != y[p-1] or y[p] != y[p+1]:
                        new_x.append(x[p])
                        new_y.append(y[p])
                new_x.append(x[-1])
                new_y.append(y[-1])
                    

            if y[-1] == 1:
                count += 1
            xs.append(new_x)
            ys.append(new_y)
    
    # pick random set of walks
    sample_size = 50
    idx = np.random.choice(len(xs), size=sample_size)

    x_plot = [x for i, x in enumerate(xs) if i in idx]
    y_plot = [y for i, y in enumerate(ys) if i in idx]
    for x, y in zip(x_plot, y_plot):
        if y[-1] == 1:
            ax.plot(x, y, color="red", linewidth=1, marker="o", markersize=1, alpha=0.5)
        else:
            ax.plot(x, y, color="black", linewidth=1, marker="o", markersize=1, alpha=0.3)
    c = 0
    for y in y_plot:
        if y[-1] == 1:
            c += 1
    print(args.input, np.round(c/sample_size, 2))
    ax.grid()

    plt.savefig(args.output, format="pdf", dpi=30)
