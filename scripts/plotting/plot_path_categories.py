#!/usr/bin/env python

import argparse
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Output files containing path " \
                        "categories", nargs="+", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()
    fig, ax = plt.subplots(figsize=(8, 5))
    
    for i, path in enumerate(args.input, start=2):
        categs = pickle.load(open(path, "rb"))

        s = sum(list(categs.values()))
        w = 0.5
        sum_height = 0
        colors = ["green", "lightgreen", "orangered", "darkred"]
        for j, label in enumerate(["successful_mono", "successful_nonmono", "unsuccessful_mono", "unsuccessful_nonmono"]):
            height = categs[label]/s
            ax.bar(i, height, w, bottom=sum_height, color=colors[j])
            sum_height += height

    p1 = mpatches.Patch(color='green', label='Succ. Mono.')
    p2 = mpatches.Patch(color='lightgreen', label='Succ. Non-mono.')
    p3 = mpatches.Patch(color='orangered', label='Unsucc. Mono.')
    p4 = mpatches.Patch(color='darkred', label='Unsucc. Non-mono.')

    plt.legend(handles=[p1, p2, p3, p4], bbox_to_anchor=(1, .3))

    ax.set_xticks(list(range(2, 12)))
    ax.set_xlabel("GP map")
    ax.set_ylabel("Fraction of all evo. paths")

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
    
