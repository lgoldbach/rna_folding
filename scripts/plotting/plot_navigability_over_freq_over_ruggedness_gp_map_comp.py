#!/usr/bin/env python

import argparse
import pickle
import RNA
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats.stats import pearsonr, linregress

from rna_folding.parsing import load_phenotype_and_metric_from_file, read_ruggedness_file



if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--navigability", help="Navigability per peak ", nargs="+", required=True)
    parser.add_argument("-n", "--nc_graph", help="NC graph", nargs="+", required=True)
    parser.add_argument("-r", "--ruggedness", help="Peak sizes from sample", nargs="+", required=True)
    parser.add_argument("-p", "--ph_dist", help="phenotype distribution", nargs="+", required=True)
    parser.add_argument("-o", "--output", help="Output file for plot (.pdf)",
                        required=True)
    
    args = parser.parse_args()
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))
    ax1.set_box_aspect(1)
    ax2.set_aspect("equal")

    print(args.nc_graph)
    labels = ["RNA12", "HP5x5", "HP20", "S_2,8"]
    for i, (nc_graph, navig, rugg, ph_dist, label) in enumerate(zip(args.nc_graph, args.navigability, args.ruggedness, args.ph_dist, labels)): 
        nc_graph = pickle.load(open(nc_graph, "rb"))

        ph_to_nc = {nc_graph.nodes[nc]["phenotype"]: nc for nc in nc_graph.nodes}        

        ph, counts = load_phenotype_and_metric_from_file(ph_dist)
        ph_to_count = dict(zip(ph, counts))

        peak_sizes = read_ruggedness_file(rugg)

        # rugg_mean = np.mean([sum(s) for s in peak_sizes])  # average of peak size 

        y_nav = []
        x = []
        global_peaks = []
        local_peaks = []
        with open(navig, "r") as f:
            for line_ in f:
                line = line_.strip().split(" ")
                global_peak_ph = line[0]
                nav = float(line[1])
                rugg = int(line[2])
                y_nav.append(nav)

                ph_count = ph_to_count[global_peak_ph]
                x.append(ph_count/rugg)

                local_peaks.append(rugg-ph_count)
                global_peaks.append(ph_count)

        
        ax1.scatter(x, y_nav, color=f"C{i}", s=10, alpha=.6, linewidths=0)

        x_m = np.mean(x)
        y_m = np.mean(y_nav)

        xq1 = [max(0, x_m-np.percentile(x, q=25))]
        xq2 = [max(0, np.percentile(x, q=75)-x_m)]
    
        yq1 = [max(0, y_m-np.percentile(y_nav, q=25))]
        yq2 = [np.percentile(y_nav, q=75)-y_m]

        print(xq1, xq2, x_m, np.percentile(x, q=25), np.percentile(x, q=25))
        ax1.errorbar(x_m, y_m, xerr=(xq1, xq2), yerr=(yq1, yq2), elinewidth=1, marker="s", markersize=7, label=label)

        p, res, l, o, k = np.polyfit(x, y_nav, 1, full=True)
        poly1d_fn = np.poly1d(p) 
        if label != "HP20":
            ax1.plot(x, poly1d_fn(x), linestyle=(0, (5, 10)), color=f"C{i}", linewidth=1, label=label)

        slope, intercept, r_value, p_value, std_err = linregress(x, y_nav)
        print(label, slope, r_value, p_value)
        p_str = "%.3g" % p_value
        # if label == "HP5x5":
        #     ax1.text(.8, .45, f'r = {np.round(r_value, 2)}\np = {p_str}', transform=ax1.transAxes, horizontalalignment='left', size=12)
        # elif label == "RNA12":
        #     ax1.text(.9, .5, f'r = {np.round(r_value, 2)}\np = {p_str}', transform=ax1.transAxes, horizontalalignment='left', size=12)
    

        ax2.scatter(local_peaks, global_peaks, label=label, color=f"C{i}")

    ax1.set_xlabel("global peak size / combined local peak size", fontsize=15)
    ax1.set_ylabel("Navigability", fontsize=20)

    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax1.tick_params(axis='both', which='minor', labelsize=15)

    ax1.plot([0, 1], [0, 1], zorder=-10, color="0.5", linewidth=1, transform=ax1.transAxes)
    # ax1.legend(loc="lower right", fancybox=False, frameon=False, fontsize=15)
    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
