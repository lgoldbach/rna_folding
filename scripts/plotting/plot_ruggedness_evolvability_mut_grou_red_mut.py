#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib.lines as mlines
import matplotlib as mpl
from scipy.stats.stats import pearsonr

from rna_folding.parsing import load_phenotype_and_metric_from_file, read_ruggedness_per_ph_file
from rna_folding.utils import count_bp


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--phenotype_distributions", help="Phenotype distribution files ", 
                        required=True, nargs='+')
    parser.add_argument("-n", "--navigability", help="Navigablity per ph per fl ", 
                        required=True, nargs='+')
    parser.add_argument("-r", "--ruggedness", help="Peak sizes files ", 
                        required=True, nargs='+')
    parser.add_argument("-c", "--nc_graphs", help="nc graph pickle files ", 
                        required=True, nargs='+')
    parser.add_argument("-k", "--rugg_sample_size", help="Ruggedness sample size ", 
                        required=True, type=int)
    parser.add_argument("-o", "--output", help="Output file name "
                        "(should end in .pdf)", required=True)
    
    mut_graph_n = {2: 2,
                3: 2,
                4: 2,
                5: 4,
                6: 1,
                7: 2,
                8: 1,
                9: 2,
                10: 1,
                11: 1}

    mut_graph_s = {2: 1,
                3: 2,
                4: 3,
                5: 1,
                6: 6,
                7: 3,
                8: 8,
                9: 4,
                10: 10,
                11: 12}

    max_bp = {2: 2,
                3: 2,
                4: 2,
                5: 2,
                6: 2,
                7: 2,
                8: 3,
                9: 3,
                10: 4,
                11: 4}

    num_red_mut = {2: 2,
                    3: 2,
                    4: 6,
                    5: 0,
                    6: 0,
                    7: 0,
                    8: 0,
                    9: 4,
                    10: 2,
                    11: 0}
    
    navig = {1: 71,
            2: 77,
            3: 83,
            4: 81,
            5: 86,
            6: 53,
            7: 56,
            8: 73,
            9: 62,
            10: 67}

    args = parser.parse_args()
    fig, axes = plt.subplots(nrows=2, ncols=10, figsize=(50, 10))

    ### ruggedness vs ruggedness
    # rug data laden
    # local peak sizes
    ph_bias_over_mut_handles = []
    mut_groups = []
    red_mut = []
    evolvs = []
    navig_vals = []
    rugged_es_all = []
    rugged_all = []
    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]
    labels = ["1", "2", "3", "natural", "5", "6", "7", "8", "9", "10"]
    for i, (idx, label) in enumerate(zip(new_order_idx, labels)):
        phenotypes, ph_count = load_phenotype_and_metric_from_file(args.phenotype_distributions[idx])
        ph_to_count = dict(zip(phenotypes, ph_count))

        sum_of_folded = sum([c for p, c in zip(phenotypes, ph_count) if p != "............"])

        peak_sizes = read_ruggedness_per_ph_file(args.ruggedness[idx], n=args.rugg_sample_size)

        rugged_av_per_ph = {}  # compute average local peak size = ruggedness
        for ph in peak_sizes:
            # sum peak sizes and take average over sums
            peak_size_sums = [sum(ps) for ps in peak_sizes[ph]]
            rugged_av_per_ph[ph] = np.mean(peak_size_sums)  # average size of local peaks
    
        # nc graph load
        nc_graph = pickle.load(open(args.nc_graphs[idx], "rb"))

        nc_sizes = []
        nc_phenos = []
        evolvab = []
        for node in nc_graph.nodes:
            nc_sizes.append(nc_graph.nodes[node]["size"])  # size
            nc_phenos.append(nc_graph.nodes[node]["phenotype"])  # size
            evolvab.append(len(np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(node)]))) # unique neighbors

        nc_sizes_sort, evo_sort, nc_phenos_sort = zip(*[(nc, e, p) for nc, e, p in sorted(zip(nc_sizes, evolvab, nc_phenos), reverse=True)])

        rug_ph = {}
        for ph in rugged_av_per_ph:
            rug_ph[ph] = 0
            for node in nc_graph.nodes:
                nc_size = nc_graph.nodes[node]["size"]  # size
                evo = len(np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(node)])) # unique neighbors
                rug_ph[ph] += nc_size/(evo+1)

        rugged_es = np.mean(list(rug_ph.values()))
        rugged = np.mean([np.mean(rugged_av_per_ph[ph])+ph_to_count[ph] for ph in rugged_av_per_ph])

        axes[0][0].scatter(rugged_es, rugged, label=label)
        axes[0][0].set_xlabel("<Ruggedness> prediction from\nNeut. comp. size and evolvability", size=15)
        axes[0][0].set_ylabel("<Ruggedness>", size=15)
        axes[0][0].legend(title="GP map")
        
        rugged_es_all.append(rugged_es)
        rugged_all.append(rugged)

        ### evolvability as func of neutral component size with horizontal line
        # axes[1][1].scatter(np.log10(nc_sizes_sort), np.log10(evo_sort), label=label)
        p, res, l, o, k = np.polyfit(np.log10(nc_sizes_sort), np.log10(evo_sort), 1, full=True)
        poly1d_fn = np.poly1d(p) 
        # axes[1][1].plot(np.log10(nc_sizes_sort), poly1d_fn(np.log10(nc_sizes_sort)), label=label)
        
        rug_terms = [siz/ev for siz, ev in zip(nc_sizes_sort, evo_sort)]
        axes[1][1].scatter(np.log10(nc_sizes_sort), np.log10(rug_terms), label=label, s=1)

        rug_truth = sum([siz/ev for siz, ev in zip(nc_sizes_sort, evo_sort)])
        rug1 = sum([siz for siz, ev in zip(nc_sizes_sort, evo_sort)])
        rug1 = sum(nc_sizes_sort)
        # axes[1][2].scatter(rug_truth/sum_of_folded, rug1, label=label)
        axes[1][2].scatter(i+1, rug_truth/sum_of_folded, label=label)
        print(rug_truth/sum_of_folded)
        # rug2 = sum([siz/max(evo_sort) for siz, ev in zip(nc_sizes_sort, evo_sort)])
        rug2 = sum([siz/(ev/max(evo_sort)) for siz, ev in zip(nc_sizes_sort, evo_sort)])

        # axes[1][3].scatter(rug_truth, rug2, label=label)
        axes[1][3].scatter(i+1, rug2/sum_of_folded, label=label)

        non_red_frac = 1-(num_red_mut[new_order[i]]/12)
        print(i+1, non_red_frac)
        # rug3 = sum([(siz*non_red_frac)/ev for siz, ev in zip(nc_sizes_sort, evo_sort)])
        rug3 = sum([siz/(ev/np.mean(evo_sort)) for siz, ev in zip(nc_sizes_sort, evo_sort)])
        # axes[1][4].scatter(rug_truth, rug3, label=label)
        axes[1][4].scatter(i+1, rug3/sum_of_folded, label=label)
        
        # axes[1][1].scatter(np.log10(nc_sizes_sort), [siz/ev for siz, ev in zip(nc_sizes_sort, evo_sort)], label=label)
        # axes[1][i].axhline(y=len(np.unique(nc_phenos)), color="orange", label="No. Phenotypes", linewidth=3)
        # axes[1][i].axhline(y=np.mean(evo_sort), color="green", label="Mean. NC evo.", linewidth=3)
        # axes[1][i].set_ylim(0, 4)
        # axes[1][i].set_xlim(1.5, 6)
        # axes[1][i].set_xlabel("Neutral component size", size=15)
        # axes[1][i].set_ylabel("Evolvability", size=15)
        # axes[1][i].legend()
        axes[1][i].set_title(f"Alphabet: " + label, size=15)
        # axes[1][1].scatter(np.mean(np.log10(nc_sizes_sort)), np.mean(np.log10(evo_sort)), s=25, marker="s", label=label)
        # axes[1][i].scatter(np.mean(np.log10(nc_sizes_sort)), np.mean(evo_sort), color="orange", s=25, marker="s", label=label)
        y = np.mean([e/nc_s for e, nc_s in zip(evo_sort, nc_sizes_sort)])
        axes[0][7].scatter(i, y)

        ### mut group vs number of phenotypes in the top X, X and X percent.
        ph_counts_sort, ph_sort = zip(*sorted(zip(ph_count, phenotypes), reverse=True))
        ph_counts_sort = ph_counts_sort[1:]
        ph_sort = ph_sort[1:]

        folded_sum = sum(ph_counts_sort)
        top90_ph_count = 0
        top75_ph_count = 0
        top50_ph_count = 0

        s = 0
        for count in ph_counts_sort:
            s += count
            frac = s/folded_sum
            if frac < .9:
                top90_ph_count += 1
            if frac < .75:
                top75_ph_count += 1
            if frac < .5:
                top50_ph_count += 1

        if i == 5:
            shift = 0.1
        elif i == 2:
            shift = -0.1
        elif i == 7:
            shift = -0.1
        elif i == 4:
            shift = 0.1
        else:
            shift = 0
        i_ = new_order[i]
        sc = axes[0][1].scatter([mut_graph_s[i_]+shift], [top90_ph_count], marker="v")
        col = sc.get_facecolors()[0].tolist()
        axes[0][1].scatter([mut_graph_s[i_]+shift], [top75_ph_count], color=col, marker="o")
    
        axes[0][1].scatter([mut_graph_s[i_]+shift], [top50_ph_count], color=col, marker="s")
        
        if label=="natural":
            vl = axes[0][1].vlines(x=mut_graph_s[i_]+shift, ymin=top90_ph_count, ymax=top50_ph_count, color=col, label = "nat.")
        else:
            vl = axes[0][1].vlines(x=mut_graph_s[i_]+shift, ymin=top90_ph_count, ymax=top50_ph_count, color=col, label = label)
        ph_bias_over_mut_handles.append(vl)

        ### mut group vs evolvability
        axes[0][2].scatter(mut_graph_s[i_], np.mean(evo_sort), label=label)

        ### evolv vs number of redundant mutations
        axes[0][3].scatter(num_red_mut[i_]/12, np.mean(evo_sort), label=label)

        # collect data for plot
        mut_groups.append(mut_graph_s[i_])
        red_mut.append(num_red_mut[i_]/12)
        evolvs.append(np.mean(evo_sort))

        navig_vals.append(navig[i+1])

    ### mut group vs redundant, evo colored and annotate
    t = [e/37 for e in evolvs]
    cm = plt.cm.get_cmap('plasma')
    axes[0][4].scatter(red_mut, mut_groups, c=navig_vals, cmap=cm)
    # axes[0][4].scatter(red_mut, mut_groups, c=evolvs, cmap=cm)

    # sc = axes[0][5].scatter(red_mut, mut_groups, c=navig_vals, cmap=cm)  # for nav
    sc = axes[0][5].scatter(red_mut, mut_groups, c=evolvs, cmap=cm)
    plt.colorbar(sc, label="Mean evolvability")

    for i, idx in enumerate(new_order):
        print(i, idx)
        l = labels[i]
        axes[0][4].annotate(l, (red_mut[i], mut_groups[i]+.2))


    # plot stuff
    top90mark = mlines.Line2D([], [], color="black", marker='v', linestyle='None',
                        markersize=5, label='Top 90')
    top75mark = mlines.Line2D([], [], color="black", marker='o', linestyle='None',
                        markersize=5, label='Top 75')
    top50mark = mlines.Line2D([], [], color="black", marker='s', linestyle='None',
                        markersize=5, label='Top 50')
    l1 = axes[0][1].legend(handles=[top90mark, top75mark, top50mark], frameon=False, bbox_to_anchor=(0.83, 1), loc="upper right")
    axes[0][1].add_artist(l1)
    l2 = axes[0][1].legend(handles=ph_bias_over_mut_handles, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper right", title="Alphabet")
    axes[0][1].add_artist(l2)
    axes[0][1].set_ylabel("No. phenotypes in the top X percenile", size=15)
    axes[0][1].set_xlabel("Mutational group size", size=15)

    axes[0][2].legend(title="Alphabet")
    axes[0][2].set_ylabel("Mean evolvability of neutral components", size=15)
    axes[0][2].set_xlabel("Mutational group size", size=15)

    axes[0][3].legend(title="Alphabet")
    axes[0][3].set_ylabel("Mean evolvability of neutral components", size=15)
    axes[0][3].set_xlabel("Fraction of redundant mutations", size=15)

    # plt.colorbar()
    axes[0][4].set_ylabel("Mutational group size", size=15)
    axes[0][4].set_xlabel("Fraction of redundant mutations", size=15)
    
    axes[0][6].scatter(navig_vals, evolvs)

    axes[1][1].legend()
    axes[1][2].legend()
    axes[1][3].legend()
    axes[1][4].legend()

    r, p = pearsonr(rugged_es_all, rugged_all)
    p_str = "%.3g" % p
    print(r, p)

    plt.tight_layout()
    plt.savefig(args.output, format="pdf", dpi=30)
