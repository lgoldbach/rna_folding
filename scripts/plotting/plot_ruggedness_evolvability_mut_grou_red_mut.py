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

    bp_e = {2: 1,
                3: 2,
                4: 3,
                5: 2,
                6: 3,
                7: 3,
                8: 4,
                9: 4,
                10: 5,
                11: 6}

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
    
    navig = {1: .7182,
            2: .7690,
            3: .8249,
            4: .8168,
            5: .8674,
            6: .5128,
            7: .5745,
            8: .7399,
            9: .6197,
            10: .6710}
    
    navig_nc0 = {1: .54,
            2: .60,
            3: .56,
            4: .62,
            5: .68,
            6: .45,
            7: .42,
            8: .61,
            9: .50,
            10: .41}

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
    navig_vals_nc0 = []
    rugged_es_all = []
    rugged_all = []
    rugged_all_normed = []
    new_order = [2, 3, 5, 7, 6, 4, 9, 8, 10, 11]
    new_order_idx = [i - 2 for i in new_order]
    labels = ["1", "2", "3", "canonical", "5", "6", "7", "8", "9", "10"]
    bp_e_all = []
    for i, (idx, label) in enumerate(zip(new_order_idx, labels)):
        # phenotypes, ph_count = load_phenotype_and_metric_from_file(args.phenotype_distributions[idx])
        # ph_to_count = dict(zip(phenotypes, ph_count))

        # sum_of_folded = sum([c for p, c in zip(phenotypes, ph_count) if p != "............"])

        # peak_sizes = read_ruggedness_per_ph_file(args.ruggedness[idx], n=args.rugg_sample_size)

        # rugged_av_per_ph = {}  # compute average local peak size = ruggedness
        # for ph in peak_sizes:
        #     # sum peak sizes and take average over sums
        #     peak_size_sums = [sum(ps) for ps in peak_sizes[ph]]
        #     rugged_av_per_ph[ph] = np.mean(peak_size_sums)  # average size of local peaks
    
        # # nc graph load
        # nc_graph = pickle.load(open(args.nc_graphs[idx], "rb"))

        # nc_sizes = []
        # nc_phenos = []
        # evolvab = []
        # for node in nc_graph.nodes:
        #     nc_sizes.append(nc_graph.nodes[node]["size"])  # size
        #     nc_phenos.append(nc_graph.nodes[node]["phenotype"])  # size
        #     evolvab.append(len(np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(node)]))) # unique neighbors

        # nc_sizes_sort, evo_sort, nc_phenos_sort = zip(*[(nc, e, p) for nc, e, p in sorted(zip(nc_sizes, evolvab, nc_phenos), reverse=True)])

        # rug_ph = {}
        # for ph in rugged_av_per_ph:
        #     rug_ph[ph] = 0
        #     for node in nc_graph.nodes:
        #         nc_size = nc_graph.nodes[node]["size"]  # size
        #         evo = len(np.unique([nc_graph.nodes[neigh]["phenotype"] for neigh in nc_graph.neighbors(node)])) # unique neighbors
        #         rug_ph[ph] += nc_size/(evo+1)

        # rugged_es = np.mean(list(rug_ph.values()))
        # rugged = np.mean([np.mean(rugged_av_per_ph[ph])+ph_to_count[ph] for ph in rugged_av_per_ph])

        # axes[0][0].scatter(rugged_es/sum_of_folded, rugged/sum_of_folded, label=label)
        # axes[0][0].set_xlabel("Ruggedness (analytical prediction)", size=15)
        # axes[0][0].set_ylabel("Ruggedness (numerical estimate)", size=15)
        # axes[0][0].legend(title="g-p map", frameon=False, fancybox=False, loc="upper left")
        # axes[0][0].plot([0, 1], [0, 1], transform=axes[0][0].transAxes, color="0.2", zorder=-10, linewidth=0.5, linestyle="--")
        # axes[0][0].tick_params(axis='y', which='major', labelsize=15)
        # axes[0][0].tick_params(axis='y', which='minor', labelsize=15)
        # axes[0][0].tick_params(axis='x', which='major', labelsize=15)
        # axes[0][0].tick_params(axis='x', which='minor', labelsize=15)

        # rugged_es_all.append(rugged_es)
        # rugged_all.append(rugged)

        # ### evolvability as func of neutral component size with horizontal line
        # # axes[1][1].scatter(np.log10(nc_sizes_sort), np.log10(evo_sort), label=label)
        # p, res, l, o, k = np.polyfit(np.log10(nc_sizes_sort), np.log10(evo_sort), 1, full=True)
        # poly1d_fn = np.poly1d(p) 
        # # axes[1][1].plot(np.log10(nc_sizes_sort), poly1d_fn(np.log10(nc_sizes_sort)), label=label)
        
        # rug_terms = [siz/ev for siz, ev in zip(nc_sizes_sort, evo_sort)]
        # axes[1][1].scatter(np.log10(nc_sizes_sort), np.log10(rug_terms), label=label, s=1)

        # rug_truth = sum([siz/ev for siz, ev in zip(nc_sizes_sort, evo_sort)])
        # rug1 = sum([siz for siz, ev in zip(nc_sizes_sort, evo_sort)])
        # rug1 = sum(nc_sizes_sort)

        # rugged_all_normed.append(rug_truth/sum_of_folded)
        # # axes[1][2].scatter(rug_truth/sum_of_folded, rug1, label=label)
        # if i == 1:
        #     axes[1][2].bar(i+1-0.2, rug_truth/sum_of_folded, width=0.4, color="grey", label="Ruggedness")
        # else:
        #     axes[1][2].bar(i+1-0.2, rug_truth/sum_of_folded, width=0.4, color="grey")
        
        # ratio = max(evo_sort)/35
        # # rug2 = sum([siz/max(evo_sort) for siz, ev in zip(nc_sizes_sort, evo_sort)])
        # rug2 = sum([siz/(ev/ratio) for siz, ev in zip(nc_sizes_sort, evo_sort)])

        # # axes[1][3].scatter(rug_truth, rug2, label=label)
        # if i == 1:
        #     axes[1][2].bar(i+1+0.2, rug2/sum_of_folded, width=0.4, color="firebrick", label="Ruggedness\nmax -normalized")
        # else:
        #     axes[1][2].bar(i+1+0.2, rug2/sum_of_folded, width=0.4, color="firebrick")

        # non_red_frac = 1-(num_red_mut[new_order[i]]/12)
        # # print(i+1, non_red_frac)
        # # rug3 = sum([(siz*non_red_frac)/ev for siz, ev in zip(nc_sizes_sort, evo_sort)])
        # rug3 = sum([siz/(ev/np.mean(evo_sort)) for siz, ev in zip(nc_sizes_sort, evo_sort)])
        # # axes[1][4].scatter(rug_truth, rug3, label=label)
        # axes[1][4].bar(i+1, rug3/sum_of_folded, label=label)
        
        # # axes[1][1].scatter(np.log10(nc_sizes_sort), [siz/ev for siz, ev in zip(nc_sizes_sort, evo_sort)], label=label)
        # # axes[1][i].axhline(y=len(np.unique(nc_phenos)), color="orange", label="No. Phenotypes", linewidth=3)
        # # axes[1][i].axhline(y=np.mean(evo_sort), color="green", label="Mean. NC evo.", linewidth=3)
        # # axes[1][i].set_ylim(0, 4)
        # # axes[1][i].set_xlim(1.5, 6)
        # # axes[1][i].set_xlabel("Neutral component size", size=15)
        # # axes[1][i].set_ylabel("Evolvability", size=15)
        # # axes[1][i].legend()
    
        # # axes[1][1].scatter(np.mean(np.log10(nc_sizes_sort)), np.mean(np.log10(evo_sort)), s=25, marker="s", label=label)
        # # axes[1][i].scatter(np.mean(np.log10(nc_sizes_sort)), np.mean(evo_sort), color="orange", s=25, marker="s", label=label)
        # y = np.mean([e/nc_s for e, nc_s in zip(evo_sort, nc_sizes_sort)])
        # axes[0][7].scatter(i, y)

        # ### mut group vs number of phenotypes in the top X, X and X percent.
        # ph_counts_sort, ph_sort = zip(*sorted(zip(ph_count, phenotypes), reverse=True))
        # ph_counts_sort = ph_counts_sort[1:]
        # ph_sort = ph_sort[1:]

        # folded_sum = sum(ph_counts_sort)
        # top90_ph_count = 0
        # top75_ph_count = 0
        # top50_ph_count = 0

        # s = 0
        # for count in ph_counts_sort:
        #     s += count
        #     frac = s/folded_sum
        #     if frac < .9:
        #         top90_ph_count += 1
        #     if frac < .75:
        #         top75_ph_count += 1
        #     if frac < .5:
        #         top50_ph_count += 1

        # if i == 5:
        #     shift = 0.1
        # elif i == 2:
        #     shift = -0.1
        # elif i == 7:
        #     shift = -0.1
        # elif i == 4:
        #     shift = 0.1
        # else:
        #     shift = 0
        i_ = new_order[i]
        # sc = axes[0][1].scatter([mut_graph_s[i_]+shift], [top90_ph_count], marker="v")
        # col = sc.get_facecolors()[0].tolist()
        # axes[0][1].scatter([mut_graph_s[i_]+shift], [top75_ph_count], color=col, marker="o")
    
        # axes[0][1].scatter([mut_graph_s[i_]+shift], [top50_ph_count], color=col, marker="s")
        
        # if label=="natural":
        #     vl = axes[0][1].vlines(x=mut_graph_s[i_]+shift, ymin=top90_ph_count, ymax=top50_ph_count, color=col, label = "nat.")
        # else:
        #     vl = axes[0][1].vlines(x=mut_graph_s[i_]+shift, ymin=top90_ph_count, ymax=top50_ph_count, color=col, label = label)
        # ph_bias_over_mut_handles.append(vl)

        # ### mut group vs evolvability
        # axes[0][2].scatter(mut_graph_s[i_], np.mean(evo_sort), label=label)

        # ### evolv vs number of redundant mutations
        # axes[0][3].scatter(num_red_mut[i_]/12, np.mean(evo_sort), label=label)


        # collect data for plot
        mut_groups.append(mut_graph_s[i_])
        if i == 4:
            red_mut.append((num_red_mut[i_]/12)-0.01)
        elif i == 3:
            red_mut.append((num_red_mut[i_]/12)+0.01)
        else:
            red_mut.append(num_red_mut[i_]/12)
        bp_e_all.append(bp_e[i_])
        # evolvs.append(np.mean(evo_sort))

        navig_vals.append(navig[i+1])
        navig_vals_nc0.append(navig_nc0[i+1])

    ### mut group vs redundant, evo colored and annotate
    t = [e/37 for e in evolvs]
    cm = plt.cm.get_cmap('Oranges')
    ##### NC 0 !!!!!!!!!!!!!!!! ################
    sc = axes[0][4].scatter(red_mut, mut_groups, c=navig_vals_nc0, cmap=cm, s=10)
    # sc = axes[0][4].scatter(red_mut, mut_groups, c=navig_vals, cmap=cm, s=10)
    plt.colorbar(sc, label="Navigability")

    sc = axes[0][8].scatter(red_mut, np.array(bp_e_all)/2, c=navig_vals, cmap=cm, s=100, edgecolor="black")
    axes[0][8].set_box_aspect(1)
    fontsize = 18
    x_tick_fontsize=18
    y_tick_fontsize=18
    axes[0][8].tick_params(axis='y', which='major', labelsize=y_tick_fontsize)
    axes[0][8].tick_params(axis='y', which='minor', labelsize=y_tick_fontsize)
    axes[0][8].tick_params(axis='x', which='major', labelsize=x_tick_fontsize)
    axes[0][8].tick_params(axis='x', which='minor', labelsize=x_tick_fontsize)
    axes[0][8].set_ylabel("Promiscuity\nof base-pairing rule", size=fontsize)
    axes[0][8].set_xlabel("Redundancy\nof base-pairing rule", size=fontsize)
    axes[0][8].set_xlim(-.05, 0.54)
    axes[0][8].set_ylim(0, 3.2)

    # axes[0][4].scatter(red_mut, mut_groups, c=evolvs, cmap=cm)

    # sc = axes[0][5].scatter(red_mut, mut_groups, c=navig_vals, cmap=cm)  # for nav
    # sc = axes[0][5].scatter(red_mut, mut_groups, c=evolvs, cmap=cm)


    for i, idx in enumerate(new_order):
        l = labels[i]
        axes[0][4].annotate(l, (red_mut[i], mut_groups[i]+.2))
        if l == "canonical":
            axes[0][8].annotate(l, (red_mut[i]+0.02, (bp_e_all[i]/2)-0.15), size=15)
        elif l == "5":
            axes[0][8].annotate(l, (red_mut[i]-0.03, (bp_e_all[i]/2)+0.05), size=15)
        else:
            axes[0][8].annotate(l, (red_mut[i]+0.01, (bp_e_all[i]/2)+0.05), size=15)

    # plot stuff
    # top90mark = mlines.Line2D([], [], color="black", marker='v', linestyle='None',
    #                     markersize=5, label='Top 90')
    # top75mark = mlines.Line2D([], [], color="black", marker='o', linestyle='None',
    #                     markersize=5, label='Top 75')
    # top50mark = mlines.Line2D([], [], color="black", marker='s', linestyle='None',
    #                     markersize=5, label='Top 50')
    # l1 = axes[0][1].legend(handles=[top90mark, top75mark, top50mark], frameon=False, bbox_to_anchor=(0.83, 1), loc="upper right")
    # axes[0][1].add_artist(l1)
    # l2 = axes[0][1].legend(handles=ph_bias_over_mut_handles, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper right", title="Alphabet")
    # axes[0][1].add_artist(l2)
    # axes[0][1].set_ylabel("No. phenotypes in the top X percenile", size=15)
    # axes[0][1].set_xlabel("Mutational group size", size=15)

    # axes[0][2].legend(title="Alphabet")
    # axes[0][2].set_ylabel("Mean evolvability of neutral components", size=15)
    # axes[0][2].set_xlabel("Mutational group size", size=15)

    # axes[0][3].legend(title="Alphabet")
    # axes[0][3].set_ylabel("Mean evolvability of neutral components", size=15)
    # axes[0][3].set_xlabel("Fraction of redundant mutations", size=15)

    # # plt.colorbar()
    # axes[0][4].set_ylabel("Mutational group size", size=15)
    # axes[0][4].set_xlabel("Fraction of redundant mutations", size=15)
    
    # axes[0][6].scatter(navig_vals, evolvs)

    # axes[1][1].legend()
    # # axes[1][2].legend()
    
    # fontsize = 18
    # x_tick_fontsize=15
    # y_tick_fontsize=15
    # axes[1][2].set_box_aspect(1)
    # axes[1][2].set_xticks(list(range(1, 11)))
    # p =[1, 2, 3, " ", 5, 6, 7, 8, 9, 10]
    # axes[1][2].set_xticklabels(p) #(list(range(1, 11)))
    # axes[1][2].tick_params(axis='both', which='major', labelsize=fontsize)
    # axes[1][2].tick_params(axis='both', which='minor', labelsize=fontsize)
    # axes[1][2].set_ylabel("Ruggedness", size=fontsize)
    # axes[1][2].set_xlabel("g-p map", size=fontsize)
    # axes[1][2].tick_params(axis='y', which='major', labelsize=y_tick_fontsize)
    # axes[1][2].tick_params(axis='y', which='minor', labelsize=y_tick_fontsize)
    # axes[1][2].tick_params(axis='x', which='major', labelsize=x_tick_fontsize)
    # axes[1][2].tick_params(axis='x', which='minor', labelsize=x_tick_fontsize)
    # axes[1][2].set_title("canon. ", size=15)

    # axes[1][3].set_box_aspect(1)
    # axes[1][3].set_xticks(list(range(1, 11)))
    # p =[1, 2, 3, " ", 5, 6, 7, 8, 9, 10]
    # axes[1][3].set_xticklabels(p) #(list(range(1, 11)))
    # axes[1][3].tick_params(axis='y', which='major', labelsize=y_tick_fontsize)
    # axes[1][3].tick_params(axis='y', which='minor', labelsize=y_tick_fontsize)
    # axes[1][3].tick_params(axis='x', which='major', labelsize=x_tick_fontsize)
    # axes[1][3].tick_params(axis='x', which='minor', labelsize=x_tick_fontsize)
    # axes[1][3].set_ylabel("Ruggedness", size=fontsize)
    # axes[1][3].set_xlabel("g-p map", size=fontsize)


    # # legenfontsize=14
    # axes[1][2].legend(loc="upper left", fontsize=8, title_fontsize=8, frameon=False)
    # # axes[1][3].legend(loc="upper left", fontsize=legenfontsize, title_fontsize=legenfontsize)
    # for row in axes:
    #     for ax in row:
    #         ax.set_box_aspect(1)
    # r, p = pearsonr(rugged_es_all, rugged_all)
    # p_str = "%.3g" % p

    # r, p = pearsonr(navig_vals, rugged_all_normed)
    # axes[1][8].scatter(rugged_all_normed, navig_vals)
    # p_str = "%.3g" % p
    # print(r, p)

    plt.tight_layout(pad=1.5)
    plt.savefig(args.output, format="pdf", dpi=30)
