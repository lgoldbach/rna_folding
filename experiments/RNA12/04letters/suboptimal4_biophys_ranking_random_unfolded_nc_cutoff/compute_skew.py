
import numpy as np
from scipy.stats import skew as skw
from scipy.stats import linregress
import matplotlib.pyplot as plt
from scipy.stats.stats import pearsonr, linregress

from rna_folding.utils import load_phenotype_and_metric_from_file

skews = []
lens = []

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


edge_numbers = np.array([1, 2, 2, 3, 3, 3, 4, 4, 5, 6])/6
fig, axes = plt.subplots(ncols=10, nrows=2, figsize=(50, 5))

axislabel_size = 10
labelsize = 10
label_size = 10
skews_all = []
mut_graph_s_all = []
for c, (i, e) in enumerate(zip([2, 3, 5, 7, 6, 4, 9, 8, 10, 11], edge_numbers), start=0):
    axes[0][c].set_box_aspect(1)
    axes[1][c].set_box_aspect(1)

    file = f"bp_graph{i}/ranking1/phenotype_distribution.txt"
    
    ph, counts = load_phenotype_and_metric_from_file(file)
    
    counts = np.array(sorted(counts)[:-1])[::-1]
    counter = 0
    for p, k in zip(ph, counts):
        if k > 0 and p != "............":
            counter+=1

    counts = [count/4**12 for count in counts if count > 0]
    
    # print(sum(counts)/4**12)
    # print(counts)
    log_counts = np.log10(counts)
    # print(log_counts)
    # x = np.log10(np.arange(stop=len(log_counts)+1, start=1))
    x = np.arange(stop=len(log_counts)+1, start=1)

    p, res, l, o, k = np.polyfit(x, log_counts, 1, full=True)
    
    poly1d_fn = np.poly1d(p) 
    m = p[0]
    b = p[1]

    slope, intercept, r_value, p_value, std_err = linregress(x, log_counts)
    print("XXXX", c+1, slope, r_value, p_value)



    skew = 3*(np.mean(log_counts) - np.median(log_counts))/np.std(log_counts)
    # print(np.mean(counts), np.median(counts), np.std(counts))

    skew2 = skw(log_counts)

    skews.append(skew2)
    lens.append(len(log_counts))

    axes[1][c].scatter(x, log_counts, marker='o', label=f"g-p map {c+1}", color=f"C{c}", s=4, zorder=10)
    axes[1][c].plot(x, poly1d_fn(x), linestyle="--", color="black", linewidth=0.75)

    p_str = "%.3g" % p_value
    r_str = np.round(r_value, 2)
    axes[1][c].text(0.05, 0.05, f'r={r_str}\np={p_str}', horizontalalignment='left', verticalalignment='bottom',transform = axes[1][c].transAxes, size=8)

    axes[1][c].set_xlim(0, 37)
    axes[1][c].set_ylim(-6.5, -1)


    
    if c == 3:
        axes[0][0].scatter(x, log_counts, marker='o', label="canon. g-p map", color='C3', s=4, zorder=10)
        axes[0][0].plot(x, poly1d_fn(x), linestyle="--", color="black", linewidth=0.75)
        for c_ in range(10):
            axes[1][c_].scatter(x, log_counts, marker='o', label="canon. g-p map", color='C3', s=4, zorder=10)
            axes[1][c_].plot(x, poly1d_fn(x), linestyle="--", color="black", linewidth=0.75)
        
    if c == 8:
        axes[0][0].scatter(x, log_counts, marker='o', label="g-p map 9", color='C8', s=4, zorder=10)
        axes[0][0].plot(x, poly1d_fn(x), linestyle="--", color="black", linewidth=0.75)
    if c == 9:
        axes[0][0].scatter(x, log_counts, marker='o', label="g-p map 10", color='C9', s=4, zorder=10)
        axes[0][0].plot(x, poly1d_fn(x), linestyle="--", color="black", linewidth=0.75)
    # if c == 6:
    #     axes[0][0].scatter(x, log_counts, marker='o', label="g-p map 7", color='C6', s=4, zorder=10)
    #     axes[0][0].plot(x, poly1d_fn(x), linestyle="--", color="black", linewidth=0.75)


    axes[0][0].set_ylabel("Phenotype\nfrequency (log10)", fontsize=axislabel_size)
    axes[0][0].set_xlabel("Rank", fontsize=axislabel_size)
    axes[0][0].tick_params(axis='both', which='major', labelsize=labelsize)
    axes[0][0].tick_params(axis='both', which='minor', labelsize=labelsize)
    axes[0][0].set_xlim(0, 38)

    y_pred = poly1d_fn(x)
    linr = linregress(log_counts, y_pred)

    # axes[c].set_title(f"Slope: {np.round(m, 3)}, r2: {np.round(linr[2], 3)}")
    # if i == 4:
    #     axes[1].plot(x, log_counts, color="black")
    #     axes[2].plot(x, log_counts, color="black")
    
    # if i < 9 and i != 4:
    #     axes[2].plot(x, log_counts)
    # elif i > 8 and i != 4:
    #     axes[1].plot(x, log_counts)

    # print("BP", c)
    # # print(skew, skew2)
    # print(m, b, res)
    # print("\n\n")

    # axes[0].scatter(len(log_counts), m, label=c)
    # axes[0].scatter(e, m, label=c)
    if c+1 == 2:
        axes[0][1].scatter(e-0.05, np.abs(m), label=f"{c+1}", s=10)
    elif c+1 == 3:
        axes[0][1].scatter(e+0.05, np.abs(m), label=f"{c+1}", s=10)
    elif c+1 == 5:
        axes[0][1].scatter(e-0.05, np.abs(m), label=f"{c+1}", s=10)
    elif c+1 == 6:
        axes[0][1].scatter(e+0.05, np.abs(m), label=f"{c+1}", s=10)
    elif c+1 == 7:
        axes[0][1].scatter(e-0.05, np.abs(m), label=f"{c+1}", s=10)
    elif c+1 == 8:
        axes[0][1].scatter(e+0.05, np.abs(m), label=f"{c+1}", s=10)
    else:
        axes[0][1].scatter(e, np.abs(m), label=f"{c+1}", s=10)
    
    # axes[0][1].set_xticks(np.arange(1, 7))
    skews_all.append(m)
    mut_graph_s_all.append(mut_graph_s[i])
    # axes[0].bar(c, m, label=c)
for c_ in range(10):
    axes[1][c_].legend(loc="upper right", frameon=False, fancybox=False, prop={'size': 8})
    axes[1][c_].set_xlabel("Rank", fontsize=axislabel_size)
    axes[1][c_].set_ylabel("Phenotype frequency (log 10)", fontsize=axislabel_size)
    axes[1][c_].tick_params(axis='both', which='major', labelsize=labelsize)
    axes[1][c_].tick_params(axis='both', which='minor', labelsize=labelsize)

axes[0][0].legend(loc="upper right", frameon=False, fancybox=False, prop={'size': 8})
axes[0][1].legend(loc="upper left", frameon=False, fancybox=False, prop={'size': 5}, title="Base-pairing rule")

axes[0][1].set_xlabel("Base-pairing rule promiscuity", fontsize=axislabel_size)
axes[0][1].set_ylabel("Phenotype bias", fontsize=axislabel_size)
axes[0][1].tick_params(axis='both', which='major', labelsize=labelsize)
axes[0][1].tick_params(axis='both', which='minor', labelsize=labelsize)

r, p = pearsonr(skews_all, mut_graph_s_all)
p_str = "%.3g" % p
# print(r, p)
plt.tight_layout()
plt.savefig("skew_vs_ph_count.pdf", format="pdf", dpi=10)