
import numpy as np
from scipy.stats import skew as skw
from scipy.stats import linregress
import matplotlib.pyplot as plt

from rna_folding.utils import load_phenotype_and_metric_from_file

skews = []
lens = []

edge_numbers = [1, 2, 2, 3, 3, 3, 4, 4, 5, 6]
fig, axes = plt.subplots(ncols=11, nrows=1, figsize=(55, 5))
for c, (i, e) in enumerate(zip([2, 3, 5, 7, 6, 4, 9, 8, 10, 11], edge_numbers), start=1):
    file = f"bp_graph{i}/ranking1/phenotype_distribution.txt"
    
    ph, counts = load_phenotype_and_metric_from_file(file)

    counts = np.array(sorted(counts)[:-1])[::-1]
    counts = [count for count in counts if count > 0]

    # print(counts)
    log_counts = np.log10(counts)
    # print(log_counts)
    # x = np.log10(np.arange(stop=len(log_counts)+1, start=1))
    x = np.arange(stop=len(log_counts)+1, start=1)

    p, res, l, o, k = np.polyfit(x, log_counts, 1, full=True)
    
    poly1d_fn = np.poly1d(p) 
    m = p[0]
    b = p[1]

    skew = 3*(np.mean(log_counts) - np.median(log_counts))/np.std(log_counts)
    # print(np.mean(counts), np.median(counts), np.std(counts))

    skew2 = skw(log_counts)

    skews.append(skew2)
    lens.append(len(log_counts))

    
    axes[c].plot(x, log_counts, 'yo', x, poly1d_fn(x), '--k', label=c)
    axes[c].legend()
    axes[c].set_xlim(0, 37)
    axes[c].set_ylim(0, 6.5)

    y_pred = poly1d_fn(x)
    linr = linregress(log_counts, y_pred)
    print(linr[2])

    axes[c].set_title(f"Slope: {np.round(m, 3)}, r2: {np.round(linr[2], 3)}")
    # if i == 4:
    #     axes[1].plot(x, log_counts, color="black")
    #     axes[2].plot(x, log_counts, color="black")
    
    # if i < 9 and i != 4:
    #     axes[2].plot(x, log_counts)
    # elif i > 8 and i != 4:
    #     axes[1].plot(x, log_counts)

    print("BP", c)
    # print(skew, skew2)
    print(m, b, res)
    print("\n\n")

    # axes[0].scatter(len(log_counts), m, label=c)
    axes[0].scatter(e, m, label=c)
    # axes[0].bar(c, m, label=c)

    axes[0].legend()

plt.savefig("skew_vs_ph_count.pdf", format="pdf", dpi=10)