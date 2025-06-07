
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.stats.stats import pearsonr
import seaborn as sns

with open("phenotype_props.csv", "r") as f:
    ph = []
    score = []
    de = []
    bul = []
    bp = []
    unp = []
    labels = ["Dangling ends", "Bulges", "Unpaired sites", "Dangling ends + Bulges + Unpaired sites"]
    for line_ in f:
        line = line_.strip().split(" ")
        ph.append(line[0])
        score.append(float(line[1]))
        de.append(int(line[2]))
        bul.append(int(line[3]))
        bp.append(int(len([s for s in ph[-1] if s == "(" or s ==")"])/2))
        unp.append(len([s for s in ph[-1] if s == "."]))

    
    fig, axes = plt.subplots(nrows=1, ncols=4, figsize=(20, 5), sharey=True)

    summ = [bul[i]+(unp[i]/2) for i in range(len(de))]
    for i, (x, ax, l) in enumerate(zip([de, bul, unp, summ], axes, labels)):
        y = score
        xy = np.vstack([x, y])
        z = gaussian_kde(xy)(xy)
        cax = ax.scatter(x, score, c="black", s=20, marker="x")
        # fig.colorbar(cax)
        if i == 0:
            ax.set_ylabel("Bradley-Terry score (a.u.)", size=15)
        ax.set_xlabel(l, size=15)
        r, p = pearsonr(x, score)
        p_str = "%.3g" % p
        ax.text(.7, .9, f'r = {np.round(r, 2)}\np = {p_str}', transform=ax.transAxes, horizontalalignment='left', size=10)
    
    plt.savefig("ranking_props.pdf", format="pdf", dpi=30)
