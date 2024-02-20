#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import numpy as np

from rna_folding.utils import ranked_ph_log_distribution

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--files", help="Files containing phenotype counts "
                        "file", required=True, nargs="+")
    parser.add_argument("-n", "--num_of_ranks", help="Number of rankings per bp rule",
                        required=True)
    parser.add_argument("-o", "--output", help="File output for plot",
                        required=True)

    args = parser.parse_args()
    rn = int(args.num_of_ranks)

    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 4))

    data = []
    # loop over files 
    for bp_rule, i in enumerate(range(0, len(args.files), rn)):
        data.append([])
        for rank, file in enumerate(args.files[i:i+rn]):
            phenotypes, distr = ranked_ph_log_distribution(ph_distr_file=file)
            data[bp_rule].append(distr)
    
    data= np.array(data)  # dims: #{bp_rules} x #{rankings} x #{phenotypes}

    ## full scatter of all data points
    for i in range(data.shape[2]):
        for j, all_rankings in enumerate(data[:, :, i]):
            if i == 0:
                label = f"Rule {j+1}"
            else:
                label = ""

            ax2.scatter(np.full_like(all_rankings, i), all_rankings, 
                       alpha=0.8, s=1,
                       label=label, color=plt.cm.tab10(j))
    for j, r in enumerate(data):
        r[r==-np.inf] = np.nan
        mean = np.nanmean(r, axis=0)
        stderr = np.nanstd(r, axis=0)
        # print(mean, list(range(data.shape[2])))
        # ax.scatter(range(data.shape[2]), mean, color=plt.cm.tab10(j), s=1, alpha=0.4)

        rang = np.arange(0, data.shape[2])

        mask = np.where(rang%(10+j)==0)[0]
        # ax.errorbar(rang[mask], mean[mask], yerr=stderr[mask], 
        #             color=plt.cm.tab10(j),
        #             label=f"Rule {j+2}")
        ax1.errorbar(rang, mean, yerr=stderr,
                    color=plt.cm.tab10(j),
                    label=f"Rule {j+2}", alpha=0.6)
        
    ax1.set_title("Phenotype distribution of different base-pairing rules", loc="left")
    ax1.set_xlabel("Rank")
    ax1.set_ylabel("Phenotype frequency (log10)")
    ax2.set_xlabel("Rank")
    ax2.set_ylabel("Phenotype frequency (log10)")
    plt.legend()
    plt.savefig(args.output, format="pdf")






            
