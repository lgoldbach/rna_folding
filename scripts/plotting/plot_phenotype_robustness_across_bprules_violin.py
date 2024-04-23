#!/usr/bin/env python

import argparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from rna_folding.utils import load_phenotype_and_metric_from_file

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

    fig, (ax1, ax2)  = plt.subplots(nrows=2, figsize=(12, 5))

    data = []
    # loop over files
    for bp_rule, i in enumerate(range(0, len(args.files), rn)):
        data.append([])
        for rank, file in enumerate(args.files[i:i+rn]):
            phenotypes, robust = load_phenotype_and_metric_from_file(file)
            robust = np.pad(robust, (0, 274-len(robust)), mode='constant')
            data[bp_rule].append(robust)
    
    data = np.array(data)  # dims: #{bp_rules} x #{rankings} x #{phenotypes}

    # use np.where to get indices along axis
    x = np.where(data>0)[0]
    hue = np.where(data>0)[-2]
    data_f = data[data>0].flatten()
    
    sns.violinplot(x=x, y=data_f, hue=hue, ax=ax1, legend=False)
    bp_rules = int(len(args.files)/int(args.num_of_ranks))
    plt.xticks(range(0, bp_rules), labels=range(2, bp_rules+2))

    ax1.set_xlabel("bp_rules")
    ax1.set_ylabel("Phenotype robustness")
    ax1.set_title("Phenotype robustness per base-pairing rule and ranking")
    ax1.legend(title="Rankings", fontsize=9, ncol=2)

    sns.violinplot(x=x, y=data_f, ax=ax2)
    ax2.set_xlabel("bp_rules")
    ax2.set_ylabel("Phenotype robustness")
    ax2.set_title("Phenotype robustness per base-pairing rule (rankings combined)")


    plt.tight_layout()

    ### full scatter of all data points
    # for i in range(data.shape[2]):
    #     for j, all_rankings in enumerate(data[:, :, i]):
    #         if i == 0:
    #             label = f"Rule {j+1}"
    #         else:
    #             label = ""

    #         ax.scatter(np.full_like(all_rankings, i), all_rankings, 
    #                    alpha=0.6, s=1,
    #                    label=label, color=plt.cm.tab10(j))
    
    # for j, r in enumerate(data):
    #     mean = np.nanmean(r, axis=0)
    #     stderr = np.nanstd(r, axis=0)
    #     # print(mean, list(range(data.shape[2])))
    #     # ax.scatter(range(data.shape[2]), mean, color=plt.cm.tab10(j), s=1, alpha=0.4)

    #     rang = np.arange(0, data.shape[2])
    #     mask = np.where(rang%(10+j)==0)[0]
    #     ax.errorbar(rang[mask], mean[mask], yerr=stderr[mask], 
    #                 color=plt.cm.tab10(j),
    #                 label=f"Rule {j+2}")

    # plt.legend()
    plt.savefig(args.output, format="pdf")






            
