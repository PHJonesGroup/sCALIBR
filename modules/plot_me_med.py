import numpy as np
import matplotlib.pyplot as plt
import os

def plot_me_med(tab_norm_T0, tab_norm_T1, output_dir):
    # all replicate columns for each condition (whatever the count)
    T0 = tab_norm_T0.iloc[:, 2:].to_numpy()
    T1 = tab_norm_T1.iloc[:, 2:].to_numpy()

    me0,  med0  = np.mean(T0, axis=0), np.median(T0, axis=0)
    me1,  med1  = np.mean(T1, axis=0), np.median(T1, axis=0)

    # concatenate across both conditions; row 0 = means, row 1 = medians
    means   = np.concatenate([me0,  me1])
    medians = np.concatenate([med0, med1])
    me_med  = np.vstack([means, medians])

    # column labels straight from the tables, so they match the data
    labels = list(tab_norm_T0.columns[2:]) + list(tab_norm_T1.columns[2:])

    # ---------- visualization (mean + median , tallest at back) ----------
    plt.figure(figsize=(max(8, 0.6 * len(labels)), 5))
    indices = np.arange(len(labels))
    bar_width = 0.6

    # For each sample, plot the larger value first (back) then the smaller (front),
    # so both remain visible in a single bar position.
    for i in range(len(labels)):
        if means[i] >= medians[i]:
            plt.bar(i, means[i],   width=bar_width, color='tab:blue')    # back
            plt.bar(i, medians[i], width=bar_width, color='tab:orange')  # front
        else:
            plt.bar(i, medians[i], width=bar_width, color='tab:orange')  # back
            plt.bar(i, means[i],   width=bar_width, color='tab:blue')    # front

    # legend proxies (loop bars aren't labelled individually)
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(color='tab:blue',   label='Mean'),
        Patch(color='tab:orange', label='Median'),
    ]

    plt.ylabel('CPM', fontsize=12)
    plt.title('Mean and Median of normalised read counts', fontsize=14)
    plt.xticks(indices, labels, rotation=45, ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.5)

    plt.legend(handles=legend_handles, loc='center left',
               bbox_to_anchor=(1.02, 0.5), borderaxespad=0., frameon=False)
    plt.tight_layout()

    out_path = os.path.join(output_dir, "normalised_counts_mean_med.png")
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    return me_med, labels