import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

from .find_columns_indiv import find_columns_indiv

def format_data(index_scheme, indiv, norm_tab, rep, output_dir):
    """
    Reformats data depending on scheme (1 or 2) and computes mean/median.
    Bars are stacked (median overlays mean).

    Parameters
    ----------
    index_scheme : int
        Scheme type (1 or 2). Determines how individual columns are selected.
    indiv : str
        Individual ID (used if index_scheme == 1).
    norm_tab : DataFrame
        Normalized read-count table with columns like:
    output_dir : str
        Directory to save plots.

    Returns
    -------
    T_norm_indiv : DataFrame
        Normalisation table
    norm_dat_values : ndarray
        Raw numeric matrix from norm_tab (without gRNA/gene columns).
    me_med_nd : ndarray
        2×N array of mean and median values for each column.
    """

    siz = norm_tab.shape

    me_med_nd = None
    norm_dat_values = None
    T_norm_indiv = None

    if siz[1] > 0:
        # ---------- compute mean & median across rows (per column) ----------
        norm_dat_values = norm_tab.iloc[:, 2:].values  # drop gRNA, gene
        me_nd  = np.mean(norm_dat_values, axis=0)
        med_nd = np.median(norm_dat_values, axis=0)
        me_med_nd = np.vstack([me_nd, med_nd])

        # ---------- visualization (stacked mean + median) ----------
        plt.figure(figsize=(8, 5))
        indices = np.arange(len(me_nd))
        bar_width = 0.6

        # Plot mean
        plt.bar(indices, me_nd, width=bar_width,
                color='tab:blue', alpha=1, label='Mean')

        # Overlay median on top of mean (semi-transparent)
        plt.bar(indices, med_nd, width=bar_width,
                color='tab:orange', alpha=1, label='Median')

        # Label axes and ticks
        plt.ylabel('CPM', fontsize=12)
        plt.title('Mean and Median of normalised read counts', fontsize=14)

        xlabels = norm_tab.columns[2:]
        plt.xticks(indices, xlabels, rotation=0)
        plt.grid(axis='y', linestyle='--', alpha=0.5)

        plt.legend(
            loc='center left',        
            bbox_to_anchor=(1.02, 0.5), 
            borderaxespad=0.,         
            frameon=False
        )
        plt.tight_layout()

        out_path = os.path.join(output_dir, "normalised_counts_mean_med.png")
        plt.savefig(out_path, dpi=300, bbox_inches="tight")
        plt.close()

        # ---------- scheme-specific table ----------
        if index_scheme == 2:
            T_norm_indiv = norm_tab.copy()

        elif index_scheme == 1:
            ind_found, num_found = find_columns_indiv.find_columns_indiv(norm_tab, indiv)
            if num_found >= 2:
                T0 = norm_tab.iloc[:, ind_found[0]]
                T1 = norm_tab.iloc[:, ind_found[1]]
                T_norm_indiv = pd.DataFrame({
                    'gRNA': norm_tab.iloc[:, 0],
                    'gene': norm_tab.iloc[:, 1],
                    'T0_F': T0,
                    'T0_M': T0,
                    'T1_F': T1,
                    'T1_M': T1
                })
                norm_dat = T_norm_indiv.iloc[:, 2:].values
            else:
                print("Not enough matching columns for selected individual.")

    return (
        T_norm_indiv,
        norm_dat_values if index_scheme == 2 else norm_dat,
        me_med_nd
    )
