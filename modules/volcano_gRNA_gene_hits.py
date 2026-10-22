import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

from .plot_volcano import plot_volcano
from .perGene_med_horiz import perGene_med_horiz
from .collapse_table import collapse_table


def volcano_gRNA_gene_hits(alf, sfdr_corr, thr_lfchz, thr_lfcdz,
                           T_vert, plot_name, output_dir,
                           norm_method="pooled"):
    """
    Draw per-gRNA and per-gene Z volcano plots side by side and return their tables.

    Panels: per-gRNA Z (left), per-gene Z (right). Replicate columns are
    collapsed to a per-gRNA median before plotting; per-gene values come from
    perGene_med_horiz.

    Parameters
    ----------
    alf : float
        q-value significance threshold for hit calling.
    sfdr_corr : float
        Pseudocount added to q-value before -log10 (avoids log10(0) on the y-axis).
    thr_lfchz, thr_lfcdz : float
        Right / left critical thresholds on the normalised lfc scale.
    T_vert : pandas.DataFrame
        Per-gRNA table: [gRNA, gene, gene_type, lfc(_<rep>)...,
        normalised_lfc(_<rep>)..., Q(_<rep>)...].
    plot_name : str
        Contrast name, used in the title and output filename.
    output_dir : str
        Directory where the figure is saved.
    norm_method : {'pooled', 'sample'}, optional
        Gene-level aggregation mode, passed to perGene_med_horiz.

    Returns
    -------
    T_gRNA : pandas.DataFrame
        Per-gRNA table [gRNA, gene, normalised_lfc, q]: replicate-median
        normalised LFC and q-value.
    T_gene : pandas.DataFrame
        Per-gene table [gene, normalised_lfc, q]: per-gene median normalised
        LFC and q-value; empty if there are no genes.
    """
    grnas = T_vert['gRNA']
    genes = T_vert['gene']

    # collapse replicates to one value per gRNA
    scoreZ = collapse_table(T_vert, prefix="normalised_lfc", exact="normalised_lfc")
    fdr    = collapse_table(T_vert, prefix="Q",              exact="Q")

    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    plt.suptitle(f"{plot_name}", fontsize=20)

    # --- Left: per-gRNA Z ---
    T_gRNA = plot_volcano(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz,
        scoreZ, fdr, genes, ax=axs[0], grna=grnas)
    axs[0].set_xlabel('Normalised LFC', fontsize=14)
    axs[0].set_title('Volcano per gRNA', fontsize=16)

    # --- per-gene aggregation ---
    T_lfc_z_q_med = perGene_med_horiz(T_vert, norm_method)

    T_gene = pd.DataFrame()
    if not T_lfc_z_q_med.empty:
        # --- Right: per-gene Z ---
        T_gene = plot_volcano(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med['median_normalised_LFC'].values,
            T_lfc_z_q_med['median_q'].values, T_lfc_z_q_med['genes'].values,
            ax=axs[1])
        axs[1].set_xlabel('Normalised LFC', fontsize=14)
        axs[1].set_title('Volcano per gene', fontsize=16)
    else:
        axs[1].axis('off')

    # one shared legend instead of two
    for ax in axs.flat:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc='center left', bbox_to_anchor=(1.0, 0.5),
               frameon=True, edgecolor='black', facecolor='white',
               framealpha=1.0, fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"volcano_{plot_name}.pdf"),
                bbox_inches="tight")
    plt.close()

    # output tables: raw q-value (without the plotting pseudocount)
    T_gRNA = (T_gRNA.drop(columns='fdr_corr')
                    .rename(columns={'score': 'normalised_lfc'})
                    .assign(q=np.asarray(fdr, dtype=float)))
    if not T_gene.empty:
        T_gene = (T_gene.drop(columns='fdr_corr')
                        .rename(columns={'score': 'normalised_lfc'})
                        .assign(q=T_lfc_z_q_med['median_q'].values))

    return T_gRNA, T_gene