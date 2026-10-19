import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

from .general_volcano import general_volcano
from .perGene_med_horiz import perGene_med_horiz
from .collapse_table import collapse_table


def volcano_gRNA_gene_hits(alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz,
                           T_vert, plot_name, control_type, output_dir,
                           norm_method="pooled"):
    """Draw a 2x2 grid of volcano plots and return the gRNA- and gene-level tables.

    Panels: per-gRNA LFC, per-gRNA Z (top row); per-gene LFC, per-gene Z
    (bottom row). Replicate columns are collapsed to a per-gRNA median before
    plotting; per-gene values come from perGene_med_horiz.

    Parameters
    ----------
    alf : float
        FDR significance threshold for hit calling.
    sfdr_corr : float
        Pseudocount added to FDR before -log10 (avoids log10(0) on the y-axis).
    thr_lfch, thr_lfcd : float
        Right / left critical thresholds on the raw LFC scale.
    thr_lfchz, thr_lfcdz : float
        Right / left critical thresholds on the Z-corrected scale.
    T_vert : pandas.DataFrame
        Per-gRNA table: [gRNA, gene, gene_type, lfc(_<rep>)...,
        normalised_lfc(_<rep>)..., Q(_<rep>)...].
    plot_name : str
        Contrast name, used in the title and output filename.
    control_type : str
        Control category the Z columns were calibrated against.
    output_dir : str
        Directory where the figure is saved.
    norm_method : {'pooled', 'sample'}, optional
        Gene-level aggregation mode, passed to perGene_med_horiz.

    Returns
    -------
    T_gRNA : pandas.DataFrame
        Per-gRNA volcano table [gRNA, gene, score, fdr_corr] from the Z panel.
    T_gene : pandas.DataFrame
        Per-gene volcano table [gene, score, fdr_corr] from the Z panel;
        empty if there are no genes.
    """
    grnas = T_vert['gRNA']
    genes = T_vert['gene']

    # collapse replicates to one value per gRNA
    scoreL = collapse_table(T_vert, prefix="lfc",            exact="lfc")
    scoreZ = collapse_table(T_vert, prefix="normalised_lfc", exact="normalised_lfc")
    fdr    = collapse_table(T_vert, prefix="Q",              exact="Q")

    fig, axs = plt.subplots(2, 2, figsize=(14, 12))
    plt.suptitle(f"{plot_name}", fontsize=20)

    # --- Row 0, Col 0: per-gRNA LFC ---
    general_volcano(alf, sfdr_corr, thr_lfch, thr_lfcd,
                    scoreL, fdr, genes, ax=axs[0, 0], grna=grnas)
    axs[0, 0].set_xlabel('LFC', fontsize=14)
    axs[0, 0].set_title('Volcano per gRNA', fontsize=16)

    # --- Row 0, Col 1: per-gRNA Z ---
    T_gRNA = general_volcano(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz,
        scoreZ, fdr, genes, ax=axs[0, 1], grna=grnas)
    axs[0, 1].set_xlabel('Z-corrected LFC', fontsize=14)
    axs[0, 1].set_title('Volcano per gRNA', fontsize=16)

    # --- per-gene aggregation ---
    (T_lfc_z_q_med, T_lfc_z_q_me,
     T_LFC, T_Q, T_Z) = perGene_med_horiz(T_vert, control_type, norm_method)

    T_gene = pd.DataFrame()
    if not T_lfc_z_q_med.empty:
        genes_gene = T_lfc_z_q_med['genes'].values
        fdr_gene   = T_lfc_z_q_med['median_q'].values

        # --- Row 1, Col 0: per-gene LFC ---
        general_volcano(alf, sfdr_corr, thr_lfch, thr_lfcd,
                        T_lfc_z_q_med['median_LFC'].values, fdr_gene, genes_gene,
                        ax=axs[1, 0])
        axs[1, 0].set_xlabel('LFC', fontsize=14)
        axs[1, 0].set_title('Volcano per gene', fontsize=16)

        # --- Row 1, Col 1: per-gene Z ---
        T_gene = general_volcano(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med['median_normalised_LFC'].values, fdr_gene, genes_gene,
            ax=axs[1, 1])
        axs[1, 1].set_xlabel('Z-corrected LFC', fontsize=14)
        axs[1, 1].set_title('Volcano per gene', fontsize=16)
    else:
        axs[1, 0].axis('off')
        axs[1, 1].axis('off')

    # one shared legend instead of four
    for ax in axs.flat:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    handles, labels = axs[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc='center left', bbox_to_anchor=(1.0, 0.5),
               frameon=True, edgecolor='black', facecolor='white',
               framealpha=1.0, fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"volcano_{plot_name}.pdf"),
                bbox_inches="tight")
    plt.close()

    return T_gRNA, T_gene