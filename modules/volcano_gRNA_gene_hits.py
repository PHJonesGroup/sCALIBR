import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
from .general_volcano import general_volcano
from .pergene_hits_med_horiz import pergene_hits_med_horiz
from .collapse_table import collapse_table

def volcano_grna_gene_hits(alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz,
                              T_vert, baseline, cond, control_type, output_dir):
    """
    Draw a 2x2 grid of volcano plots and return hit/depleted indices.

    Panels: per-gRNA LFC, per-gRNA Z (top row); per-gene LFC, per-gene Z (bottom
    row).

    Parameters
    ----------
    alf : float
        Significance level for hit calling.
    sfdr_corr : float
        Pseudocount added to FDR before -log10 (avoids log10(0) on the y-axis).
    thr_lfch, thr_lfcd : float
        Right / left critical LFC thresholds (enrichment / depletion), raw LFC scale.
    thr_lfchz, thr_lfcdz : float
        Right / left critical thresholds on the Z-corrected scale.
    T_vert : pandas.DataFrame
        Per-gRNA table: [gRNA, gene, lfc(_<rep>)..., Z_<control_type>_lfc(_<rep>)...,
        Q(_<rep>)...].
    baseline, cond : str
        Condition labels, used in titles and the output filename.
    control_type : str
        Control category used for the Z columns (e.g. 'Zero-expressed gene'),
        needed to locate the 'Z_<control_type>_lfc' columns.
    output_dir : str
        Directory where the figure is saved.

    Returns
    -------
    T_gRNA : pandas.DataFrame
        Per-gRNA volcano table (from the per-gRNA Z panel).
    T_lfc_z_q_med : pandas.DataFrame
        Per-gene median LFC / Z / q table.
    T_lfc_z_q_me : pandas.DataFrame
        Per-gene mean LFC / Z / q table.
    T_LFC, T_Q, T_Z : pandas.DataFrame
        Per-gene LFC, Q, and Z summary tables.
    indha, indda : array-like
        Indices of enriched / depleted hits, per-gRNA on the Z scale.
    indhaz, inddaz : array-like
        Indices of enriched / depleted hits, per-gene on the Z scale.
    """
    genes  = T_vert['gene']
    # collapse reps to median per gRNA (or use the single column if no reps)
    scoreL = collapse_table(T_vert, prefix="lfc",                     exact="lfc")
    scoreZ = collapse_table(T_vert, prefix=f"Z_{control_type}_lfc",   exact=f"Z_{control_type}_lfc")
    fdr    = collapse_table(T_vert, prefix="Q",                       exact="Q")

    # 2x2 grid: row 0 = per gRNA (LFC, Z), row 1 = per gene (LFC, Z)
    fig, axs = plt.subplots(2, 2, figsize=(14, 12))
    plt.suptitle(f"{cond} vs {baseline}", fontsize=20)

    # --- Row 0, Col 0: per-gRNA LFC ---
    general_volcano(
        alf, sfdr_corr, thr_lfch, thr_lfcd,
        scoreL, fdr, cond, genes, ax=axs[0, 0]
    )
    axs[0, 0].set_xlabel('LFC', fontsize=14)
    axs[0, 0].set_title('Volcano per gRNA', fontsize=16)

    # --- Row 0, Col 1: per-gRNA Z ---
    LPV, indha, indda, T_ds, T_hs, T_gRNA = general_volcano(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz,
        scoreZ, fdr, cond, genes, ax=axs[0, 1]
    )
    axs[0, 1].set_xlabel('Z-corrected LFC', fontsize=14)
    axs[0, 1].set_title('Volcano per gRNA', fontsize=16)

    # --- Per-gene stats ---
    (
        T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z,
        indha_gene, indda_gene, indhaz, inddaz
    ) = pergene_hits_med_horiz(
        alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond, control_type
    )

    if not T_lfc_z_q_med.empty:
        fdr_gene   = T_lfc_z_q_med.iloc[:, 3]
        genes_gene = T_lfc_z_q_med.iloc[:, 0]

        # --- Row 1, Col 0: per-gene LFC ---
        general_volcano(
            alf, sfdr_corr, thr_lfch, thr_lfcd,
            T_lfc_z_q_med.iloc[:, 1].values, fdr_gene.values, cond, genes_gene.values,
            ax=axs[1, 0]
        )
        axs[1, 0].set_xlabel('LFC', fontsize=14)
        axs[1, 0].set_title('Volcano per gene', fontsize=16)

        # --- Row 1, Col 1: per-gene Z ---
        general_volcano(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med.iloc[:, 2].values, fdr_gene.values, cond, genes_gene.values,
            ax=axs[1, 1]
        )
        axs[1, 1].set_xlabel('Z-corrected LFC', fontsize=14)
        axs[1, 1].set_title('Volcano per gene', fontsize=16)
    else:
        axs[1, 0].axis('off')
        axs[1, 1].axis('off')

    # remove the per-subplot legends that general_volcano added
    for ax in axs.flat:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    handles, labels = axs[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc='center left', bbox_to_anchor=(1.0, 0.5),
               frameon=True,
               edgecolor='black',        # box border colour
               facecolor='white',        # box background
               framealpha=1.0,           # opaque
               fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"volcano_gRNA_{cond}_vs_{baseline}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    return T_gRNA, T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z, indha, indda, indhaz, inddaz