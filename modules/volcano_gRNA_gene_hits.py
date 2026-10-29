import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

from .plot_volcano import plot_volcano
from .gene_fisher import gene_fisher
from .collapse_table import collapse_table


def volcano_gRNA_gene_hits(alf, sfdr_corr, thr_lfchz, thr_lfcdz,
                           T_vert, plot_name, output_dir,
                           norm_method="pooled", fisher_scale_factor=1/3):
    """
    Draw per-gRNA and per-gene Z volcano plots side by side and return their tables.

    Panels: per-gRNA Z (left), per-gene Z (right). Replicate columns are
    collapsed to a per-gRNA median before plotting; per-gene values come from
    gene_fisher (median normalised LFC, Fisher-combined p, BH across genes).
    In 'sample' mode one extra figure per sample, volcano_<plot_name>_<rep>.pdf,
    shows that sample's gRNA values and per-sample gene results.

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
        Gene-level aggregation mode, passed to gene_fisher.
    fisher_scale_factor : float, optional
        Scale factor for Fisher's method, passed to gene_fisher.

    Returns
    -------
    T_gRNA : pandas.DataFrame
        Per-gRNA table [gRNA, gene, normalised_lfc, q]: replicate-median
        normalised LFC and q-value.
    T_gene : pandas.DataFrame
        Per-gene table from gene_fisher: [gene, normalised_lfc, fisher_p, q, ...]
        (per-replicate columns too in 'sample' mode); empty if there are no genes.
    """
    grnas = T_vert['gRNA']
    genes = T_vert['gene']

    # collapse replicates to one value per gRNA
    scoreZ = collapse_table(T_vert, prefix="normalised_lfc", exact="normalised_lfc")
    fdr    = collapse_table(T_vert, prefix="Q",              exact="Q")

    # --- per-gene aggregation ---
    T_gene = gene_fisher(T_vert, norm_method, fisher_scale_factor)

    # final volcano (all samples)
    T_gRNA = _volcano_figure(alf, sfdr_corr, thr_lfchz, thr_lfcdz,
                             scoreZ, fdr, genes, grnas, T_gene,
                             'normalised_lfc', 'q', plot_name,
                             os.path.join(output_dir, f"volcano_{plot_name}.pdf"))

    # sample mode: one volcano per sample, from that sample's gRNA and gene values
    if norm_method == "sample" and not T_gene.empty:
        for sfx in [c[len('Q'):] for c in T_vert.columns if c.startswith('Q_')]:
            _volcano_figure(alf, sfdr_corr, thr_lfchz, thr_lfcdz,
                            T_vert[f'normalised_lfc{sfx}'], T_vert[f'Q{sfx}'], genes, grnas,
                            T_gene, f'normalised_lfc{sfx}', f'q{sfx}',
                            f"{plot_name} ({sfx.lstrip('_')})",
                            os.path.join(output_dir, f"volcano_{plot_name}{sfx}.pdf"))

    # output tables: raw q-value (without the plotting pseudocount)
    T_gRNA = (T_gRNA.drop(columns='fdr_corr')
                    .rename(columns={'score': 'normalised_lfc'})
                    .assign(q=np.asarray(fdr, dtype=float)))

    return T_gRNA, T_gene


def _volcano_figure(alf, sfdr_corr, thr_lfchz, thr_lfcdz, scoreZ, fdr, genes, grnas,
                    T_gene, gene_lfc_col, gene_q_col, title, path):
    """Two-panel volcano (per gRNA, per gene) saved to path; returns the gRNA table."""
    fig, axs = plt.subplots(1, 2, figsize=(14, 6))
    plt.suptitle(title, fontsize=20)

    # --- Left: per-gRNA Z ---
    T_gRNA = plot_volcano(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz,
        scoreZ, fdr, genes, ax=axs[0], grna=grnas)
    axs[0].set_xlabel('Normalised LFC', fontsize=14)
    axs[0].set_title('Volcano per gRNA', fontsize=16)

    if not T_gene.empty:
        # --- Right: per-gene Z ---
        plot_volcano(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_gene[gene_lfc_col].values,
            T_gene[gene_q_col].values, T_gene['gene'].values,
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
    plt.savefig(path, bbox_inches="tight")
    plt.close()
    return T_gRNA
