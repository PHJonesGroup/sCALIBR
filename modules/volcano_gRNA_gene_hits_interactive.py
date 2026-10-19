import plotly.graph_objs as go
import os
import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
from .general_volcano_interactive import general_volcano_interactive
from .perGene_med_horiz import perGene_med_horiz
from .collapse_table import collapse_table


def volcano_gRNA_gene_hits_interactive(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz,
    T_vert, plot_name, control_type, output_dir, norm_method="pooled"
):
    """Draw an interactive 2x2 grid of volcano plots and write it as HTML.

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
        Directory where the HTML figure is saved.
    norm_method : {'pooled', 'sample'}, optional
        Gene-level aggregation mode, passed to perGene_med_horiz.

    Returns
    -------
    None
    """
    grnas = T_vert['gRNA']
    genes = T_vert['gene']

    scoreL = collapse_table(T_vert, prefix="lfc",            exact="lfc")
    scoreZ = collapse_table(T_vert, prefix="normalised_lfc", exact="normalised_lfc")
    fdr    = collapse_table(T_vert, prefix="Q",              exact="Q")

    fig = make_subplots(rows=2, cols=2, subplot_titles=(
        "Volcano per gRNA (LFC)",
        "Volcano per gRNA (Z)",
        "Volcano per gene (LFC)",
        "Volcano per gene (Z)",
    ))

    def add_volcano(T_all, T_hs, T_ds, row, col, label, show_legend):
        fig.add_trace(go.Scattergl(
            x=T_all["score"], y=T_all["LPV"], mode="markers",
            marker=dict(color="lightgray", size=6),
            text=T_all[label], name="", hoverinfo="text+x+y", showlegend=False
        ), row=row, col=col)
        if not T_hs.empty:
            fig.add_trace(go.Scattergl(
                x=T_hs["score"], y=T_hs["LPV"], mode="markers",
                marker=dict(color="red", size=7),
                text=T_hs[label], name="Hits", hoverinfo="text+x+y",
                showlegend=show_legend
            ), row=row, col=col)
        if not T_ds.empty:
            fig.add_trace(go.Scattergl(
                x=T_ds["score"], y=T_ds["LPV"], mode="markers",
                marker=dict(color="blue", size=7),
                text=T_ds[label], name="Depleted", hoverinfo="text+x+y",
                showlegend=show_legend
            ), row=row, col=col)

    # --- Row 1: per gRNA ---
    _, _, _, T_dsL, T_hsL, T_gRNA_L = general_volcano_interactive(
        alf, sfdr_corr, thr_lfch, thr_lfcd, scoreL, fdr, genes, grnas)
    add_volcano(T_gRNA_L, T_hsL, T_dsL, row=1, col=1, label="gRNA", show_legend=True)

    _, _, _, T_ds, T_hs, T_gRNA = general_volcano_interactive(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz, scoreZ, fdr, genes, grnas)
    add_volcano(T_gRNA, T_hs, T_ds, row=1, col=2, label="gRNA", show_legend=False)

    # --- per-gene aggregation ---
    (T_lfc_z_q_med, T_lfc_z_q_me,
     T_LFC, T_Q, T_Z) = perGene_med_horiz(T_vert, control_type, norm_method)

    if not T_lfc_z_q_med.empty:
        genes_gene = T_lfc_z_q_med['genes'].values
        fdr_gene   = T_lfc_z_q_med['median_q'].values

        _, _, _, T_ds2, T_hs2, T_gene2 = general_volcano_interactive(
            alf, sfdr_corr, thr_lfch, thr_lfcd,
            T_lfc_z_q_med['median_LFC'].values, fdr_gene, genes_gene)
        add_volcano(T_gene2, T_hs2, T_ds2, row=2, col=1, label="gene", show_legend=False)

        _, _, _, T_ds3, T_hs3, T_gene3 = general_volcano_interactive(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med['median_normalised_LFC'].values, fdr_gene, genes_gene)
        add_volcano(T_gene3, T_hs3, T_ds3, row=2, col=2, label="gene", show_legend=False)
    else:
        print("No per-gene data available, skipping gene volcano plots")

    fig.update_xaxes(title_text="LFC", row=1, col=1)
    fig.update_xaxes(title_text="Z-corrected LFC", row=1, col=2)
    fig.update_xaxes(title_text="LFC", row=2, col=1)
    fig.update_xaxes(title_text="Z-corrected LFC", row=2, col=2)
    for r in (1, 2):
        for c in (1, 2):
            fig.update_yaxes(title_text="-log10(FDR)", row=r, col=c)

    fig.update_layout(title=f"{plot_name}", showlegend=True)
    fig.write_html(
        os.path.join(output_dir, f"volcano_plot_interactive_{plot_name}.html"),
        include_plotlyjs="cdn")