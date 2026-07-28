import plotly.graph_objs as go
import os
import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
from .general_volcano_interactive import general_volcano_interactive
from .pergene_hits_med_horiz import pergene_hits_med_horiz
from .collapse_table import collapse_table

def volcano_grna_gene_hits_interactive(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, baseline, cond, control_type, output_dir
):
    """
    Draw an interactive 2x2 grid of volcano plots and return hit/depleted indices.

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
    None
    """
    grnas = T_vert['gRNA']
    genes  = T_vert['gene']
    # collapse reps to median per gRNA (or use the single column if no reps)
    scoreL = collapse_table(T_vert, prefix="lfc",                     exact="lfc")
    scoreZ = collapse_table(T_vert, prefix=f"Z_{control_type}_lfc",   exact=f"Z_{control_type}_lfc")
    fdr    = collapse_table(T_vert, prefix="Q",                       exact="Q")


    # 2x2 grid: row1 = per gRNA (LFC, Z), row2 = per gene (LFC, Z)
    fig = make_subplots(rows=2, cols=2, subplot_titles=(
        f"Volcano per gRNA",
        f"Volcano per gRNA",
        f"Volcano per gene",
        f"Volcano per gene",
    ))

    # helper to add the three trace layers (all / hits / depleted)
    def add_volcano(T_all, T_hs, T_ds, row, col, label, show_legend):
        fig.add_trace(go.Scatter(
            x=T_all["score"], y=T_all["LPV"], mode="markers",
            marker=dict(color="lightgray", size=6),
            text=T_all[label], name="", hoverinfo="text+x+y", showlegend=False
        ), row=row, col=col)
        if not T_hs.empty:
            fig.add_trace(go.Scatter(
                x=T_hs["score"], y=T_hs["LPV"], mode="markers",
                marker=dict(color="red", size=7),
                text=T_hs[label], name="Hits", hoverinfo="text+x+y",
                showlegend=show_legend
            ), row=row, col=col)
        if not T_ds.empty:
            fig.add_trace(go.Scatter(
                x=T_ds["score"], y=T_ds["LPV"], mode="markers",
                marker=dict(color="blue", size=7),
                text=T_ds[label], name="Depleted", hoverinfo="text+x+y",
                showlegend=show_legend
            ), row=row, col=col)

    # --- Row 1, Col 1: per-gRNA LFC ---
    _, _, _, T_dsL, T_hsL, T_gRNA_L = general_volcano_interactive(
        alf, sfdr_corr, thr_lfch, thr_lfcd, scoreL, fdr, cond, genes, grnas
    )
    add_volcano(T_gRNA_L, T_hsL, T_dsL, row=1, col=1, label="gRNA", show_legend=True)  # legend from here only

    # --- Row 1, Col 2: per-gRNA Z ---
    LPV, indha, indda, T_ds, T_hs, T_gRNA = general_volcano_interactive(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz, scoreZ, fdr, cond, genes, grnas
    )
    add_volcano(T_gRNA, T_hs, T_ds, row=1, col=2, label="gRNA", show_legend=False)

    # --- Per-gene stats ---
    (
        T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z,
        indha_gene, indda_gene, indhaz, inddaz
    ) = pergene_hits_med_horiz(
        alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond, control_type
    )

    if not T_lfc_z_q_med.empty:
        fdr_gene   = T_lfc_z_q_med.iloc[:, 3].values
        genes_gene = T_lfc_z_q_med.iloc[:, 0].values

        # --- Row 2, Col 1: per-gene LFC ---
        _, _, _, T_ds2, T_hs2, T_gene2 = general_volcano_interactive(
            alf, sfdr_corr, thr_lfch, thr_lfcd,
            T_lfc_z_q_med.iloc[:, 1].values, fdr_gene, cond, genes_gene
        )
        add_volcano(T_gene2, T_hs2, T_ds2, row=2, col=1, label="gene", show_legend=False)

        # --- Row 2, Col 2: per-gene Z ---
        _, _, _, T_ds3, T_hs3, T_gene3 = general_volcano_interactive(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med.iloc[:, 2].values, fdr_gene, cond, genes_gene
        )
        add_volcano(T_gene3, T_hs3, T_ds3, row=2, col=2, label="gene", show_legend=False)
    else:
        print("No per-gene data available, skipping gene volcano plots")

    # axis titles
    fig.update_xaxes(title_text="LFC", row=1, col=1)
    fig.update_xaxes(title_text="Z-corrected LFC", row=1, col=2)
    fig.update_xaxes(title_text="LFC", row=2, col=1)
    fig.update_xaxes(title_text="Z-corrected LFC", row=2, col=2)
    for r in (1, 2):
        for c in (1, 2):
            fig.update_yaxes(title_text="-log10(FDR)", row=r, col=c)

    fig.update_layout(title=f"{cond} vs {baseline}", showlegend=True)
    fig.write_html(os.path.join(output_dir, f"volcano_plot_interactive_{cond}_vs_{baseline}.html"))

    return 1