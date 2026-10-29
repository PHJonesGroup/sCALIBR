import plotly.graph_objs as go
import os
from plotly.subplots import make_subplots
from .plot_volcano_interactive import plot_volcano_interactive
from .gene_fisher import gene_fisher
from .collapse_table import collapse_table


def volcano_gRNA_gene_hits_interactive(
    alf, sfdr_corr, thr_lfchz, thr_lfcdz,
    T_vert, plot_name, output_dir, norm_method="pooled", fisher_scale_factor=1/3
):
    """
    Draw interactive per-gRNA and per-gene Z volcano plots and write them as HTML.

    Panels: per-gRNA Z (left), per-gene Z (right). Replicate columns are collapsed to a per-gRNA median before
    plotting; per-gene values come from gene_fisher (median normalised LFC,
    Fisher-combined p, BH across genes).

    Parameters
    ----------
    alf : float
        Significance threshold for hit calling: on the raw gRNA p-value in the
        gRNA panel, on the BH q-value in the gene panel.
    sfdr_corr : float
        Pseudocount added to p / q before -log10 (avoids log10(0)).
    thr_lfchz, thr_lfcdz : float
        Right / left critical thresholds on the normalised LFC scale.
    T_vert : pandas.DataFrame
        Per-gRNA table: [gRNA, gene, gene_type, lfc(_<rep>)...,
        normalised_lfc(_<rep>)..., p(_<rep>)...]; p = raw empirical gRNA p-value.
    plot_name : str
        Contrast name, used in the title and output filename.
    output_dir : str
        Directory where the HTML figure is saved.
    norm_method : {'pooled', 'sample'}, optional
        Gene-level aggregation mode, passed to gene_fisher.
    fisher_scale_factor : float, optional
        Scale factor for Fisher's method, passed to gene_fisher.

    Returns
    -------
    None
    """
    grnas = T_vert['gRNA']
    genes = T_vert['gene']

    scoreZ = collapse_table(T_vert, prefix="normalised_lfc", exact="normalised_lfc")
    pval   = collapse_table(T_vert, prefix="p_",             exact="p")

    fig = make_subplots(rows=1, cols=2, subplot_titles=(
        "Volcano per gRNA",
        "Volcano per gene",
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

    # --- Left: per-gRNA Z ---
    _, _, _, T_ds, T_hs, T_gRNA = plot_volcano_interactive(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz, scoreZ, pval, genes, grnas)
    add_volcano(T_gRNA, T_hs, T_ds, row=1, col=1, label="gRNA", show_legend=True)

    # --- per-gene aggregation ---
    T_gene = gene_fisher(T_vert, norm_method, fisher_scale_factor)

    if not T_gene.empty:
        # --- Right: per-gene Z ---
        _, _, _, T_ds2, T_hs2, T_gene2 = plot_volcano_interactive(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_gene['normalised_lfc'].values,
            T_gene['q'].values, T_gene['gene'].values)
        add_volcano(T_gene2, T_hs2, T_ds2, row=1, col=2, label="gene", show_legend=False)
    else:
        print("No per-gene data available, skipping gene volcano plot")

    for c, y_title in ((1, "-log10(p-value)"), (2, "-log10(q-value)")):
        fig.update_xaxes(title_text="Normalised LFC", row=1, col=c)
        fig.update_yaxes(title_text=y_title, row=1, col=c)

    fig.update_layout(title=f"{plot_name}", showlegend=True)
    fig.write_html(
        os.path.join(output_dir, f"volcano_plot_interactive_{plot_name}.html"),
        include_plotlyjs="cdn")