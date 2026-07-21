import plotly.graph_objs as go
import os
from plotly.subplots import make_subplots
from .general_volcano_interactive import general_volcano_interactive
from .perGene_4_hits_med_horiz import perGene_4_hits_med_horiz

def volcano_gRNA_gene_hits_wt_interactive(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond, output_dir
):
    genes = T_vert['gene']
    scoreZ = T_vert['Z_zGE_LFC']
    fdr = T_vert['Q']

    # --- create subplot layout with 3 plots ---
    fig = make_subplots(rows=1, cols=3, subplot_titles=(
        f"Volcano per gRNA (Z): {cond}",
        f"Volcano per gene (LFC): {cond}",
        f"Volcano per gene (Z): {cond}"
    ))

    # 1. Volcano per gRNA
    LPV, indha, indda, T_ds, T_hs, T_gRNA = general_volcano_interactive(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz, scoreZ, fdr, cond, genes
    )
    # All points
    fig.add_trace(go.Scatter(
        x=T_gRNA["score"], y=T_gRNA["LPV"],
        mode="markers", marker=dict(color="lightgray", size=6),
        text=T_gRNA["gene"], name="", hoverinfo="text+x+y",showlegend=False
    ), row=1, col=1)
    # Hits
    if not T_hs.empty:
        fig.add_trace(go.Scatter(
            x=T_hs["score"], y=T_hs["LPV"],
            mode="markers", marker=dict(color="red", size=7),
            text=T_hs["gene"], name="Hits", hoverinfo="text+x+y"
        ), row=1, col=1)
    # Depleted
    if not T_ds.empty:
        fig.add_trace(go.Scatter(
            x=T_ds["score"], y=T_ds["LPV"],
            mode="markers", marker=dict(color="blue", size=7),
            text=T_ds["gene"], name="Depleted", hoverinfo="text+x+y"
        ), row=1, col=1)

    # 2 & 3. Per gene (LFC and Z)
    (
        T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z,
        indha_gene, indda_gene, indhaz, inddaz
    ) = perGene_4_hits_med_horiz(
        alf, sfdr_corr, thr_lfch, thr_lfcd,
        thr_lfchz, thr_lfcdz, T_vert, cond, plot=False
    )

    if not T_lfc_z_q_med.empty:
        fdr_gene = T_lfc_z_q_med.iloc[:, 3].values
        genes_gene = T_lfc_z_q_med.iloc[:, 0].values

        # --- Volcano per gene LFC ---
        LPV2, indh2, indd2, T_ds2, T_hs2, T_gene2 = general_volcano_interactive(
            alf, sfdr_corr, thr_lfch, thr_lfcd,
            T_lfc_z_q_med.iloc[:, 1].values, fdr_gene, cond, genes_gene
        )
        fig.add_trace(go.Scatter(
            x=T_gene2["score"], y=T_gene2["LPV"],
            mode="markers", marker=dict(color="lightgray", size=6),
            text=T_gene2["gene"], name="", hoverinfo="text+x+y",showlegend=False
        ), row=1, col=2)
        if not T_hs2.empty:
            fig.add_trace(go.Scatter(
                x=T_hs2["score"], y=T_hs2["LPV"],
                mode="markers", marker=dict(color="red", size=7),
                text=T_hs2["gene"], name="", hoverinfo="text+x+y",showlegend=False
            ), row=1, col=2)
        if not T_ds2.empty:
            fig.add_trace(go.Scatter(
                x=T_ds2["score"], y=T_ds2["LPV"],
                mode="markers", marker=dict(color="blue", size=7),
                text=T_ds2["gene"], name="", hoverinfo="text+x+y",showlegend=False
            ), row=1, col=2)

        # --- Volcano per gene Z ---
        LPV3, indh3, indd3, T_ds3, T_hs3, T_gene3 = general_volcano_interactive(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med.iloc[:, 2].values, fdr_gene, cond, genes_gene
        )
        fig.add_trace(go.Scatter(
            x=T_gene3["score"], y=T_gene3["LPV"],
            mode="markers", marker=dict(color="lightgray", size=6),
            text=T_gene3["gene"], name="", hoverinfo="text+x+y",showlegend=False
        ), row=1, col=3)
        if not T_hs3.empty:
            fig.add_trace(go.Scatter(
                x=T_hs3["score"], y=T_hs3["LPV"],
                mode="markers", marker=dict(color="red", size=7),
                text=T_hs3["gene"], name="", hoverinfo="text+x+y",showlegend=False
            ), row=1, col=3)
        if not T_ds3.empty:
            fig.add_trace(go.Scatter(
                x=T_ds3["score"], y=T_ds3["LPV"],
                mode="markers", marker=dict(color="blue", size=7),
                text=T_ds3["gene"], name="", hoverinfo="text+x+y",showlegend=False
            ), row=1, col=3)
    else:
        print("No per-gene data available, skipping gene volcano plots")

    # Final layout
    fig.update_xaxes(title_text="Score", row=1, col=1)
    fig.update_yaxes(title_text="-log10(FDR)", row=1, col=1)
    fig.update_layout(title=f"Interactive Volcano Plots: {cond}", showlegend=True)

    fig.write_html(os.path.join(output_dir, f"volcano_plot_interactive_{cond}.html"))

    return T_gRNA, T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z, indha, indda, indhaz, inddaz
