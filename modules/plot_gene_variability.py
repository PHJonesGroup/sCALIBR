import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os

def plot_gene_variability(T_norm, output_dir, baseline, cond1,
                          target_type="GOI", rep_pairs=None):
    """
    Interactive per-gene CPM boxplots: two boxes per gene (baseline vs cond1),
    coloured by condition. Points coloured by replicate name (consistent across
    both boxes). Genes sorted by overall median CPM.

    T_norm : [sgRNA_name, gene, gene_type, <cond>_<rep> count columns...]
    """
    gt = T_norm["gene_type"].astype(str).str.strip()
    goi = T_norm[gt == target_type].copy()
    if goi.empty:
        raise ValueError(f"No rows with gene_type == '{target_type}'")

    id_cols = list(goi.columns[:2]) + ["gene_type"]           # sgRNA_name, gene, gene_type

    # map each count column -> (condition, replicate)
    def parse_col(c):
        for cond in (baseline, cond1):
            if c == cond:                       # no-rep case: column IS the condition
                return cond, "rep1"
            if c.startswith(f"{cond}_"):        # rep case: <cond>_<rep>
                return cond, c[len(cond) + 1:]
        return None, None                       # not one of our two conditions

    col_info = {c: parse_col(c) for c in goi.columns if c not in id_cols}
    count_cols = [c for c, (cond, _) in col_info.items() if cond is not None]
    if not count_cols:
        raise ValueError(f"No columns for {baseline} or {cond1}. Cols: {goi.columns.tolist()}")

    # long format with condition + rep columns
    long = goi.melt(id_vars=id_cols, value_vars=count_cols,
                    var_name="sample", value_name="CPM")
    long["CPM"] = pd.to_numeric(long["CPM"], errors="coerce")
    long = long.dropna(subset=["CPM"])
    long["condition"] = long["sample"].map(lambda c: col_info[c][0])
    long["rep"]       = long["sample"].map(lambda c: col_info[c][1])

    # order genes by overall median CPM
    order = long.groupby("gene")["CPM"].median().sort_values().index.tolist()

    gene_to_x = {g: i for i, g in enumerate(order)}
    long["x"] = long["gene"].map(gene_to_x)

    fig = go.Figure()

    # box offset: baseline shifted left, cond1 right (matches boxmode='group' spacing)
    offset = 0.2
    cond_offset = {baseline: -offset, cond1: +offset}
    cond_colors = {baseline: "rgba(31,119,180,0.5)", cond1: "rgba(255,127,14,0.5)"}

    # --- one box per condition, placed at the offset numeric x ---
    for cond in (baseline, cond1):
        sub = long[long["condition"] == cond]
        fig.add_trace(go.Box(
            x=sub["x"] + cond_offset[cond], y=sub["CPM"],
            name=cond, marker_color=cond_colors[cond],
            width=0.35, boxpoints=False, showlegend=True
        ))

    # --- points coloured by rep, offset to sit inside their condition's box ---
    reps = sorted(long["rep"].unique())
    palette = ["#111111", "#e41a1c", "#377eb8", "#4daf4a",
               "#984ea3", "#ff7f00", "#a65628", "#f781bf"]
    rep_color = {r: palette[i % len(palette)] for i, r in enumerate(reps)}

    rng = np.random.default_rng(0)                 # small jitter so points don't stack
    for r in reps:
        sub = long[long["rep"] == r]
        xpos = (sub["x"]
                + sub["condition"].map(cond_offset)
                + rng.uniform(-0.05, 0.05, len(sub)))   # jitter within the box
        fig.add_trace(go.Scatter(
                    x=xpos, y=sub["CPM"], mode="markers",
                    marker=dict(color=rep_color[r], size=10, opacity=0.85),
                    name=f"rep {r}", text=sub["sgRNA_name"] + " | " + sub["condition"],
                    hoverinfo="text+y"
                ))
    # restore gene names on the x-axis (we used numeric positions)
    fig.update_layout(
        title=f"Per-gene CPM variability (GOI): {baseline} vs {cond1}",
        xaxis_title="Gene", yaxis_title="Normalised counts (CPM)",
        xaxis=dict(tickmode="array",
                   tickvals=list(gene_to_x.values()),
                   ticktext=list(gene_to_x.keys())),
        yaxis_title_standoff=10,
        legend_title="Condition / Replicate",
    )

    n_box = 2                              # baseline + cond1 boxes
    n_pts = len(reps)                      # one scatter trace per rep
    # visibility patterns: boxes always on; points toggled
    boxes_on_points_on  = [True]*n_box + [True]*n_pts
    boxes_on_points_off = [True]*n_box + [False]*n_pts

    fig.update_layout(
        title=f"Per-gene CPM variability (GOI): {baseline} vs {cond1}",
        xaxis_title="Gene", yaxis_title="Normalised counts (CPM)",
        xaxis=dict(tickmode="array",
                   tickvals=list(gene_to_x.values()),
                   ticktext=list(gene_to_x.keys())),
        legend_title="Condition / Replicate",
        updatemenus=[dict(
            type="buttons", direction="right",
            x=1.0, y=1.12, xanchor="right",
            buttons=[
                dict(label="Points: on",  method="update",
                     args=[{"visible": boxes_on_points_on}]),
                dict(label="Points: off", method="update",
                     args=[{"visible": boxes_on_points_off}]),
            ],
        )],
    )
    out = os.path.join(output_dir, "gene_variability_boxplots.html")
    fig.write_html(out)
    return out