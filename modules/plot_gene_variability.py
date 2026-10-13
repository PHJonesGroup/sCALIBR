import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os

def plot_gene_variability(T_norm, output_dir, baseline, cond1,
                          target_type="GOI"):
    """
    Interactive per-gene CPM boxplots: two boxes per gene (baseline vs cond1),
    coloured by condition. Points coloured by replicate name (consistent across
    both boxes). Genes sorted by overall median CPM.

    Parameters
    ----------
    T_norm : pandas.DataFrame
        Normalised counts; ID columns first, then '<condition>_<rep>' count columns.
    output_dir : str
        Directory where the plot 'gene_variability_boxplots.html' is saved.
    baseline : str
        Baseline condition name / column prefix
    cond1: str
        Treatment condition name / column prefix
    target_type : str, optional
        Value in the 'gene_type' column selecting which genes to plot
        (default 'GOI').
    
    Returns
    -------
    None
    """
    if isinstance(baseline, str):
        baseline = [baseline]
    if isinstance(cond1, str):
        cond1 = [cond1]

    gt = T_norm["gene_type"].astype(str).str.strip()
    goi = T_norm[gt == target_type].copy()
    if goi.empty:
        raise ValueError(f"No rows with gene_type == '{target_type}'")

    id_cols = list(goi.columns[:2]) + ["gene_type"]      # sgRNA_name, gene, gene_type

    # map each count column -> (group, sub_condition)
    # group is 'baseline'/'cond1' (for the two boxes); sub is the column itself (F/M)
    col_info = {}
    for c in baseline:
        if c in goi.columns:
            col_info[c] = ("baseline", c)
    for c in cond1:
        if c in goi.columns:
            col_info[c] = ("cond1", c)

    count_cols = list(col_info.keys())
    if not count_cols:
        raise ValueError(f"No columns for {baseline} or {cond1}. Cols: {goi.columns.tolist()}")

    # long format
    long = goi.melt(id_vars=id_cols, value_vars=count_cols,
                    var_name="sample", value_name="CPM")
    long["CPM"] = pd.to_numeric(long["CPM"], errors="coerce")
    long = long.dropna(subset=["CPM"])
    long["group"] = long["sample"].map(lambda c: col_info[c][0])   # baseline / cond1
    long["sub"]   = long["sample"].map(lambda c: col_info[c][1])   # the column (F / M)

    # order genes by overall median CPM
    order = long.groupby("gene")["CPM"].median().sort_values().index.tolist()
    gene_to_x = {g: i for i, g in enumerate(order)}
    long["x"] = long["gene"].map(gene_to_x)

    fig = go.Figure()

    offset = 0.2
    group_offset = {"baseline": -offset, "cond1": +offset}
    group_colors = {"baseline": "rgba(31,119,180,0.5)", "cond1": "rgba(255,127,14,0.5)"}
    group_label  = {"baseline": "baseline", "cond1": "cond1"}

    # one box per group (baseline / cond1), pooling its sub-condition columns
    for grp in ("baseline", "cond1"):
        sub = long[long["group"] == grp]
        if sub.empty:
            continue
        fig.add_trace(go.Box(
            x=sub["x"] + group_offset[grp], y=sub["CPM"],
            name=group_label[grp], marker_color=group_colors[grp],
            width=0.35, boxpoints=False, showlegend=True
        ))

    # points coloured by sub-condition (F / M), offset into their group's box
    subs = sorted(long["sub"].unique())
    palette = ["#111111", "#e41a1c", "#377eb8", "#4daf4a",
               "#984ea3", "#ff7f00", "#a65628", "#f781bf"]
    sub_color = {s: palette[i % len(palette)] for i, s in enumerate(subs)}

    rng = np.random.default_rng(0)
    for s in subs:
        ss = long[long["sub"] == s]
        xpos = (ss["x"]
                + ss["group"].map(group_offset)
                + rng.uniform(-0.05, 0.05, len(ss)))
        fig.add_trace(go.Scattergl(
            x=xpos, y=ss["CPM"], mode="markers",
            marker=dict(color=sub_color[s], size=10, opacity=0.85),
            name=str(s), text=ss["sgRNA_name"] + " | " + ss["sub"],
            hoverinfo="text+y"
        ))

    n_box = 2
    n_pts = len(subs)
    boxes_on_points_on  = [True]*n_box + [True]*n_pts
    boxes_on_points_off = [True]*n_box + [False]*n_pts

    fig.update_layout(
        title=f"Per-gene CPM variability ({target_type})",
        xaxis_title="Gene", yaxis_title="Normalised counts (CPM)",
        xaxis=dict(tickmode="array",
                   tickvals=list(gene_to_x.values()),
                   ticktext=list(gene_to_x.keys())),
        legend_title="Condition / Sub-condition",
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
    fig.write_html(out, include_plotlyjs="cdn")