import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os
import plotly.io as pio
import json
import re

def _sample_group(name):
    """Strip a trailing _R<N> replicate tag: 'T0_F_R2' -> 'T0_F'."""
    m = re.search(r"_R\d+", name)
    return name[:m.start()] if m else name

def _replicate_tag(name):
    """Extract the trailing _R<N> replicate tag: 'T0_F_R2' -> 'R2'."""
    m = re.search(r"(R\d+)$", name)
    return m.group(1) if m else name

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
    baseline = [_sample_group(b) for b in baseline]
    cond1 = [_sample_group(c) for c in cond1]

    baseline_label = "_".join(baseline)
    cond1_label = "_".join(cond1)

    group_label = {"baseline": baseline_label, "cond1": cond1_label}

    gt = T_norm["gene_type"].astype(str).str.strip()
    goi = T_norm[gt == target_type].copy()
    if goi.empty:
        raise ValueError(f"No rows with gene_type == '{target_type}'")

    id_cols = list(goi.columns[:2]) + ["gene_type"]      # sgRNA_name, gene, gene_type

    # Determine which columns belong to baseline/cond1 based on their
    # sample group (i.e. column name with any trailing _R<N> stripped off),
    # rather than requiring the caller to enumerate every replicate column.
    candidate_cols = [c for c in goi.columns if c not in id_cols]

    col_info = {}
    for c in candidate_cols:
        grp = _sample_group(c)
        if grp in baseline:
            col_info[c] = ("baseline", _replicate_tag(c))
        elif grp in cond1:
            col_info[c] = ("cond1", _replicate_tag(c))

    count_cols = list(col_info.keys())
    if not count_cols:
        raise ValueError(
            f"No columns matched group(s) {baseline} or {cond1}. "
            f"Available columns: {candidate_cols}"
        )

    # long format
    long = goi.melt(id_vars=id_cols, value_vars=count_cols,
                    var_name="sample", value_name="CPM")
    long["CPM"] = pd.to_numeric(long["CPM"], errors="coerce")
    long = long.dropna(subset=["CPM"])
    long["group"] = long["sample"].map(lambda c: col_info[c][0])   # baseline / cond1
    long["sub"]   = long["sample"].map(lambda c: col_info[c][1])   # replicate tag, e.g. R1/R2

    # keep only the top-N most variable genes
    top_n = 100
    sd_by_group = (long.groupby(["gene", "group"])["CPM"]
                       .std()
                       .unstack("group"))                      # columns: baseline, cond1
    gene_variability = sd_by_group.max(axis=1)                 # largest SD across the two groups
    top_genes = gene_variability.sort_values(ascending=False).head(top_n).index
    long = long[long["gene"].isin(top_genes)].copy()

    # order genes by largest SD (most variable first)
    order = gene_variability.loc[top_genes].sort_values(ascending=False).index.tolist()
    gene_to_x = {g: i for i, g in enumerate(order)}
    long["x"] = long["gene"].map(gene_to_x)

    fig = go.Figure()

    offset = 0.2
    group_offset = {"baseline": -offset, "cond1": +offset}
    group_colors = {"baseline": "rgba(31,119,180,0.5)", "cond1": "rgba(255,127,14,0.5)"}
    baseline_label = baseline[0]
    cond1_label = cond1[0]
    group_label  = {"baseline": baseline_label, "cond1": cond1_label}

    # one box per group (baseline / cond1), pooling its replicate columns
    for grp in ("baseline", "cond1"):
        sub = long[long["group"] == grp]
        if sub.empty:
            continue
        fig.add_trace(go.Box(
            x=sub["x"] + group_offset[grp], y=sub["CPM"],
            name=group_label[grp], marker_color=group_colors[grp],
            width=0.35, boxpoints=False, showlegend=True
        ))

    # points coloured by replicate tag, consistent across both boxes
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
            name=str(s), text=ss["sgRNA_name"] + " | " + ss["sample"],
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
    plot_div = pio.to_html(fig, include_plotlyjs="cdn", full_html=False,
                           div_id="genePlot")

    # a gene -> x-position lookup the JS can use to locate/zoom a gene
    import json
    gene_x_json = json.dumps(gene_to_x)

    html = f"""<!DOCTYPE html>
<html>
<head><meta charset="utf-8"></head>
<body>
  <div style="margin:10px 0; font-family:sans-serif;">
    <label>Search gene: </label>
    <input type="text" id="geneSearch" placeholder="e.g. Trp53"
           style="padding:4px; width:200px;">
    <button onclick="clearSearch()">Reset</button>
    <span id="searchMsg" style="margin-left:10px; color:#666;"></span>
  </div>
  {plot_div}

  <script>
    const geneX = {gene_x_json};
    const plot = document.getElementById('genePlot');

    function searchGene() {{
      const q = document.getElementById('geneSearch').value.trim().toLowerCase();
      const msg = document.getElementById('searchMsg');
      if (!q) {{ clearSearch(); return; }}

      // find matching gene(s)
      const matches = Object.keys(geneX).filter(g => g.toLowerCase().includes(q));
      if (matches.length === 0) {{
        msg.textContent = 'no match';
        return;
      }}
      // zoom the x-axis to the matching gene(s)
      const xs = matches.map(g => geneX[g]);
      const lo = Math.min(...xs) - 1;
      const hi = Math.max(...xs) + 1;
      Plotly.relayout(plot, {{'xaxis.range': [lo, hi]}});
      msg.textContent = matches.length + ' match' + (matches.length>1?'es':'') + ': ' + matches.join(', ');
    }}

    function clearSearch() {{
      document.getElementById('geneSearch').value = '';
      document.getElementById('searchMsg').textContent = '';
      Plotly.relayout(plot, {{'xaxis.autorange': true}});
    }}

    // search as you type (debounced lightly)
    document.getElementById('geneSearch').addEventListener('input', searchGene);
  </script>
</body>
</html>"""

    out = os.path.join(output_dir, "gene_variability_boxplots.html")
    with open(out, "w") as f:
        f.write(html)