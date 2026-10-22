import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
import seaborn as sns
import os

def plot_heatmap(T_norm, output_dir, id_cols=None,
                           method="pearson", target_type=None,
                           cluster=True):
    """
    Sample-by-sample similarity heatmap based on normalised counts, saved as PDF.

    Parameters
    ----------
    T_norm : pandas.DataFrame
        Normalised counts; ID/metadata columns first, then count columns
        named '<group>_R<N>'.
    output_dir : str
        Directory where the plot 'sample_similarity_heatmap.pdf' is saved.
    method : str, optional
        Correlation method passed to pandas .corr(): 'pearson', 'spearman',
        or 'kendall' (default 'pearson').
    target_type : str, optional
        If given, restricts the similarity calculation to rows where
        'gene_type' equals this value (requires a 'gene_type' column).
        If None, all rows are used.
    cluster : bool, optional
        If True, reorders samples by hierarchical clustering so similar
        samples are adjacent (default True).

    Returns
    -------
    corr : pandas.DataFrame
        The sample x sample similarity matrix, in the order plotted.
    """
    df = T_norm.copy()

    if target_type is not None:
        if "gene_type" not in df.columns:
            raise ValueError("target_type given but no 'gene_type' column present")
        gt = df["gene_type"].astype(str).str.strip()
        df = df[gt == target_type].copy()
        if df.empty:
            raise ValueError(f"No rows with gene_type == '{target_type}'")

    if id_cols is None:
        sample_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    else:
        sample_cols = [c for c in df.columns if c not in id_cols]

    if len(sample_cols) < 2:
        raise ValueError(f"Need at least 2 sample columns, found: {sample_cols}")

    counts = df[sample_cols].apply(pd.to_numeric, errors="coerce").dropna()

    corr = counts.corr(method=method)

    sample_order = list(corr.columns)
    if cluster:
        dist = 1 - corr
        np.fill_diagonal(dist.values, 0)
        condensed = squareform(dist.values, checks=False)
        link = linkage(condensed, method="average")
        dendro = dendrogram(link, labels=corr.columns.tolist(), no_plot=True)
        sample_order = dendro["ivl"]
        corr = corr.loc[sample_order, sample_order]

    n = len(sample_order)
    fig_size = max(6, 0.4 * n)
    _, ax = plt.subplots(figsize=(fig_size, fig_size))

    sns.heatmap(
        corr,
        ax=ax,
        cmap="RdBu_r",
        vmin=-1, vmax=1,
        square=True,
        cbar_kws={"label": method.capitalize()},
        xticklabels=sample_order,
        yticklabels=sample_order,
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=8)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Sample")
    ax.set_title(f"Sample similarity ({method}{', ' + target_type if target_type else ''})")

    os.makedirs(output_dir, exist_ok=True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "sample_similarity_heatmap.pdf"), bbox_inches="tight")
    plt.close()

    return corr