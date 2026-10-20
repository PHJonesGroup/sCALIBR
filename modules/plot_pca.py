import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from .sample_group import sample_group

def _label_points(ax, xy, labels, colors, min_sep_frac=0.12):
    """
    Annotate points; draw a headless leader only when the label is offset.
    """
    xs, ys = xy[:, 0], xy[:, 1]
    span = max(np.ptp(xs), np.ptp(ys)) or 1.0
    close_thresh = span * min_sep_frac
    cx, cy = xs.mean(), ys.mean()

    label_xy = []
    for i, ((x, y), lab, col) in enumerate(zip(xy, labels, colors)):
        # is any other point within the crowding threshold?
        d = np.hypot(xs - x, ys - y)
        d[i] = np.inf
        crowded = d.min() < close_thresh

        if crowded:
            dx, dy = x - cx, y - cy
            norm = np.hypot(dx, dy) or 1.0
            ox, oy = dx / norm * span * 0.16, dy / norm * span * 0.16
            ax.annotate(lab, xy=(x, y), xytext=(x + ox, y + oy),
                        fontsize=8, color=col,
                        ha="left" if ox >= 0 else "right",
                        va="bottom" if oy >= 0 else "top",
                        arrowprops=dict(arrowstyle="-", lw=0.5, color="grey",
                                        shrinkA=0, shrinkB=3))
            label_xy.append((x + ox, y + oy))
        else:
            ax.annotate(lab, xy=(x, y), fontsize=8, color=col,
                        textcoords="offset points", xytext=(5, 4),
                        ha="left", va="bottom")
            label_xy.append((x, y))

    # expand limits so offset labels and their text stay inside the axes
    allx = np.concatenate([xs, [p[0] for p in label_xy]])
    ally = np.concatenate([ys, [p[1] for p in label_xy]])
    padx = (np.ptp(allx) or 1.0) * 0.22
    pady = (np.ptp(ally) or 1.0) * 0.18
    ax.set_xlim(allx.min() - padx, allx.max() + padx)
    ax.set_ylim(ally.min() - pady, ally.max() + pady)


def plot_pca(norm_dat, output_dir, log_transform=True):
    """Plot sample PCA side by side: all gRNAs, and within each gene_type.

    Samples are the observations and gRNAs the features, so each point is one
    sample. Points are coloured by sample group — the sample name with any
    trailing '_R<N>' replicate tag removed, so replicates of the same group
    share a colour and a single legend entry. The left panel uses all gRNAs;
    the right panel runs a separate PCA within each gene_type, keeping the
    sample-group colours and encoding gene_type as marker shape.

    Parameters
    ----------
    norm_dat : pandas.DataFrame
        Normalised counts: three ID columns (sgRNA_name, gene, gene_type),
        then one count column per sample.
    output_dir : str
        Directory where plot is saved.
    log_transform : bool, optional
        Apply log2(x + 1) before PCA (default True), so a few high-count gRNAs
        don't dominate the components.

    Returns
    -------
    None
    """
    count_cols = list(norm_dat.columns[3:])
    if not count_cols:
        raise ValueError(f"No count columns after the first 3. "
                         f"Columns: {norm_dat.columns.tolist()}")

    # sample -> group (replicate tag stripped); one colour per group
    groups = [sample_group(s) for s in count_cols]
    uniq_groups = sorted(set(groups))
    palette = ["#377EB8", "#E41A1C", "#1D9E75", "#7F77DD",
               "#888780", "#BA7517", "#D4537E", "#D85A30"]
    group_color = {g: palette[i % len(palette)] for i, g in enumerate(uniq_groups)}
    sample_color = [group_color[g] for g in groups]

    def pca_of(sub):
        """PCA on a gRNA x sample block -> (coords, explained variance %)."""
        X = sub[count_cols].astype(float).to_numpy()
        X = np.log2(X + 1) if log_transform else X
        X = X.T                                   # samples as rows
        X = X - X.mean(axis=0, keepdims=True)     # centre each gRNA
        n_comp = min(2, X.shape[0], X.shape[1])
        p = PCA(n_components=n_comp)
        coords = p.fit_transform(X)
        if coords.shape[1] == 1:
            coords = np.column_stack([coords, np.zeros(len(coords))])
        var = list(p.explained_variance_ratio_ * 100) + [0.0]
        return coords, var

    fig, axs = plt.subplots(1, 2, figsize=(13, 5.5))

    # ---- panel 1: all gRNAs ----
    coords, var = pca_of(norm_dat)
    axs[0].scatter(coords[:, 0], coords[:, 1], color=sample_color, s=40, zorder=3)
    _label_points(axs[0], coords, count_cols, sample_color)
    axs[0].set_title("All gRNAs", fontsize=13)
    axs[0].set_xlabel(f"PC1 ({var[0]:.1f}%)")
    axs[0].set_ylabel(f"PC2 ({var[1]:.1f}%)")

    # ---- panel 2: separate PCA per gene_type, shape = gene_type ----
    gt = norm_dat["gene_type"].astype(str).str.strip()
    cats = sorted(gt.unique())
    markers = ["o", "s", "^", "D", "v", "P", "X"]
    cat_marker = {c: markers[i % len(markers)] for i, c in enumerate(cats)}

    for cat in cats:
        sub = norm_dat[gt == cat]
        if len(sub) < 2:
            continue
        coords_c, var_c = pca_of(sub)
        axs[1].scatter(coords_c[:, 0], coords_c[:, 1],
                       color=sample_color, marker=cat_marker[cat],
                       s=40, zorder=3)
    axs[1].set_title("gRNAs within each gene type", fontsize=13)
    axs[1].set_xlabel(f"PC1 ({var_c[0]:.1f}%)")
    axs[1].set_ylabel(f"PC2 ({var_c[1]:.1f}%)")
    axs[1].legend(handles=[Line2D([], [], marker=cat_marker[c], linestyle="",
                                  color="#444441", label=c, markersize=5)
                           for c in cats if (gt == c).sum() >= 2],
                  fontsize=9, frameon=True, edgecolor="black", title="gene type")

    sample_leg = axs[1].legend(
        handles=[Line2D([], [], marker="o", linestyle="",
                        color=group_color[g], label=g, markersize=7)
                 for g in uniq_groups],
        title="sample", fontsize=9, frameon=True, edgecolor="black",
        loc="upper left", bbox_to_anchor=(1.02, 1.0))
    axs[1].add_artist(sample_leg)          # keep it when the second legend is added

    axs[1].legend(
        handles=[Line2D([], [], marker=cat_marker[c], linestyle="",
                        color="#444441", label=c, markersize=7)
                 for c in cats if (gt == c).sum() >= 2],
        title="gene type", fontsize=9, frameon=True, edgecolor="black",
        loc="lower left", bbox_to_anchor=(1.02, 0.0))
    
    for ax in axs:
        ax.axhline(0, color="grey", linewidth=0.5, zorder=1)
        ax.axvline(0, color="grey", linewidth=0.5, zorder=1)
        ax.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "pca_qc.pdf"), bbox_inches="tight")
    plt.close()