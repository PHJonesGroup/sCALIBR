import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42


def normalise_total(ssc, raw_ind, output_dir):
    """Normalise counts to the mean library size across samples.

    Each column is rescaled so its total equals the average total across all
    samples: norm_ij = round(ssc + mean(column sums) * count_ij / sum(column j)).
    This corrects for sequencing-depth differences while keeping values on the
    original count scale, so normalised numbers stay comparable in magnitude to
    the raw counts.

    Parameters
    ----------
    ssc : int
        Small scalar offset added after scaling, to avoid zeros.
    raw_ind : pandas.DataFrame
        3 ID columns (gRNA, gene, gene_type) followed by count columns.
    output_dir : str
        Directory where the QC plot and summary CSV are saved.

    Returns
    -------
    norm_df : pandas.DataFrame
        [gRNA, gene, gene_type, <normalised count columns>].
    """
    data_cols = list(raw_ind.columns[3:])
    if not data_cols:
        raise ValueError("No count columns found after the 3 ID columns")

    datc = raw_ind[data_cols].astype(float).to_numpy()

    sc = datc.sum(axis=0)                       # raw column sums (library sizes)
    if np.any(sc == 0):
        raise ValueError(f"Zero total counts in column(s): "
                         f"{[c for c, s in zip(data_cols, sc) if s == 0]}")

    factor_av = sc.mean()                       # average library size
    norm_dat = np.round(ssc + factor_av * datc / sc)   # scale each column to factor_av

    sn = norm_dat.sum(axis=0)                   # normalised column sums

    count_table = pd.DataFrame({
        "sample_name": data_cols,
        "sum_raw_counts": sc,
        "sum_normalised_counts": sn})
    count_table.to_csv(
        os.path.join(output_dir, "gRNA_counts_normalisation_data.csv"), index=False)

    plt.figure(figsize=(max(8, 0.6 * len(data_cols)), 6))

    plt.subplot(2, 1, 1)
    plt.bar(range(len(sc)), sc, width=0.3)
    plt.title('Raw gRNA counts per sample')
    plt.xticks([])
    plt.ylabel("Sum of counts")

    plt.subplot(2, 1, 2)
    plt.bar(range(len(sn)), sn, width=0.3)
    plt.title('Normalised gRNA counts per sample')
    plt.xticks(range(len(sn)), data_cols, rotation=45, ha='right')
    plt.ylabel("Sum of counts")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gRNA_counts_normalisation.pdf"),
                bbox_inches="tight")
    plt.close()

    norm_df = raw_ind.iloc[:, :3].copy()
    norm_df[data_cols] = norm_dat

    return norm_df