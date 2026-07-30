import matplotlib.pyplot as plt
import pandas as pd
import os

def normalise_prop(ssc, raw_ind, output_dir):
    """
    Normalise counts per column to CPM: (count_j / sum(counts)) * 1e6 + ssc offset.

    Parameters
    ----------
    ssc: int
        Scalar offset added after scaling, to avoid zeros
    raw_ind: pandas.DataFrame
        DataFrame with 3 ID columns (gRNA, gene, gene_type) followed by count columns
    output_dir : str
        Directory where the plot 'gRNA_counts_normalisation.png' is saved.

    Returns
    -------
    norm_df: pandas.DataFrame
        DataFrame [gRNA, gene, <normalised count columns>]
    """
    # all data columns = everything after the 3 ID columns
    data_cols = list(raw_ind.columns[3:])
    if not data_cols:
        raise ValueError("No count columns found after the 3 ID columns")

    datc = raw_ind[data_cols].astype(float).to_numpy()

    sn1 = datc.sum(axis=0)                          # raw column sums
    sc  = datc.sum(axis=0)
    norm_dat = ssc + 1_000_000 * datc / sc          # CPM + offset
    sn  = norm_dat.sum(axis=0)                       # normalised column sums

    count_table = pd.DataFrame({
        "sample_name": data_cols,
        "sum_raw_counts": sn1,
        "sum_normalised_counts": sn})
    count_table.to_csv(
        os.path.join(output_dir, "gRNA_counts_normalisation_data.csv"), index=False)

    # QC plot: raw vs normalised column totals
    plt.figure(figsize=(max(8, 0.6 * len(data_cols)), 6))

    plt.subplot(2, 1, 1)
    plt.bar(range(len(sn1)), sn1, width=0.3)
    plt.title('Raw gRNA counts')
    plt.xticks([])
    plt.ylabel("Sum of counts")

    plt.subplot(2, 1, 2)
    plt.bar(range(len(sn)), sn, width=0.3)
    plt.title('Normalised gRNA counts')
    plt.xticks(range(len(sn)), data_cols, rotation=45, ha='right')
    plt.ylabel("Normalised proportion")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gRNA_counts_normalisation.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    norm_df = raw_ind.iloc[:, :3].copy()
    norm_df[data_cols] = norm_dat

    return norm_df