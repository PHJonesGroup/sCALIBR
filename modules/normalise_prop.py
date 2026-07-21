import matplotlib.pyplot as plt
import os 

def normalise_prop(ssc, raw_ind, rep, output_dir):
    """
    Normalise counts per column as [ count_j / sum(counts) * 1,000,000 ] + ssc offset.
    ssc is a small offset added to avoid zeros (e.g., 1 or 0.01).

    Parameters:
    - ssc: scalar offset added to all values after scaling
    - raw_ind: pandas DataFrame, with numeric count columns at positions 3:6 

    Returns:
    - norm_dat: numpy array with normalised data
    """

    # Extract columns 3 to 6 
    datc = raw_ind.iloc[:, 2:(2 + 2*rep)].to_numpy() 

    sn1 = datc.sum(axis=0)  # sum per column (raw counts)

    # Normalisation: (count / sum_counts) * 1,000,000 + ssc offset
    sc = datc.sum(axis=0)
    norm_dat = ssc + 1_000_000 * datc / sc

    sn = norm_dat.sum(axis=0)  # sum per column (normalised)
    
    labels = raw_ind.columns[2:(2 + 2*rep)]

    # Plot bar charts
    plt.figure(figsize=(8, 6))

    plt.subplot(2, 1, 1)
    plt.bar(range(len(sn1)), sn1, width=0.3)
    plt.title('Raw gRNA counts')
    plt.xticks([])
    plt.ylabel("Sum of counts")

    plt.subplot(2, 1, 2)
    plt.bar(range(len(sn)), sn, width=0.3)
    plt.title('Normalised gRNA counts')
    plt.xticks(range(len(sn)), labels)
    plt.ylabel("Normalised proportion")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "gRNA_counts_normalisation"), dpi=300, bbox_inches="tight")
    plt.close()

    data_cols = raw_ind.columns[2:2 + norm_dat.shape[1]]
    norm_df = raw_ind.iloc[:, :2].copy()              # sgRNA_name, gene
    norm_df[list(data_cols)] = norm_dat               # attach normalized data under original names
    return norm_df