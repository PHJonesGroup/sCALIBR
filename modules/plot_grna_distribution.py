import matplotlib.pyplot as plt
import os
from collections import Counter

def plot_grna_distribution(gg, output_dir):
    """
    Plot the distribution of gRNA counts per gene as a bar chart.
    
    Parameters
    ----------
    gg : list of int
        Number of gRNAs for each gene (one entry per gene).
    output_dir : str
        Directory where the plot 'gene_distribution_gRNA.png' is saved.
    
    Returns
    -------
    None
    """

    counts = Counter(gg)                       # {gRNA_count: number_of_genes}
    d = max(gg)
    x = list(range(1, d + 1))
    heights = [counts[k] for k in x]

    plt.figure()
    plt.bar(x, heights, edgecolor='black')
    plt.xticks(x)
    plt.xlabel('number of gRNA per gene')
    plt.ylabel('number of genes')
    plt.title('number of gRNA per gene distribution')
    plt.savefig(os.path.join(output_dir, "gene_distribution_gRNA.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    