import numpy as np

def separate_genes_by_grna_count(T_vert, ggenes, gene_names, gg, ind_gn, rep):
    """
    Buckets every gene into a dictionary keyed by its gRNA count.

    Parameters
    ----------
    T_vert : pandas.DataFrame
        Table with ID columns followed by count columns
    ggenes : list of str
        List of unique gene names
    gene_names : lsit of str
        Gene name for every gRNA row
    gg : list of int
        Number of gRNAs for each gene in ``ggenes``
    ind_gn : list of int
        List of start indices (1-based) for each group in ``gene_names`` / T_vert.
    rep : int
        Number of replicates per condition (used to slice the count columns).

    Returns
    -------
    groups : dict
        Maps gRNA-count -> dict with keys:
        - 'genes'      : list of str, gene names with that gRNA count
        - 'gene_names' : list of str, per-gRNA gene names for those genes
        - 'ind'        : list of int, 0-based row indices into T_vert
        - 'wt'         : ndarray (n_gRNA x n_count_cols), stacked count data
    """
    wt = T_vert.iloc[:, 3:(3 + 2*rep)].to_numpy()

    groups = {}   # count -> dict with genes, gene_names, indices, data
    for i in range(len(gg)):
        c = gg[i]
        start = ind_gn[i] - 1
        g = groups.setdefault(c, {"genes": [], "gene_names": [], "ind": [], "wt": []})
        g["genes"].append(ggenes[i])
        g["gene_names"].extend(gene_names[start:start + c])
        g["ind"].extend(range(start, start + c))
        g["wt"].append(wt[start:start + c, :])

    for c, g in groups.items():
        g["wt"] = np.vstack(g["wt"]) if g["wt"] else np.empty((0, wt.shape[1]))

    return groups