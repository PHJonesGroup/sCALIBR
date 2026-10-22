import numpy as np

def separate_genes_by_grna_count(T_vert, ggenes, gene_names, gg, ind_gn):
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

    Returns
    -------
    groups : dict
        Maps gRNA-count -> dict with key:
        - 'ind' : list of int, 0-based row indices into T_vert
    """
    groups = {}
    for i in range(len(gg)):
        c = gg[i]
        start = ind_gn[i] - 1
        groups.setdefault(c, {"ind": []})["ind"].extend(range(start, start + c))

    return groups
