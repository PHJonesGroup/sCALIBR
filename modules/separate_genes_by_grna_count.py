def separate_genes_by_grna_count(gg, ind_gn):
    """
    Buckets every gene into a dictionary keyed by its gRNA count.

    Parameters
    ----------
    gg : list of int
        Number of gRNAs for each gene (one entry per run of identical gene names).
    ind_gn : list of int
        List of start indices (1-based) for each gene in the count table.

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
