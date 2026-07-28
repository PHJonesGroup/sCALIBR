def collapse_table(T_vert, prefix, exact=None):
    """
    Collapse a per-replicate LFC to a single value per gRNA.
    - one replicate only: use LFC column directly.
    - more than one replicate: median across all columns starting with `prefix`.
    
    Parameters
    ----------
    T_vert : pandas.DataFrame
        Table containing lfc values, either as one column (``exact``) or as
        several replicate columns sharing ``prefix``.
    prefix : str
        Column-name prefix identifying the replicate columns (e.g. 'lfc', 'Q').
    exact : str, optional
        Exact column name to use if present (the no-replicate case). If None,
        only the prefix match is attempted.

    Returns
    -------
    pandas.Series
        One value per gRNA (float): the single column, or the row-wise median
        across the replicate columns.

    Raises
    ------
    KeyError
        If no column matches ``exact`` and none start with ``prefix``.
    """
    if exact is not None and exact in T_vert.columns:
        return T_vert[exact].astype(float)
    cols = [c for c in T_vert.columns if c.startswith(prefix)]
    if not cols:
        raise KeyError(f"No columns for '{exact or prefix}'. Columns: {T_vert.columns.tolist()}")
    return T_vert[cols].astype(float).median(axis=1)