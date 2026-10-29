import pandas as pd
import numpy as np

def make_tables_Z(target, Z_t):
    """Build a per-replicate gRNA table with LFC, control-calibrated Z, and Q.

    Parameters
    ----------
    target : pandas.DataFrame
        Target gRNA table [gRNA, gene, gene_type, lfc_<r>, ..., p_<r>, ...].
    Z_t : numpy.ndarray
        (n_gRNA, n_reps) per-replicate Z-scores, column order matching the
        lfc_* columns of ``target``.

    Returns
    -------
    T_wt : pandas.DataFrame
        [gRNA, gene, gene_type, lfc_<r>, normalised_lfc_<r>, p_<r>, ...] (p = raw
        empirical p-value),
        one set of columns per replicate.

    Raises
    ------
    ValueError
        If the number of columns in ``Z_t`` does not match the number of
        lfc columns in ``target``.
    """
    Z_t = np.asarray(Z_t, dtype=float)
    if Z_t.ndim == 1:
        Z_t = Z_t[:, None]

    lfc_cols = [c for c in target.columns if c.startswith('lfc')]

    if Z_t.shape[1] != len(lfc_cols):
        raise ValueError(f"Z_t has {Z_t.shape[1]} cols but found "
                         f"{len(lfc_cols)} lfc columns: {lfc_cols}")

    T_wt = pd.DataFrame({
        'gRNA':      target.iloc[:, 0].values,
        'gene':      target.iloc[:, 1].values,
        'gene_type': target.iloc[:, 2].values,
    })

    for k, lc in enumerate(lfc_cols):
        rep_suffix = lc[len('lfc'):]          # '' for bare 'lfc', '_F' for 'lfc_F'

        T_wt[f'lfc{rep_suffix}']                  = target[lc].values
        T_wt[f'normalised_lfc{rep_suffix}'] = Z_t[:, k]

        qc = f'p{rep_suffix}'
        if qc in target.columns:
            T_wt[qc] = target[qc].values

    return T_wt