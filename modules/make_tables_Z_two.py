import pandas as pd
import numpy as np

def make_tables_Z_two(target, Z_t):
    """
    Build a per-replicate gRNA table with LFC, control-calibrated Z, and Q,
    one set of columns per replicate.

    Parameters
    ----------
    target : pd.DataFrame with [gRNA, gene, lfc_<r1>, ..., Q_<r1>, ...]
    Z_t    : 2-D ndarray (n_gRNA x n_reps), per-rep Z-scores aligned to the
             lfc_* columns (same column order)

    Returns
    -------
    T_wt : DataFrame [gRNA, gene, lfc_<r>, Z_zGE_<r>, Q_<r>, ...] per replicate
    """
    Z_t = np.asarray(Z_t, dtype=float)
    if Z_t.ndim == 1:
        Z_t = Z_t[:, None]

    lfc_cols = [c for c in target.columns if c.startswith('lfc')]
    q_cols   = [c for c in target.columns if c.startswith('Q')]

    if Z_t.shape[1] != len(lfc_cols):
        raise ValueError(f"Z_t has {Z_t.shape[1]} cols but found {len(lfc_cols)} lfc columns")

    T_wt = pd.DataFrame({
        'gRNA': target.iloc[:, 0].values,
        'gene': target.iloc[:, 1].values,
    })

    for k, lc in enumerate(lfc_cols):
        rep = lc.replace('lfc_', '')
        T_wt[f'lfc_{rep}']    = target[lc].values
        T_wt[f'Z_zGE_{rep}']  = Z_t[:, k]
        # match the Q column for this rep if it exists
        qc = f'Q_{rep}'
        if qc in target.columns:
            T_wt[qc] = target[qc].values

    return T_wt