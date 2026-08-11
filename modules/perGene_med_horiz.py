import numpy as np
import pandas as pd
from .collapse_table import collapse_table

def pergene_med_horiz(T_vert, control_type):
    """
    Per-gene median/mean of LFC, Z, Q, computed from a per-gRNA table that may
    have any number of gRNAs per gene and any number of replicate columns.

    Parameters
    ----------
    T_vert : pandas.DataFrame
        Per-gRNA table: [gRNA, gene, lfc(_<rep>)..., Z_<control_type>_lfc(_<rep>)...,
        Q(_<rep>)...].
    control_type : str
        Control category label (e.g. 'Zero-expressed gene'), used in labels.

    Returns
    -------
    T_lfc_z_q_med : pandas.DataFrame
        [genes, median_LFC, median_Z, median_q] — per-gene medians.
    T_lfc_z_q_me : pandas.DataFrame
        [genes, mean_LFC, mean_Z, mean_q] — per-gene means.
    T_LFC : pandas.DataFrame
        [genes, median_LFC, mean_LFC, std_LFC, n_gRNA] — LFC summary + gRNA count.
    T_Q : pandas.DataFrame
        [genes, median_q, mean_q, std_q] — q-value summary.
    T_Z : pandas.DataFrame
        [genes, median_Z, mean_Z, std_Z] — Z summary.
    """
    # 1. collapse replicates -> one value per gRNA
    per_grna = pd.DataFrame({
        'gene': T_vert['gene'].values,
        'LFC':  collapse_table(T_vert, 'lfc', 'lfc').values,
        'Z':    collapse_table(T_vert, f'Z_{control_type}_lfc', f'Z_{control_type}_lfc').values,
        'Q':    collapse_table(T_vert, 'Q', 'Q').values,
    })

    # 2. aggregate gRNAs -> one row per gene (median + mean + std)
    g = per_grna.groupby('gene')
    med  = g.median(numeric_only=True)
    mean = g.mean(numeric_only=True)
    std  = g.std(ddof=1, numeric_only=True)
    n_grna = g.size().rename('n_gRNA')

    eps = 1e-6
    genes = med.index.values

    T_lfc_z_q_med = pd.DataFrame({
        'genes':     genes,
        'median_LFC': med['LFC'].values,
        'median_Z':   med['Z'].values,
        'median_q':   med['Q'].values + eps,
    })
    T_lfc_z_q_me = pd.DataFrame({
        'genes':    genes,
        'mean_LFC': mean['LFC'].values,
        'mean_Z':   mean['Z'].values,
        'mean_q':   mean['Q'].values + eps,
    })
    T_LFC = pd.DataFrame({'genes': genes, 'median_LFC': med['LFC'].values,
                          'mean_LFC': mean['LFC'].values, 'std_LFC': std['LFC'].values,
                          'n_gRNA': n_grna.values})
    T_Z   = pd.DataFrame({'genes': genes, 'median_Z': med['Z'].values,
                          'mean_Z': mean['Z'].values, 'std_Z': std['Z'].values})
    T_Q   = pd.DataFrame({'genes': genes, 'median_q': med['Q'].values + eps,
                          'mean_q': mean['Q'].values + eps, 'std_q': std['Q'].values})

    return T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z