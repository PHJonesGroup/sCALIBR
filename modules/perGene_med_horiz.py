import numpy as np
import pandas as pd

def perGene_med_horiz(T_vert, control_type):
    """
    Per-gene median/mean of LFC, Z, Q, computed from a per-gRNA table that may
    have any number of gRNAs per gene and any number of replicate columns.

    T_vert : [gRNA, gene, lfc(_<r>)..., Z_<control>_lfc(_<r>)..., Q(_<r>)...]
    """
    def collapse(prefix, exact):
        if exact in T_vert.columns:
            return T_vert[exact].astype(float)
        cols = [c for c in T_vert.columns if c.startswith(prefix)]
        if not cols:
            raise KeyError(f"No columns for {exact or prefix}. Cols: {T_vert.columns.tolist()}")
        return T_vert[cols].astype(float).median(axis=1)

    # 1. collapse replicates -> one value per gRNA
    per_grna = pd.DataFrame({
        'gene': T_vert['gene'].values,
        'LFC':  collapse('lfc', 'lfc').values,
        'Z':    collapse(f'Z_{control_type}_lfc', f'Z_{control_type}_lfc').values,
        'Q':    collapse('Q', 'Q').values,
    })

    # 2. aggregate gRNAs -> one row per gene (median + mean + std), ANY gRNA count
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