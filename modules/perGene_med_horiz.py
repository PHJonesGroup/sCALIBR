import pandas as pd
from .summarise import summarise


def _value_cols(T_vert, prefix):
    """Column names for one quantity: the exact name if present, else all
    replicate columns starting with `prefix`."""
    if prefix in T_vert.columns:
        return [prefix]
    cols = [c for c in T_vert.columns if c.startswith(prefix)]
    if not cols:
        raise KeyError(f"No columns for '{prefix}'. Columns: {T_vert.columns.tolist()}")
    return cols

def perGene_med_horiz(T_vert, norm_method="pooled"):
    """Per-gene medians of LFC, normalised LFC, and Q.

    Two aggregation modes:
      'pooled' — one median over all (guides x replicates) values for the gene.
      'sample' — median across the gene's guides within each replicate, then
                 the median of those per-replicate medians.

    Parameters
    ----------
    T_vert : pandas.DataFrame
        Per-gRNA table: [gRNA, gene, gene_type, lfc(_<rep>)...,
        normalised_lfc(_<rep>)..., Q(_<rep>)...].
    norm_method : {'pooled', 'sample'}
        Aggregation mode, as above.

    Returns
    -------
    T_lfc_normalised_q_med : pandas.DataFrame
        Per-gene medians [genes, median_LFC, median_normalised_LFC, median_q].
    """
    blocks = {
        'LFC':            _value_cols(T_vert, 'lfc'),
        'normalised_LFC': _value_cols(T_vert, 'normalised_lfc'),
        'Q':              _value_cols(T_vert, 'Q'),
    }

    med = {key: summarise(T_vert, cs, norm_method)[0] for key, cs in blocks.items()}

    eps = 1e-6
    T_lfc_normalised_q_med = pd.DataFrame({
        'genes':      med['LFC'].index.values,
        'median_LFC': med['LFC'].values,
        'median_normalised_LFC':   med['normalised_LFC'].values,
        'median_q':   med['Q'].values + eps,
    })

    return T_lfc_normalised_q_med
