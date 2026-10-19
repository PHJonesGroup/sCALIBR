import numpy as np
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

def perGene_med_horiz(T_vert, control_type, norm_method="pooled"):
    """Per-gene summaries of LFC, normalised LFC, and Q.

    Two aggregation modes:
      'pooled' — one median over all (guides x replicates) values for the gene.
      'sample' — median across the gene's guides within each replicate, then
                 the median of those per-replicate medians.

    Parameters
    ----------
    T_vert : pandas.DataFrame
        Per-gRNA table: [gRNA, gene, gene_type, lfc(_<rep>)...,
        normalised_lfc(_<rep>)..., Q(_<rep>)...].
    control_type : str
        Control category the Z columns were calibrated against (labels only).
    norm_method : {'pooled', 'sample'}
        Aggregation mode, as above.

    Returns
    -------
    T_lfc_z_q_med, T_lfc_z_q_me : pandas.DataFrame
        Per-gene medians / means [genes, *_LFC, *_Z, *_q].
    T_LFC, T_Q, T_Z : pandas.DataFrame
        Per-gene summary tables with median, mean, std (and n_gRNA for LFC).
    """
    blocks = {
        'LFC':            _value_cols(T_vert, 'lfc'),
        'normalised_LFC': _value_cols(T_vert, 'normalised_lfc'),
        'Q':              _value_cols(T_vert, 'Q'),
    }

    med, mean, std = {}, {}, {}
    for key, cs in blocks.items():
        med[key], mean[key], std[key] = summarise(T_vert, cs, norm_method)

    genes = med['LFC'].index.values
    n_grna = T_vert.groupby('gene')['gRNA'].nunique().reindex(genes).values
    eps = 1e-6

    T_lfc_normalised_q_med = pd.DataFrame({
        'genes':      genes,
        'median_LFC': med['LFC'].values,
        'median_normalised_LFC':   med['normalised_LFC'].values,
        'median_q':   med['Q'].values + eps,
    })
    T_lfc_normalised_q_me = pd.DataFrame({
        'genes':    genes,
        'mean_LFC': mean['LFC'].values,
        'mean_normalised_LFC':   mean['normalised_LFC'].values,
        'mean_q':   mean['Q'].values + eps,
    })
    T_LFC = pd.DataFrame({'genes': genes,
                          'median_LFC': med['LFC'].values,
                          'mean_LFC': mean['LFC'].values,
                          'std_LFC': std['LFC'].values,
                          'n_gRNA': n_grna})
    T_normalised   = pd.DataFrame({'genes': genes,
                          'median_normalised_LFC': med['normalised_LFC'].values,
                          'mean_normalised_LFC': mean['normalised_LFC'].values,
                          'std_normalised_LFC': std['normalised_LFC'].values})
    T_Q   = pd.DataFrame({'genes': genes,
                          'median_q': med['Q'].values + eps,
                          'mean_q': mean['Q'].values + eps,
                          'std_q': std['Q'].values})

    return T_lfc_normalised_q_med, T_lfc_normalised_q_me, T_LFC, T_Q, T_normalised