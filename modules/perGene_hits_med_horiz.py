import numpy as np
from .pergene_med_horiz import pergene_med_horiz
from .general_volcano import general_volcano

def pergene_hits_med_horiz(alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz,
                             T_vert, cond, control_type):
    """
    Aggregate per-gRNA data to per-gene medians and identify gene-level hits.

    Parameters
    ----------
    alf : float
        FDR significance threshold for hit calling.
    sfdr_corr : float
        Pseudocount added to FDR before -log10 (avoids log10(0)).
    thr_lfch, thr_lfcd : float
        Right / left critical LFC thresholds (enrichment / depletion), raw LFC scale.
    thr_lfchz, thr_lfcdz : float
        Right / left critical thresholds on the Z-corrected scale.
    T_vert : pandas.DataFrame
        Per-gRNA table: [gRNA, gene, lfc(_<rep>)..., Z_<control_type>_lfc(_<rep>)...,
        Q(_<rep>)...].
    cond : str
        Condition label (passed through to the hit-calling step).
    control_type : str
        Control category used in the Z column names (e.g. 'Zero-expressed gene'),
        needed to locate the 'Z_<control_type>_lfc' columns.

    Returns
    -------
    T_lfc_z_q_med : pandas.DataFrame
        Per-gene medians [genes, median_LFC, median_Z, median_q].
    T_lfc_z_q_me : pandas.DataFrame
        Per-gene means [genes, mean_LFC, mean_Z, mean_q].
    T_LFC, T_Q, T_Z : pandas.DataFrame
        Per-gene LFC, Q, and Z summary tables.
    indha, indda : list of int
        Indices of enriched / depleted genes on the median-LFC scale.
    indhaz, inddaz : list of int
        Indices of enriched / depleted genes on the median-Z scale.
    """
    T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z = pergene_med_horiz(T_vert, control_type)

    if not T_lfc_z_q_med.empty:
        fdr   = T_lfc_z_q_med['median_q'].values
        genes = T_lfc_z_q_med['genes'].values

    else:
        indha = indda = indhaz = inddaz = []

    return T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z
