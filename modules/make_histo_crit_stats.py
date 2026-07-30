import pandas as pd
import numpy as np
from .make_histo_LFC import make_histo_LFC
from .compute_p_critLFC import compute_p_critLFC
from .med_mad_MZNP import med_mad_MZNP
from .plot_histograms import plot_histograms

def make_histo_crit_stats(alf, st, en, step, T_control, baseline, cond1, control, output_dir):
    """
    Histograms, critical LFC thresholds, and stats for control genes
    (baseline vs cond1), pooling all replicate LFC columns.
    
    Parameters
    ----------
    alf : float
        Significance level for the critical-value / p-curve calculation.
    st, en : float
        Lower / upper bounds of the histogram (LFC) range.
    step : float
        Histogram bin width.
    T_control : pandas.DataFrame
        Control table: [gRNA, gene, lfc_<rep>, ...] (per-replicate LFC columns).
    baseline, cond1 : str
        Condition labels.
    control : str
        Control category label (e.g. 'Zero-expressed gene'), used in labels.
    output_dir : str
        Directory where the plot is saved.

    Returns
    -------
    bin_ : numpy.ndarray
        Histogram bin edges (left edges).
    his : numpy.ndarray
        Counts per bin of the pooled control LFCs.
    perc : numpy.ndarray
        Percentages per bin.
    crit_LR : numpy.ndarray
        [left, right] critical LFC thresholds.
    bin_p : numpy.ndarray
        [bin, p] array (bin edges with the control p-curve).
    med_mad : numpy.ndarray
        [median, MAD] of the pooled control LFC.
    me_sd : numpy.ndarray
        [mean, SD] of the pooled control LFC.
    mod : float
        Mode (or mode-like statistic) of the pooled control LFC.
    MZ : numpy.ndarray
        Modified (median/MAD-based) Z-scores, one per pooled LFC value.
    Z : numpy.ndarray
        Standard (mean/SD-based) Z-scores, one per pooled LFC value.
    n : int
        Number of finite pooled LFC values used.
    """
    # 1. Pool all replicate LFC columns
    lfc_cols = list(T_control.columns[3:])
    if not lfc_cols:
        raise ValueError("No LFC columns found after gRNA, gene")
    LFC = T_control[lfc_cols].astype(float).values.ravel()   # flatten reps into one vector
    LFC = LFC[np.isfinite(LFC)]                           # drop NaN/inf

    # 2. Histogram of the pooled control LFCs
    bin_, his, perc = make_histo_LFC(step, LFC, st, en)

    # 3. Critical LFC thresholds & p-curve
    bin_p, crit_LR = compute_p_critLFC(alf, bin_, his, baseline, cond1, control, output_dir)

    # 4. Robust stats, Z / MZ on the pooled control LFCs
    med_mad, MZZ, MZ, Z, me_sd, mod = med_mad_MZNP(LFC)

    # 5. Diagnostic plot
    _, _, perc_z  = make_histo_LFC(step, Z,  st, en)
    _, _, perc_mz = make_histo_LFC(step, MZ, st, en)
    plot_histograms(
        bin_, perc, perc_z, perc_mz,
        title=f'Distribution of {control} gRNAs: {cond1} vs {baseline}',
        color='b', xlab='LFC',
        save_name=f'distri_LFC_{control}_{cond1}_vs_{baseline}', output_dir=output_dir)

    return bin_, his, perc, crit_LR, bin_p, med_mad, me_sd, mod, MZ, Z, len(LFC)