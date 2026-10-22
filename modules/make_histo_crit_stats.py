import pandas as pd
import numpy as np
from .make_histo_LFC import make_histo_LFC
from .compute_p_critLFC import compute_p_critLFC
from .med_mad_MZNP import med_mad_MZNP
from .plot_histograms import plot_histograms

def make_histo_crit_stats(alf, st, en, step, T_control, baseline, cond1, plot_name, control, output_dir):
    """
    Histogram, p-curve, and mean/SD for control genes
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
    bin_p : numpy.ndarray
        [bin, p] array (bin edges with the control p-curve).
    me_sd : numpy.ndarray
        [mean, SD] of the pooled control LFC.
    """
    # 1. Pool all replicate LFC columns
    lfc_cols = list(T_control.columns[3:])
    if not lfc_cols:
        raise ValueError("No LFC columns found after gRNA, gene")
    LFC = T_control[lfc_cols].astype(float).values.ravel()      # flatten reps into one vector
    LFC = LFC[np.isfinite(LFC)]                                 # drop NaN/inf

    # 2. Histogram of the pooled control LFCs
    bin_, his, _ = make_histo_LFC(step, LFC, st, en)

    # 3. p-curve
    bin_p = compute_p_critLFC(alf, bin_, his)

    # 4. Mean / SD and Z on the pooled control LFCs
    Z, me_sd = med_mad_MZNP(LFC)

    # 5. Diagnostic plot, on [st, en] widened to fit the normalised LFCs
    #    (bin_/his/perc above stay on [st, en]: they define the p-curve bins)
    z_fin = Z[np.isfinite(Z)]
    st_p, en_p = st, en
    if z_fin.size:
        st_p = min(st, float(np.round((np.floor(z_fin.min() / step) - 1) * step, 10)))
        en_p = max(en, float(np.round((np.ceil(z_fin.max() / step) + 1) * step, 10)))
    bin_plot, _, perc_plot = make_histo_LFC(step, LFC, st_p, en_p)
    _,        _, perc_z    = make_histo_LFC(step, Z,   st_p, en_p)
    plot_histograms(
        bin_plot, perc_plot, perc_z,
        title=f'Distribution of {control} gRNAs: {plot_name}',
        color='b', xlab='LFC',
        save_name=f'distri_LFC_{control}_{plot_name}', output_dir=output_dir)

    return bin_, bin_p, me_sd