import numpy as np
import matplotlib.pyplot as plt
import os
from .q_val_frequentist_critical import q_val_frequentist_critical

def compute_p_critLFC(alf, binn, hiss, baseline, cond, control, output_dir):
    """
    Compute the p-curve and critical LFC values for a single LFC distribution
    (contrast cond vs baseline).

    Parameters
    ----------
    alf   : float
        Significance level for the critical-value / p-curve calculation.
    binn : numpy.ndarray
        Bin edges (left edges) for the LFC histogram.
    hiss : numpy.ndarray
        Histogram counts (single series), same length as ``binn``.
    baseline, cond : str
        Condition labels.
    control : str
        Control category label (e.g. 'Zero-expressed gene'), used in labels.
    output_dir : str
        Directory where the plot is saved.

    Returns
    -------
    bin_p : numpy.ndarray
        Two-column array [bin, p]: bin edges with the corresponding p-curve.
    crit_LR : numpy.ndarray
        [left_critical, right_critical] critical LFC thresholds.
    """
    p, cl, cr, bin_pz, med_LFCp, hiss4p = q_val_frequentist_critical(alf, binn, hiss)
    crit_LR = np.array([cl, cr])

    bin_p = np.column_stack((binn, p))

    return bin_p, crit_LR