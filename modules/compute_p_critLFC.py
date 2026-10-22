import numpy as np
import matplotlib.pyplot as plt
import os
from .q_val_frequentist_critical import q_val_frequentist_critical

def compute_p_critLFC(alf, binn, hiss):
    """
    Compute the p-curve for the control LFC distribution

    Parameters
    ----------
    alf   : float
        Significance level for the critical-value / p-curve calculation.
    binn : numpy.ndarray
        Bin edges (left edges) for the LFC histogram.
    hiss : numpy.ndarray
        Histogram counts (single series), same length as ``binn``.

    Returns
    -------
    bin_p : numpy.ndarray
        Two-column array [bin, p]: bin edges with the corresponding p-curve.
    """
    p = q_val_frequentist_critical(alf, binn, hiss)
    bin_p = np.column_stack((binn, p))

    return bin_p