import numpy as np
import matplotlib.pyplot as plt
from .make_histo_LFC import make_histo_LFC

def computeZ(st, en, step, LFC_t, me_sd_z):
    """
    Compute per-replicate Z-scores for targets against the pooled-control null.
    
    Parameters
    ----------
    st, en : float
        Lower / upper bounds of the histogram (LFC) range.
    step : float
        Histogram bin width.
    LFC_t : numpy.ndarray
        Target LFCs, shape (n_gRNA, n_reps); a 1-D array is treated as one column.
    me_sd_z : sequence of float
        [mean, SD] of the pooled control LFC (used as the standardisation null).

    Returns
    -------
    Z_t : numpy.ndarray
        Z-scored LFCs, shape (n_gRNA, n_reps), one column per replicate.
    binn : numpy.ndarray
        Histogram bin edges.
    perc_t : numpy.ndarray
        Pooled percentage histogram of the raw LFC values (all reps combined).
    perc_zt : numpy.ndarray
        Pooled percentage histogram of the Z-scored values (all reps combined).

    Raises
    ------
    ValueError
        If the control SD is zero or non-finite (cannot standardise).
    """
    LFC_t = np.asarray(LFC_t, dtype=float)
    if LFC_t.ndim == 1:
        LFC_t = LFC_t[:, None]                 # treat single series as one column

    mu, sd = me_sd_z[0], me_sd_z[1]

    if sd == 0 or not np.isfinite(sd):
        raise ValueError("Control SD is zero or invalid — cannot compute Z-scores")

    print(LFC_t)
    # per-replicate Z: each rep's LFC standardized against the pooled-control null
    Z_t = (LFC_t - mu) / sd                     # (n_gRNA x n_reps)

    # pooled histograms (all reps' values combined) for the diagnostic plot
    lfc_flat = LFC_t[np.isfinite(LFC_t)]
    z_flat   = Z_t[np.isfinite(Z_t)]
    binn, _, perc_t  = make_histo_LFC(step, lfc_flat, st, en)
    _,    _, perc_zt = make_histo_LFC(step, z_flat,   st, en)

    return Z_t, binn, perc_t, perc_zt
