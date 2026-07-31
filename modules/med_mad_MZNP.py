import numpy as np
from scipy import stats

def med_mad_MZNP(LFC_t):
    """
    Compute median, MAD, modified Z-scores, mean, std, and mode for a vector.
    Implements Leys 2013 robust Z-score method.

    Parameters
    ----------
    LFC_t : array_like
        all replicate LFC columns.

    Returns
    -------
    med_mad : list
        [median, adjusted MAD]
    MZ2 : ndarray
        Modified Z-scores (normalized by adjusted MAD).
    MZ : ndarray
        Modified Z-scores (Leys version, scaled by MAD).
    Z : ndarray
        Standard Z-scores (mean/std normalization).
    me_sd : list
        [mean, standard deviation]
    mod : float
        Mode of the input vector.
    """

    LFC_t = np.array(LFC_t)
    med = np.median(LFC_t)
    
    mode_result = stats.mode(LFC_t, nan_policy='omit')
    mod = mode_result.mode.item() 
    
    b = 1.4826  # scale factor for MAD assuming normal distribution

    mad = np.median(np.abs(LFC_t - med))
    mad2 = b * mad

    med_mad = [med, mad2]

    me_sd = [np.mean(LFC_t), np.std(LFC_t, ddof=1)]

    # Modified Z-scores
    if mad == 0:
        MZ = np.zeros_like(LFC_t)
    else:
        MZ = 0.6745 * (LFC_t - med) / mad

    if mad2 == 0:
        MZ2 = np.zeros_like(LFC_t)
    else:
        MZ2 = (LFC_t - med) / mad2

    # Standard Z-scores
    mean_val = me_sd[0]
    std_val = me_sd[1]
    if std_val == 0:
        Z = np.zeros_like(LFC_t)
    else:
        Z = (LFC_t - mean_val) / std_val

    return med_mad, MZ2, MZ, Z, me_sd, mod
