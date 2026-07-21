import numpy as np
from scipy import stats

def med_mad_MZNP_2(GE_vec):
    """
    Compute median, MAD, modified Z-scores, mean, std, and mode for a vector.
    Implements Leys 2013 robust Z-score method.

    Parameters
    ----------
    GE_vec : array_like
        Input data vector.

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

    GE_vec = np.array(GE_vec)
    med = np.median(GE_vec)
    
    # scipy.stats.mode returns mode and count; mode is an array
    mode_result = stats.mode(GE_vec, nan_policy='omit')
    mod = mode_result.mode.item()  # safely get the scalar value from numpy array
    
    b = 1.4826  # scale factor for MAD assuming normal distribution

    mad = np.median(np.abs(GE_vec - med))
    mad2 = b * mad

    med_mad = [med, mad2]

    me_sd = [np.mean(GE_vec), np.std(GE_vec, ddof=1)]

    # Modified Z-scores
    # Avoid division by zero if mad or mad2 is zero:
    if mad == 0:
        MZ = np.zeros_like(GE_vec)
    else:
        MZ = 0.6745 * (GE_vec - med) / mad

    if mad2 == 0:
        MZ2 = np.zeros_like(GE_vec)
    else:
        MZ2 = (GE_vec - med) / mad2

    # Standard Z-scores
    mean_val = me_sd[0]
    std_val = me_sd[1]
    if std_val == 0:
        Z = np.zeros_like(GE_vec)
    else:
        Z = (GE_vec - mean_val) / std_val

    return med_mad, MZ2, MZ, Z, me_sd, mod
