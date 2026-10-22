import numpy as np

def med_mad_MZNP(LFC_t):
    """
    Compute mean, standard deviation and standard Z-scores for a vector.

    Parameters
    ----------
    LFC_t : array_like
        all replicate LFC columns.

    Returns
    -------
    Z : ndarray
        Standard Z-scores (mean/std normalization).
    me_sd : list
        [mean, standard deviation]
    """

    LFC_t = np.array(LFC_t)

    me_sd = [np.mean(LFC_t), np.std(LFC_t, ddof=1)]

    # Standard Z-scores
    mean_val = me_sd[0]
    std_val = me_sd[1]
    if std_val == 0:
        Z = np.zeros_like(LFC_t)
    else:
        Z = (LFC_t - mean_val) / std_val

    return Z, me_sd
