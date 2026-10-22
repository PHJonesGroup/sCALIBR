import numpy as np

def q_val_frequentist_critical(his: np.ndarray):
    """
    Given a histogram of LFC values (counts per bin), compute the combined
    two-tailed p-curve: the left-tail cumulative fraction up to the point
    where the left and right curves cross, the right-tail fraction beyond it.

    Parameters
    ----------
    his : numpy.ndarray
        1-D array of counts per bin.

    Returns
    -------
    p : numpy.ndarray
        Combined p-curve (left tail then right tail); length == len(his).
    """
    S = his.sum()         # total gRNA count

    # Cumulative fractions (right & left)
    cum_R  = np.cumsum(his[::-1])[::-1]   # cumulative counts from right
    cum_L  = np.cumsum(his)               # cumulative counts from left
    p_right = cum_R / S
    p_left  = cum_L / S

    #intersection (mode‑like) where |p_R - p_L| is minimal
    ind_min = np.argmin(np.abs(p_right - p_left))

    # combined p vector
    p_combined = np.concatenate([p_left[:ind_min+1],
                                 p_right[ind_min+1:]])

    return p_combined
