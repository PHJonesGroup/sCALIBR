import numpy as np

def q_val_frequentist_critical(his: np.ndarray):
    """
    Given a histogram of LFC values (counts per bin), compute the one-sided
    empirical p-curves for depletion (left tail) and enrichment (right tail).

    Each tail uses p = (k + 1) / (n + 1), where k is the number of values in
    that bin or further out and n is the total. Counting the tested value as
    one more draw from the null keeps p above zero: a value beyond every
    control gets 1 / (n + 1), the smallest p n controls can support.

    Each tail on its own runs from about 0 to 1 under the null, which Fisher's
    method needs; taking the smaller of the two would cap p at about 0.5.

    Parameters
    ----------
    his : numpy.ndarray
        1-D array of counts per bin.

    Returns
    -------
    p_left : numpy.ndarray
        Depletion p per bin: share of values in that bin or further left.
    p_right : numpy.ndarray
        Enrichment p per bin: share of values in that bin or further right.
    """
    S = his.sum()         # total gRNA count

    # Cumulative counts (right & left) -> empirical p, never 0
    cum_R  = np.cumsum(his[::-1])[::-1]   # cumulative counts from right
    cum_L  = np.cumsum(his)               # cumulative counts from left
    p_right = (cum_R + 1) / (S + 1)
    p_left  = (cum_L + 1) / (S + 1)

    return p_left, p_right


def two_sided(p_left, p_right):
    """Two-sided p: twice the smaller tail, capped at 1."""
    return np.minimum(1.0, 2 * np.minimum(p_left, p_right))
