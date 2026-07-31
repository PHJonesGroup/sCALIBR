import numpy as np

def q_val_frequentist_critical(alf: float,
                               bin_edges: np.ndarray,
                               his: np.ndarray):
    """
    Given a histogram of LFC values (counts per bin), compute   \
    • two‑tailed p‑curves (left & right),
    • critical LFC thresholds cL / cR at level `alf`,
    • median‑like intersection point,
    • and a diagnostics matrix `bin_pi`.

    Parameters
    ----------
    alf : float
        Significance level (e.g. 0.05).
    bin_edges : numpy.ndarray
        1-D array of left-edge bin positions (same length as ``his``).
    his : numpy.ndarray
        1-D array of counts per bin (same length as ``bin_edges``).

    Returns
    -------
    p : numpy.ndarray
        Combined p-curve (left tail then right tail); length == len(his).
    cL : float
        Left-tail critical LFC value.
    cR : float
        Right-tail critical LFC value.
    bin_pi : numpy.ndarray
        (n_bins, 6) diagnostics matrix: [bin, p_right, p_left, his, cum_R, cum_L].
    med_LFCp : float
        LFC value where the left and right p-curves intersect.
    his4p : numpy.ndarray
        Combined cumulative counts (left then right); length == len(his).
    """

    thrLowR = -4.5        # hard‑coded right‑tail lower bound
    S = his.sum()         # total gRNA count
    N = len(his)

    # Cumulative fractions (right & left)
    cum_R  = np.cumsum(his[::-1])[::-1]   # cumulative counts from right
    cum_L  = np.cumsum(his)               # cumulative counts from left
    cum_fracR = cum_R / S
    cum_fracL = cum_L / S

    # p‑curves
    p_right = cum_fracR
    p_left  = cum_fracL

    #critical right‑tail threshold 
    crit_right_ind = np.where(p_right <= alf)[0][0]          # first ≤ alf
    lfc_crit_right = bin_edges[crit_right_ind]
    cR             = lfc_crit_right

    # If right tail is “too low”, tighten criterion (alf -> alf/2)
    delta = 0.0
    if lfc_crit_right < thrLowR:
        delta = alf / 2.0
        crit_right_ind = np.where(p_right <= (alf - delta))[0][0]
        lfc_crit_right = bin_edges[crit_right_ind]
        cR = lfc_crit_right

    
    #critical left‑tail threshold
    crit_left_ind = np.where(p_left >= (alf + delta))[0][0]
    lfc_crit_left = bin_edges[crit_left_ind]
    cL            = lfc_crit_left

    #intersection (mode‑like) where |p_R - p_L| is minimal
    ind_min = np.argmin(np.abs(p_right - p_left))
    med_LFCp = bin_edges[ind_min]

    # combined p vector 
    p_combined = np.concatenate([p_left[:ind_min+1],
                                 p_right[ind_min+1:]])

    his4p = np.concatenate([cum_L[:ind_min+1],
                            cum_R[ind_min+1:]])

    # Diagnostics matrix
    bin_pi = np.column_stack([bin_edges,
                              p_right,
                              p_left,
                              his,
                              cum_R,
                              cum_L])

    return (
        p_combined,        # p
        cL,
        cR,
        bin_pi,
        med_LFCp,
        his4p
    )

