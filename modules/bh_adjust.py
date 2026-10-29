import numpy as np

def bh_adjust(p):
    """
    Benjamini-Hochberg adjusted p-values (q-values controlling the FDR).

    Each finite p is ranked among the m finite p-values; q = p * m / rank,
    then made monotone (cumulative minimum from the largest p down) and
    capped at 1. NaN entries are ignored for m and stay NaN.

    Parameters
    ----------
    p : array-like
        p-values for one family of tests (e.g. all target gRNAs in one replicate).

    Returns
    -------
    q : numpy.ndarray
        BH-adjusted values, same shape and order as ``p``.
    """
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    m = ok.sum()
    if m == 0:
        return q

    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]   # monotone from the top
    q_ok = np.empty(m)
    q_ok[order] = np.minimum(ranked, 1.0)
    q[ok] = q_ok
    return q
