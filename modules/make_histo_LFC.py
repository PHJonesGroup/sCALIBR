import numpy as np

def make_histo_LFC(step, vec_num, st, en):
    """
    Histogram counts and percentages on a fixed bin grid (st..en, width step).
    
    Parameters
    ----------
    step : float
        Histogram bin width.
    vec_num : array-like
        Values to histogram (e.g. LFC values); NaN/inf entries are ignored.
    st, en : float
        Lower / upper bounds of the histogram range.

    Returns
    -------
    bin_return : numpy.ndarray
        Left edges of the bins (length = number of bins).
    val : numpy.ndarray of int
        Count of values falling in each bin.
    perc : numpy.ndarray of float
        Percentage of values in each bin (sums to 100 over finite values,
        or all zeros if there are no finite values).
    """
    vec_num = np.asarray(vec_num, dtype=float)
    vec_num = vec_num[np.isfinite(vec_num)]        # drop NaN/inf

    bins = np.arange(st, en + step, step)
    val, _ = np.histogram(vec_num, bins=bins)

    total = val.sum()
    perc = 100 * val / total if total > 0 else np.zeros_like(val, dtype=float)

    return bins[:-1], val, perc
