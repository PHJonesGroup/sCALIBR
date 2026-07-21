import numpy as np
def make_histo_LFC(step, vec_num, st, en):
    """Histogram counts and percentages on a fixed bin grid (st..en, width step)."""
    vec_num = np.asarray(vec_num, dtype=float)
    vec_num = vec_num[np.isfinite(vec_num)]        # drop NaN/inf

    bins = np.arange(st, en + step, step)
    val, _ = np.histogram(vec_num, bins=bins)

    total = val.sum()
    perc = 100 * val / total if total > 0 else np.zeros_like(val, dtype=float)

    return bins[:-1], val, perc
