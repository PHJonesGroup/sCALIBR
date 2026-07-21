
import numpy as np
from .make_histo_LFC import make_histo_LFC

def make_histo_vec_rep12(step: float,
                         vec_1: np.ndarray,
                         vec_2: np.ndarray,
                         st: float,
                         en: float):
    """
    Build histograms (counts & %s) for two replicate vectors and for the
    pooled vector.

    Parameters
    ----------
    step : float
        Bin width.
    vec_1, vec_2 : 1‑D arrays
        Data vectors (e.g. LFC values) for replicate‑1 and replicate‑2.
    st, en : float
        Lower / upper bounds for binning.

    Returns
    -------
    bin        : 1‑D array of bin left‑edges (len == #bins)
    perc_1     : % histogram for vec_1
    perc_2     : % histogram for vec_2
    perc_12    : % histogram for pooled vec_1 + vec_2
    his_1      : counts per bin for vec_1
    his_2      : counts per bin for vec_2
    his_12     : counts per bin for pooled data
    n1_n2_n12  : tuple  (len(vec_1), len(vec_2), len(vec_1)+len(vec_2))
    """

    # Histogram for replicate‑1
    bin_edges, his_1, perc_1 = make_histo_LFC(step, np.asarray(vec_1), st, en)

    # Histogram for replicate‑2
    _,        his_2, perc_2 = make_histo_LFC(step, np.asarray(vec_2), st, en)

    # Histogram for pooled replicates
    vec_12 = np.concatenate([vec_1, vec_2])
    _, his_12, perc_12 = make_histo_LFC(step, vec_12, st, en)

    n1_n2_n12 = (len(vec_1), len(vec_2), len(vec_12))

    return (
        bin_edges,
        perc_1,
        perc_2,
        perc_12,
        his_1,
        his_2,
        his_12,
        n1_n2_n12,
    )
