import numpy as np
import pandas as pd
from .make_histo_LFC import make_histo_LFC

import numpy as np
import pandas as pd
from .make_histo_LFC import make_histo_LFC

def compute_hiss_LFC(chr_lfc, st, en, step):
    """
    Compute histograms for individual gRNAs across one or more LFC columns.
    
    Parameters
    ----------
    chr_lfc : pandas.DataFrame
        [gRNA, gene, lfc]                 if no replicates, or
        [gRNA, gene, lfc_<rep>, ...]      one LFC column per replicate.
    st, en : float
        Lower / upper bounds of the histogram (LFC) range.
    step : float
        Histogram bin width.

    Returns
    -------
    binn : numpy.ndarray
        Bin edges (left edges), shared across all LFC columns.
    hiss : numpy.ndarray, shape (n_bins, n_lfc)
        Counts per bin, one column per LFC column.
    perc : numpy.ndarray, shape (n_bins, n_lfc)
        Percentages per bin, one column per LFC column.
    LFC : numpy.ndarray, shape (n_rows, n_lfc)
        The raw LFC values used.
    labels : list of str
        The LFC column names, in the order of the hiss/perc columns.
    st1, en1 : float
        Data-driven LFC bounds (min-1, max+1 across all columns); reported only,
        not used for the actual binning.
    
    Raises
    ------
    ValueError
        If LFC columns are missing from ``chr_lfc``
    """
    # all LFC columns (everything after gRNA, gene)
    lfc_cols = list(chr_lfc.columns[3:])
    if not lfc_cols:
        raise ValueError("No LFC columns found (expected columns after gRNA, gene)")

    LFC = chr_lfc[lfc_cols].astype(float).values      # (rows, n_lfc)

    # data-driven bounds across all columns (finite values only)
    finite = LFC[np.isfinite(LFC)]
    if finite.size:
        st1, en1 = finite.min() - 1, finite.max() + 1
    else:
        st1, en1 = st, en

    # histogram each LFC column on the same fixed bins
    hiss_list, perc_list = [], []
    binn = None
    for k in range(LFC.shape[1]):
        b, hiss_k, perc_k = make_histo_LFC(step, LFC[:, k], st, en)
        binn = b
        hiss_list.append(hiss_k)
        perc_list.append(perc_k)

    hiss = np.column_stack(hiss_list)
    perc = np.column_stack(perc_list)

    return binn, hiss, perc, LFC, lfc_cols, st1, en1