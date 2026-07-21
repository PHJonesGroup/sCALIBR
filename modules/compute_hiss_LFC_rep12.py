import numpy as np
import pandas as pd
from .make_histo_LFC import make_histo_LFC

import numpy as np
import pandas as pd
from .make_histo_LFC import make_histo_LFC

def compute_hiss_LFC_rep12(chr_lfc, st, en, step):
    """
    Compute histograms for individual gRNAs across one or more LFC columns.

    Parameters:
        chr_lfc (pd.DataFrame): [gRNA, gene, lfc_<r1>, lfc_<r2>, ...]
                                (any number of LFC columns from col index 2 on)
        st, en, step: histogram range and bin width

    Returns:
        binn : bin edges (shared across all columns)
        hiss : (#bins, n_lfc) counts, one column per LFC column
        perc : (#bins, n_lfc) percentages, one column per LFC column
        LFC  : (#rows, n_lfc) the LFC values
        labels : list of the LFC column names
        st1, en1 : data-driven bounds across all LFC columns
    """
    # all LFC columns (everything after gRNA, gene)
    lfc_cols = list(chr_lfc.columns[2:])
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