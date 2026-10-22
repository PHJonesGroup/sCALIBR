import numpy as np
from .make_histo_LFC import make_histo_LFC

def computeZ(st, en, step, LFC_t, me_sd_z):
    """Compute per-replicate Z-scores for targets against the control null.

    Accepts either a single pooled null or one null per replicate:
      - me_sd_z of shape (2,)        -> every replicate is standardised against
                                        the same pooled-control mean and SD.
      - me_sd_z of shape (n_rep, 2)  -> column j of LFC_t is standardised against
                                        row j's mean and SD (per-sample calibration).

    Parameters
    ----------
    st, en : float
        Lower / upper bounds of the histogram (LFC) range; widened if needed
        so the Z-scored values also fit.
    step : float
        Histogram bin width.
    LFC_t : numpy.ndarray
        Target LFCs, shape (n_gRNA, n_reps); a 1-D array is treated as one column.
    me_sd_z : numpy.ndarray
        [mean, SD] of the pooled control LFC, or an (n_rep, 2) array of
        [mean, SD] per replicate, ordered to match the columns of LFC_t.

    Returns
    -------
    Z_t : numpy.ndarray
        Z-scored LFCs, shape (n_gRNA, n_reps), one column per replicate.
    binn : numpy.ndarray
        Histogram bin edges.
    perc_t : numpy.ndarray
        Pooled percentage histogram of the raw LFC values (all reps combined).
    perc_zt : numpy.ndarray
        Pooled percentage histogram of the Z-scored values (all reps combined).

    Raises
    ------
    ValueError
        If a control SD is zero or non-finite, or if a per-replicate me_sd_z
        has a different number of rows than LFC_t has columns.
    """
    LFC_t = np.asarray(LFC_t, dtype=float)
    if LFC_t.ndim == 1:
        LFC_t = LFC_t[:, None]

    me_sd_z = np.asarray(me_sd_z, dtype=float)

    if me_sd_z.ndim == 1:                         # pooled: one null for all reps
        mu, sd = me_sd_z[0], me_sd_z[1]
        bad = (sd == 0) or (not np.isfinite(sd))
        if bad:
            raise ValueError("Control SD is zero or invalid — cannot compute Z-scores")
    else:                                         # per-replicate nulls
        if me_sd_z.shape[0] != LFC_t.shape[1]:
            raise ValueError(
                f"me_sd_z has {me_sd_z.shape[0]} rows but LFC_t has "
                f"{LFC_t.shape[1]} replicate columns — they must match in order")
        mu = me_sd_z[:, 0][None, :]               # (1, n_rep), broadcasts over gRNAs
        sd = me_sd_z[:, 1][None, :]
        bad = np.any(sd == 0) | np.any(~np.isfinite(sd))
        if bad:
            raise ValueError("A control SD is zero or invalid — cannot compute Z-scores")

    # each replicate's LFC standardised against its own (or the pooled) null
    Z_t = (LFC_t - mu) / sd                       # (n_gRNA x n_reps)

    # pooled histograms for the diagnostic plot; widen [st, en] so the
    # normalised LFCs (which can exceed the raw LFC range) are all shown
    lfc_flat = LFC_t[np.isfinite(LFC_t)]
    z_flat   = Z_t[np.isfinite(Z_t)]
    if z_flat.size:
        st = min(st, float(np.round((np.floor(z_flat.min() / step) - 1) * step, 10)))
        en = max(en, float(np.round((np.ceil(z_flat.max() / step) + 1) * step, 10)))
    binn, _, perc_t  = make_histo_LFC(step, lfc_flat, st, en)
    _,    _, perc_zt = make_histo_LFC(step, z_flat,   st, en)

    return Z_t, binn, perc_t, perc_zt