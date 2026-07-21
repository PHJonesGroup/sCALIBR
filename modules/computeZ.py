import numpy as np
import matplotlib.pyplot as plt
from .make_histo_LFC import make_histo_LFC

def computeZ(st, en, step, LFC_t, me_sd_z):
    """
    Compute per-replicate Z-scores for targets, using the mean/SD from the
    pooled control (zGE / intergenic).

    Inputs:
    - LFC_t: 2-D array (n_gRNA x n_reps) of target LFCs, one column per rep
    - me_sd_z: [mean, SD] of the pooled control LFC
    - st, en, step: histogram bin parameters

    Outputs:
    - Z_t: 2-D array (n_gRNA x n_reps) of Z-scored LFCs, one column per rep
    - binn: bin edges
    - perc_t:  pooled % histogram of LFC  (all reps combined)
    - perc_zt: pooled % histogram of Z    (all reps combined)
    """
    LFC_t = np.asarray(LFC_t, dtype=float)
    if LFC_t.ndim == 1:
        LFC_t = LFC_t[:, None]                 # treat single series as one column

    mu, sd = me_sd_z[0], me_sd_z[1]
    if sd == 0 or not np.isfinite(sd):
        raise ValueError("Control SD is zero or invalid — cannot compute Z-scores")

    # per-replicate Z: each rep's LFC standardized against the pooled-control null
    Z_t = (LFC_t - mu) / sd                     # (n_gRNA x n_reps)

    # pooled histograms (all reps' values combined) for the diagnostic plot
    lfc_flat = LFC_t[np.isfinite(LFC_t)]
    z_flat   = Z_t[np.isfinite(Z_t)]
    binn, _, perc_t  = make_histo_LFC(step, lfc_flat, st, en)
    _,    _, perc_zt = make_histo_LFC(step, z_flat,   st, en)

    return Z_t, binn, perc_t, perc_zt
