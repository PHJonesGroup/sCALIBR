import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from .make_histo_crit_stats import make_histo_crit_stats
from .q_val_frequentist_critical import q_val_frequentist_critical

def CTR_stats(alf, st, en, step, T_control, hist, baseline, cond, control, output_dir):
    """
    Calibrate a control null and score targets against it.

    Parameters
    ----------
    alf : float
        Significance level for the critical-value / p-curve calculation.
    st, en : float
        Lower / upper bounds of the histogram (LFC) range.
    step : float
        Histogram bin width.
    T_control : pandas.DataFrame
        Control table: [gRNA, gene, lfc_<rep>, ...] (per-replicate LFC columns).
    hist : numpy.ndarray
        Target histogram; 1-D counts, or 2-D (n_bins x n_reps) which is pooled
        (summed across reps) internally.
    baseline, cond : str
        Condition labels.
    control : str
        Control category label (e.g. 'Zero-expressed gene'), used in labels.
    output_dir : str
        Directory where the target-vs-control p-curve plot is saved.

    Returns
    -------
    crit_LR : numpy.ndarray
        [left, right] critical LFC thresholds, from the pooled control.
    me_sd : numpy.ndarray
        [mean, SD] of the pooled control LFC.
    med_mad : numpy.ndarray
        [median, MAD] of the pooled control LFC.
    binn : numpy.ndarray
        Histogram bin edges.
    p_cont : numpy.ndarray
        Control p-curve (per-bin probability).
    hiss_cont : numpy.ndarray
        Control histogram counts.
    p_targ : numpy.ndarray
        Target p-curve (per-bin probability).
    """
    # 1. Pooled-control histogram & stats
    (binn, hiss_z, perc, crit_LR, bin_pz, med_mad, me_sd, mod, MZ, Z, n_z
    ) = make_histo_crit_stats(alf, st, en, step, T_control, baseline, cond, control, output_dir)

    # 3. Control p-curve & histogram
    p_cont    = bin_pz[:, 1]
    hiss_cont = hiss_z

    # 4. Target p-curve (pool reps if hist is 2-D)
    hist = np.asarray(hist)
    if hist.ndim == 2:
        hist = hist.sum(axis=1)          # pooled counts across reps -> 1-D
    p_targ, cL, cR, bin_pi, med_LFCp, his4p = q_val_frequentist_critical(alf, binn, hist)

    p_table = pd.DataFrame({
        "bin" : np.round(binn,1),
        f"control_{control}": np.round(p_cont,1),
        f"target": np.round(p_targ,1)})
    p_table.to_csv(
        os.path.join(output_dir, f"p_distri_targ_cont_{cond}_vs_{baseline}_data.csv"), index=False)

    # 5. Plot: target vs. control p-curves
    plt.figure()
    plt.plot(binn, p_cont, color="0.5", linestyle="--", label=f"control ({control})")
    plt.plot(binn, p_targ, color="C3", linestyle="-",  label="target")
    plt.title(f"Target vs. {control} gRNAs ({cond} vs {baseline})")
    plt.xlabel("LFC bin")
    plt.ylabel("Probability")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"p_distri_targ_cont_{cond}_vs_{baseline}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    return crit_LR, me_sd, med_mad, binn, p_cont, hiss_cont, p_targ