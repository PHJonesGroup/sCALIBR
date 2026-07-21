import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from .make_histo_crit_stats import make_histo_crit_stats
from .q_val_frequentist_critical import q_val_frequentist_critical

def CTR_stats_zGE(alf, st, en, step, T_zGE, hist, cond1, cond2, condz, control, output_dir):
    """
    Control vs. target analysis for a single contrast (cond1 vs cond2),
    preserving multiple replicate LFC columns.

    T_zGE : [gRNA, gene, lfc_<r1>, lfc_<r2>, ...]  (control table)
    hist  : target histogram (1-D, or 2-D bins x reps -> pooled here)

    Returns
    -------
    crit_LR   : [left, right] critical LFC (from pooled control)
    me_sd     : [mean, SD]    of pooled control LFC
    med_mad   : [median, MAD] of pooled control LFC
    binn      : bin edges
    p_cont    : control p-curve
    hiss_cont : control histogram
    p_targ    : target p-curve
    T_zGE_out : control table, unchanged: [gRNA, gene, lfc_<r1>, ...]
    """

    # 1. Pooled-control histogram & stats
    (binn, hiss_z, perc, crit_LR, bin_pz, med_mad, me_sd, mod, MZ, Z, n_z
    ) = make_histo_crit_stats(alf, st, en, step, T_zGE, cond1, cond2, control, output_dir)

    # 2. Control table passed through unchanged (gRNA, gene, per-rep LFCs)
    T_zGE_out = T_zGE.copy()

    # 3. Control p-curve & histogram
    p_cont    = bin_pz[:, 1]
    hiss_cont = hiss_z

    # 4. Target p-curve (pool reps if hist is 2-D)
    hist = np.asarray(hist)
    if hist.ndim == 2:
        hist = hist.sum(axis=1)          # pooled counts across reps -> 1-D
    p_targ, cL, cR, bin_pi, med_LFCp, his4p = q_val_frequentist_critical(alf, binn, hist)

    # 5. Plot: target vs. control p-curves
    plt.figure()
    plt.plot(binn, p_cont, color="0.5", linestyle="--", label=f"control ({control})")
    plt.plot(binn, p_targ, color="C3", linestyle="-",  label="target")
    plt.title(f"Target vs. {control} gRNAs ({cond1} vs {cond2})")
    plt.xlabel("LFC bin")
    plt.ylabel("Probability")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"p_distri_targ_cont_{cond1}_vs_{cond2}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    return crit_LR, me_sd, med_mad, binn, p_cont, hiss_cont, p_targ, T_zGE_out