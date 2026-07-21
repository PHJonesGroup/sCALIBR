import pandas as pd
import numpy as np
from .make_histo_LFC import make_histo_LFC
from .compute_p_critLFC import compute_p_critLFC
from .med_mad_MZNP_2 import med_mad_MZNP_2
from .plot_three_panels import plot_three_panels

def make_histo_crit_stats(alf, st, en, step, T_zGE, cond1, cond2, control, output_dir):
    """
    Histograms, critical LFC thresholds, and robust stats for control genes
    (cond1 vs cond2), pooling all replicate LFC columns.

    T_zGE : DataFrame [gRNA, gene, lfc_<rep1>, lfc_<rep2>, ...]
    """
    # 1. Pool all replicate LFC columns (everything after gRNA, gene)
    lfc_cols = list(T_zGE.columns[2:])
    if not lfc_cols:
        raise ValueError("No LFC columns found after gRNA, gene")
    LFC = T_zGE[lfc_cols].astype(float).values.ravel()   # flatten reps into one vector
    LFC = LFC[np.isfinite(LFC)]                           # drop NaN/inf

    # 2. Histogram of the pooled control LFCs
    bin_, his, perc = make_histo_LFC(step, LFC, st, en)

    # 3. Critical LFC thresholds & p-curve
    bin_p, crit_LR = compute_p_critLFC(alf, bin_, his, cond1, cond2, control, output_dir)

    # 4. Robust stats, Z / MZ on the pooled control LFCs
    med_mad, MZZ, MZ, Z, me_sd, mod = med_mad_MZNP_2(LFC)

    # 5. Diagnostic plot
    _, _, perc_z  = make_histo_LFC(step, Z,  st, en)
    _, _, perc_mz = make_histo_LFC(step, MZ, st, en)
    plot_three_panels(
        bin_, perc, perc_z, perc_mz,
        title=f'Distribution of {control} genes: {cond1} vs {cond2}',
        color='b', xlab='LFC',
        save_name=f'distri_{control}_{cond1}_vs_{cond2}', output_dir=output_dir)

    return bin_, his, perc, crit_LR, bin_p, med_mad, me_sd, mod, MZ, Z, len(LFC)