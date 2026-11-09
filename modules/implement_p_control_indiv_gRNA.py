import numpy as np
import pandas as pd
from .q_val_frequentist_critical import two_sided

def implement_p_control_indiv_gRNA(T_gRNA_LFC, binn, p_cont):
    """
    Implement p-controls by LFC targets (instead of p-targets) at gRNA level.

    Parameters
    ----------
    T_gRNA_LFC : pandas.DataFrame
        Target gRNA LFC table: [gRNA, gene, gene_type, lfc(_<rep>)...].
    binn : array-like
        Histogram bin edges.
    p_cont : array-like
        Per-bin control p-curves, last axis (depletion, enrichment):
        (n_bins, 2) shared by all replicates, or (n_bins, n_rep, 2).

    Returns
    -------
    out : pandas.DataFrame
            [gRNA, gene, gene_type, lfc(_<rep>), p(_<rep>), dep_p(_<rep>),
            enr_p(_<rep>), ...]: per gRNA, the control-calibrated empirical
            p-values for depletion (left tail) and enrichment (right tail),
            each running 0 to 1, and p = two-sided (2 x the smaller, capped
            at 1). None is corrected for multiple testing.
    """
    lfc_cols = [c for c in T_gRNA_LFC.columns if c.startswith('lfc')]

    out = pd.DataFrame({
        'gRNA': T_gRNA_LFC.iloc[:, 0].values,
        'gene': T_gRNA_LFC.iloc[:, 1].values,
        'gene_type': T_gRNA_LFC.iloc[:, 2].values,
    })

    # pooled: one pair of curves shared by all reps; sample: one pair per rep
    p_cont = np.asarray(p_cont, dtype=float)
    if p_cont.ndim == 3 and p_cont.shape[1] != len(lfc_cols):
        raise ValueError(f"p_cont has {p_cont.shape[1]} replicates but there are "
                         f"{len(lfc_cols)} lfc columns: {lfc_cols}")

    for j, col in enumerate(lfc_cols):
        suffix = col[len('lfc'):]
        LFC = T_gRNA_LFC[col].astype(float).to_numpy()
        p_col = p_cont[:, j, :] if p_cont.ndim == 3 else p_cont    # (n_bins, 2)

        dep = np.full(len(LFC), np.nan)
        enr = np.full(len(LFC), np.nan)
        for i in range(len(binn) - 1):
            idx = (LFC > binn[i]) & (LFC <= binn[i + 1])
            dep[idx] = p_col[i, 0]
            enr[idx] = p_col[i, 1]

        out[f'lfc{suffix}']   = LFC
        out[f'p{suffix}']     = two_sided(dep, enr)
        out[f'dep_p{suffix}'] = dep
        out[f'enr_p{suffix}'] = enr

    return out
