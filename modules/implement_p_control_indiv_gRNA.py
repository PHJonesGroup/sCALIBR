import numpy as np
import pandas as pd

def implement_p_control_indiv_grna(T_gRNA_LFC, binn, p_cont):
    """
    Implement p-controls by LFC targets (instead of p-targets) at gRNA level.
    
    Parameters
    ----------
    T_target : pandas.DataFrame
        Target gRNA LFC table: [gRNA, gene, lfc(_<rep>)...].
    binn : array-like
        Histogram bin edges.
    p_cont : array-like
        Per-bin control p-values (1-D), aligned to ``binn``.    
    Returns
    -------
    out : pandas.DataFrame
            [gRNA, gene, lfc(_<rep>)..., Q(_<rep>)...] with control-calibrated q-values.
    """
    #implement p-controls by LFC targets (instead of p-tagets)'
    lfc_cols = [c for c in T_gRNA_LFC.columns if c.startswith('lfc')]

    out = pd.DataFrame({
        'gRNA': T_gRNA_LFC.iloc[:, 0].values,
        'gene': T_gRNA_LFC.iloc[:, 1].values,
    })

    for col in lfc_cols:
        suffix = col[len('lfc'):]                
        LFC = T_gRNA_LFC[col].astype(float).to_numpy()

        qq = np.full(len(LFC), np.nan)
        for i in range(len(binn) - 1):
            idx = (LFC > binn[i]) & (LFC <= binn[i + 1])
            qq[idx] = p_cont[i]

        out[f'lfc{suffix}'] = LFC
        out[f'Q{suffix}']   = qq

    return out
