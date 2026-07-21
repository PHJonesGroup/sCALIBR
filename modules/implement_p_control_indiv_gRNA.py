import numpy as np
import pandas as pd

def implement_p_control_indiv_gRNA(T_gRNA_LFC, binn, p_contr):
    """
    Implement p-controls by LFC targets (instead of p-targets) at gRNA level.
    
    Parameters:
    - T_gRNA_LFC: pd.DataFrame with columns ['gRNA', 'gene', 'LFCM', 'LFCF']
    - binn: array-like, bin edges
    - p_contr: array-like, p/q values corresponding to binn for LFC
    
    Returns:
    - T_q_vertical: pd.DataFrame with columns ['gRNA', 'gene', 'lfc', 'Q']
    """
    #implement p-controls by LFC targets (instead of p-tagets)'
    lfc_cols = list(T_gRNA_LFC.columns[2:])
    out = pd.DataFrame({
        'gRNA': T_gRNA_LFC.iloc[:, 0].values,
        'gene': T_gRNA_LFC.iloc[:, 1].values,
    })
    for col in lfc_cols:
        LFC = T_gRNA_LFC[col].astype(float).to_numpy()
        qq = np.full(len(LFC), np.nan)
        for i in range(len(binn) - 1):
            idx = (LFC > binn[i]) & (LFC <= binn[i + 1])
            qq[idx] = p_contr[i]
        rep = col.replace('lfc_', '')
        out[f'lfc_{rep}'] = LFC
        out[f'Q_{rep}'] = qq
    return out

