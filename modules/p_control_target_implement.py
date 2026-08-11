import pandas as pd
from .implement_p_control_indiv_grna import implement_p_control_indiv_grna

def p_control_target_implement(T_target, T_control, binn, p_cont):
    """
    Assign control-calibrated p/q values to gRNAs based on their LFC bin.

    Parameters
    ----------
    T_target : pandas.DataFrame
        Target gRNA LFC table: [gRNA, gene, lfc(_<rep>)...].
    T_control : pandas.DataFrame
        Control gRNA LFC table: [gRNA, gene, lfc(_<rep>)...].
    binn : array-like
        Histogram bin edges.
    p_cont : array-like
        Per-bin control p-values (1-D), aligned to ``binn``.

    Returns
    -------
    T_q_verti : pandas.DataFrame
            [gRNA, gene, lfc(_<rep>)..., Q(_<rep>)...] with control-calibrated q-values.
    """
    # Concatenate tables vertically 
    T_target_T_control = pd.concat([T_target, T_control], ignore_index=True)

    # Assign p/q values based on LFC bins
    T_q_verti = implement_p_control_indiv_grna(T_target_T_control, binn, p_cont)
    
    return T_q_verti