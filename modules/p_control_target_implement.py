import pandas as pd
from .implement_p_control_indiv_gRNA import implement_p_control_indiv_gRNA

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
        Per-bin control p-curves (depletion, enrichment) from CTR_stats.

    Returns
    -------
    T_q_verti : pandas.DataFrame
            [gRNA, gene, gene_type, lfc(_<rep>), p(_<rep>), dep_p(_<rep>),
            enr_p(_<rep>), ...]: control-calibrated empirical p-values, two-sided,
            depletion and enrichment (not corrected for multiple testing).
    """
    # Concatenate tables vertically 
    T_target_T_control = pd.concat([T_target, T_control], ignore_index=True)

    # Assign p/q values based on LFC bins
    T_q_verti = implement_p_control_indiv_gRNA(T_target_T_control, binn, p_cont)
    
    return T_q_verti