import pandas as pd
from .implement_p_control_indiv_gRNA import implement_p_control_indiv_gRNA

def p_control_target_implement(T_target, T_zGE, binn, p_cont):
    """
    Given p per bin, infer p/q values for [T_target; T_zGE] from LFC values.

    Parameters:
    - T_target: pd.DataFrame with gRNA target data
    - T_zGE: pd.DataFrame with zGE gene data
    - binn: array-like, bin edges
    - p_cont: 2D array-like, calibration p/q values per bin

    Returns:
    - T_q_verti: pd.DataFrame with gRNAs, genes, lfc, Q1
    """

    #FDR  FRR (false rejection)-correct p-values for both tails of LFC NNT distribution
    #implement p-controls into LFC (instead of p-targets)

    # Concatenate tables vertically 
    T_target_zGE = pd.concat([T_target, T_zGE], ignore_index=True)

    # Extract p/q control columns
    p_contr = p_cont

    # Call the function that assigns p/q values based on LFC bins
    T_q_verti = implement_p_control_indiv_gRNA(T_target_zGE, binn, p_contr)
    
    return T_q_verti