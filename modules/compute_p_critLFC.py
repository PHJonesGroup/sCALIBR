import numpy as np
import matplotlib.pyplot as plt
import os
from .q_val_frequentist_critical import q_val_frequentist_critical

def compute_p_critLFC(alf, binn, hiss, cond1, cond2, control, output_dir):
    """
    Compute the p-curve and critical LFC values for a single LFC distribution
    (contrast cond1 vs cond2).

    Parameters
    ----------
    alf   : significance level
    binn  : bin edges/centers for the LFC histogram
    hiss  : histogram counts (single series)
    cond1, cond2 : condition labels for the plot title

    Returns
    -------
    bin_p   : ndarray, columns [bin, p]
    crit_LR : ndarray, [left_critical, right_critical]
    """
    p, cl, cr, bin_pz, med_LFCp, hiss4p = q_val_frequentist_critical(alf, binn, hiss)
    crit_LR = np.array([cl, cr])

    bin_p = np.column_stack((binn, p))

    # Plot the single p-curve
    plt.figure(figsize=(8, 6))
    plt.plot(binn, p, 'b', label=f'{cond1} vs {cond2}')
    plt.grid(True)
    plt.xlabel('LFC')
    plt.ylabel('Probability')
    plt.title('Zero Gene Expression gRNAs', fontsize=14)
    plt.legend()
    plt.savefig(os.path.join(output_dir, f"p_controls_{control}_{cond1}_vs_{cond2}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    return bin_p, crit_LR