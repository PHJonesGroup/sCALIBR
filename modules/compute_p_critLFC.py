import numpy as np
import matplotlib.pyplot as plt
import os
from .q_val_frequentist_critical import q_val_frequentist_critical

def compute_p_critLFC(alf, binn, hiss, baseline, cond, control, output_dir):
    """
    Compute the p-curve and critical LFC values for a single LFC distribution
    (contrast cond vs baseline).

    Parameters
    ----------
    alf   : float
        Significance level for the critical-value / p-curve calculation.
    binn : numpy.ndarray
        Bin edges (left edges) for the LFC histogram.
    hiss : numpy.ndarray
        Histogram counts (single series), same length as ``binn``.
    baseline, cond : str
        Condition labels.
    control : str
        Control category label (e.g. 'Zero-expressed gene'), used in labels.
    output_dir : str
        Directory where the plot is saved.

    Returns
    -------
    bin_p : numpy.ndarray
        Two-column array [bin, p]: bin edges with the corresponding p-curve.
    crit_LR : numpy.ndarray
        [left_critical, right_critical] critical LFC thresholds.
    """
    p, cl, cr, bin_pz, med_LFCp, hiss4p = q_val_frequentist_critical(alf, binn, hiss)
    crit_LR = np.array([cl, cr])

    bin_p = np.column_stack((binn, p))

    # Plot the single p-curve
    plt.figure(figsize=(8, 6))
    plt.plot(binn, p, 'b', label=f'{cond} vs {baseline}')
    plt.grid(True)
    plt.xlabel('LFC')
    plt.ylabel('Probability')
    plt.title('Zero Gene Expression gRNAs', fontsize=14)
    plt.legend()
    plt.savefig(os.path.join(output_dir, f"p_controls_{control}_{cond}_vs_{baseline}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    return bin_p, crit_LR