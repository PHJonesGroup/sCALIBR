import matplotlib.pyplot as plt
import numpy as np
import os

def plot_three_panels(bin_edges, perc_lfc, perc_z=None, perc_mz=None,
                       title='', color='b', xlab='LFC', save_name ='', output_dir=''):
    """
    Plot up to three stacked histogram panels: raw LFC, Z-LFC, and MZ-LFC.
    
    Parameters
    ----------
    bin_edges : numpy.ndarray
        Bin left-edges (x positions for the bars); must have length >= 2.
    perc_lfc : numpy.ndarray
        Percentage histogram of the raw LFC values (top panel).
    perc_z : numpy.ndarray, optional
        Percentage histogram of the Z-scored values (middle panel); skipped if None.
    perc_mz : numpy.ndarray, optional
        Percentage histogram of the MZ-scored values (bottom panel); skipped if None.
    title : str, optional
        Figure-level title (suptitle).
    color : str, optional
        Bar colour (default 'b').
    xlab : str, optional
        X-axis label for the top (raw LFC) panel (default 'LFC').
    save_name : str, optional
        Filename for the saved figure (within output_dir).
    output_dir : str, optional
        Directory where the figure is saved.

    Returns
    -------
    None
    """

    plt.figure()
    plt.suptitle(title, fontsize=14)
    plt.subplot(3, 1, 1)
    plt.bar(bin_edges, perc_lfc, width=np.diff(bin_edges)[0], color=color)
    plt.grid(True)
    plt.xlabel(xlab)

    if perc_z is not None:
        plt.subplot(3, 1, 2)
        plt.bar(bin_edges, perc_z, width=np.diff(bin_edges)[0], color=color)
        plt.grid(True)
        plt.xlabel('Z LFC')

    if perc_mz is not None:
        plt.subplot(3, 1, 3)
        plt.bar(bin_edges, perc_mz, width=np.diff(bin_edges)[0], color=color)
        plt.grid(True)
        plt.xlabel('MZ LFC')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, save_name), dpi=300, bbox_inches="tight")
    plt.close()