import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

def plot_histograms(bin_edges, perc_lfc, perc_z=None,
                    title='', color='b', xlab='LFC', save_name='', output_dir=''):
    """Plot stacked histogram panels: raw LFC and normalised LFC.

    The normalised panel is drawn only if perc_z is given, so the figure has
    one or two rows. The plotted values are also written to a CSV alongside
    the figure.

    Parameters
    ----------
    bin_edges : numpy.ndarray
        Bin left-edges (x positions for the bars); must have length >= 2.
    perc_lfc : numpy.ndarray
        Percentage histogram of the raw LFC values (top panel).
    perc_z : numpy.ndarray, optional
        Percentage histogram of the normalised LFC values; panel skipped if None.
    title : str, optional
        Figure-level title.
    color : str, optional
        Bar colour (default 'b').
    xlab : str, optional
        X-axis label for the top (raw LFC) panel (default 'LFC').
    save_name : str, optional
        Base filename (without extension) for the saved figure and CSV.
    output_dir : str, optional
        Directory where the figure and CSV are saved.

    Returns
    -------
    None
    """
    bin_edges = np.asarray(bin_edges, dtype=float)
    n_bins = len(perc_lfc)

    # --- CSV of the plotted values: only columns that exist ---
    data = {"bin": np.round(bin_edges[:n_bins], 2),
            "LFC": np.round(np.asarray(perc_lfc, dtype=float), 3)}
    if perc_z is not None:
        data["normalised_LFC"] = np.round(np.asarray(perc_z, dtype=float), 3)

    base = os.path.splitext(save_name)[0]
    pd.DataFrame(data).to_csv(
        os.path.join(output_dir, f"{base}_data.csv"), index=False)

    # --- figure: one row per panel actually supplied ---
    panels = [(perc_lfc, xlab)]
    if perc_z is not None:
        panels.append((perc_z, 'Normalised LFC'))

    width = np.diff(bin_edges)[0] if len(bin_edges) > 1 else 0.05

    fig, axs = plt.subplots(len(panels), 1, figsize=(8, 2.6 * len(panels)),
                            squeeze=False)
    fig.suptitle(title, fontsize=14)

    for ax, (vals, lab) in zip(axs[:, 0], panels):
        ax.bar(bin_edges[:n_bins], vals, width=width, color=color)
        ax.grid(True)
        ax.set_xlabel(lab)
        ax.set_ylabel('%')

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(output_dir, f"{base}.pdf"), bbox_inches="tight")
    plt.close()
