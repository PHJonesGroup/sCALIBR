import matplotlib.pyplot as plt
import numpy as np
import os

def plot_three_panels(bin_edges, perc_lfc, perc_z=None, perc_mz=None,
                       title='', color='b', xlab='LFC', save_name ='', output_dir=''):
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