import matplotlib.pyplot as plt
import os

def plot_histograms(binn, perc_t1, perc_zt1, baseline, cond1, gene_type, st, en, output_dir):

    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)
    plt.bar(binn, perc_t1, color='k')
    plt.grid(True)
    plt.xlim([st, en])
    plt.xlabel('LFC')

    plt.subplot(2, 1, 2)
    plt.bar(binn, perc_zt1, color='g')
    plt.grid(True)
    plt.xlim([st, en])
    plt.xlabel('Z-corrected LFC')


    plt.suptitle(f'gRNA distribution: {cond1} vs {baseline}', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_dir, f"target_distri_LFC_Z_corrected_{cond1}_vs_{baseline}"), dpi=300, bbox_inches="tight")
    plt.close()