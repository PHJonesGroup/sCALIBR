import matplotlib.pyplot as plt
import os

def plot_histograms(binn, perc_t1, perc_zt1, cond1, cond2, condition_label, st, en, output_dir):

    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)
    plt.bar(binn, perc_t1, color='k')
    plt.grid(True)
    plt.title(f'LFC distri {cond1} vs {cond2}')
    plt.xlim([st, en])
    plt.xlabel('LFC')

    plt.subplot(2, 1, 2)
    plt.bar(binn, perc_zt1, color='g')
    plt.grid(True)
    plt.title(f'Z LFC distri {cond1} vs {cond2}')
    plt.xlim([st, en])
    plt.xlabel('LFC Z-normalised')


    plt.suptitle(f'gRNA distribution: {condition_label}', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_dir, f"target_distri_LFC_Z_corr_{condition_label}"), dpi=300, bbox_inches="tight")
    plt.close()