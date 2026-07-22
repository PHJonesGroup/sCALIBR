import numpy as np
import matplotlib.pyplot as plt
import os

def _as_pooled_1d(perc):
    if perc is None:
        return None
    perc = np.asarray(perc, dtype=float)
    return np.nanmean(perc, axis=1) if perc.ndim == 2 else perc

def distri_target_contr_plots_all(binn, category_percs, highlight, title, output_dir):
    """
    Plot pooled LFC distributions for any number of categories.

    category_percs : dict {category_name -> 1-D (or 2-D) perc array}
    highlight      : category name to emphasise (e.g. the target)
    """
    # keep only categories that have data
    items = [(name, _as_pooled_1d(p)) for name, p in category_percs.items()]
    items = [(n, p) for n, p in items if p is not None]
    if not items:
        return 1

    width = np.diff(binn)[0]
    cmap = plt.get_cmap('tab10')                       # distinct colours for N categories
    colors = {name: cmap(i % 10) for i, (name, _) in enumerate(items)}

    # ---- Figure 1: one panel per category ----
    n = len(items)
    plt.figure(figsize=(10, 2.2 * n))
    plt.suptitle(title, fontsize=16)
    for k, (name, perc) in enumerate(items, start=1):
        plt.subplot(n, 1, k)
        plt.bar(binn, perc, width=width, color=colors[name])
        plt.grid(True)
        plt.title(name + (" (target)" if name == highlight else ""), fontsize=12)
        plt.ylabel('%', fontsize=12)
        if k == n:
            plt.xlabel('LFC', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"distri_separate_{title}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    # ---- Figure 2: all categories overlaid ----
    plt.figure(figsize=(10, 5))
    for name, perc in items:
        lw = 2.5 if name == highlight else 1.2
        plt.step(binn, perc, where='mid', color=colors[name], linewidth=lw, label=name)
    plt.grid(True)
    plt.title(title, fontsize=14)
    plt.xlabel('LFC', fontsize=12)
    plt.ylabel('%', fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"distri_overlaid_{title}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()
    return 1