from .compute_hiss_LFC import compute_hiss_LFC
from .distri_target_contr_plots_all import distri_target_contr_plots_all
from .lfc_table import lfc_table

import numpy as np
import pandas as pd

def separate_target_control(step, T_norm_indiv,
                            target_type, control_type,
                            baseline, cond1, plot_name, output_dir,
                            gene_type_column="gene_type"):
    """
    Split gRNAs by the gene_type column, compute per-rep LFC distributions per category, and plot them.

    The histogram range is set from the raw LFCs themselves: it spans the
    smallest to the largest finite LFC across all categories and replicates,
    aligned to ``step`` with one bin of margin on each side, so that every
    gRNA falls inside a bin and later receives a Q value.

    Parameters
    ----------
    step : float
        Histogram bin width.
    T_norm_indiv : pandas.DataFrame
        DataFrame [gRNA, gene, <normalised count columns>]
    target_type : str
        gene_type value denoting targets
    control_type : str
        gene_type value to calibrate on
    baseline : str
        Baseline condition name / column prefix
    cond1 : str
        Treatment condition name / column prefix
    output_dir : str
        Directory where the distribution plot is saved.
    gene_type_column : str, optional
        Name of the category column (default 'gene_type').

    Returns
    -------
    T_target : pandas.DataFrame
        Per-replicate LFC table for targets: [gRNA, gene, lfc_<rep>, ...].
    T_control : pandas.DataFrame
        Per-replicate LFC table for control_type (empty if that category is
        absent): [gRNA, gene, lfc_<rep>, ...].
    hist : numpy.ndarray
        Target histogram counts, (n_bins, n_reps), one column per replicate.
    binn : numpy.ndarray
        Histogram bin edges (shared across all categories)
    st, en : float
        Lower / upper bounds of the histogram (LFC) range, derived from the data;
        pass these on so every later histogram uses the same bins.
    """
    if gene_type_column not in T_norm_indiv.columns:
        raise KeyError(f"'{gene_type_column}' not found. Columns: {T_norm_indiv.columns.tolist()}")

    gt = T_norm_indiv[gene_type_column].astype(str).str.strip()
    categories = sorted(gt.unique())

    if target_type not in categories:
        raise ValueError(f"target_type '{target_type}' not in {gene_type_column}. "
                         f"Available: {categories}")

    control_present = control_type in categories
    if not control_present:
        print(f"[warning] control_type '{control_type}' not found in data "
              f"(available: {categories}). No control group built for it.")

    # one per-rep LFC table per discovered category
    tables = {cat: lfc_table(T_norm_indiv[gt == cat], baseline, cond1)
              for cat in categories}

    # histogram range from the raw LFCs: min..max over all categories and reps,
    # aligned to step, plus one bin of margin so no gRNA falls on/over an edge
    all_lfc = np.concatenate([t.iloc[:, 3:].to_numpy(dtype=float).ravel()
                              for t in tables.values()])
    all_lfc = all_lfc[np.isfinite(all_lfc)]
    if all_lfc.size == 0:
        raise ValueError(f"No finite LFC values for contrast '{plot_name}'")
    st = float(np.round((np.floor(all_lfc.min() / step) - 1) * step, 10))
    en = float(np.round((np.ceil(all_lfc.max() / step) + 1) * step, 10))

    # pooled histogram per category on the shared bins
    groups = {}
    binn = None
    for cat, tab in tables.items():
        binn, hiss, perc = compute_hiss_LFC(tab, st, en, step)
        pooled_perc = np.nanmean(perc, axis=1) if perc.ndim == 2 else perc
        groups[cat] = {'table': tab, 'perc': pooled_perc, 'hiss' : hiss}

    T_target  = groups[target_type]['table']
    hist = groups[target_type]['hiss']
    T_control = groups[control_type]['table'] if control_present else pd.DataFrame(columns=['gRNA', 'gene'])

    # plot: target + control + any other categories (skips None automatically)
    plot_categories = {cat: groups[cat]['perc'] for cat in categories}

    distri_target_contr_plots_all(
        binn, plot_categories,
        highlight=target_type,                
        title=f"{plot_name}",
        output_dir=output_dir,
    )

    return T_target, T_control, hist, binn, st, en