from .compute_hiss_LFC import compute_hiss_LFC
from .distri_target_contr_plots_all import distri_target_contr_plots_all
from .lfc_table import lfc_table

import numpy as np
import pandas as pd

def separate_target_control(st, en, step, T_norm_indiv,
                            target_type, control_type,
                            baseline, cond1, rep_pairs, output_dir,
                            gene_type_column="gene_type"):
    """
    Split gRNAs by the gene_type column, compute per-rep LFC distributions per category, and plot them.

    Parameters
    ----------
    st, en : float
        Lower / upper bounds of the histogram (LFC) range.
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
    rep_pairs : list of str or None
        Name of replicates
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
    groups : dict
        Maps each gene_type value -> {'table': DataFrame, 'perc': 1-D ndarray}
        for every category present in the data.
    binn : numpy.ndarray
        Histogram bin edges (shared across all categories)
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

    # one per-rep LFC table + pooled histogram per discovered category
    groups = {}
    binn = None

    for cat in categories:
        rows = T_norm_indiv[gt == cat]
        tab = lfc_table(rows, baseline, cond1, rep_pairs)
        # pooled histogram: flatten all rep LFC columns into one distribution
        binn, hiss, perc, LFC, lfc_cols, st1, en1 = compute_hiss_LFC(tab, st, en, step)
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
        title=f"{cond1}_vs_{baseline}",
        output_dir=output_dir,
    )

    return T_target, T_control, hist, groups, binn