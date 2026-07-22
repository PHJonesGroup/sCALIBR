from .compute_hiss_LFC_rep12 import compute_hiss_LFC_rep12
from .distri_target_contr_plots_all import distri_target_contr_plots_all
import numpy as np
import pandas as pd


def _lfc_table(tab, baseline, cond1, rep_pairs):
    """
    Per-rep LFC table [gRNA, gene, lfc_<r>...] for a subset of rows.
    If rep_pairs is None/empty, columns are used directly (no rep suffix):
    one LFC column from baseline vs cond1.
    """
    eps = 1e-6
    out = pd.DataFrame({
        'gRNA': tab.iloc[:, 0].values,
        'gene': tab.iloc[:, 1].astype(str).values,
    })

    if not rep_pairs:                      # None or empty -> no replicate suffixes
        c1, c2 = baseline, cond1
        if c1 not in tab.columns or c2 not in tab.columns:
            raise ValueError(f"Missing {c1} or {c2}. Columns: {tab.columns.tolist()}")
        out['lfc'] = np.log2(
            (tab[c2].astype(float) + eps) / (tab[c1].astype(float) + eps)
        ).values
    else:                                  # replicate suffixes present
        for r in rep_pairs:
            c1, c2 = f"{baseline}_{r}", f"{cond1}_{r}"
            if c1 not in tab.columns or c2 not in tab.columns:
                raise ValueError(f"Missing {c1} or {c2}")
            out[f"lfc_{r}"] = np.log2(
                (tab[c2].astype(float) + eps) / (tab[c1].astype(float) + eps)
            ).values
    return out


def separate_target_control(st, en, step, T_norm_indiv,
                            target_type, control_type,
                            baseline, cond1, rep_pairs, output_dir,
                            gene_type_column="gene_type"):
    """
    Split gRNAs by the gene_type column (categories discovered from the data),
    compute per-rep LFC distributions per category, and plot them.

    target_type  : gene_type value denoting targets  (config, e.g. 'GOI')
    control_type : gene_type value to calibrate on    (config, e.g. 'Zero-expressed gene')
    Works with any number of replicates (rep_pairs) and any set of categories.

    Returns
    -------
    T_target  : per-rep LFC table for targets       [gRNA, gene, lfc_<r>...]
    T_control : per-rep LFC table for control_type   [gRNA, gene, lfc_<r>...]
    groups    : dict {gene_type value -> {'table', 'perc'}} for every category present
    binn      : histogram bin edges
    """
    if gene_type_column not in T_norm_indiv.columns:
        raise KeyError(f"'{gene_type_column}' not found. Columns: {T_norm_indiv.columns.tolist()}")

    gt = T_norm_indiv[gene_type_column].astype(str).str.strip()
    categories = sorted(gt.unique())

    # target must exist; control must exist ONLY if you intend to calibrate on it
    if target_type not in categories:
        raise ValueError(f"target_type '{target_type}' not in {gene_type_column}. "
                         f"Available: {categories}")

    control_present = control_type in categories
    if not control_present:
        # e.g. dataset has no zero-expressed genes — warn, don't crash
        print(f"[warning] control_type '{control_type}' not found in data "
              f"(available: {categories}). No control group built for it.")

    # one per-rep LFC table + pooled histogram per discovered category
    groups = {}
    binn = None
    for cat in categories:
        rows = T_norm_indiv[gt == cat]
        tab = _lfc_table(rows, baseline, cond1, rep_pairs)
        # pooled histogram: flatten all rep LFC columns into one distribution
        binn, hiss, perc, LFC, lfc_cols, st1, en1 = compute_hiss_LFC_rep12(tab, st, en, step)
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