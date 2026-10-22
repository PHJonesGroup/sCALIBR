import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from .make_histo_crit_stats import make_histo_crit_stats
from .q_val_frequentist_critical import q_val_frequentist_critical
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

def CTR_stats(alf, st, en, step, T_control, hist, baseline, cond, plot_name,
              control, norm_method, output_dir):
    """Calibrate a control null and score targets against it.

    Two calibration modes:
      'pooled' — all control LFC values across replicates are pooled into one
                 distribution, giving a single null (one p-curve, one mean/SD)
                 shared by every replicate.
      'sample' — each replicate column is calibrated against its own control
                 LFC distribution, giving a separate null per replicate.

    Parameters
    ----------
    alf : float
        Significance level, passed to the p-curve calculation.
    st, en : float
        Lower / upper bounds of the histogram (LFC) range.
    step : float
        Histogram bin width.
    T_control : pandas.DataFrame
        Control table: [gRNA, gene, gene_type, lfc_<rep>, ...].
    hist : numpy.ndarray
        Target histogram; 1-D counts, or 2-D (n_bins x n_reps).
    baseline, cond : str
        Condition labels.
    plot_name : str
        Contrast name, used in titles and filenames.
    control : str
        Control category label (e.g. 'Zero-expressed gene').
    norm_method : {'pooled', 'sample'}
        Whether to build one null from all replicates, or one null per replicate.
    output_dir : str
        Directory where plots and tables are saved.

    Returns
    -------
    me_sd : numpy.ndarray
        'pooled': [mean, SD]. 'sample': (n_rep, 2), one row per replicate.
    binn : numpy.ndarray
        Histogram bin edges (shared across modes and replicates).
    p_cont : numpy.ndarray
        'pooled': (n_bins,). 'sample': (n_bins, n_rep), one column per replicate.
    hiss_cont : numpy.ndarray
        Control histogram counts, same convention as p_cont.
    p_targ : numpy.ndarray
        Target p-curve (per-bin probability). 'pooled': (n_bins,).
        'sample': (n_bins, n_rep), one column per replicate.
    lfc_cols : list of str
        The replicate LFC column names, in the order used for the per-replicate
        axis of the arrays above (empty-meaning for 'pooled', but returned for
        reference).
    """
    lfc_cols = [c for c in T_control.columns if c.startswith('lfc')]
    if not lfc_cols:
        raise ValueError(f"No lfc columns in control table: "
                         f"{T_control.columns.tolist()}")

    if norm_method == "pooled":
        (binn, hiss_z, perc, bin_pz, me_sd, Z, n_z
         ) = make_histo_crit_stats(alf, st, en, step, T_control,
                                   baseline, cond, plot_name, control, output_dir)
        p_cont = bin_pz[:, 1]
        hiss_cont = hiss_z

    elif norm_method == "sample":
        p_list, hiss_list, mesd_list = [], [], []
        binn = None
        for col in lfc_cols:
            sub = T_control[list(T_control.columns[:3]) + [col]]
            (binn, hiss_z, perc, bin_pz, me_sd_c,
             Z, n_z) = make_histo_crit_stats(
                alf, st, en, step, sub, baseline, cond,
                f"{plot_name}_{col}", control, output_dir)
            p_list.append(bin_pz[:, 1])
            hiss_list.append(hiss_z)
            mesd_list.append(me_sd_c)

        p_cont    = np.column_stack(p_list)       # (n_bins, n_rep)
        hiss_cont = np.column_stack(hiss_list)
        me_sd     = np.vstack(mesd_list)          # (n_rep, 2)

    else:
        raise ValueError(f"norm_method must be 'pooled' or 'sample', got '{norm_method}'")

    # target p-curve: pooled -> one curve from all reps; sample -> one per rep
    hist = np.asarray(hist)
    if norm_method == "pooled":
        if hist.ndim == 2:
            hist = hist.sum(axis=1)
        p_targ = q_val_frequentist_critical(alf, binn, hist)
    else:
        if hist.ndim != 2 or hist.shape[1] != len(lfc_cols):
            raise ValueError(f"target hist has shape {hist.shape}; expected "
                             f"(n_bins, {len(lfc_cols)}) for lfc columns {lfc_cols}")
        p_targ = np.column_stack([q_val_frequentist_critical(alf, binn, hist[:, j])
                                  for j in range(hist.shape[1])])   # (n_bins, n_rep)

    # save the p-curves
    cols = {"bin": np.round(binn, 2)}
    if norm_method == "pooled":
        cols["target"] = np.round(p_targ, 4)
        cols[f"control_{control}"] = np.round(p_cont, 4)
    else:
        for j, col in enumerate(lfc_cols):
            cols[f"target_{col}"] = np.round(p_targ[:, j], 4)
            cols[f"control_{control}_{col}"] = np.round(p_cont[:, j], 4)
    pd.DataFrame(cols).to_csv(
        os.path.join(output_dir, f"p_distri_targ_cont_{plot_name}_data.csv"), index=False)

    plt.figure()
    if norm_method == "pooled":
        plt.plot(binn, p_cont, color="0.5", linestyle="--", label=f"control ({control})")
        plt.plot(binn, p_targ, color="C3", linestyle="-", label="target")
    else:
        # same colour per replicate: dashed = control, solid = target
        for j, col in enumerate(lfc_cols):
            plt.plot(binn, p_cont[:, j], color=f"C{j}", linestyle="--", label=f"control {col}")
            plt.plot(binn, p_targ[:, j], color=f"C{j}", linestyle="-",  label=f"target {col}")
    plt.title(f"Target vs. {control} gRNAs ({plot_name}, {norm_method})")
    plt.xlabel("LFC bin")
    plt.ylabel("Probability")
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(output_dir, f"p_distri_targ_cont_{plot_name}.pdf"),
                bbox_inches="tight")
    plt.close()

    return me_sd, binn, p_cont, hiss_cont, p_targ, lfc_cols