import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from scipy.stats import spearmanr, mannwhitneyu
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

MAX_ZERO_FRAC_BASELINE = 0.01    # zero-count gRNAs in baseline (T0/plasmid) samples
MAX_GINI_RATIO         = 1.5     # Gini vs median Gini of the other samples in its group
MIN_DEPTH_RATIO        = 0.5     # reads per gRNA vs median of the other samples in its group
MAX_CTRL_SD_RATIO      = 1.5     # control SD vs median control SD of the other replicates
MAX_CTRL_TARGET_SD     = 0.8     # control SD / target SD: controls nearly as wide as targets
MIN_SIGNAL_GAP         = 0.1     # LOO Spearman, targets minus controls


def gini(x):
    """Gini index of non-negative counts: 0 = perfectly even, 1 = all in one gRNA."""
    x = np.sort(np.asarray(x, dtype=float))
    n, s = x.size, x.sum()
    if n == 0 or s == 0:
        return np.nan
    return 2 * np.sum(np.arange(1, n + 1) * x) / (n * s) - (n + 1) / n


def sample_metrics(raw, baseline_samples, output_dir):
    """
    Count-level QC per sample, on raw counts after the keep_counts filter.

    Parameters
    ----------
    raw : pandas.DataFrame
        3 ID columns (gRNA, gene, gene_type) followed by raw count columns.
    baseline_samples : collection of str
        Samples used as a baseline in any contrast (stricter thresholds).
    output_dir : str
        Where qc_sample_metrics.csv is written.

    Returns
    -------
    pandas.DataFrame
        One row per sample: total reads, reads per gRNA, zero-count fraction,
        Gini index, and a 'flags' column ('' when the sample passes).
    """
    rows = []
    for s in raw.columns[3:]:
        c = raw[s].astype(float).to_numpy()
        rows.append({"sample": s, "baseline": s in baseline_samples,
                     "total_reads": int(c.sum()), "reads_per_gRNA": c.sum() / len(c),
                     "zero_count_frac": float(np.mean(c == 0)), "gini": gini(c)})
    out = pd.DataFrame(rows)

    flags = []
    for i, r in out.iterrows():
        peers = out[(out["baseline"] == r["baseline"]) & (out.index != i)]
        f = []
        if r["baseline"] and r["zero_count_frac"] > MAX_ZERO_FRAC_BASELINE:
            f.append(f"zero-count gRNAs {r['zero_count_frac']:.1%} in a baseline sample")
        if len(peers) and r["gini"] > MAX_GINI_RATIO * peers["gini"].median():
            f.append(f"Gini {r['gini']:.2f}, {r['gini'] / peers['gini'].median():.1f}x similar samples")
        if len(peers) and r["reads_per_gRNA"] < MIN_DEPTH_RATIO * peers["reads_per_gRNA"].median():
            f.append(f"depth {r['reads_per_gRNA'] / peers['reads_per_gRNA'].median():.1f}x similar samples")
        flags.append("; ".join(f))
    out["flags"] = flags
    out.to_csv(os.path.join(output_dir, "qc_sample_metrics.csv"), index=False)
    return out


def _ssmd_auroc(pos, ctrl):
    """SSMD and AUROC of positive controls vs controls (depletion: pos expected lower)."""
    pos, ctrl = pos[np.isfinite(pos)], ctrl[np.isfinite(ctrl)]
    if pos.size < 2 or ctrl.size < 2:
        return np.nan, np.nan
    ssmd = (pos.mean() - ctrl.mean()) / np.sqrt(pos.var(ddof=1) + ctrl.var(ddof=1))
    u = mannwhitneyu(ctrl, pos, alternative="two-sided").statistic
    return ssmd, u / (pos.size * ctrl.size)      # P(control > positive)


def replicate_metrics(lfc_all, target_type, control_types, positive_type,
                      name, output_dir):
    """
    Replicate-level QC for one contrast, on per-replicate LFCs of all gene types.

    Controls need not be centred on 0 (targets are measured relative to them);
    what matters is that they are tight, so flags are based on their spread.

    Per replicate:
      - control median (reported only) and SD (noise), and SD relative to the
        other replicates;
      - control SD / target SD: controls should be much tighter than targets;
      - leave-one-out Spearman correlation with the mean of the other replicates,
        for targets and for controls; the gap is the shared signal (the control
        correlation alone is reported, not flagged);
      - SSMD and AUROC of positive controls vs controls, if positive_type is set.

    Parameters
    ----------
    lfc_all : pandas.DataFrame
        [gRNA, gene, gene_type, lfc_<rep>, ...] for all gene types.
    target_type : str
    control_types : list of str
    positive_type : str or None
        gene_type of positive controls (e.g. essential genes); None to skip.
    name : str
        Contrast name, used in file names and titles.
    output_dir : str
        Where qc_replicate_metrics_<name>.csv and .pdf are written.

    Returns
    -------
    pandas.DataFrame
        One row per replicate, with a 'flags' column ('' when it passes).
    """
    gt = lfc_all["gene_type"].astype(str).str.strip()
    is_t = (gt == target_type).to_numpy()
    is_c = gt.isin(control_types).to_numpy()
    is_p = (gt == positive_type).to_numpy() if positive_type else np.zeros(len(gt), bool)

    reps = list(lfc_all.columns[3:])
    L = lfc_all[reps].to_numpy(dtype=float)
    ctrl_sd = np.array([np.nanstd(L[is_c, j], ddof=1) for j in range(len(reps))])

    rows = []
    for j, r in enumerate(reps):
        others = [k for k in range(len(reps)) if k != j]
        t_sd = np.nanstd(L[is_t, j], ddof=1)
        row = {"replicate": r,
               "n_control": int(is_c.sum()),
               "control_median": np.nanmedian(L[is_c, j]),
               "control_sd": ctrl_sd[j],
               "control_sd_ratio": ctrl_sd[j] / np.median(ctrl_sd[others]) if others else np.nan,
               "control_target_sd": ctrl_sd[j] / t_sd if t_sd > 0 else np.nan,
               "loo_spearman_targets": np.nan, "loo_spearman_controls": np.nan}
        if others:
            mean_others = np.nanmean(L[:, others], axis=1)
            for key, m in [("loo_spearman_targets", is_t), ("loo_spearman_controls", is_c)]:
                if m.sum() > 2:
                    row[key] = spearmanr(L[m, j], mean_others[m], nan_policy="omit").correlation
        row["signal_gap"] = row["loo_spearman_targets"] - row["loo_spearman_controls"]
        if is_p.any():
            row["ssmd_pos_vs_ctrl"], row["auroc_pos_vs_ctrl"] = _ssmd_auroc(L[is_p, j], L[is_c, j])

        flags = []
        if row["control_sd_ratio"] > MAX_CTRL_SD_RATIO:
            flags.append(f"control SD {row['control_sd_ratio']:.1f}x other replicates")
        if row["control_target_sd"] > MAX_CTRL_TARGET_SD:
            flags.append(f"controls nearly as dispersed as targets "
                         f"(SD ratio {row['control_target_sd']:.2f})")
        if row["signal_gap"] < MIN_SIGNAL_GAP:
            flags.append(f"no shared signal over controls (gap {row['signal_gap']:.2f})")
        row["flags"] = "; ".join(flags)
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(output_dir, f"qc_replicate_metrics_{name}.csv"), index=False)
    _plot(out, name, output_dir)
    return out


def _plot(m, name, output_dir):
    """One panel per replicate check; each bar red when beyond its dashed threshold."""
    x = np.arange(len(m))
    ok, bad = "#1D9E75", "#D85A30"
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2))

    v = m["control_sd_ratio"]
    axs[0].bar(x, v, color=np.where(v > MAX_CTRL_SD_RATIO, bad, ok))
    axs[0].axhline(MAX_CTRL_SD_RATIO, ls="--", c="grey", lw=1)
    axs[0].axhline(1, c="grey", lw=0.5)
    axs[0].set_title("Control SD vs other replicates")

    v = m["control_target_sd"]
    axs[1].bar(x, v, color=np.where(v > MAX_CTRL_TARGET_SD, bad, ok))
    axs[1].axhline(MAX_CTRL_TARGET_SD, ls="--", c="grey", lw=1)
    axs[1].set_title("Control SD / target SD")
    axs[1].set_ylim([0,1])

    w = 0.38
    gap = (m["signal_gap"] < MIN_SIGNAL_GAP).to_numpy()
    axs[2].bar(x - w / 2, m["loo_spearman_targets"], w, color=np.where(gap, bad, ok))
    axs[2].bar(x + w / 2, m["loo_spearman_controls"], w, color="#888780")
    axs[2].axhline(0, c="grey", lw=0.5)
    axs[2].set_title("Leave-one-out Spearman ρ")

    # legends under each panel
    threshold = lambda t: Line2D([], [], ls="--", c="grey", lw=1, label=f"flag threshold ({t})")
    legends = [
        [Patch(color=ok, label="pass"), Patch(color=bad, label="flagged"),
         threshold(f"{MAX_CTRL_SD_RATIO}×"), Line2D([], [], c="grey", lw=0.5, label="equal to others (1)")],
        [Patch(color=ok, label="pass"), Patch(color=bad, label="flagged"),
         threshold(MAX_CTRL_TARGET_SD)],
        [Patch(color=ok, label="targets, pass"),
         Patch(color=bad, label=f"targets, flagged (< {MIN_SIGNAL_GAP} above controls)"),
         Patch(color="#888780", label="controls")],
    ]
    for ax, handles, ncol in zip(axs, legends, [2, 2, 1]):
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.28),
                  ncol=ncol, fontsize=9, frameon=False)

    for ax in axs:
        ax.set_xticks(x)
        ax.set_xticklabels(m["replicate"], rotation=30, ha="right")
        ax.spines[["top", "right"]].set_visible(False)
    plt.suptitle(f"Replicate QC: {name}")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"qc_replicate_{name}.pdf"), bbox_inches="tight")
    plt.close()
