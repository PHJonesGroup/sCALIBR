import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

def find_columns_indiv(norm_tab, indiv):
    """
    Find columns containing individual identifiers like 'F' or 'M' for scheme 1.
    Returns column indices.
    """
    col_labels = norm_tab.columns.tolist()
    matches = [i for i, col in enumerate(col_labels) if indiv in col]
    if len(matches) < 2:
        print(f"Warning: Could not find two matching columns for '{indiv}'")
    return matches, len(matches)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def module1_format_data(index_scheme, indiv, norm_tab):
    """
    Reformats data depending on scheme (1 or 2) and computes mean/median
    to visualize T0 vs T1 read‑count changes with overlaid bars.
    """
    siz = norm_tab.shape
    print(f"Shape of norm_tab: {siz}")

    norm_dat       = None
    me_med_nd      = None
    T_norm_indiv   = None
    norm_dat_values = None      # keep for return

    if siz[1] > 0:
        # ---------- compute mean & median across rows (per column) ----------
        norm_dat_values = norm_tab.iloc[:, 2:].values  # drop gRNA, gene
        me_nd  = np.mean   (norm_dat_values, axis=0)
        med_nd = np.median (norm_dat_values, axis=0)
        me_med_nd = np.vstack([me_nd, med_nd])

        # ---------- visualisation ----------
        indices   = np.arange(len(me_nd))
        bar_width = 0.6                         # single (shared) bar width

        plt.figure(figsize=(10, 4))
        plt.bar( indices,  me_nd,
                 width     = bar_width,
                 color     = 'tab:blue',
                 alpha     = 0.7,
                 label='Mean'   )

        plt.bar( indices,  med_nd,
                 width     = bar_width,
                 color     = 'tab:orange',
                 alpha     = 0.7,
                 label='Median' )

        plt.ylabel('RPKM')
        plt.title('Mean *vs* Median of T0 and T1 (May 2 mice)')
        plt.xlabel('T0 replicates                                 T1 replicates')
        plt.xticks([])               # hide tick labels (purely aesthetic here)
        plt.legend()
        plt.tight_layout()
        plt.show()

        # ---------- scheme‑specific table ----------
        if index_scheme == 2:        # use full table
            T_norm_indiv = norm_tab.copy()

        elif index_scheme == 1:      # pick columns for one individual
            print('Finding columns for this individual …')
            ind_found, num_found = find_columns_indiv(norm_tab, indiv)
            if num_found >= 2:
                T0   = norm_tab.iloc[:, ind_found[0]]
                T1   = norm_tab.iloc[:, ind_found[1]]
                T_norm_indiv = pd.DataFrame({
                    'gRNA' : norm_tab.iloc[:, 0],
                    'gene' : norm_tab.iloc[:, 1],
                    'T0_1' : T0,
                    'T0_2' : T0,
                    'T1_1' : T1,
                    'T1_2' : T1
                })
                norm_dat = T_norm_indiv.iloc[:, 2:].values
            else:
                print("Error: Not enough matching columns for selected individual")

    # Return exactly what the old MATLAB version produced
    return (
        T_norm_indiv,
        norm_dat_values if index_scheme == 2 else norm_dat,
        me_med_nd
    )


def filter_pattern_distri(raw_ind, pat, n, st, en, step):
    """
    Python equivalent of MATLAB's filter_pattern_distri.m
    """

    # 1. Extract genes and gRNAs
    genes = raw_ind.iloc[:, 1].astype(str)
    gRNA = raw_ind.iloc[:, 0]

    # 2. Apply pattern filter using only the first n characters
    def matches_pattern(g):
        return g[:n].lower() == pat.lower()
    
    matches = genes.apply(matches_pattern)
    ind_out = matches[matches].index  # pattern-matching rows
    ind_in = matches[~matches].index  # rest

    indiv_Chr = raw_ind.loc[ind_out]
    indiv_noChr = raw_ind.loc[ind_in]
    num_out_in = len(ind_out)

    # 3. Compute LFCs from columns 3 to 6
    # Assume: col3 = rep1_T1, col4 = rep1_T2, col5 = rep2_T1, col6 = rep2_T2
    # Correct LFCs: T1_F / T0_F and T1_M / T0_M
    counts = indiv_Chr[['T0_F', 'T0_M', 'T1_F', 'T1_M']].astype(float)

    LFC_1 = np.log2((counts['T1_F'] + 1e-6) / (counts['T0_F'] + 1e-6))
    LFC_2 = np.log2((counts['T1_M'] + 1e-6) / (counts['T0_M'] + 1e-6))

    LFC_sum = LFC_1 + LFC_2

    bins = np.arange(st, en + step, step)
    hissFM, _ = np.histogram(LFC_sum, bins)
    percFM = 100 * hissFM / hissFM.sum() if hissFM.sum() > 0 else np.zeros_like(hissFM)

    # Optional: total raw counts (for reference, not returned)
    s_chr = counts.sum().values
    st1, en1 = st, en

    # 4. Construct LFC table
    T_lfc_pat = pd.DataFrame({
        'gRNA_chr': gRNA.loc[ind_out].values,
        'genes_chr': genes.loc[ind_out].values,
        'lfc1': LFC_1.values,
        'lfc2': LFC_2.values
    })

    return (
        indiv_Chr, indiv_noChr,
        ind_out.tolist(), ind_in.tolist(), num_out_in,
        bins, hissFM, percFM,
        LFC_1.values, LFC_2.values,
        s_chr, st1, en1, T_lfc_pat
    )

def zGE_target_distri(d, genes_nnt, zGE_genes, T_norm, st, en, step):
    """
    Separate zGE and non-zGE target genes and compute LFC distributions.
    """
    zGE_set = set(zGE_genes)
    nzGE_set = set(T_norm.iloc[:, 1]) - zGE_set

    # Select rows for zGE genes
    T_zGE = T_norm[T_norm.iloc[:, 1].isin(zGE_set)].copy()
    # Select rows for non-zGE genes
    T_nzGE = T_norm[T_norm.iloc[:, 1].isin(nzGE_set)].copy()

    # Compute LFCs for zGE genes
    counts_zGE = T_zGE.iloc[:, 2:6].astype(float)
    LFC_1 = np.log2((counts_zGE['T1_F'] + 1e-6) / (counts_zGE['T0_F'] + 1e-6))
    LFC_2 = np.log2((counts_zGE['T1_M'] + 1e-6) / (counts_zGE['T0_M'] + 1e-6))
    LFC_sum = LFC_1 + LFC_2

    bins = np.arange(st, en + step, step)
    hissFM, _ = np.histogram(LFC_sum, bins)
    percFM = 100 * hissFM / hissFM.sum() if hissFM.sum() > 0 else np.zeros_like(hissFM)

    hiss1, _ = np.histogram(LFC_1, bins)
    hiss2, _ = np.histogram(LFC_2, bins)
    perc1 = 100 * hiss1 / hiss1.sum() if hiss1.sum() > 0 else np.zeros_like(hiss1)
    perc2 = 100 * hiss2 / hiss2.sum() if hiss2.sum() > 0 else np.zeros_like(hiss2)

    T_lfc_zGE = pd.DataFrame({
        'gRNA': T_zGE.iloc[:, 0].values,
        'gene': T_zGE.iloc[:, 1].values,
        'lfc1': LFC_1.values,
        'lfc2': LFC_2.values
    })

    # Compute LFCs for non-zGE genes
    counts_nzGE = T_nzGE.iloc[:, 2:6].astype(float)
    LFC_1t = np.log2((counts_nzGE['T1_F'] + 1e-6) / (counts_nzGE['T0_F'] + 1e-6))
    LFC_2t = np.log2((counts_nzGE['T1_M'] + 1e-6) / (counts_nzGE['T0_M'] + 1e-6))
    LFC_sum_t = LFC_1t + LFC_2t

    hissFMt, _ = np.histogram(LFC_sum_t, bins)
    percFMt = 100 * hissFMt / hissFMt.sum() if hissFMt.sum() > 0 else np.zeros_like(hissFMt)

    hiss1t, _ = np.histogram(LFC_1t, bins)
    hiss2t, _ = np.histogram(LFC_2t, bins)
    perc1t = 100 * hiss1t / hiss1t.sum() if hiss1t.sum() > 0 else np.zeros_like(hiss1t)
    perc2t = 100 * hiss2t / hiss2t.sum() if hiss2t.sum() > 0 else np.zeros_like(hiss2t)

    T_lfc_nzGE = pd.DataFrame({
        'gRNA': T_nzGE.iloc[:, 0].values,
        'gene': T_nzGE.iloc[:, 1].values,
        'lfc1': LFC_1t.values,
        'lfc2': LFC_2t.values
    })

    gzn = [len(zGE_set), len(nzGE_set)]

    return (
        T_lfc_zGE, T_lfc_nzGE,
        bins,
        hiss1, perc1, hiss2, perc2, hissFM, percFM,
        hiss1t, perc1t, hiss2t, perc2t, hissFMt, percFMt,
        gzn
    )


def distri_target_contr_plots_all(binn, perc_z, perc_i, perc_t, perc_n, cond1):
    """
    Plot LFC distribution for zGE, intergenic, target, and NT control.
    """
    # Plot 1: Separate subplots
    plt.figure(figsize=(10, 8))
    
    plt.subplot(3, 1, 1)
    plt.bar(binn, perc_z, width=np.diff(binn)[0], color='c')
    plt.grid(True)
    plt.title(f'LFC distribution zGE, one replica: {cond1}', fontsize=14)
    
    plt.subplot(3, 1, 2)
    plt.bar(binn, perc_i, width=np.diff(binn)[0], color='m')
    plt.grid(True)
    plt.title('LFC distribution intergenic', fontsize=12)
    
    plt.subplot(3, 1, 3)
    plt.bar(binn, perc_t, width=np.diff(binn)[0], color='k')
    plt.grid(True)
    plt.title('LFC distribution target genes', fontsize=12)
    plt.xlabel('LFC', fontsize=12)
    
    plt.tight_layout()
    plt.show()

    # Plot 2: Overlaid plot
    
    plt.figure(figsize=(10, 5))
    plt.bar(binn, perc_z, width=np.diff(binn)[0], color='c', label='zGE')
    plt.bar(binn, perc_i, width=np.diff(binn)[0], color='m', alpha=0.7, label='intergenic')
    plt.bar(binn, perc_t, width=np.diff(binn)[0], color='k', alpha=0.5, label='target')
    plt.bar(binn, perc_n, width=np.diff(binn)[0], color='g', alpha=0.3, label='NT control')

    plt.grid(True)
    plt.title(f'LFC distribution - Target Controls (one replica): {cond1}', fontsize=14)
    plt.xlabel('LFC', fontsize=12)
    plt.ylabel('Percentage (%)', fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return 1


def separate_target_control(d, st, en, step, T_norm_indiv, zGE, pat1, pat2, n):
    """
    Python equivalent of MATLAB's separate_target_control.m

    Parameters:
    - d: dataset or config used in zGE_target_distri
    - st, en, step: for histogram binning
    - T_norm_indiv: full normalized gRNA table
    - zGE: list (table-like) of zGE genes (first column = gene names)
    - pat1: pattern for intergenic control (e.g., 'chr')
    - pat2: pattern for NT control (e.g., 'Non')
    - n: number of characters to match pattern

    Returns:
    - T_target: target genes (excluding controls + zGE)
    - T_target_zGE: target + zGE genes
    - T_lfc_chr: intergenic control genes
    - T_lfc_nt: NT control genes
    - T_zGE: genes from zGE list (sorted)
    - bin: histogram bins
    - All histograms (his*, perc*) and gzn table
    """

    # Initialize outputs
    T_target = None
    T_target_zGE = None
    T_lfc_chr = None
    T_lfc_nt = None
    T_zGE = None
    gzn = None

    # -----------------------
    # 1. Separate intergenic (e.g., 'chr')
    pat = pat1
    (
        indiv_Chr, indiv_noChr,
        ind_out1, ind_in1, num_out_in1,
        bin, hisFMi, percFMi,
        LFC_1i, LFC_2i,
        s_chr, st1, en1, T_lfc_chr
    ) = filter_pattern_distri(T_norm_indiv, pat, n, st, en, step)
    
    num_chr_rest = num_out_in1
    print("Remaining after removing 'chr' genes:", num_chr_rest)

    # -----------------------
    # 2. Separate NT (e.g., 'Non')
    pat = pat2
    (
        indiv_NT, indiv_noChr_noNT,
        ind_out2, ind_in2, num_out_in2,
        bin, hisFMnt, percFMnt,
        LFC_1nt, LFC_2nt,
        s_nt, st1, en1, T_lfc_nt
    ) = filter_pattern_distri(indiv_noChr, pat, n, st, en, step)

    num_NT_rest = num_out_in2
    print("Remaining after removing NT genes:", num_NT_rest)

    # -----------------------
    # 3. Separate zGE and target genes
    zGE_genes = zGE.iloc[:, 0].tolist()
    genes_nnt = indiv_noChr_noNT.iloc[:, 1].tolist()  # gene column

    (
        T_zGE, T_nzGE, bin,
        his1z, perc1z, his2z, perc2z, hisFMz, percFMz,
        his1t, perc1t, his2t, perc2t, hisFMt, percFMt,
        gzn
    ) = zGE_target_distri(d, genes_nnt, zGE_genes, indiv_noChr_noNT, st, en, step)

    # -----------------------
    # 4. Compute distributions for chr and NT controls
    bin, his1i,perc1i,his2i,perc2i,hisFMi,percFMi,LFC_1i,LFC_2i, st1i, en1i = compute_hiss_LFC_rep12(T_lfc_chr, st, en, step)
    bin, his1n,perc1n,his2n,perc2n,hisFMn,percFMn,LFC_1n,LFC_2n, st1n, en1n = compute_hiss_LFC_rep12(T_lfc_nt, st, en, step)

    # -----------------------
    # 5. Plot histograms (female and male)
    cond1 = "F"
    cond2 = "M"
    distri_target_contr_plots_all(bin, perc1z, perc1i, perc1t, perc1n, cond1)
    distri_target_contr_plots_all(bin, perc2z, perc2i, perc2t, perc2n, cond2)

    # -----------------------
    # 6. Final Outputs
    T_target = T_nzGE  # targets only
    T_target_zGE = indiv_noChr_noNT  # target + zGE
    T_zGE = T_zGE.sort_values(by="gene")  # sort alphabetically

    return (
        T_target, T_target_zGE, T_lfc_chr, T_lfc_nt, T_zGE,
        bin, his1z, perc1z, his2z, perc2z, hisFMz, percFMz,
        his1t, perc1t, his2t, perc2t, hisFMt, percFMt, gzn
    )




import numpy as np
import pandas as pd

def compute_hiss_LFC_rep12(chr_lfc_12, st, en, step):
    """
    Compute histograms (distributions) for individual gRNAs, for both replicates and the average.
    
    Parameters:
        chr_lfc_12 (pd.DataFrame): DataFrame with columns:
            - gRNA_chr (str)
            - genes_chr (str)
            - lfc1 (float): replicate 1
            - lfc2 (float): replicate 2
        st (float): start of histogram range
        en (float): end of histogram range
        step (float): bin width

    Returns:
        binn: histogram bin edges
        hiss1: histogram counts for replicate 1
        perc1: histogram percentages for replicate 1
        hiss2: histogram counts for replicate 2
        perc2: histogram percentages for replicate 2
        hissFM: histogram counts for average
        percFM: histogram percentages for average
        LFC_1: replicate 1 LFCs (numpy array)
        LFC_2: replicate 2 LFCs (numpy array)
        st1: min(LFCs) - 1
        en1: max(LFCs) + 1
    """
    
    # Extract values
    LFC_1 = chr_lfc_12.iloc[:, 2].astype(float).values
    LFC_2 = chr_lfc_12.iloc[:, 3].astype(float).values

    # Compute new bounds
    st1 = min(LFC_1.min(), LFC_2.min()) - 1
    en1 = max(LFC_1.max(), LFC_2.max()) + 1

    # Compute average LFC (mean of both replicates)
    meLFC = np.mean(np.vstack((LFC_1, LFC_2)), axis=0)

    # Histograms
    binn, hissFM, percFM = make_histo_LFC(step, meLFC, st, en)
    _, hiss1, perc1 = make_histo_LFC(step, LFC_1, st, en)
    _, hiss2, perc2 = make_histo_LFC(step, LFC_2, st, en)

    return binn, hiss1, perc1, hiss2, perc2, hissFM, percFM, LFC_1, LFC_2, st1, en1



def make_histo_LFC(step, vec_num, st, en):
    """
    Compute histogram counts and percentages with fixed bin size.

    Parameters:
    - step: bin width
    - vec_num: 1D array-like numeric values (e.g., LFC values)
    - st: minimum value (start of bins)
    - en: maximum value (end of bins)

    Returns:
    - bin: array of bin edges excluding the last edge (like MATLAB's bin[:-1])
    - his: counts of values in each bin
    - perc: percentage of counts in each bin
    """

    # Create bins from st to en with step size
    bins = np.arange(st, en + step, step)

    # Initialize counts array (length = number of bins - 1)
    val = np.zeros(len(bins) - 1, dtype=int)

    # Count values in each bin (vectorized approach)
    for i in range(len(bins) - 1):
        val[i] = np.sum((vec_num >= bins[i]) & (vec_num < bins[i + 1]))

    # Compute percentage
    total = val.sum()
    perc = 100 * val / total if total > 0 else np.zeros_like(val)

    # Return bins excluding the last edge (to match MATLAB behavior)
    bin_return = bins[:-1]

    return bin_return, val, perc



def make_histo_crit_stats(
    alf: float,
    st: float,
    en: float,
    step: float,
    T_zGE_12: pd.DataFrame,
    cond1: str,
    cond2: str,
    cond: str,
):
    """
    Compute histograms, critical LFC thresholds, and robust stats for zGE
    controls in two replicates plus the pooled set.

    Parameters
    ----------
    alf   : significance level (e.g. 0.06)
    st/en : histogram lower / upper bounds for LFC
    step  : histogram bin width
    T_zGE_12 : DataFrame with columns [gRNA, gene, lfc1, lfc2]
    cond1, cond2, cond : strings for figure titles / legends

    Returns (same order/meaning as MATLAB)
    -------
    bin           : 1‑D ndarray of bin left‑edges
    his           : (#bins, 3) counts for rep1, rep2, pooled
    perc          : (#bins, 3) percentages for rep1, rep2, pooled
    crit_1_2_12   : (3, ?) matrix of critical LFC thresholds
    bin_p12b      : (#bins, 4)  [bin  p_rep1  p_rep2  p_pooled]
    med_mad       : (3, 2)   [median, MAD] for rep1, rep2, pooled
    me_sd         : (3, 2)   [mean,   SD]  for rep1, rep2, pooled
    mod           : length‑3 array of modes
    MZ, Z         : (#rows, 2) MZ / Z for rep1 & rep2
    MZ_12, Z_12   : 1‑D arrays for pooled replicate
    n1_n2_n12     : (#bins, 3) cumulative counts (from make_histo_vec_rep12)
    """

    # ------------------------------------------------------------------
    # --- 1.  Pull LFC vectors for two replicates -----------------------
    # ------------------------------------------------------------------
    LFC_1 = T_zGE_12.iloc[:, 2].astype(float).values
    LFC_2 = T_zGE_12.iloc[:, 3].astype(float).values
    LFC_12 = np.concatenate([LFC_1, LFC_2])

    # ------------------------------------------------------------------
    # --- 2.  Histograms for rep1 / rep2 / pooled ----------------------
    # ------------------------------------------------------------------
    (
        bin_edges,
        perc_1,
        perc_2,
        perc_12,
        his_1,
        his_2,
        his_12,
        n1_n2_n12,
    ) = make_histo_vec_rep12(step, LFC_1, LFC_2, st, en)

    perc = np.column_stack([perc_1, perc_2, perc_12])
    his  = np.column_stack([his_1,  his_2,  his_12])
    bin_ = bin_edges  # rename to match MATLAB name

    # ------------------------------------------------------------------
    # --- 3.  Critical LFC thresholds & p‑curves ------------------------
    # ------------------------------------------------------------------
    bin_p12b, crit_1_2_12 = compute_p_critLFC(
        alf, bin_, his_1, his_2, his_12, cond
    )

    # ------------------------------------------------------------------
    # --- 4.  Med‑MAD, Mean‑SD, Mode, Z‑scores, MZ‑scores --------------
    # ------------------------------------------------------------------
    med_mad_1, MZZ1, MZ_1, Z_1, me_sd_1, mod_1 = med_mad_MZNP_2(LFC_1)
    med_mad_2, MZZ2, MZ_2, Z_2, me_sd_2, mod_2 = med_mad_MZNP_2(LFC_2)
    med_mad_12, MZZ12, MZ_12, Z_12, me_sd_12, mod_12 = med_mad_MZNP_2(LFC_12)

    med_mad = np.vstack([med_mad_1, med_mad_2, med_mad_12])
    me_sd   = np.vstack([me_sd_1,   me_sd_2,   me_sd_12])
    mod     = np.array([mod_1, mod_2, mod_12])

    # Stack MZ/Z for two reps (pooled goes out separately)
    MZ = np.column_stack([MZ_1, MZ_2])
    Z  = np.column_stack([Z_1,  Z_2])

    # ------------------------------------------------------------------
    # --- 5.  Diagnostic plots (like MATLAB) ---------------------------
    # ------------------------------------------------------------------
    # --- Compute Z and MZ histograms just like MATLAB ---
    _, perc_z1, perc_z2, perc_z12, _, _, _, _ = make_histo_vec_rep12(step, Z_1, Z_2, st, en)
    _, perc_mz1, perc_mz2, perc_mz12, _, _, _, _ = make_histo_vec_rep12(step, MZ_1, MZ_2, st, en)

    # LFC / Z / MZ – replicate‑1
    _plot_three_panels(bin_, perc_1, perc_z1, perc_mz1,
                   title=f'distri LFC rep i, {cond1}', color='b', xlab='LFC')

    # replicate‑2
    _plot_three_panels(bin_, perc_2, perc_z2, perc_mz2,
                   title=f'distri LFC rep i+1, {cond2}', color='r', xlab='LFC')

    # ------------------------------------------------------------------
    return (
        bin_,
        his,
        perc,
        crit_1_2_12,
        bin_p12b,
        med_mad,
        me_sd,
        mod,
        MZ,
        Z,
        MZ_12,
        Z_12,
        n1_n2_n12,
    )


# ----------------------------------------------------------------------
# --- helper: replicate the 3‑row bar‑plot MATLAB produced -------------
# ----------------------------------------------------------------------
def _plot_three_panels(bin_edges, perc_lfc, perc_z=None, perc_mz=None,
                       title='', color='b', xlab='LFC'):
    """Quick 3‑row barplot helper (LFC / Z / MZ)."""
    plt.figure()
    plt.subplot(3, 1, 1)
    plt.bar(bin_edges, perc_lfc, width=np.diff(bin_edges)[0], color=color)
    plt.grid(True)
    plt.title(title, fontsize=14)
    plt.xlabel(xlab)

    if perc_z is not None:
        plt.subplot(3, 1, 2)
        plt.bar(bin_edges, perc_z, width=np.diff(bin_edges)[0], color=color)
        plt.grid(True)
        plt.title('distri Z LFC', fontsize=10)
        plt.xlabel('Z LFC')

    if perc_mz is not None:
        plt.subplot(3, 1, 3)
        plt.bar(bin_edges, perc_mz, width=np.diff(bin_edges)[0], color=color)
        plt.grid(True)
        plt.title('distri MZ LFC', fontsize=10)
        plt.xlabel('MZ LFC')

    plt.tight_layout()



def make_LFC_Z_MZ_tables_two(zGE_12: pd.DataFrame,
                             MZ12: np.ndarray,
                             Z12: np.ndarray):
    """
    Build per‑replicate tables that append Z‑score and MZ‑score columns
    to the zGE LFC data.

    Parameters
    ----------
    zGE_12 : pd.DataFrame
        Columns order: 0‑gRNA, 1‑gene, 2‑lfc1, 3‑lfc2
    MZ12   : ndarray, shape (n_rows, 2)
        MZ scores for replicate‑1 (col‑0) and replicate‑2 (col‑1)
    Z12    : ndarray, shape (n_rows, 2)
        Z  scores for replicate‑1 (col‑0) and replicate‑2 (col‑1)

    Returns
    -------
    T_zGE_1 : DataFrame  (replicate‑1)  columns: gRNA, gene, LFC, Z, MZ
    T_zGE_2 : DataFrame  (replicate‑2)  columns: gRNA, gene, LFC, Z, MZ
    """

    # shared identifiers
    gRNA = zGE_12.iloc[:, 0].values
    gene = zGE_12.iloc[:, 1].values

    # ---------- replicate‑1 ----------
    T_zGE_1 = pd.DataFrame({
        "gRNA": gRNA,
        "gene": gene,
        "LFC":  zGE_12.iloc[:, 2].astype(float).values,   # lfc1
        "Z":    Z12[:, 0].astype(float),
        "MZ":   MZ12[:, 0].astype(float),
    })

    # ---------- replicate‑2 ----------
    T_zGE_2 = pd.DataFrame({
        "gRNA": gRNA,
        "gene": gene,
        "LFC":  zGE_12.iloc[:, 3].astype(float).values,   # lfc2
        "Z":    Z12[:, 1].astype(float),
        "MZ":   MZ12[:, 1].astype(float),
    })

    return T_zGE_1, T_zGE_2



def q_val_frequentist_critical(alf: float,
                               bin_edges: np.ndarray,
                               his: np.ndarray):
    """
    Given a histogram of LFC values (counts per bin), compute   \
    • two‑tailed p‑curves (left & right),
    • critical LFC thresholds cL / cR at level `alf`,
    • median‑like intersection point,
    • and a diagnostics matrix `bin_pi`.

    Parameters
    ----------
    alf       : significance level (e.g. 0.05)
    bin_edges : 1‑D array of left‑edge bin positions (same length as `his`)
    his       : 1‑D array of counts per bin (same length as `bin_edges`)

    Returns
    -------
    p         : combined p‑curve (left then right), length == len(his)
    cL, cR    : critical LFC values (left/right tails)
    bin_pi    : (#bins, 6) matrix  [bin, p_right, p_left, his, cum_R, cum_L]
    med_LFCp  : LFC value at which left & right p‑curves intersect (~mode)
    his4p     : combined cumulative counts (left then right), len == len(his)
    """

    # ------------------------------------------------------------------
    thrLowR = -4.5        # hard‑coded right‑tail lower bound
    S = his.sum()         # total gRNA count
    N = len(his)
    step = bin_edges[1] - bin_edges[0]

    # Cumulative fractions (right & left)
    cum_R  = np.cumsum(his[::-1])[::-1]   # cumulative counts from right
    cum_L  = np.cumsum(his)               # cumulative counts from left
    cum_fracR = cum_R / S
    cum_fracL = cum_L / S

    # p‑curves
    p_right = cum_fracR
    p_left  = cum_fracL

    # --- critical right‑tail threshold --------------------------------
    crit_right_ind = np.where(p_right <= alf)[0][0]          # first ≤ alf
    lfc_crit_right = bin_edges[crit_right_ind]
    cR             = lfc_crit_right

    # If right tail is “too low”, tighten criterion (alf -> alf/2)
    delta = 0.0
    if lfc_crit_right < thrLowR:
        delta = alf / 2.0
        crit_right_ind = np.where(p_right <= (alf - delta))[0][0]
        lfc_crit_right = bin_edges[crit_right_ind]
        cR = lfc_crit_right

    # --- critical left‑tail threshold ---------------------------------
    crit_left_ind = np.where(p_left >= (alf + delta))[0][0]
    lfc_crit_left = bin_edges[crit_left_ind]
    cL            = lfc_crit_left

    # --- intersection (mode‑like) where |p_R - p_L| is minimal --------
    ind_min = np.argmin(np.abs(p_right - p_left))
    med_LFCp = bin_edges[ind_min]

    # combined p vector (left part then right part, matching MATLAB logic)
    p_combined = np.concatenate([p_left[:ind_min+1],
                                 p_right[ind_min+1:]])

    his4p = np.concatenate([cum_L[:ind_min+1],
                            cum_R[ind_min+1:]])

    # Diagnostics matrix (same layout as MATLAB)
    bin_pi = np.column_stack([bin_edges,
                              p_right,
                              p_left,
                              his,
                              cum_R,
                              cum_L])

    return (
        p_combined,        # p
        cL,
        cR,
        bin_pi,
        med_LFCp,
        his4p
    )


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def CTR_stats_zGE(
    alf: float,
    st: float,
    en: float,
    step: float,
    T_zGE_12: pd.DataFrame,
    his12t: np.ndarray,
    condz1: str,
    condz2: str,
    cond1: str,
    cond2: str,
    condz: str,
):
    """
    Analyse zGE controls vs. targets for two replicates.

    Returns
    -------
    crit_LR12   : (2, ?)  critical LFC limits (rep1, rep2)
    me_sd12     : (2, 2)  mean & SD per replicate
    med_mad12   : (2, 2)  median & MAD per replicate
    binn        : 1‑D bin edges
    p_cont12    : (#bins, 2) control p‑curves (rep1, rep2)
    hiss_cont12 : (#bins, 2) control histograms (rep1, rep2)
    p_targ12    : (#bins, 2) target‑gene p‑curves (rep1, rep2)
    T_zGE_1/2   : DataFrames with LFC, Z, MZ for zGE (rep1 / rep2)
    """

    # ------------------------------------------------------------------
    # 1. Histogram & critical stats for zGE controls
    # ------------------------------------------------------------------
    (
        binn,
        hiss_z12,
        perc_z12,
        crit_1_2_z12,
        bin_p12bz,
        med_mad_z12,
        me_sd_z12,
        mod_z12,
        MZ12,
        Z12,
        MZ_12,
        Z_12,
        n1_n2_n12_z,
    ) = make_histo_crit_stats(
        alf, st, en, step, T_zGE_12, condz1, condz2, condz
    )

    me_sd12   = me_sd_z12[:2, :]      # rows 0 & 1  (rep1, rep2)
    med_mad12 = med_mad_z12[:2, :]
    crit_LR12 = crit_1_2_z12[:2, :]

    # ------------------------------------------------------------------
    # 2. zGE tables with LFC / Z / MZ
    # ------------------------------------------------------------------
    T_zGE_1, T_zGE_2 = make_LFC_Z_MZ_tables_two(T_zGE_12, MZ12, Z12)

    # ------------------------------------------------------------------
    # 3. Control p‑curves & histograms (rep1, rep2)
    # ------------------------------------------------------------------
    p_cont12    = bin_p12bz[:, 1:3]    # columns 1 & 2 => rep1 / rep2
    p_cont1, p_cont2 = p_cont12.T
    hiss_cont12 = hiss_z12[:, :2]

    # ------------------------------------------------------------------
    # 4. Target‑gene p‑curves for each replicate
    # ------------------------------------------------------------------
    p_targ1, cL1, cR1, bin_pi1, med_LFCp1, his4p1 = q_val_frequentist_critical(
        alf, binn, his12t[:, 0]
    )
    p_targ2, cL2, cR2, bin_pi2, med_LFCp2, his4p2 = q_val_frequentist_critical(
        alf, binn, his12t[:, 1]
    )
    p_targ12 = np.column_stack([p_targ1, p_targ2])

    # ------------------------------------------------------------------
    # 5. Diagnostic plots (optional, comment out for batch runs)
    # ------------------------------------------------------------------
    plt.figure()
    plt.plot(binn, p_cont1, "k--", label="control")
    plt.plot(binn, p_targ1, "k",  label="target")
    plt.title(f"p‑control (--), p‑target (k), indiv gRNA, {cond1}")
    plt.xlabel("LFC bin"); plt.ylabel("probability to be extreme"); plt.grid(True)

    plt.figure()
    plt.plot(binn, p_cont2, "k--", label="control")
    plt.plot(binn, p_targ2, "k",  label="target")
    plt.title(f"p‑control (--), p‑target (k), indiv gRNA, {cond2}")
    plt.xlabel("LFC bin"); plt.ylabel("probability to be extreme"); plt.grid(True)

    # ------------------------------------------------------------------
    # 6. Return in the same order MATLAB provided
    # ------------------------------------------------------------------
    return (
        crit_LR12,
        me_sd12,
        med_mad12,
        binn,
        p_cont12,
        hiss_cont12,
        p_targ12,
        T_zGE_1,
        T_zGE_2,
    )


import numpy as np

def make_histo_vec_rep12(step: float,
                         vec_1: np.ndarray,
                         vec_2: np.ndarray,
                         st: float,
                         en: float):
    """
    Build histograms (counts & %s) for two replicate vectors and for the
    pooled vector.

    Parameters
    ----------
    step : float
        Bin width.
    vec_1, vec_2 : 1‑D arrays
        Data vectors (e.g. LFC values) for replicate‑1 and replicate‑2.
    st, en : float
        Lower / upper bounds for binning.

    Returns
    -------
    bin        : 1‑D array of bin left‑edges (len == #bins)
    perc_1     : % histogram for vec_1
    perc_2     : % histogram for vec_2
    perc_12    : % histogram for pooled vec_1 + vec_2
    his_1      : counts per bin for vec_1
    his_2      : counts per bin for vec_2
    his_12     : counts per bin for pooled data
    n1_n2_n12  : tuple  (len(vec_1), len(vec_2), len(vec_1)+len(vec_2))
    """

    # Histogram for replicate‑1
    bin_edges, his_1, perc_1 = make_histo_LFC(step, np.asarray(vec_1), st, en)

    # Histogram for replicate‑2
    _,        his_2, perc_2 = make_histo_LFC(step, np.asarray(vec_2), st, en)

    # Histogram for pooled replicates
    vec_12 = np.concatenate([vec_1, vec_2])
    _, his_12, perc_12 = make_histo_LFC(step, vec_12, st, en)

    n1_n2_n12 = (len(vec_1), len(vec_2), len(vec_12))

    return (
        bin_edges,
        perc_1,
        perc_2,
        perc_12,
        his_1,
        his_2,
        his_12,
        n1_n2_n12,
    )


import numpy as np
import matplotlib.pyplot as plt

def compute_p_critLFC(alf, binn, hiss_1, hiss_2, hiss_12, cond):
    """
    Compute p-values and critical values for LFC distributions of replicates and pooled.

    Parameters
    ----------
    alf : float
        Significance level.
    binn : array_like
        Bin edges or centers for LFC histogram.
    hiss_1, hiss_2, hiss_12 : array_like
        Histogram counts for replicate 1, replicate 2, and pooled data.
    cond : str
        Condition label for plot title.

    Returns
    -------
    bin_pzit : ndarray
        Array with columns: [bin, p_rep1, p_rep2, p_pooled]
    crit_12_targ : ndarray
        Critical values for rep1, rep2, and pooled (shape 3x2)
    """

    # Replicate 1
    p1, clz, crz, bin_pz, med_LFCpz, hiss4pz = q_val_frequentist_critical(alf, binn, hiss_1)
    critz = [clz, crz]

    # Replicate 2 (intergenic)
    p2, cli, cri, bin_pi, med_LFCpi, hiss4pi = q_val_frequentist_critical(alf, binn, hiss_2)
    criti = [cli, cri]

    # Both replicates pooled
    p12, clT, crT, bin_pT, med_LFCpT, hiss4pT = q_val_frequentist_critical(alf, binn, hiss_12)
    critT = [clT, crT]

    crit_12_targ = np.array([critz, criti, critT])

    # Assemble bin_pzit matrix with bin and p-values for rep1, rep2, pooled
    bin_pzit = np.column_stack((binn, p1, p2, p12))

    # Plotting
    plt.figure()
    plt.plot(binn, p1, 'b', label='rep1')
    plt.plot(binn, p12, 'k:', label='repBoth')
    plt.plot(binn, p2, 'r', label='rep2')
    plt.grid(True)
    plt.xlabel('LFC')
    plt.title(f'control rep1 rep2 repBoth(:), {cond}', fontsize=14)
    plt.legend()
    plt.show()

    return bin_pzit, crit_12_targ


import numpy as np
from scipy import stats

def med_mad_MZNP_2(GE_vec):
    """
    Compute median, MAD, modified Z-scores, mean, std, and mode for a vector.
    Implements Leys 2013 robust Z-score method.

    Parameters
    ----------
    GE_vec : array_like
        Input data vector.

    Returns
    -------
    med_mad : list
        [median, adjusted MAD]
    MZ2 : ndarray
        Modified Z-scores (normalized by adjusted MAD).
    MZ : ndarray
        Modified Z-scores (Leys version, scaled by MAD).
    Z : ndarray
        Standard Z-scores (mean/std normalization).
    me_sd : list
        [mean, standard deviation]
    mod : float
        Mode of the input vector.
    """

    GE_vec = np.array(GE_vec)
    med = np.median(GE_vec)
    
    # scipy.stats.mode returns mode and count; mode is an array
    mode_result = stats.mode(GE_vec, nan_policy='omit')
    mod = mode_result.mode.item()  # safely get the scalar value from numpy array
    
    b = 1.4826  # scale factor for MAD assuming normal distribution

    mad = np.median(np.abs(GE_vec - med))
    mad2 = b * mad

    med_mad = [med, mad2]

    me_sd = [np.mean(GE_vec), np.std(GE_vec, ddof=1)]

    # Modified Z-scores
    # Avoid division by zero if mad or mad2 is zero:
    if mad == 0:
        MZ = np.zeros_like(GE_vec)
    else:
        MZ = 0.6745 * (GE_vec - med) / mad

    if mad2 == 0:
        MZ2 = np.zeros_like(GE_vec)
    else:
        MZ2 = (GE_vec - med) / mad2

    # Standard Z-scores
    mean_val = me_sd[0]
    std_val = me_sd[1]
    if std_val == 0:
        Z = np.zeros_like(GE_vec)
    else:
        Z = (GE_vec - mean_val) / std_val

    return med_mad, MZ2, MZ, Z, me_sd, mod

import numpy as np
import pandas as pd

def implement_p_control_indiv_gRNA(T_gRNA_LFC, binn, p_contr1, p_contr2):
    """
    Implement p-controls by LFC targets (instead of p-targets) at gRNA level.
    
    Parameters:
    - T_gRNA_LFC: pd.DataFrame with columns ['gRNA', 'gene', 'LFCM', 'LFCF']
    - binn: array-like, bin edges
    - p_contr1: array-like, p/q values corresponding to binn for LFC_1
    - p_contr2: array-like, p/q values corresponding to binn for LFC_2
    
    Returns:
    - T_12_q_vertical: pd.DataFrame with columns ['gRNA', 'gene', 'lfc_1', 'lfc_2', 'Q1', 'Q2']
    """
    print('----------implement p-controls by LFC targets (instead of p-tagets)')

    LFC_1 = T_gRNA_LFC.iloc[:, 2].to_numpy()  # LFCM
    LFC_2 = T_gRNA_LFC.iloc[:, 3].to_numpy()  # LFCF

    qq1 = np.zeros(len(LFC_1))
    qq2 = np.zeros(len(LFC_2))

    for i in range(len(binn) - 1):
        bin_start = binn[i]
        bin_end = binn[i + 1]
        
        # For LFC_1
        idx_1 = (LFC_1 > bin_start) & (LFC_1 <= bin_end)
        qq1[idx_1] = p_contr1[i]

        # For LFC_2
        idx_2 = (LFC_2 > bin_start) & (LFC_2 <= bin_end)
        qq2[idx_2] = p_contr2[i]

    # Assemble new DataFrame
    T_12_q_vertical = pd.DataFrame({
        'gRNA': T_gRNA_LFC.iloc[:, 0],
        'gene': T_gRNA_LFC.iloc[:, 1],
        'lfc_1': LFC_1,
        'lfc_2': LFC_2,
        'Q1': qq1,
        'Q2': qq2,
    })

    return T_12_q_vertical


import pandas as pd

def module_p_control_target_implement(T_target, T_zGE, binn, p_cont12):
    """
    Given p per bin, infer p/q values for [T_target; T_zGE] from LFC values.

    Parameters:
    - T_target: pd.DataFrame with gRNA target data
    - T_zGE: pd.DataFrame with zGE gene data
    - binn: array-like, bin edges
    - p_cont12: 2D array-like, calibration p/q values per bin (columns for rep1 and rep2)

    Returns:
    - T_ij_q_verti: pd.DataFrame with gRNAs, genes, lfc_1, lfc_2, Q1, Q2
    """

    print('module------FDR  FRR (false rejection)-correct p-values for both tails of LFC NNT distribution')
    print('%-------chosen controls ')
    print('%-------------computed p-vals from Control distribution, frequentist')

    print('module ----------implement p-controls into LFC (instead of p-targets)')
    print('T for gRNAs (verti) with re-calibrated (fdr and frr) Qs, Target with zGE genes')

    # Concatenate tables vertically (like MATLAB's [T_target; T_zGE])
    T_target_zGE = pd.concat([T_target, T_zGE], ignore_index=True)

    # Extract p/q control columns
    p_contr1 = p_cont12[:, 0]
    p_contr2 = p_cont12[:, 1]

    # Call the function that assigns p/q values based on LFC bins
    T_ij_q_verti = implement_p_control_indiv_gRNA(T_target_zGE, binn, p_contr1, p_contr2)

    print(T_ij_q_verti.head(5))

    return T_ij_q_verti


import numpy as np
import matplotlib.pyplot as plt

def computeZ(st, en, step, LFC_t1, LFC_t2, me_sd_z12):
    """
    Compute Z scores for targets based on mu and sd from control (CTR).

    Inputs:
    - LFC_t1, LFC_t2: numpy arrays of LFC values for replicate 1 and 2
    - me_sd_z12: 2x2 array with mean and sd for replicate 1 and 2 (rows)
    - st, en, step: float, histogram bin parameters

    Outputs:
    - Z_t1, Z_t2: Z-normalized LFC vectors for replicate 1 and 2
    - binn: bin centers for histograms
    - perc_t1, perc_zt1, perc_t2, perc_zt2: percentage histograms for LFC and Z-LFC
    """
    mu1, sd1 = me_sd_z12[0, 0], me_sd_z12[0, 1]
    mu2, sd2 = me_sd_z12[1, 0], me_sd_z12[1, 1]

    binn, _, perc_t1 = make_histo_LFC(step, LFC_t1, st, en)
    Z_t1 = (LFC_t1 - mu1) / sd1
    _, _, perc_zt1 = make_histo_LFC(step, Z_t1, st, en)

    _, _, perc_t2 = make_histo_LFC(step, LFC_t2, st, en)
    Z_t2 = (LFC_t2 - mu2) / sd2
    _, _, perc_zt2 = make_histo_LFC(step, Z_t2, st, en)

    return Z_t1, Z_t2, binn, perc_t1, perc_zt1, perc_t2, perc_zt2

import pandas as pd

def make_tables_Z_two(target_12, Z_t1, Z_t2, cond1, cond2):
    """
    Takes a vertical gRNA gene table with LFC1 & LFC2, Q1 & Q2,
    adds Z1, Z2 columns, and makes two separate tables with Z for rep1 and rep2.
    
    Parameters:
    - target_12: pd.DataFrame with columns [gRNA, gene, lfc_1, lfc_2, Q1, Q2]
    - Z_t1, Z_t2: numpy arrays or lists of Z scores for replicate 1 and 2
    - cond1, cond2: string labels for conditions for rep1 and rep2
    
    Returns:
    - T_t1_wt, T_t2_wt: pd.DataFrames with columns [gRNA, gene, LFC, Z_zGE_LFC, Q, condition]
    """

    # Extract columns
    gRNA = target_12.iloc[:, 0]
    gene = target_12.iloc[:, 1]
    LFC_t1 = target_12.iloc[:, 2]
    LFC_t2 = target_12.iloc[:, 3]
    Q1 = target_12.iloc[:, 4]
    Q2 = target_12.iloc[:, 5]

    # Replicate 1 table
    condition1 = [cond1] * len(LFC_t1)
    T_t1_wt = pd.DataFrame({
        'gRNA': gRNA,
        'gene': gene,
        'LFC': LFC_t1,
        'Z_zGE_LFC': Z_t1,
        'Q': Q1,
        'condition': condition1
    })

    # Replicate 2 table
    condition2 = [cond2] * len(LFC_t2)
    T_t2_wt = pd.DataFrame({
        'gRNA': gRNA,
        'gene': gene,
        'LFC': LFC_t2,
        'Z_zGE_LFC': Z_t2,
        'Q': Q2,
        'condition': condition2
    })

    return T_t1_wt, T_t2_wt

import pandas as pd
import numpy as np

def z_p_CTR_any(
    st: float,
    en: float,
    step: float,
    T_any: pd.DataFrame,
    binn: np.ndarray,
    p_cont12: np.ndarray,
    me_sd_12: np.ndarray,
    cond1: str,
    cond2: str,
):
    """
    Calibrate any control set (e.g. intergenic or NT) against zGE controls,
    assign p/q values, compute Z‑scores, and return per‑replicate tables + histos.
    """
    T_ij_q_verti_any = module_p_control_any_implement(T_any, binn, p_cont12)

    LFC_t1 = T_ij_q_verti_any.iloc[:, 2].to_numpy()
    LFC_t2 = T_ij_q_verti_any.iloc[:, 3].to_numpy()

    Z_t1, Z_t2, binn, perc_t1, perc_zt1, perc_t2, perc_zt2 = computeZ(
        st, en, step, LFC_t1, LFC_t2, me_sd_12
    )

    T_z_q_any1, T_z_q_any2 = make_tables_Z_two(
        T_ij_q_verti_any, Z_t1, Z_t2, cond1, cond2
    )

    return T_z_q_any1, T_z_q_any2, binn, perc_t1, perc_zt1, perc_t2, perc_zt2

def plot_histograms(binn, perc_t1, perc_zt1, perc_t2, perc_zt2, condition_label):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))

    plt.subplot(2, 2, 1)
    plt.bar(binn, perc_t1, color='k')
    plt.grid(True)
    plt.title('LFC distri rep1')
    plt.xlim([-16, 8])
    plt.xlabel('LFC')

    plt.subplot(2, 2, 3)
    plt.bar(binn, perc_zt1, color='g')
    plt.grid(True)
    plt.title('Z LFC distri rep1')
    plt.xlim([-16, 8])
    plt.xlabel('Z LFC')

    plt.subplot(2, 2, 2)
    plt.bar(binn, perc_t2, color='k')
    plt.grid(True)
    plt.title('LFC distri rep2')
    plt.xlim([-16, 8])
    plt.xlabel('LFC')

    plt.subplot(2, 2, 4)
    plt.bar(binn, perc_zt2, color='g')
    plt.grid(True)
    plt.title('Z LFC distri rep2')
    plt.xlim([-16, 8])
    plt.xlabel('Z LFC')

    plt.suptitle(f'LFC gRNA distri: {condition_label}', fontsize=14)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

import numpy as np
import pandas as pd
from typing import Union

def module_p_control_target_implement(
    T_target: pd.DataFrame,
    T_zGE: pd.DataFrame,
    binn: np.ndarray,
    p_cont12: Union[np.ndarray, list],
):
    """
    Calibrate p/Q‑values for target + zGE gRNAs using control p‑curves.

    Parameters
    ----------
    T_target : pd.DataFrame
        Target table (columns: gRNA, gene, lfc_1, lfc_2, …).
    T_zGE    : pd.DataFrame
        zGE control table (same column structure).
    binn     : 1‑D array
        Bin edges (length N) used to map LFC→Q.
    p_cont12 : (#bins, 2) np.ndarray or list‑of‑lists
        Calibration p/Q per bin: col‑0 = replicate 1, col‑1 = replicate 2.

    Returns
    -------
    T_ij_q_verti : pd.DataFrame
        Combined table with new Q1/Q2 columns appended.
    """

    print(
        "module —— FDR / FRR recalibration: applying control p‑values "
        "to both tails of the LFC distribution"
    )
    print("%------- chosen controls / calibration curves loaded")

    # ------------------------------------------------------------------
    # 1.  Merge target + zGE vertically (row‑wise)
    # ------------------------------------------------------------------
    T_target_zGE = pd.concat([T_target, T_zGE], ignore_index=True)

    # ------------------------------------------------------------------
    # 2.  Split calibration matrix into two replicate vectors
    # ------------------------------------------------------------------
    p_cont12 = np.asarray(p_cont12)
    p_contr1 = p_cont12[:, 0]    # replicate 1
    p_contr2 = p_cont12[:, 1]    # replicate 2

    # ------------------------------------------------------------------
    # 3.  Delegate p‑assignment to helper
    # ------------------------------------------------------------------
    T_ij_q_verti = implement_p_control_indiv_gRNA(
        T_target_zGE, binn, p_contr1, p_contr2
    )

    # Preview first 5 rows (mirrors MATLAB's T_ij_q_verti_head)
    print("\nHead of recalibrated table (first 5 rows):")
    print(T_ij_q_verti.head(5))

    return T_ij_q_verti


import pandas as pd

def module_p_control_target_implement(T_target: pd.DataFrame, T_zGE: pd.DataFrame, binn, p_cont12):
    """
    Given p-values per bin, infer recalibrated Q-values for combined T_target and T_zGE based on LFC values.

    Parameters:
    - T_target: pd.DataFrame with columns ['gRNA', 'gene', 'lfc_1', 'lfc_2']
    - T_zGE: pd.DataFrame with same structure as T_target
    - binn: bin edges or categories used for calibration
    - p_cont12: 2D array or DataFrame, calibration p-values per bin (shape Nx2)

    Returns:
    - T_ij_q_verti: pd.DataFrame with columns
      ['gRNA', 'gene', 'lfc_1', 'lfc_2', 'Q1', 'Q2']
    """

    print("module------FDR  FRR (false rejection)-correct p-values for both tails of LFC NNT distribution")
    print("%-------chosen controls")
    print("%-------------computed p-vals from Control distribution, frequentist")

    print("module ----------implement p-controls into LFC (instead of p-targets)")
    print("T for gRNAs (verti) with re-calibrated (fdr and frr) Qs, Target with zGE genes")

    # Concatenate target and zGE DataFrames vertically
    T_target_zGE = pd.concat([T_target, T_zGE], ignore_index=True)

    # Extract calibration p-values for each replicate
    p_contr1 = p_cont12[:, 0]  # first column
    p_contr2 = p_cont12[:, 1]  # second column

    # Call the helper function to implement p-control corrections
    T_ij_q_verti = implement_p_control_indiv_gRNA(T_target_zGE, binn, p_contr1, p_contr2)

    print("Top 5 rows of recalibrated table:")
    print(T_ij_q_verti.head(5))

    return T_ij_q_verti


import pandas as pd

def module_p_control_any_implement(T_any: pd.DataFrame, binn, p_cont12):
    """
    Given p-values per bin, infer recalibrated Q-values for T_any based on its LFC values.

    Parameters:
    - T_any: pd.DataFrame with columns ['gRNA', 'gene', 'lfc_1', 'lfc_2']
    - binn: bin edges or categories used for calibration
    - p_cont12: 2D array or DataFrame, calibration p-values per bin (shape Nx2)

    Returns:
    - T_ij_q_verti: pd.DataFrame with columns
      ['gRNA', 'gene', 'lfc_1', 'lfc_2', 'Q1', 'Q2']
    """

    print("module------FDR  FRR (false rejection)-correct p-values for both tails of LFC NNT distribution")
    print("%-------chosen controls")
    print("%-------------computed p-vals from Control distribution, frequentist")

    print("module ----------implement p-controls into LFC (instead of p-targets)")
    print("T for gRNAs (verti) with re-calibrated (fdr and frr) Qs, Any set with lfc1 lfc2")

    # Extract calibration p-values for each replicate
    p_contr1 = p_cont12[:, 0]  # first column
    p_contr2 = p_cont12[:, 1]  # second column

    # Call the helper function to implement p-control corrections
    T_ij_q_verti = implement_p_control_indiv_gRNA(T_any, binn, p_contr1, p_contr2)

    print("Top 5 rows of recalibrated table:")
    print(T_ij_q_verti.head(5))

    return T_ij_q_verti


import matplotlib.pyplot as plt

def volcano_gRNA_gene_hits_wt(alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond):
    genes = T_vert['gene']
    scoreZ = T_vert['Z_zGE_LFC']
    fdr = T_vert['Q']

    # Prepare figure with 3 subplots side by side
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))

    # 1. Volcano per gRNA (Z-normalized)
    LPV, indha, indda, T_ds, T_hs, T_gRNA = general_volcano(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz, scoreZ, fdr, cond, genes, plot=True, ax=axs[0]
    )
    axs[0].set_xlabel('LFC standardized', fontsize=14)
    axs[0].set_title(f'Volcano per gRNA: {cond}', fontsize=16)

    # 2 & 3. Per gene (LFC and Z-normalized LFC)
    (
        T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z,
        indha_gene, indda_gene, indhaz, inddaz
    ) = perGene_4_hits_med_horiz(
        alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond, plot=False
    )

    if not T_lfc_z_q_med.empty:
        fdr_gene = T_lfc_z_q_med.iloc[:, 3]
        genes_gene = T_lfc_z_q_med.iloc[:, 0]

        # Per gene LFC volcano
        general_volcano(
            alf, sfdr_corr, thr_lfch, thr_lfcd,
            T_lfc_z_q_med.iloc[:, 1].values, fdr_gene.values, cond, genes_gene.values,
            plot=True, ax=axs[1]
        )
        axs[1].set_xlabel('LFC', fontsize=14)
        axs[1].set_title(f'Volcano per gene LFC: {cond}', fontsize=16)

        # Per gene Z volcano
        general_volcano(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med.iloc[:, 2].values, fdr_gene.values, cond, genes_gene.values,
            plot=True, ax=axs[2]
        )
        axs[2].set_xlabel('Z-normalised LFC', fontsize=14)
        axs[2].set_title(f'Volcano per gene Z: {cond}', fontsize=16)
    else:
        # If no gene data, hide those axes
        axs[1].axis('off')
        axs[2].axis('off')

    plt.tight_layout()
    plt.show()

    return T_gRNA, T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z, indha, indda, indhaz, inddaz


import plotly.graph_objects as go
import numpy as np
import pandas as pd

def general_volcano_interactive(alf, sfdr_corr, thr_scoreh, thr_scored, score, fdr, cond, genes, plot=True, ax=None):
    thr_fdr = alf
    fdr_corr = fdr + sfdr_corr  # Avoid log10(0)
    LPV = -np.log10(fdr_corr)
    
    T_gene = pd.DataFrame({'gene': genes, 'score': score, 'fdr_corr': fdr_corr})
    
    # Hits and depleted indices
    indh = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s >= thr_scoreh)]
    indd = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s <= thr_scored)]

    # Default color for all points
    colors = ['black'] * len(genes)
    for i in indh:
        colors[i] = 'red'  # hits
    for i in indd:
        colors[i] = 'cyan'  # depleted

    fig = go.Figure()

    # All genes scatter plot
    fig.add_trace(go.Scatter(
        x=score,
        y=LPV,
        mode='markers',
        marker=dict(color=colors, size=7, line=dict(width=0.5, color='DarkSlateGrey')),
        text=genes,  # tooltip text
        hovertemplate='<b>%{text}</b><br>Score: %{x:.2f}<br>-Log10 FDR: %{y:.2f}<extra></extra>'
    ))

    fig.update_layout(
        title=f'Volcano plot: {cond}',
        xaxis_title='Score',
        yaxis_title='-Log10 FDR',
        template='plotly_white',
        xaxis=dict(range=[min(score)-1, max(score)+1]),
        yaxis=dict(range=[0, max(LPV) + 1]),
        hovermode='closest',
        width=900,
        height=600,
    )
    
    fig.show()

def general_volcano(alf, sfdr_corr, thr_scoreh, thr_scored, score, fdr, cond, genes, plot=True, ax=None):
    thr_fdr = alf
    fdr_corr = fdr + sfdr_corr  # Avoid log10(0)
    LPV = -np.log10(fdr_corr)

    T_gene = pd.DataFrame({'gene': genes, 'score': score, 'fdr_corr': fdr_corr})

    indh = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s >= thr_scoreh)]
    indd = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s <= thr_scored)]

    T_hs = T_gene.loc[indh].copy() if indh else pd.DataFrame()
    T_ds = T_gene.loc[indd].copy() if indd else pd.DataFrame()

    if plot and ax is not None:
        ma = max(5, np.max(LPV) + 1 if LPV.size > 0 else 5)
        ax.plot(score, LPV, 'pk', markersize=4, label='All genes')
        ax.grid(True)
        ax.set_ylabel('-Log10 FDR', fontsize=14)
        ax.set_xlabel('score', fontsize=14)
        ax.set_xlim([-20, 10])
        ax.set_ylim([0, ma])

        if not T_hs.empty:
            ax.plot(T_hs['score'], -np.log10(T_hs['fdr_corr']), 'pr', linewidth=2, label='Hits')
            for i, gene_name in enumerate(T_hs['gene'][:20]):
                ax.text(T_hs['score'].iloc[i], -np.log10(T_hs['fdr_corr']).iloc[i], gene_name, va='bottom', ha='right')

        if not T_ds.empty:
            ax.plot(T_ds['score'], -np.log10(T_ds['fdr_corr']), 'pc', linewidth=2, label='Depleted')
            for i, gene_name in enumerate(T_ds['gene'][:20]):
                ax.text(T_ds['score'].iloc[i], -np.log10(T_ds['fdr_corr']).iloc[i], gene_name, va='bottom', ha='right')

        ax.legend()

    return LPV, indh, indd, T_ds, T_hs, T_gene


def perGene_4_hits_med_horiz(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond, plot=False):
    """
    Analyze gene hits using median of 4-gRNA targets. Generates volcano plots
    only when plot=True.

    Returns:
    - T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z
    - indha, indda, indhaz, inddaz: hit/depleted indices for LFC and Z
    """

    print('check how many four three gRNAs we have')

    gene_names = T_vert['gene'].values

    ggenes, gg, ind_gn = count_vertical_names(gene_names)
    ind_num = np.column_stack((ind_gn, gg))

    print('Unique gene names count:', len(set(T_vert['gene'])))
    print('Gene counts sample:')
    gene_counts = T_vert['gene'].value_counts()
    print(gene_counts.head(10))
    print("gRNA count per gene distribution:")
    print(gene_counts.value_counts().sort_index())
    (Genes_2, Genes_3, Genes_4, gn_2, gn_3, gn_4,
     wt_2, wt_3, wt_4, d4, nums4) = separate_fours_threes_twos_genes(
        T_vert, ggenes, gene_names, ind_gn, gg
    )
    print('Genes_4 length:', len(Genes_4))
    print('Sample Genes_4:', Genes_4[:10])
    sum_nums = np.sum(nums4)

    if sum_nums > 0:
        T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z = perGene_4_med_horiz(Genes_4, wt_4)

        fdr = T_lfc_z_q_med.iloc[:, 3]
        genes = T_lfc_z_q_med.iloc[:, 0]
        score = T_lfc_z_q_med.iloc[:, 1]  
        # median_LFC
        print('sum_nums:', sum_nums)
        print('T_lfc_z_q_med shape:', T_lfc_z_q_med.shape)
        print('T_lfc_z_q_med head:\n', T_lfc_z_q_med.head())
        print('Genes used for plotting:', genes.values if not T_lfc_z_q_med.empty else 'None')
        print('Score values:', score.values if not T_lfc_z_q_med.empty else 'None')
        print('FDR values:', fdr.values if not T_lfc_z_q_med.empty else 'None')
        # Just return the data; do NOT plot here

        _, indha, indda, _, _, _  = general_volcano(
            alf, sfdr_corr, thr_lfch, thr_lfcd,
            T_lfc_z_q_med.iloc[:, 1].values, fdr.values, cond, genes.values, plot=False
        )
        _, indhaz, inddaz, _, _, _ = general_volcano(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med.iloc[:, 2].values, fdr.values, cond, genes.values, plot=False
        )
    else:
        # empty returns as before
        indha = indda = indhaz = inddaz = []

    return T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z, indha, indda, indhaz, inddaz


def count_vertical_names(gene_names):
    """
    Given a list of gene names (strings), this function groups consecutive identical names,
    returning:
    - ggenes: list of unique gene names for each consecutive group
    - gg: list of counts of repeats per unique gene name group
    - ind_gn: list of start indices (1-based) for each group in the original list
    
    Example:
    gene_names = ['Aars2', 'Aars2', 'Aars2', 'Aars2', 'Aasdhppt', 'Aasdhppt', 'Aasdhppt', 'Aasdhppt', ...]
    ggenes = ['Aars2', 'Aasdhppt', ...]
    gg = [4, 4, ...]  # counts
    ind_gn = [1, 5, ...]  # start positions (1-based indexing for MATLAB compatibility)
    """
    ggenes = []
    gg = []
    ind_gn = []

    j = 0  # group index in Python (0-based)
    k = 1  # count for current group
    
    ggenes.append(gene_names[0])
    ind_gn.append(1)  # MATLAB style 1-based indexing
    
    for i in range(len(gene_names) - 1):
        s1 = gene_names[i]
        s2 = gene_names[i + 1]
        
        if s1 == s2:
            k += 1
        else:
            gg.append(k)
            # add the unique gene for the completed group
            ggenes[j] = gene_names[i]  # update last added to be explicit (optional)
            j += 1
            k = 1
            ggenes.append(gene_names[i + 1])  # start new group
            ind_gn.append(i + 2)  # next group start index (1-based)
    
    # Handle the last group
    # After the loop ends, 'k' is the count for the last group
    gg.append(k)
    ggenes[j] = gene_names[-1]

    return ggenes, gg, ind_gn


def separate_fours_threes_twos_genes(T_vert, ggenes, gene_names, ind_gn, gg):
    """
    Separates genes by the number of gRNAs per gene (1, 2, 3, or 4).
    
    Inputs:
        T_vert: pandas DataFrame with columns ['gRNA_targ', 'gene_targ', 'LFC', 'Z_zGE_LFC', 'Q', 'condition']
        ggenes: list of unique gene names (compressed)
        gene_names: list of all gene names, sorted/grouped
        ind_gn: list of start indices (1-based) for groups in gene_names
        gg: list of counts of gRNA per gene
    
    Outputs:
        Genes_2, Genes_3, Genes_4: lists of gene names with 2, 3, or 4 gRNAs respectively
        gn_2, gn_3, gn_4: lists of gene names repeated per gRNA counts (detailed)
        wt_2, wt_3, wt_4: numpy arrays of [LFC, Z, Q] values for each group
        d: max number of gRNAs found per gene if all equal, otherwise variability
        nums: counts of how many genes have 4, 3, 2, or 1 gRNA(s)
    """
    
    import numpy as np
    
    k4 = k3 = k2 = k1 = 0
    gn_3, wt_3 = [], []
    gn_4, wt_4 = [], []
    gn_2, wt_2 = [], []
    Genes_2, Genes_3, Genes_4, Genes_1 = [], [], [], []
    ind1 = []

    # Extract numeric columns: LFC, Z, Q as numpy array for slicing
    wt = T_vert.iloc[:, 2:5].to_numpy()
    
    for i in range(len(gg)):
        start_idx = ind_gn[i] - 1  # Convert MATLAB 1-based to Python 0-based indexing
        
        if gg[i] == 3:
            k3 += 1
            Genes_3.append(ggenes[i])
            gn = gene_names[start_idx:start_idx+3]
            gn_3.extend(gn)
            wtt = wt[start_idx:start_idx+3, :]
            wt_3.append(wtt)
            
        elif gg[i] == 4:
            k4 += 1
            Genes_4.append(ggenes[i])
            gn = gene_names[start_idx:start_idx+4]
            gn_4.extend(gn)
            wtt = wt[start_idx:start_idx+4, :]
            wt_4.append(wtt)
            
        elif gg[i] == 2:
            k2 += 1
            Genes_2.append(ggenes[i])
            gn = gene_names[start_idx:start_idx+2]
            gn_2.extend(gn)
            wtt = wt[start_idx:start_idx+2, :]
            wt_2.append(wtt)
            
        elif gg[i] == 1:
            k1 += 1
            Genes_1.append(ggenes[i])
            ind1.append(i)
    
    # Convert wt lists to numpy arrays (if not empty)
    wt_2 = np.vstack(wt_2) if wt_2 else np.empty((0,3))
    wt_3 = np.vstack(wt_3) if wt_3 else np.empty((0,3))
    wt_4 = np.vstack(wt_4) if wt_4 else np.empty((0,3))
    
    max_num_gRNA = max(gg)
    min_num_gRNA = min(gg)
    
    if max_num_gRNA == min_num_gRNA:
        print("All genes have the same number of gRNAs")
        d = min_num_gRNA
    else:
        print("Variability of gRNAs per gene: make decision what to do")
        d = max_num_gRNA
        
    nums = [k4, k3, k2, k1]
    
    return Genes_2, Genes_3, Genes_4, gn_2, gn_3, gn_4, wt_2, wt_3, wt_4, d, nums

import numpy as np
import pandas as pd

def perGene_4_med_horiz(Genes_4, wt_4):
    """
    Processes sets of 4 gRNAs per gene, calculates statistics (mean, median, std)
    for LFC, Z, and Q values and returns pandas DataFrames similar to MATLAB tables.
    
    Inputs:
        Genes_4: list of gene names (length = number of genes)
        wt_4: numpy array or 2D list of shape (4 * number_of_genes, 3) 
              columns correspond to [LFC, Z, Q] per gRNA, stacked vertically
    
    Outputs:
        T_lfc_z_q_med: DataFrame with genes and median statistics for LFC, Z, Q
        T_lfc_z_q_me: DataFrame with genes and mean statistics for LFC, Z, Q
        T_LFC: DataFrame with LFC values per gRNA and summary stats
        T_Q: DataFrame with Q values per gRNA and summary stats
        T_Z: DataFrame with Z values per gRNA and summary stats
    """
    
    d = 4  # number of gRNAs per gene
    n_genes = len(Genes_4)
    
    # Initialize arrays to hold the values
    LFC = np.zeros((n_genes, d))
    Z = np.zeros((n_genes, d))
    Q = np.zeros((n_genes, d))
    
    # Compute statistics for each gene
    mean_LFC = np.zeros(n_genes)
    median_LFC = np.zeros(n_genes)
    std_LFC = np.zeros(n_genes)
    
    mean_Z = np.zeros(n_genes)
    median_Z = np.zeros(n_genes)
    std_Z = np.zeros(n_genes)
    
    mean_Q = np.zeros(n_genes)
    median_Q = np.zeros(n_genes)
    std_Q = np.zeros(n_genes)
    
    for n in range(n_genes):
        n1 = n * d
        n2 = n1 + d
        
        wt_slice = wt_4[n1:n2, :]  # shape (4,3) [LFC,Z,Q]
        
        lfc = wt_slice[:, 0]
        z = wt_slice[:, 1]
        q = wt_slice[:, 2]
        
        LFC[n, :] = lfc
        Z[n, :] = z
        Q[n, :] = q
        
        mean_LFC[n] = np.mean(lfc)
        median_LFC[n] = np.median(lfc)
        std_LFC[n] = np.std(lfc, ddof=1)  # sample std dev
        
        mean_Z[n] = np.mean(z)
        median_Z[n] = np.median(z)
        std_Z[n] = np.std(z, ddof=1)
        
        mean_Q[n] = np.mean(q)
        median_Q[n] = np.median(q)
        std_Q[n] = np.std(q, ddof=1)
    
    # Add a tiny constant to Q means/medians to avoid zero
    insteadZ = 1e-6
    mean_Q_adj = mean_Q + insteadZ
    median_Q_adj = median_Q + insteadZ
    
    # Create pandas DataFrames resembling MATLAB tables
    
    # DataFrames for LFC, Z, Q with per gRNA and summary stats
    T_LFC = pd.DataFrame({
        'genes': Genes_4,
        'lfc_1': LFC[:, 0],
        'lfc_2': LFC[:, 1],
        'lfc_3': LFC[:, 2],
        'lfc_4': LFC[:, 3],
        'mean_LFC': mean_LFC,
        'median_LFC': median_LFC,
        'std_LFC': std_LFC,
        'mean_q': mean_Q_adj,
        'median_q': median_Q_adj
    })
    
    T_Z = pd.DataFrame({
        'genes': Genes_4,
        'z_1': Z[:, 0],
        'z_2': Z[:, 1],
        'z_3': Z[:, 2],
        'z_4': Z[:, 3],
        'mean_Z': mean_Z,
        'median_Z': median_Z,
        'std_Z': std_Z,
        'mean_q': mean_Q_adj,
        'median_q': median_Q_adj
    })
    
    T_Q = pd.DataFrame({
        'genes': Genes_4,
        'q_1': Q[:, 0],
        'q_2': Q[:, 1],
        'q_3': Q[:, 2],
        'q_4': Q[:, 3],
        'mean_q': mean_Q_adj,
        'median_q': median_Q_adj,
        'std_q': std_Q
    })
    
    T_lfc_z_q_med = pd.DataFrame({
        'genes': Genes_4,
        'median_LFC': median_LFC,
        'median_Z': median_Z,
        'median_q': median_Q_adj
    })
    
    T_lfc_z_q_me = pd.DataFrame({
        'genes': Genes_4,
        'mean_LFC': mean_LFC,
        'mean_Z': mean_Z,
        'mean_q': mean_Q_adj
    })
    
    return T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z
