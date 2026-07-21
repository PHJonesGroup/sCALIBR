from .filter_pattern_distri import filter_pattern_distri
from .zGE_target_distri import zGE_target_distri
from .compute_hiss_LFC_rep12 import compute_hiss_LFC_rep12
from .distri_target_contr_plots_all import distri_target_contr_plots_all
import numpy as np
import pandas as pd 
def separate_target_control(st, en, step, T_norm_indiv, zGE, pat1, pat2, cond1, cond2, rep_pairs, output_dir):
    """
    Parameters:
    - st, en, step: for histogram binning
    - T_norm_indiv: full normalized gRNA table
    - zGE: list (table-like) of zGE genes (first column = gene names)
    - pat1: pattern for intergenic control (e.g., 'chr')
    - pat2: pattern for NT control (e.g., 'Non')

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
    n = len(pat)

    (
        indiv_Chr, indiv_noChr,
        ind_out1, ind_in1, num_out_in1,
        bin, hisFMi, percFMi,
        LFC_i,
        s_chr, st1, en1, T_lfc_chr
    ) = filter_pattern_distri(T_norm_indiv, pat, n, st, en, step, cond1, cond2, rep_pairs)
    
    num_chr_rest = num_out_in1

    # -----------------------
    # 2. Separate NT (e.g., 'Non')
    pat = pat2
    n = len(pat)
    (
        indiv_NT, indiv_noChr_noNT,
        ind_out2, ind_in2, num_out_in2,
        bin, hisFMnt, percFMnt,
        LFC_nt,
        s_nt, st1, en1, T_lfc_nt
    ) = filter_pattern_distri(indiv_noChr, pat, n, st, en, step, cond1, cond2, rep_pairs)

    num_NT_rest = num_out_in2

    # -----------------------
    # 3. Separate zGE and target genes
    zGE_genes = zGE.iloc[:, 0].tolist() if zGE is not None and len(zGE) > 0 else []
    genes_nnt = indiv_noChr_noNT.iloc[:, 1].tolist()
    has_zGE = len(zGE_genes) > 0

    if has_zGE:
        (
            T_zGE, T_nzGE, bin,
            his1z, perc1z, his2z, perc2z, hisFMz, percFMz,
            his1t, perc1t, his2t, perc2t, hisFMt, percFMt,
            gzn
        ) = zGE_target_distri(genes_nnt, zGE_genes, indiv_noChr_noNT, st, en, step)

        # collapse two-series -> single series (single averaged LFC)
        hisz, percz = his1z, perc1z
        hist, perct = his1t, perc1t
        # hisFMz, percFMz, hisFMt, percFMt already defined above
    else:
        T_zGE = pd.DataFrame(columns=["gRNA", "gene", "lfc"])
        T_nzGE = indiv_noChr_noNT
        gzn = None
        hisz = percz = hisFMz = percFMz = None

        genes = indiv_noChr_noNT.iloc[:, 1].astype(str)
        gRNA  = indiv_noChr_noNT.iloc[:, 0]

        eps = 1e-6
        # pair cond1_<rep> with cond2_<rep> by shared replicate id
        T_tar = pd.DataFrame({'gRNA': gRNA.values, 'gene': genes.values})
        for r in rep_pairs:                       # e.g. ['60159','60160','60161','60162']
            c1 = f"{cond1}_{r}"
            c2 = f"{cond2}_{r}"
            if c1 not in indiv_noChr_noNT.columns or c2 not in indiv_noChr_noNT.columns:
                raise ValueError(f"Missing {c1} or {c2}")
            T_tar[f"lfc_{r}"] = np.log2(
                (indiv_noChr_noNT[c1].astype(float) + eps) /
                (indiv_noChr_noNT[c2].astype(float) + eps)
            ).values
        bin, hist, perct, LFC_t, lfc_cols_t, st1t, en1t = compute_hiss_LFC_rep12(T_tar, st, en, step)
        hisFMt = percFMt = None

    # -----------------------
    # 4. Compute distributions for chr and NT controls

    bin, hisi, perci, LFC_i, lfc_cols_i, st1i, en1i = compute_hiss_LFC_rep12(T_lfc_chr, st, en, step)
    bin, hisn, percn, LFC_n, lfc_cols_n, st1n, en1n = compute_hiss_LFC_rep12(T_lfc_nt, st, en, step)
    # -----------------------
    # 5. Plot histograms by condition
    distri_target_contr_plots_all(bin, percz, perci, perct, percn,
                                  f"{cond1}_vs_{cond2}", output_dir)

    # -----------------------
    # 6. Final Outputs
    T_target = T_nzGE                 # targets only
    T_target_zGE = indiv_noChr_noNT   # target + zGE
    if T_zGE is not None and len(T_zGE) > 0:
        T_zGE = T_zGE.sort_values(by="gene")

    return (
        T_tar, T_target_zGE, T_lfc_chr, T_lfc_nt, T_zGE,
        bin, hisz, percz, hisFMz, percFMz,
        hist, perct, hisFMt, percFMt, gzn
    )