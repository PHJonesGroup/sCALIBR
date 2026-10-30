import numpy as np
import pandas as pd

def plot_volcano_interactive(significance_threshold, sfdr_corr, thr_scoreh, thr_scored, score, fdr, genes, grna=None):
    """
    Draw a single interactive volcano panel and identify enriched / depleted hits.

    Parameters
    ----------
    significance_threshold : float
        Hit threshold on the y-axis values (raw p for gRNAs, BH q for genes).
    sfdr_corr : float
        Pseudocount added to qvalue before -log10 (avoids log10(0) on the y-axis).
    thr_scoreh : float
        Upper score threshold; points at or above it (and significant) are enriched hits.
    thr_scored : float
        Lower score threshold; points at or below it (and significant) are depleted hits.
    score : array-like
        Per-point score on the x-axis (LFC or Z-corrected LFC, depending on caller).
    fdr : array-like
        Per-poinyt q-values, aligned to ``score`` and ``genes``.
    genes : array-like
        Gene names, aligned to ``score`` and ``fdr``, used to label hit points.
    
    Returns
    -------
    LPV : numpy.ndarray
        -log10 of the pseudocount-adjusted q-value, per point.
    indh : list of int
        Indices of enriched hits (significant, score >= thr_scoreh).
    indd : list of int
        Indices of depleted hits (significant, score <= thr_scored).
    T_ds : pandas.DataFrame
        Table of depleted-hit rows [gene, score, fdr_corr] (empty if none).
    T_hs : pandas.DataFrame
        Table of enriched-hit rows [gene, score, fdr_corr] (empty if none).
    T_gene : pandas.DataFrame
        All points as [gene, score, fdr_corr].
    """
    thr_fdr = significance_threshold
    fdr_corr = fdr + sfdr_corr  # avoid log10(0)
    LPV = -np.log10(fdr_corr)

    if grna is not None:
        T_gene = pd.DataFrame({
            'gRNA': np.asarray(grna),
            'gene': np.asarray(genes),        
            'score': np.asarray(score),
            'fdr_corr': np.asarray(fdr_corr),
            'LPV': np.asarray(LPV),            
        })
    else:
        T_gene = pd.DataFrame({
            'gene': np.asarray(genes),
            'score': np.asarray(score),
            'fdr_corr': np.asarray(fdr_corr),
            'LPV': np.asarray(LPV),           
        })

    indh = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s >= thr_scoreh)]
    indd = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s <= thr_scored)]

    T_hs = T_gene.loc[indh].copy() if indh else pd.DataFrame()
    T_ds = T_gene.loc[indd].copy() if indd else pd.DataFrame()

    return LPV, indh, indd, T_ds, T_hs, T_gene
