import numpy as np
import pandas as pd

def general_volcano(alf, sfdr_corr, thr_scoreh, thr_scored, score, fdr, cond, genes, ax=None, grna=None):
    """
    Draw a single volcano panel and identify enriched / depleted hits.

    Parameters
    ----------
    alf : float
        FDR significance threshold for hit calling.
    sfdr_corr : float
        Pseudocount added to FDR before -log10 (avoids log10(0) on the y-axis).
    thr_scoreh : float
        Upper score threshold; points at or above it (and significant) are enriched hits.
    thr_scored : float
        Lower score threshold; points at or below it (and significant) are depleted hits.
    score : array-like
        Per-point score on the x-axis (LFC or Z-corrected LFC, depending on caller).
    fdr : array-like
        Per-point FDR / q-values, aligned to ``score`` and ``genes``.
    cond : str
        Condition label (for context; not required for the computation).
    genes : array-like
        Gene names, aligned to ``score`` and ``fdr``, used to label hit points.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on. If None, the function computes hits without plotting.
    grna : array-like
         gRNA names, aligned to ``score`` and ``fdr``, used to label hit points.
    
    Returns
    -------
    LPV : numpy.ndarray
        -log10 of the pseudocount-adjusted FDR, per point.
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
    thr_fdr = alf
    fdr_corr = fdr + sfdr_corr  # Avoid log10(0)
    LPV = -np.log10(fdr_corr)

    if grna is not None:
        T_gene = pd.DataFrame({
            'gRNA': np.asarray(grna),
            'gene': np.asarray(genes),
            'score': np.asarray(score),
            'fdr_corr': np.asarray(fdr_corr),
        })
    else:
        T_gene = pd.DataFrame({
            'gene': np.asarray(genes),
            'score': np.asarray(score),
            'fdr_corr': np.asarray(fdr_corr),
        })

    indh = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s >= thr_scoreh)]
    indd = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s <= thr_scored)]

    T_hs = T_gene.iloc[indh].copy() if indh else pd.DataFrame()
    T_ds = T_gene.iloc[indd].copy() if indd else pd.DataFrame()
    
    if ax is not None:
        ma = max(5, np.max(LPV) + 1 if LPV.size > 0 else 5)
        ax.plot(score, LPV, 'pk', markersize=4, label='Not significant')
        ax.grid(True)
        ax.set_ylabel('-log10 FDR', fontsize=14)
        ax.set_ylim([0, ma])

        if not T_hs.empty:
            ax.plot(T_hs['score'], -np.log10(T_hs['fdr_corr']), 'pr', linewidth=2, label='Hits')
            for i in range(min(20, len(T_hs))):
                ax.text(T_hs['score'].iloc[i], -np.log10(T_hs['fdr_corr']).iloc[i], T_hs['gene'].iloc[i], va='bottom', ha='right')

        if not T_ds.empty:
            ax.plot(T_ds['score'], -np.log10(T_ds['fdr_corr']), 'pc', label='Depleted')
            for i in range(min(20, len(T_ds))):
                ax.text(T_ds['score'].iloc[i], -np.log10(T_ds['fdr_corr']).iloc[i], T_ds['gene'].iloc[i], va='bottom', ha='right')

        ax.legend()

    return LPV, indh, indd, T_ds, T_hs, T_gene
