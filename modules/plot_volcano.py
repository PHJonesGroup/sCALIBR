import numpy as np
import pandas as pd

def plot_volcano(significance_threshold, thr_scoreh, thr_scored, score, fdr, genes, ax=None, grna=None):
    """
    Draw a single volcano panel and identify enriched / depleted hits.

    Parameters
    ----------
    significance_threshold : float
        Hit threshold on the y-axis values (raw p for gRNAs, BH q for genes).
    thr_scoreh : float
        Upper score threshold; points at or above it (and significant) are enriched hits.
    thr_scored : float
        Lower score threshold; points at or below it (and significant) are depleted hits.
    score : array-like
        Per-point score on the x-axis (LFC or Z-corrected LFC, depending on caller).
    fdr : array-like
        Per-point p- or q-values, aligned to ``score`` and ``genes``; never 0
        (empirical p >= 1 / (n_controls + 1)), so -log10 needs no pseudocount.
    genes : array-like
        Gene names, aligned to ``score`` and ``fdr``, used to label hit points.
    ax : matplotlib.axes.Axes, optional
        Axis to draw on. If None, the function computes hits without plotting.
    grna : array-like
         gRNA names, aligned to ``score`` and ``fdr``, used to label hit points.
    
    Returns
    -------
    T_gene : pandas.DataFrame
        All points as [(gRNA,) gene, score, fdr].
    """
    thr_fdr = significance_threshold
    fdr = np.asarray(fdr, dtype=float)
    LPV = -np.log10(fdr)

    if grna is not None:
        T_gene = pd.DataFrame({
            'gRNA': np.asarray(grna),
            'gene': np.asarray(genes),
            'score': np.asarray(score),
            'fdr': fdr,
        })
    else:
        T_gene = pd.DataFrame({
            'gene': np.asarray(genes),
            'score': np.asarray(score),
            'fdr': fdr,
        })

    indh = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s >= thr_scoreh)]
    indd = [i for i, (f, s) in enumerate(zip(fdr, score)) if (f <= thr_fdr and s <= thr_scored)]

    T_hs = T_gene.iloc[indh].copy() if indh else pd.DataFrame()
    T_ds = T_gene.iloc[indd].copy() if indd else pd.DataFrame()
    
    if ax is not None:
        fin = LPV[np.isfinite(LPV)]
        ma = max(5, fin.max() + 1 if fin.size else 5)
        ax.plot(score, LPV, 'pk', markersize=4, label='Not significant')
        ax.grid(True)
        ax.set_ylabel('-log10(q-value)', fontsize=14)
        ax.set_ylim([0, ma])

        if not T_hs.empty:
            ax.plot(T_hs['score'], -np.log10(T_hs['fdr']), 'pr', linewidth=2, label='Enriched')
            for i in range(min(20, len(T_hs))):
                ax.text(T_hs['score'].iloc[i], -np.log10(T_hs['fdr']).iloc[i], T_hs['gene'].iloc[i], va='bottom', ha='right')

        if not T_ds.empty:
            ax.plot(T_ds['score'], -np.log10(T_ds['fdr']), 'pc', label='Depleted')
            for i in range(min(20, len(T_ds))):
                ax.text(T_ds['score'].iloc[i], -np.log10(T_ds['fdr']).iloc[i], T_ds['gene'].iloc[i], va='bottom', ha='right')

        ax.legend()

    return T_gene
