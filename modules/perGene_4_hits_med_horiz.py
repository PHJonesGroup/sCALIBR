import numpy as np
from .count_vertical_names import count_vertical_names
from .separate_fours_threes_twos_genes import separate_fours_threes_twos_genes
from .perGene_4_med_horiz import perGene_4_med_horiz
from .general_volcano import general_volcano

def perGene_4_hits_med_horiz(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond, plot=False):
    """
    Analyze gene hits using median of 4-gRNA targets. Generates volcano plots
    only when plot=True.

    Returns:
    - T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z
    - indha, indda, indhaz, inddaz: hit/depleted indices for LFC and Z
    """
    gene_names = T_vert['gene'].values

    ggenes, gg, ind_gn = count_vertical_names(gene_names)
    ind_num = np.column_stack((ind_gn, gg))

    gene_counts = T_vert['gene'].value_counts()
    (Genes_2, Genes_3, Genes_4, gn_2, gn_3, gn_4,
     wt_2, wt_3, wt_4, d4, nums4) = separate_fours_threes_twos_genes(
        T_vert, ggenes, gene_names, ind_gn, gg
    )
    sum_nums = np.sum(nums4)

    if sum_nums > 0:
        T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z = perGene_4_med_horiz(Genes_4, wt_4)

        fdr = T_lfc_z_q_med.iloc[:, 3]
        genes = T_lfc_z_q_med.iloc[:, 0]
        score = T_lfc_z_q_med.iloc[:, 1]  
        # median_LFC
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

