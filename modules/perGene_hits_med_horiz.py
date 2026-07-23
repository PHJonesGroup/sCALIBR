import numpy as np
from .pergene_med_horiz import pergene_med_horiz
from .general_volcano import general_volcano

def pergene_hits_med_horiz(alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz,
                             T_vert, cond, control_type, plot=False):
    T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z = pergene_med_horiz(T_vert, control_type)

    if not T_lfc_z_q_med.empty:
        fdr   = T_lfc_z_q_med['median_q'].values
        genes = T_lfc_z_q_med['genes'].values
        _, indha,  indda,  _, _, _ = general_volcano(
            alf, sfdr_corr, thr_lfch, thr_lfcd,
            T_lfc_z_q_med['median_LFC'].values, fdr, cond, genes, plot=False)
        _, indhaz, inddaz, _, _, _ = general_volcano(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med['median_Z'].values, fdr, cond, genes, plot=False)
    else:
        indha = indda = indhaz = inddaz = []

    return T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z, indha, indda, indhaz, inddaz
