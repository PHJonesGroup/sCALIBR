import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def separate_fours_threes_twos_ones_genes_gRNAs(T_vert, ggenes, gene_names, ind_gn, gg):
    """
    Given a table of genes and their data, separates genes based on the number of gRNAs (1,2,3,4).
    
    Inputs:
        T_vert: pandas DataFrame, with gene data (columns assumed: gene, lfc_1, lfc_2, dif_1, dif_2)
        ggenes: list of unique gene names corresponding to gg groups
        gene_names: list of gene names for all entries (sorted/grouped)
        ind_gn: list of start indices (1-based) for each gene group in gene_names
        gg: list of counts of gRNAs per gene (length matches ggenes)
    
    Outputs:
        Genes_1, Genes_2, Genes_3, Genes_4: lists of gene names by number of gRNAs
        gn_1, gn_2, gn_3, gn_4: lists of gene names repeated per gRNA counts
        wt_1, wt_2, wt_3, wt_4: numpy arrays of gene data per group
        d: max number of gRNAs if uniform, else max (variability warning)
        nums: list of counts of genes with 4,3,2,1 gRNAs
        ind4, ind3, ind2, ind1: indices in T_vert for each group
    """

    k4 = k3 = k2 = k1 = 0
    
    gn_4 = []
    wt_4 = []
    ind4 = []
    Genes_4 = []
    
    gn_3 = []
    wt_3 = []
    ind3 = []
    Genes_3 = []
    
    gn_2 = []
    wt_2 = []
    ind2 = []
    Genes_2 = []
    
    gn_1 = []
    wt_1 = []
    ind1 = []
    Genes_1 = []

    wt = T_vert.iloc[:, 2:6].to_numpy()

    for i in range(len(gg)):
        start_idx = ind_gn[i] - 1

        if gg[i] == 4:
            k4 += 1
            Genes_4.append(ggenes[i])
            gn = gene_names[start_idx:start_idx + 4]
            gn_4.extend(gn)
            indd4 = list(range(start_idx, start_idx + 4))
            ind4.extend(indd4)
            wtt = wt[start_idx:start_idx + 4, :]
            wt_4.append(wtt)

        elif gg[i] == 3:
            k3 += 1
            Genes_3.append(ggenes[i])
            gn = gene_names[start_idx:start_idx + 3]
            gn_3.extend(gn)
            indd3 = list(range(start_idx, start_idx + 3))
            ind3.extend(indd3)
            wtt = wt[start_idx:start_idx + 3, :]
            wt_3.append(wtt)

        elif gg[i] == 2:
            k2 += 1
            Genes_2.append(ggenes[i])
            gn = gene_names[start_idx:start_idx + 2]
            gn_2.extend(gn)
            indd2 = list(range(start_idx, start_idx + 2))
            ind2.extend(indd2)
            wtt = wt[start_idx:start_idx + 2, :]
            wt_2.append(wtt)

        elif gg[i] == 1:
            k1 += 1
            Genes_1.append(ggenes[i])
            gn = gene_names[start_idx:start_idx + 1]
            gn_1.extend(gn)
            indd1 = [start_idx]
            ind1.extend(indd1)
            wtt = wt[start_idx:start_idx + 1, :]
            wt_1.append(wtt)

    # Stack arrays vertically or create empty arrays if none
    wt_4 = np.vstack(wt_4) if wt_4 else np.empty((0, wt.shape[1]))
    wt_3 = np.vstack(wt_3) if wt_3 else np.empty((0, wt.shape[1]))
    wt_2 = np.vstack(wt_2) if wt_2 else np.empty((0, wt.shape[1]))
    wt_1 = np.vstack(wt_1) if wt_1 else np.empty((0, wt.shape[1]))

    max_num_gRNA = max(gg)
    min_num_gRNA = min(gg)

    if max_num_gRNA == min_num_gRNA:
        d = min_num_gRNA
    else:
        d = max_num_gRNA

    nums = [k4, k3, k2, k1]

    return (Genes_1, Genes_2, Genes_3, Genes_4,
            gn_1, gn_2, gn_3, gn_4,
            wt_1, wt_2, wt_3, wt_4,
            d, nums,
            ind4, ind3, ind2, ind1)
