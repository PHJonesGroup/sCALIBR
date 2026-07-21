
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
        start_idx = ind_gn[i] - 1  
        
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
        d = min_num_gRNA
    else:
        d = max_num_gRNA
        
    nums = [k4, k3, k2, k1]
    
    return Genes_2, Genes_3, Genes_4, gn_2, gn_3, gn_4, wt_2, wt_3, wt_4, d, nums
