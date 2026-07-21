import numpy as np
import pandas as pd

def perGene_4_med_horiz(Genes_4, wt_4):
    """
    Processes sets of 4 gRNAs per gene, calculates statistics (mean, median, std)
    for LFC, Z, and Q values and returns pandas DataFrames

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
    
    # Create pandas DataFrames
    
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
