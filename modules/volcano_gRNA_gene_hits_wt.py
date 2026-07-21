import matplotlib.pyplot as plt
import os
from .general_volcano import general_volcano
from .perGene_4_hits_med_horiz import perGene_4_hits_med_horiz

def volcano_gRNA_gene_hits_wt(alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond, output_dir):
    genes = T_vert['gene']
    scoreZ = T_vert['Z_zGE_LFC']
    fdr = T_vert['Q']

    # Prepare figure with 3 subplots side by side
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    plt.suptitle(f"Replica: {cond}", fontsize=20)
    # 1. Volcano per gRNA (Z-normalized)
    LPV, indha, indda, T_ds, T_hs, T_gRNA = general_volcano(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz, scoreZ, fdr, cond, genes, plot=True, ax=axs[0]
    )
    axs[0].set_xlabel('LFC standardized', fontsize=14)
    axs[0].set_title(f'Volcano per gRNA', fontsize=16)

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
        axs[1].set_title(f'Volcano per gene LFC', fontsize=16)

        # Per gene Z volcano
        general_volcano(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med.iloc[:, 2].values, fdr_gene.values, cond, genes_gene.values,
            plot=True, ax=axs[2]
        )
        axs[2].set_xlabel('Z-normalised LFC', fontsize=14)
        axs[2].set_title(f'Volcano per gene Z', fontsize=16)
    else:
        # If no gene data, hide those axes
        axs[1].axis('off')
        axs[2].axis('off')

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"volcano_gRNA_{cond}"), dpi=300, bbox_inches="tight")
    plt.close()


    return T_gRNA, T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z, indha, indda, indhaz, inddaz

