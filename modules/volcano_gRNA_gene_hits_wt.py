import matplotlib.pyplot as plt
import os
from .general_volcano import general_volcano
from .perGene_4_hits_med_horiz import perGene_4_hits_med_horiz

def volcano_gRNA_gene_hits_wt(alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz,
                              T_vert, cond, control_type, output_dir):
    genes  = T_vert['gene']
    scoreL = T_vert['lfc']                          # raw per-gRNA LFC
    scoreZ = T_vert[f'Z_{control_type}_lfc']        # standardized per-gRNA LFC
    fdr    = T_vert['Q']

    # 2x2 grid: row 0 = per gRNA (LFC, Z), row 1 = per gene (LFC, Z)
    fig, axs = plt.subplots(2, 2, figsize=(14, 12))
    plt.suptitle(f"Replica: {cond}", fontsize=20)

    # --- Row 0, Col 0: per-gRNA LFC ---
    general_volcano(
        alf, sfdr_corr, thr_lfch, thr_lfcd,
        scoreL, fdr, cond, genes, plot=True, ax=axs[0, 0]
    )
    axs[0, 0].set_xlabel('LFC', fontsize=14)
    axs[0, 0].set_title('Volcano per gRNA (LFC)', fontsize=16)

    # --- Row 0, Col 1: per-gRNA Z ---
    LPV, indha, indda, T_ds, T_hs, T_gRNA = general_volcano(
        alf, sfdr_corr, thr_lfchz, thr_lfcdz,
        scoreZ, fdr, cond, genes, plot=True, ax=axs[0, 1]
    )
    axs[0, 1].set_xlabel('LFC standardized', fontsize=14)
    axs[0, 1].set_title('Volcano per gRNA (Z)', fontsize=16)

    # --- Per-gene stats ---
    (
        T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z,
        indha_gene, indda_gene, indhaz, inddaz
    ) = perGene_4_hits_med_horiz(
        alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond, plot=False
    )

    if not T_lfc_z_q_med.empty:
        fdr_gene   = T_lfc_z_q_med.iloc[:, 3]
        genes_gene = T_lfc_z_q_med.iloc[:, 0]

        # --- Row 1, Col 0: per-gene LFC ---
        general_volcano(
            alf, sfdr_corr, thr_lfch, thr_lfcd,
            T_lfc_z_q_med.iloc[:, 1].values, fdr_gene.values, cond, genes_gene.values,
            plot=True, ax=axs[1, 0]
        )
        axs[1, 0].set_xlabel('LFC', fontsize=14)
        axs[1, 0].set_title('Volcano per gene (LFC)', fontsize=16)

        # --- Row 1, Col 1: per-gene Z ---
        general_volcano(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz,
            T_lfc_z_q_med.iloc[:, 2].values, fdr_gene.values, cond, genes_gene.values,
            plot=True, ax=axs[1, 1]
        )
        axs[1, 1].set_xlabel('Z-normalised LFC', fontsize=14)
        axs[1, 1].set_title('Volcano per gene (Z)', fontsize=16)
    else:
        axs[1, 0].axis('off')
        axs[1, 1].axis('off')

    # remove the per-subplot legends that general_volcano added
    for ax in axs.flat:
        leg = ax.get_legend()
        if leg is not None:
            leg.remove()

    handles, labels = axs[0, 1].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc='center left', bbox_to_anchor=(1.0, 0.5),
               frameon=True,
               edgecolor='black',        # box border colour
               facecolor='white',        # box background
               framealpha=1.0,           # opaque
               fontsize=12)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"volcano_gRNA_{cond}.png"),
                dpi=300, bbox_inches="tight")
    plt.close()

    return T_gRNA, T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z, indha, indda, indhaz, indda