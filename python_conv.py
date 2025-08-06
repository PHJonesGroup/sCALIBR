import pandas as pd
import numpy as np
from format_utils import module1_format_data, find_columns_indiv, filter_pattern_distri, zGE_target_distri, distri_target_contr_plots_all, separate_target_control, compute_hiss_LFC_rep12, make_histo_LFC, make_histo_crit_stats
from format_utils import make_LFC_Z_MZ_tables_two, q_val_frequentist_critical, CTR_stats_zGE, make_histo_vec_rep12, compute_p_critLFC, med_mad_MZNP_2
from format_utils import module_p_control_target_implement, implement_p_control_indiv_gRNA, computeZ, z_p_CTR_any, make_tables_Z_two, module_p_control_any_implement
from format_utils import volcano_gRNA_gene_hits_wt, general_volcano, perGene_4_hits_med_horiz, count_vertical_names, separate_fours_threes_twos_genes, perGene_4_med_horiz, plot_histograms

print("0. INPUTS: (1)zGE and (2) normalised counts both replicas")

# ==== Load possible zGE datasets (only one active)
zGE = pd.read_csv('data_norm/genes_zGE_official.csv')  # official
zGE = zGE.drop_duplicates() 
# Load WT data (choose one)
norm_wt_may = pd.read_csv('data_norm/T_norm_T0T1_indiv_WT_May_2.csv')
norm_tab = norm_wt_may.copy()

print(f"Shape of normalized table: {norm_tab.shape}")
column_labels = norm_tab.columns.tolist()
print("Column Labels:", column_labels)
print("Head of normalized table:")
print(norm_tab.head())

# Displayed for reference in MATLAB, we'll use print here
print("module 1: already normalised as RPKM, ordered T0 m1 m2... T1 m1 m2 ..")
print("1.1 define scheme (1 or 2): usually 2, for both replicas, F and M")

print("parameters and inputs module 1")
d = 4  # number of gRNAs per gene biologists planned to apply
num_gRNA = d
index_scheme = 2  # both replicas if 2; 1 if only one replica

indiv = 'FM'  # used in column labeling or processing
cond_id = 'WT_may2_F_M'  # used for saving

T_norm_indiv, norm_dat, me_med_nd = module1_format_data(index_scheme, indiv, norm_tab)
print("mean and median of T0 and T1, normalised raw data, all two WT valid mice, May2")
print("Same format as norm_tab_head if scheme=2; with double columns if scheme=1")

# Show head of T_norm_indiv
print(T_norm_indiv.head())


# ============================ MODULE 2 ============================
print('=' * 30)
print('MODULE 2: Distribution of:')
print(' (1) Controls: intergenic enzyme cuttings')
print(' (2) NT = non-targeted')
print(' (3) zGE = zero GE genes')
print(' (4) Filtered targeted genes')
print('Distributions are in the same bins [-20, 10] with step 0.05 for all histograms')

# -------------------- HARD CODED PARAMETERS --------------------
alf = 0.06  # significance level of right-left tail of distribution

st = -20
en = 10
step = 0.05
start_end = [st, en]

name = 'WT_val_May2'
cond1 = 'WTval May2 F'
cond2 = 'WTval May2 M'
condz = 'WT zGE FM'
condi = 'WT interg FM'

print('\nKasumi scheme (further) for:')
print(' (i) individual mice (twice T0 T0 T1 T1) or')
print(' (ii) M vs F comparisons')
print('\nMODULE 2.1: Separate and show controls, NTs, zGE, and normalised counts')

# Patterns for filtering (should be known in advance)
pat1 = 'chr'
pat2 = 'Non'
n = 3

# Check full duplicate rows
dup_rows = T_norm_indiv[T_norm_indiv.duplicated()]
print(f"\n T_norm_indiv: {len(dup_rows)} duplicate rows found")
if not dup_rows.empty:
    print(dup_rows.head())

# Assuming gRNA is in the first column
gRNA_col = T_norm_indiv.columns[0]
dup_gRNAs = T_norm_indiv[gRNA_col][T_norm_indiv[gRNA_col].duplicated()]
print(f"\n T_norm_indiv: {dup_gRNAs.nunique()} duplicated gRNAs")
if not dup_gRNAs.empty:
    print(dup_gRNAs.value_counts().head())


# -------------------- Run Separation --------------------
(
    T_target, T_target_zGE_counts, T_lfc_chr, T_lfc_nt, T_zGE,
    bin, his1z, perc1z, his2z, perc2z, hisFMz, percFMz,
    his1t, perc1t, his2t, perc2t, hisFMt, percFMt, num_gzt
) = separate_target_control(
    d=4, st=st, en=en, step=step,
    T_norm_indiv=T_norm_indiv, zGE=zGE,
    pat1=pat1, pat2=pat2, n=n
)


# -------------------- Inspect Results --------------------
print('\nNumber of zGE and Target genes:')
print(num_gzt)

print('\nNumber of gRNAs per data set:')
print(f'Target + zGE set: {T_target_zGE_counts.shape}')
print(f'Target only:      {T_target.shape}')
print(f'Chr controls:     {T_lfc_chr.shape}')
print(f'NT controls:      {T_lfc_nt.shape}')
print(f'zGE genes:        {T_zGE.shape}')

print('\nExample Target gRNAs:')
print(T_target)

print('\nExample Chr Controls:')
print(T_lfc_chr.head(6))

print('\nExample NT Controls:')
print(T_lfc_nt.head(6))

print('\nExample zGE Genes (sorted):')
print(T_zGE.head(6))


# ==========================================================
# MODULE 3  – choose zGE controls, compute Z/MZ/crit‑LR
# ==========================================================
print("=" * 25)
print("Module 3: choose control, implement q, adjust controls for Z")

print("module 3.1 — (i) choose between LFC distributions of zGE")
print("            (ii) compute true Z / MZ / critLR by zGE genes")
print("------------------------- rep 1 & 2 pooled")

# ----------------------------------------------------------
# Build 2‑column matrix of target histograms (rep1, rep2)
# ----------------------------------------------------------
# his1t / his2t come from Module 2
his12t = np.column_stack([his1t, his2t])        # (bins × 2)

condz1 = "WT zGE F"
condz2 = "WT zGE M"
condz  = "WT zGE pooled"

# If you ever need to compare against intergenic controls instead, keep:
# condi1 = "WT chr F"
# condi2 = "WT chr M"

# ----------------------------------------------------------
# Main zGE stats routine  ➜  returns stats + zGE tables
# ----------------------------------------------------------
(
    crit_LR12,
    me_sd12,
    med_mad12,
    binn,
    p_cont12,
    hiss_cont12,
    p_targ12,
    T_zGE_1,
    T_zGE_2,
) = CTR_stats_zGE(
    alf,
    st,
    en,
    step,
    T_zGE,          # DataFrame with gRNA / gene / lfc1 / lfc2
    his12t,         # 2‑col array of target histograms
    condz1,
    condz2,
    cond1,          # e.g. "WT val May2 F"
    cond2,          # e.g. "WT val May2 M"
    condz,
)

# ----------------------------------------------------------
# Show example output and key statistics
# ----------------------------------------------------------
print("\nFirst five zGE rows (rep‑1 view):")
print(T_zGE_1.head())

print("\ncrit_LR12 (rep1, rep2):")
print(crit_LR12)

print("\nmean / SD per replicate (me_sd12):")
print(me_sd12)

print("\nmedian / MAD per replicate (med_mad12):")
print(med_mad12)


print('module3.2 : p values to each gRNA=T_vert: currently Q12 is assigned from zGE')

print("T_target shape:", T_target.shape)
print("T_target columns:", T_target.columns.tolist())

print("T_zGE shape:", T_zGE.shape)
print("T_zGE columns:", T_zGE.columns.tolist())

print("T_zGE group counts (top):\n",
      T_zGE.groupby(['gRNA', 'gene']).size().sort_values(ascending=False).head())

print("bin shape:", bin.shape if hasattr(bin, 'shape') else type(bin))
print("p_cont12 shape:", p_cont12.shape if hasattr(p_cont12, 'shape') else type(p_cont12))

# Run the function to get recalibrated p/q values for gRNAs
T_vert_q = module_p_control_target_implement(T_target, T_zGE, bin, p_cont12)

# Display the first 5 rows similar to MATLAB's T_vert(1:5,:)
print(T_vert_q.head(5))
print("T_vert_q shape:", T_vert_q.shape)
print("T_vert_q head:\n", T_vert_q.head())

# ------------------------- scheme 2 (commented code in MATLAB, so here too)

# print('choose2: adjust me sd from chosen me med critLR')

# mu_1, sd_1, mu_2, sd_2 = compare_choose_Z_MZ_crit_LR(me_sd12, med_mad12, crit_LR12)

# print('new me sd via crit 1 2')
# me_sd_12 = np.array([[mu_1, sd_1], [mu_2, sd_2]])

# print('new me sd via me sd')
# print('new mu via me sd')

# me_sd_12[0, 0] = me_sd12[0, 0]  # Overwrite with original mean if needed
# me_sd_12[1, 0] = me_sd12[1, 0]

print('final me sd')
me_sd_12 = me_sd12


print('module3.3------Compute Z LFC for two replicas rep1 rep2, for targets only')

print(T_vert_q.columns)

LFC_t1 = T_vert_q.iloc[:, 2].to_numpy()
LFC_t2 = T_vert_q.iloc[:, 3].to_numpy()

Z_t1, Z_t2, bin, perc_t1, perc_zt1, perc_t2, perc_zt2 = computeZ(
    st, en, step, LFC_t1, LFC_t2, me_sd_12
)

plot_histograms(bin, perc_t1, perc_zt1, perc_t2, perc_zt2, "Target genes")

print('module3.4---------------make T_vert tables gRNA-based: gene4 LFC Z_zGE_LFC')
print(T_vert_q.head(10))

T_t1, T_t2 = make_tables_Z_two(T_vert_q, Z_t1, Z_t2, cond1, cond2)
print(T_t1.head(10))

print('NEW 25 June: additional module for Z any set: intergenic, NT, etc')

# ---- Intergenic ----
T_any = T_lfc_chr.copy()
cond11 = 'Interg WTval May2 F'
T_z_q_chr1, T_z_q_chr2, binn_chr, perc_t1_chr, perc_zt1_chr, perc_t2_chr, perc_zt2_chr = z_p_CTR_any(
    st, en, step, T_any, binn, p_cont12, me_sd_12, cond1, cond2
)

plot_histograms(binn_chr, perc_t1_chr, perc_zt1_chr, perc_t2_chr, perc_zt2_chr, cond11)

# ---- NT control ----
T_any = T_lfc_nt.copy()
cond12 = 'NT WTval May2 F'
T_z_q_nt1, T_z_q_nt2, binn_nt, perc_t1_nt, perc_zt1_nt, perc_t2_nt, perc_zt2_nt = z_p_CTR_any(
    st, en, step, T_any, binn, p_cont12, me_sd_12, cond1, cond2
)

plot_histograms(binn_nt, perc_t1_nt, perc_zt1_nt, perc_t2_nt, perc_zt2_nt, cond12)

print('module4.1-----find extreme sets and volcano for gRNA and gene')

alf = 0.06

sfdr_corr = 0.0001
# ------------------------- for LFC
thr_lfchz = 1.7
thr_lfcdz = -1.7

# ---------------- rep1
thr_lfch = crit_LR12[0, 1]  # MATLAB indexing (1,2) is Python (0,1)
thr_lfcd = crit_LR12[0, 0]  # MATLAB (1,1) is Python (0,0)

thrLFC_d_h = [thr_lfcd, thr_lfch]

import os

print('module 4.2 used median values for LFC Z and Q')

#-----------------------rep1
cond = cond1  # 'WT May F'
T_vert = T_t1
print(T_vert.columns)
dupes = T_vert.groupby(['gRNA', 'gene']).size()
print("Max occurrences of any (gRNA, gene):", dupes.max())
print("Any duplicates:", (dupes > 1).any())
print(T_vert['condition'].value_counts())

# Call your function (make sure it's defined/imported in Python)
(T_gRNA_1, T_lfc_z_q_med_1, T_lfc_z_q_me_1, T_LFC_1, T_Q_1, T_Z_1, 
 indha_1, indda_1, indhaz_1, inddaz_1) = volcano_gRNA_gene_hits_wt(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond)

num_hd_LFC_Z_1 = [len(indha_1), len(indda_1), len(indhaz_1), len(inddaz_1)]

#-----------------------rep2
cond = cond2  # 'WT May M'
T_vert = T_t2

numeric_cols = ['LFC', 'Z_zGE_LFC', 'Q']
T_gene = T_vert.groupby('gene')[numeric_cols].median().reset_index()

gene_counts = T_vert['gene'].value_counts()
print("Number of genes with at least 4 gRNAs:", (gene_counts >= 4).sum())
print("Genes with exactly 4 gRNAs:", (gene_counts == 4).sum())

thr_lfch = crit_LR12[1, 1]  # MATLAB (2,2) -> Python (1,1)
thr_lfcd = crit_LR12[1, 0]  # MATLAB (2,1) -> Python (1,0)

(T_gRNA_2, T_lfc_z_q_med_2, T_lfc_z_q_me_2, T_LFC_2, T_Q_2, T_Z_2,
 indha_2, indda_2, indhaz_2, inddaz_2) = volcano_gRNA_gene_hits_wt(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_vert, cond)

num_hd_LFC_Z_2 = [len(indha_2), len(indda_2), len(indhaz_2), len(inddaz_2)]

#---------------------------savings

# Save DataFrames to CSV - adjust if your tables are not pandas DataFrames

T_lfc_z_q_med_1.head(6)  # Just displaying first 6 rows (similar to MATLAB)

# Save intermediate and result files
# Directory to save files
save_dir = '/Users/sb80/Library/CloudStorage/GoogleDrive-sb80@sanger.ac.uk/My Drive/CRISPR_pipeline'
T_z_q_chr1.to_csv(os.path.join(save_dir, 'T_intergenic_F_wt_gRNA_May2.csv'), index=False)
T_z_q_chr2.to_csv(os.path.join(save_dir, 'T_intergenic_M_wt_gRNA_May2.csv'), index=False)
T_z_q_nt1.to_csv(os.path.join(save_dir, 'T_NT_F_wt_gRNA_May2.csv'), index=False)
T_z_q_nt2.to_csv(os.path.join(save_dir, 'T_NT_M_wt_gRNA_May2.csv'), index=False)
T_t2.to_csv(os.path.join(save_dir, 'T_targetM_wt_gRNA_May2.csv'), index=False)
