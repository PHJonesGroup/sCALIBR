import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import yaml
from pathlib import Path
import argparse

from modules import count_vertical_names, plot_grna_distribution
from modules import separate_fours_threes_twos_ones_genes_gRNAs
from modules import normalise_prop
from modules import norm_table_individ, plot_me_med
from modules import find_columns_indiv, format_data

from modules import compute_hiss_LFC_rep12, distri_target_contr_plots_all, filter_pattern_distri, make_histo_LFC, separate_target_control, zGE_target_distri
from modules import compute_p_critLFC, CTR_stats_zGE, make_histo_crit_stats, make_histo_vec_rep12, make_LFC_Z_MZ_tables_two, med_mad_MZNP_2, q_val_frequentist_critical
from modules import implement_p_control_indiv_gRNA, p_control_any_implement, p_control_target_implement, plot_histograms
from modules import computeZ, make_tables_Z_two, z_p_CTR_any
from modules import perGene_4_hits_med_horiz, perGene_4_med_horiz, separate_fours_threes_twos_genes, volcano_gRNA_gene_hits_wt_interactive, volcano_gRNA_gene_hits_wt
from modules import separate_genes_by_grna_count

def parse_args():
    parser = argparse.ArgumentParser(description="Run the pipeline.")
    parser.add_argument(
        "--config",
        default="config.yaml",
        help="Path to the YAML config file (default: config.yaml)",
    )
    return parser.parse_args()

def load_config(path):
    if not os.path.exists(path):
        sys.exit(f"Config error: file not found: {path}")
    with open(path, "r") as f:
        try:
            return yaml.safe_load(f)
        except yaml.YAMLError as e:
            sys.exit(f"Config error: could not parse {path}\n{e}")

args = parse_args()
config = load_config(args.config)

input_counts = Path(config["input_counts"])
input_controls = config.get("input_controls")
output_dir = Path(config["output_dir"])
rep = config["rep"]
keep_counts = config["keep_counts"]
control = config["control"]
alf = config["alf"]

pat1 = config["pat1"]
pat2 = config["pat2"]
baseline = config["baseline"]
cond1 = config["cond1"]
rep_pairs = config["rep_pairs"]

st = config["st"]
en = config["en"]
step = config["step"]
thr_lfchz = config["thr_lfchz"]
thr_lfcdz = -thr_lfchz

# condz   = f'zGE_{cond1}{cond2}'
# cond11  = 'Intergenic'
# cond12  = 'Non_Targetting'
# indiv   = f'{cond1}{cond2}'  
# cond_id = f'{cond1}_{cond2}'
sfdr_corr = 0.0001 # constant to avoid log10(0)
ssc     = 1 # small sample correction count


# ============================ MODULE 1 ============================
raw_ind = pd.read_csv(input_counts)
si_input = raw_ind.shape

T_vert = raw_ind
gRNA = raw_ind.iloc[:, 0].values  # first column (gRNA names)
gene = raw_ind.iloc[:, 1].values  # second column (gene names)

# filter dataframe for count
keep = list(T_vert.columns[:2]) + [
    c for c in T_vert.columns[2:]
    if c == baseline or c == cond1
]
T_vert = T_vert[keep]

print(T_vert)

# -------------------- 1.1 --------------------
# Counts how many times each gene name/intergenic name occurs in the raw count file & gives number of gRNA per gene/ intergenic
gene_names = gene
ggenes, gg, ind_gn = count_vertical_names.count_vertical_names(gene_names)
plot_grna_distribution.plot_grna_distribution(gg, output_dir)

# -------------------- 1.2 --------------------
# Removes genes with incorrect number of gRNAs
groups, nums = separate_genes_by_grna_count.separate_genes_by_grna_count(T_vert, ggenes, gene_names, ind_gn, gg, rep)
ind_keep = sum((groups[c]["ind"] for c in keep_counts if c in groups), [])
raw_ind_f = T_vert.iloc[ind_keep, :].reset_index(drop=True)

# -------------------- 1.3 --------------------
# Normalise (CPM) counts as percentage within a column
norm_dat = normalise_prop.normalise_prop(ssc, raw_ind_f, output_dir) 

# -------------------- 1.4 --------------------
# Format normalised counts, compute mean and median of normalised counts
tab_norm_T0, tab_norm_T1 = norm_table_individ.norm_table_individ(norm_dat, raw_ind_f, baseline, cond1, rep)

me_med_T0_T1 = plot_me_med.plot_me_med(tab_norm_T0, tab_norm_T1, output_dir)

T_norm_WT = pd.concat([tab_norm_T0, tab_norm_T1.iloc[:, 2:(2 + 2*rep)]], axis=1)
T_norm_WT.to_csv(os.path.join(output_dir, 'normalised_counts.csv'), index=False)

norm_tab = T_norm_WT.copy()
column_labels = norm_tab.columns.tolist()
