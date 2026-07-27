import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import yaml
from pathlib import Path
import argparse

from modules import count_vertical_names, plot_grna_distribution
from modules import separate_genes_by_grna_count
from modules import normalise_prop, plot_gene_variability
from modules import norm_table_individ, plot_me_med
from modules import separate_target_control
from modules import CTR_stats
from modules import p_control_target_implement
from modules import computeZ, make_tables_Z, plot_histograms
from modules import volcano_grna_gene_hits, volcano_grna_gene_hits_interactive

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
output_dir = Path(config["output_dir"])
rep = config["rep"]
keep_counts = config["keep_counts"]
alf = config["alf"]

target_type = config["target_type"]
control_type = config["control_type"]
baseline = config["baseline"]
cond1 = config["cond1"]
rep_pairs = config["rep_pairs"]

st = config["st"]
en = config["en"]
step = config["step"]
thr_lfchz = config["thr_lfchz"]
thr_lfcdz = -thr_lfchz

sfdr_corr = 0.0001 # constant to avoid log10(0)
ssc     = 1 # small sample correction count

# ============================ MODULE 1 ============================
raw_ind = pd.read_csv(input_counts)
si_input = raw_ind.shape

T_vert = raw_ind
gRNA = raw_ind.iloc[:, 0].values  # first column (gRNA names)
gene = raw_ind.iloc[:, 1].values  # second column (gene names)

# filter dataframe for count
keep = list(T_vert.columns[:3]) + [
    c for c in T_vert.columns[3:]
    if c.startswith(f"{baseline}") or c.startswith(f"{cond1}")
]
T_vert = T_vert[keep]

# -------------------- 1.1 --------------------
# Counts how many times each gene name/intergenic name occurs in the raw count file & gives number of gRNA per gene/ intergenic
gene_names = gene
ggenes, gg, ind_gn = count_vertical_names.count_vertical_names(gene_names)
plot_grna_distribution.plot_grna_distribution(gg, output_dir)

# -------------------- 1.2 --------------------
# Removes genes with incorrect number of gRNAs
groups = separate_genes_by_grna_count.separate_genes_by_grna_count(T_vert, ggenes, gene_names, gg, ind_gn, rep)
ind_keep = sum((groups[c]["ind"] for c in keep_counts if c in groups), [])
raw_ind_f = T_vert.iloc[ind_keep, :].reset_index(drop=True)

# -------------------- 1.3 --------------------
# Normalise (CPM) counts as percentage within a column
norm_dat = normalise_prop.normalise_prop(ssc, raw_ind_f, output_dir) 
plot_gene_variability.plot_gene_variability(norm_dat, output_dir, baseline, cond1, target_type)

# -------------------- 1.4 --------------------
# Format normalised counts, compute mean and median of normalised counts
tab_norm_T0, tab_norm_T1 = norm_table_individ.norm_table_individ(norm_dat, baseline, cond1)

me_med_T0_T1 = plot_me_med.plot_me_med(tab_norm_T0, tab_norm_T1, output_dir)

T_norm = pd.concat([tab_norm_T0, tab_norm_T1.iloc[:, 3:(3 + 2*rep)]], axis=1)
T_norm.to_csv(os.path.join(output_dir, 'normalised_counts.csv'), index=False)

norm_tab = T_norm.copy()
column_labels = norm_tab.columns.tolist()

# ============================ MODULE 2 ============================
T_norm_indiv = T_norm.copy()

# -------------------- 2.1 --------------------
# Separate and show controls, NTs, zGE, and normalised counts
# Distribution of:
# (1) Controls: intergenic genes
# (2) NT = non-targeted genes
# (3) zGE = zero expressed genes
# (4) Normalised targeted genes

T_target, T_control, hist, groups, bin = separate_target_control.separate_target_control(
    st, en, step, T_norm_indiv,
    target_type, control_type,
    baseline, cond1, rep_pairs, output_dir
)

# -------------------- 2.2 --------------------
# Choose control, compute Z/MZ/crit‑LR, implement q, adjust controls for Z
(
    crit_LR, 
    me_sd, 
    med_mad, 
    binn, 
    p_cont, 
    hiss_cont, 
    p_targ, 
    T_control
) = CTR_stats.CTR_stats(
    alf,
    st,
    en,
    step,
    T_control,          
    hist,         
    baseline,          
    cond1,          
    control_type,
    output_dir
)

# -------------------- 2.3 --------------------
# Get recalibrated p/q values for gRNAs
T_vert_q = p_control_target_implement.p_control_target_implement(T_target, T_control, bin, p_cont)

# -------------------- 2.4 --------------------
# Compute Z LFC for targets only
lfc_cols = [c for c in T_vert_q.columns if c.startswith('lfc')]
LFC_t = T_vert_q[lfc_cols].to_numpy(dtype=float)     # (n_gRNA x 4)
Z_t, bin, perc_t, perc_zt = computeZ.computeZ(st, en, step, LFC_t, me_sd)
plot_histograms.plot_histograms(bin, perc_t, perc_zt, baseline, cond1, "Target genes", st, en, output_dir)

T_t = make_tables_Z.make_tables_Z(T_vert_q, Z_t, control_type)

T_t.to_csv(os.path.join(output_dir, f'target_{cond1}_{baseline}_gRNA.csv'), index=False)

# -------------------- 2.5 --------------------
# Find extreme sets (enriched/ depleted) and volcano for gRNA and gene
thr_lfch = crit_LR[1]     # right-tail critical LFC (enrichment)
thr_lfcd = crit_LR[0]     # left-tail critical LFC (depletion)
thrLFC_d_h = [thr_lfcd, thr_lfch]

(T_gRNA, T_lfc_z_q_med, T_lfc_z_q_me, T_LFC, T_Q, T_Z,
 indha, indda, indhaz, inddaz) = volcano_grna_gene_hits.volcano_grna_gene_hits(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_t, baseline, cond1, control_type, output_dir)

volcano_grna_gene_hits_interactive.volcano_grna_gene_hits_interactive(
    alf, sfdr_corr, thr_lfch, thr_lfcd, thr_lfchz, thr_lfcdz, T_t, baseline, cond1, control_type, output_dir)

num_hd_LFC_Z = [len(indha), len(indda), len(indhaz), len(inddaz)]