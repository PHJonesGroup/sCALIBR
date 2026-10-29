import pandas as pd
import numpy as np
import os
import yaml
from pathlib import Path
import argparse
import sys
from importlib.metadata import version, PackageNotFoundError

from modules import run_log, plot_grna_distribution
from modules import normalise_total, plot_gene_variability, plot_pca, plot_heatmap
from modules import norm_table_individ, plot_mean_median
from modules import separate_target_control
from modules import CTR_stats
from modules import p_control_target_implement
from modules import computeZ, make_tables_Z, plot_histograms
from modules import volcano_gRNA_gene_hits, volcano_gRNA_gene_hits_interactive

try:
    __version__ = version("scalibr")
except PackageNotFoundError:
    __version__ = "unknown"        

def parse_args():
    parser = argparse.ArgumentParser(
        prog="scalibr",
        description="sCALIBR: Small CRISPR Analysis via Locally calIBrated Reference. "
                    "Identifies significant gRNA enrichment/depletion in targeted CRISPR "
                    "screens by calibrating against true-neutral control gRNAs.",
        epilog="Example: scalibr --config config_WT.yaml",
    )
    parser.add_argument(
        "--config",
        default="config.yaml",
        metavar="FILE",
        help="Path to the YAML config file (default: config.yaml)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
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

def main():
    args = parse_args()
    config = load_config(args.config)

    input_counts = Path(config["input_counts"])
    output_dir = Path(config["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True) 

    keep_counts = config["keep_counts"]
    alf = config["alf"]
    norm_method = config["normalisation_method"]
    valid_methods = {"sample", "pooled"}
    if norm_method not in valid_methods:
        sys.exit(f"Config error: normalisation_method must be one of "
                 f"{sorted(valid_methods)}, got '{norm_method}'")
    target_type = config["target_type"]
    control_type = config["control_type"]       # one gene_type, or a list to pool
    if isinstance(control_type, str):
        control_type = [control_type]
    if (not isinstance(control_type, list) or not control_type
            or not all(isinstance(c, str) for c in control_type)):
        sys.exit(f"Config error: control_type must be a gene_type name or a list of them, "
                 f"got {config['control_type']!r}")
    control_label = " + ".join(control_type)    # used in plot titles, file names, stats
    step = config["step"]       # histogram range (st, en) is set per contrast from the LFCs
    thr_lfchz = config["thr_lfchz"]
    fisher_scale_factor = config.get("fisher_scale_factor", 1/3)   # gene-level Fisher's method
    thr_lfcdz = -thr_lfchz

    sfdr_corr = 0.0001 # constant to avoid log10(0)
    ssc     = 1 # small sample correction count

    # Record this run's settings in <output_dir>/run_log_<timestamp>.txt
    log = run_log.RunLog(output_dir, args.config, __version__)
    log.section("Settings")
    log.field("normalisation_method", norm_method)
    log.field("target_type", target_type)
    log.field("control_type", f"{control_label} ({len(control_type)} type(s) pooled)")
    log.field("keep_counts", keep_counts)
    log.field("alf", alf)
    log.field("thr_lfchz / thr_lfcdz", f"{thr_lfchz} / {thr_lfcdz}")
    log.field("step", step)
    log.field("Histogram range", "set per contrast from the data (see Contrasts)")
    log.field("Fixed constants", f"sfdr_corr = {sfdr_corr}, ssc = {ssc}")
    log.field("Gene-level test", "median normalised LFC; Fisher's method on gRNA p-values, "
              "BH across target genes" + (" (per replicate, then replicates pooled by Fisher)"
                                          if norm_method == "sample" else ""))
    log.field("fisher_scale_factor", fisher_scale_factor)

    # ============================ MODULE 1 ============================
    T_vert = pd.read_csv(input_counts)
    log.input_file(input_counts, T_vert)
    gene_names = T_vert.iloc[:, 1]  # second column (gene names)

    # Number of gRNAs per gene/intergenic name, counted over the whole file
    grna_per_gene = gene_names.value_counts()
    plot_grna_distribution.plot_grna_distribution(grna_per_gene.tolist(), output_dir)

    # Removes genes with incorrect number of gRNAs
    keep_rows = gene_names.map(grna_per_gene).isin(keep_counts)
    T_vert_f = T_vert[keep_rows].reset_index(drop=True)
    log.field("Kept by keep_counts", f"{len(T_vert_f)} of {len(T_vert)} gRNAs "
                                     f"({grna_per_gene.isin(keep_counts).sum()} genes)")
    log.config_text()
    log.section("Contrasts")
    plot_pca.plot_pca(T_vert_f, output_dir, name="raw")

    # Normalise counts within a sample
    norm_dat = normalise_total.normalise_total(ssc, T_vert_f, output_dir) 
    norm_dat.to_csv(os.path.join(output_dir, 'normalised_counts.csv'), index=False)
    plot_pca.plot_pca(norm_dat, output_dir, name="normalised" )
    plot_heatmap.plot_heatmap(norm_dat, output_dir)

    # loops over contrasts
    for contrast in config["contrasts"]:
        name     = contrast["name"]
        baseline = contrast["baseline"]
        cond1    = contrast["treatment"]

        if len(baseline) != len(cond1):
            log.stop(f"Contrast {name}: baseline and cond1 must have equal length")

        c_outdir = output_dir / name
        c_outdir.mkdir(parents=True, exist_ok=True)

        # filter dataframe for count
        wanted_prefixes = baseline + cond1          # e.g. ['T0_F','T0_M','T1_F+','T1_M+']
        keep = list(norm_dat.columns[:3]) + [
            c for c in norm_dat.columns[3:]
            if any(c == p for p in wanted_prefixes)
        ]
        norm_dat_filtered = norm_dat[keep]
        plot_pca.plot_pca(norm_dat_filtered, c_outdir, name="normalised" )

        plot_gene_variability.plot_gene_variability(norm_dat_filtered, c_outdir, baseline, cond1, target_type)
        
        # Format normalised counts, compute mean and median of normalised counts
        tab_norm_T0, tab_norm_T1 = norm_table_individ.norm_table_individ(norm_dat_filtered, baseline, cond1)
        plot_mean_median.plot_mean_median(tab_norm_T0, tab_norm_T1, c_outdir, mode='replicate')
        plot_mean_median.plot_mean_median(tab_norm_T0, tab_norm_T1, c_outdir,mode='pooled')

        # ============================ MODULE 2 ============================
        # Separate and show controls, NTs, zGE, and normalised counts distribution
        T_target, T_control, hist, bin_edges, st, en = separate_target_control.separate_target_control(
            step, norm_dat_filtered, target_type, control_type, baseline, cond1, name, c_outdir)

        # Control null: mean/SD and p-curve (pooled, or one per replicate)
        me_sd, p_cont = CTR_stats.CTR_stats(
            st, en, step, T_control, hist, name, control_label, norm_method, c_outdir)
        log.write(f"{name}: {', '.join(baseline)} -> {', '.join(cond1)} | "
                  f"{len(T_target)} target, {len(T_control)} control gRNAs | "
                  f"histogram range [{st}, {en}]")

        # pooled: me_sd is (2,); sample: (n_rep, 2) -> comma-separated per replicate
        me_sd_2d = np.atleast_2d(me_sd)
        mean_str = ",".join(f"{v:.6f}" for v in me_sd_2d[:, 0])
        sd_str   = ",".join(f"{v:.6f}" for v in me_sd_2d[:, 1])

        summary_path = os.path.join(output_dir, "control_stats_by_contrast.txt")
        write_header = not os.path.exists(summary_path)
        with open(summary_path, "a") as fh:
            if write_header:
                fh.write("contrast\tbaseline\ttreatment\tcontrol_type\tnorm_method\t"
                        "n_control_gRNA\tmean\tsd\n")
            fh.write(f"{name}\t{','.join(baseline)}\t{','.join(cond1)}\t"
                    f"{control_label}\t{norm_method}\t{len(T_control)}\t"
                    f"{mean_str}\t{sd_str}\n")
            
        # Get recalibrated p/q values for gRNAs
        T_vert_q = p_control_target_implement.p_control_target_implement(T_target, T_control, bin_edges, p_cont)

        # keep targets only: controls were needed to calibrate Q / Z, not for
        # the target histogram, target table or volcano plots
        is_target = T_vert_q['gene_type'].astype(str).str.strip() == target_type
        T_vert_q = T_vert_q[is_target].reset_index(drop=True)

        # Compute Z LFC for targets only
        lfc_cols = [c for c in T_vert_q.columns if c.startswith('lfc')]
        LFC_t = T_vert_q[lfc_cols].to_numpy(dtype=float)

        Z_t, plot_bins, perc_t, perc_zt = computeZ.computeZ(st, en, step, LFC_t, me_sd)

        plot_histograms.plot_histograms(
            bin_edges=plot_bins, perc_lfc=perc_t, perc_z=perc_zt, 
            title=f'Distribution of {target_type} gRNAs: {name}',
            save_name=f'distri_LFC_{target_type}_{name}', 
            output_dir=c_outdir)

        T_t = make_tables_Z.make_tables_Z(T_vert_q, Z_t)

        T_t.to_csv(os.path.join(c_outdir, f'target_{name}_gRNA.csv'), index=False)

        # -------------------- 2.5 --------------------
        # Find extreme sets (enriched/ depleted) and volcano for gRNA and gene
        T_gRNA, T_gene = volcano_gRNA_gene_hits.volcano_gRNA_gene_hits(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz, T_t, name, c_outdir, norm_method, fisher_scale_factor)

        volcano_gRNA_gene_hits_interactive.volcano_gRNA_gene_hits_interactive(
            alf, sfdr_corr, thr_lfchz, thr_lfcdz, T_t, name, c_outdir, norm_method, fisher_scale_factor)

        T_gRNA.to_csv(os.path.join(c_outdir, f"volcano_{name}_grna.csv"), index=False)
        T_gene.to_csv(os.path.join(c_outdir, f"volcano_{name}_gene.csv"), index=False)

    log.finish()

if __name__ == "__main__":
    main()