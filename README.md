# small CRISPR Analysis via Locally calIBrated Reference (sCALIBR)

This pipeline is designed to find significant gRNA enrichment/depletion accounting for both FDR and FNR rates. This is achieved by calculating the likelihood of the observed gRNA fold change for each gene differing from the distribution of the set of true neutral control gRNAs in the screen.

## System Requirements

### Operating Systems
This package has been tested on:
- macOS: Tahoe (26.5.1)
- Linux: Ubuntu 22.04 LTS
- Windows: 11

### Software Dependencies
- Python 3.13.5

See required packages and versions in requirement.txt

### Hardware Requirements
Runs on a standard computer with enough RAM for in-memory operations (e.g., 8 GB). <br />
No non-standard hardware required.

## Installation

Clone the repository:
```
git clone https://github.com/PHJonesGroup/sCALIBR.git
cd sCALIBR/
```

Create and activate a virtual environment:
```
python -m venv myvenv
source myvenv/bin/activate      # On Windows: myvenv\Scripts\activate
```

Install dependencies:
```
pip install -r requirements.txt
```

**Typical install time:** ~10 minutes on a standard desktop with a normal internet connection.

## Running the python script
To run on your own data:

1. Copy the example config: `cp config.yaml my_config.yaml`
2. Edit the parameters (see the table below).
3. Run: `python3 scalibr.py --config my_config.yaml`

### Configuration parameters

| Parameter         | Description                                               | Example                               |
|-------------------|-----------------------------------------------------------|---------------------------------------|
| `input_counts`    | Path to the counts CSV                                    | `demo_data/WT_....csv`                |
| `output_dir`      | Directory where results are written                       | `output`                              |
| `keep_counts`     | Number of gRNAs per gene                                  | `[1, 4, 10]`                          |
| `target_type`     | Choice of target genes                                    | `GOI`                                 |
| `control_type`    | Control gene type(s) ; a list (e.g. `['Neutral expressed gene', 'Zero-expressed gene']`) pools them into one control set | `Zero-expressed genes`|
| `contrasts`       | List of contrasts to run; each has a `name`, and `baseline`/`treatment` with lists of replicates (if present) | see below          |
| `step`            | Bin sizes for histograms                                  | `0.05`                                |
| `thr_lfchz`       | LFC threshold for volcano plot                            | `2`                                 |
| `significance_threshold` | Hit threshold: raw p-value for gRNAs, BH-adjusted Fisher p (q) for genes | `0.05` |

#### Contrasts

`contrasts` defines one or more comparisons to run. The normalisation steps run
once over all samples; each contrast is then scored and plotted independently,
with results written to a subfolder named after its `name`.

Each contrast has three fields:

- **`name`** — a label for the contrast, used for output filenames and subfolders.
- **`baseline`** — list of samples treated as the reference.
- **`treatment`** — list of samples treated as the treatment.

Example:

```yaml
contrasts:
  - name: T1_vs_T0             # two replicates contrast
    baseline: [T0_F, T0_M]
    cond1:    [T1_F, T1_M]
  - name: T1_vs_T0_F_only      # one replicate contrast
    baseline: [T0_F]
    cond1:    [T1_F]
```

### Format of count file

For the pipeline to run correctly, the input count file must follow this exact column order and naming convention:

| sgRNA_name     | gene         | gene_type        | `<baseline_name>_<replicateN>`  | `<cond_name>_<replicateN>`      |
|----------------|--------------|------------------|---------------------------------|---------------------------------|
| X_sgRNA1       | X            | GOI              |                                 |                                 |
| X_sgRNA2       | X            | GOI              |                                 |                                 |
| X_sgRNA3       | X            | GOI              |                                 |                                 |
| X_sgRNA4       | X            | GOI              |                                 |                                 |
| Y_sgRNA1       | Y            | Non-targetting   |                                 |                                 |
| Y_sgRNA2       | Y            | Non-targetting   |                                 |                                 |

**Column descriptions:**

- **`sgRNA_name`** — unique identifier for each guide RNA
- **`gene`** — the target gene for that guide
- **`gene_type`** — guide/gene category (e.g. `GOI`, `Essential gene`, `Non-targeting`, `Intergenic`, `Zero-expressed gene`)
- **`<baseline_name>_<replicateN>`, `<cond_name>_<replicateN>`** — raw read counts per sample, one column per replicate per condition

**Notes:**

- Replicates are optional. If your experiment has no replicates, include a single count column per condition (e.g. `<baseline_name>`, `<cond_name>`) instead of per-replicate columns.
- Column order must be preserved exactly as shown — `sgRNA_name`, `gene`, `gene_type`, followed by all baseline/condition count columns.

### Outputs
| File Name                                             | Description                                                                                                                       |
|-------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `distri_gene_types_{contrast_name}.png`               | Individual LFC distributions for all the different gene types (eg. genes of interests, zero expressed genes, non-targeting etc)   | 
| `distri_gene_types_overlaid_{contrast_name}.png`      | Overlaid LFC distributions for all the different gene types (eg. genes of interests, zero expressed genes, non-targeting etc)     |    
| `distri_LFC_{control_type}_{contrast_name}.png`       | Distribution of LFC, Z-corrected LFC and median Z-corrected LFC of chosen control genes                                           |    
| `distri_LFC_{target_type}_{contrast_name}.png`        | LFC and Z-corrected LFC gRNA distribution for target genes                                                                        | 
| `gene_distribution_gRNA.png`                          | Number of gRNAs attributed to each gene                                                                                           | 
| `gene_variability_boxplots.html`                      | Interactive box plots of variability between gRNAs for each gene                                                                  | 
| `gRNA_counts_normalisation.png`                       | Total raw (top) and normalised (bottom) read count summed over all gRNAs in each condition                                        | 
| `normalised_counts_mean_med.png`                      | Mean versus median of the normalised gRNA counts for each of condition                                                            | 
| `normalised_counts.csv`                               | Normalised counts in csv format                                                                                                   | 
| `p_distri_targ_cont_{contrast_name}.png`              | Control-calibrated p-value curves of target vs control gRNAs                                                                      | 
| `run_log_{YYYYmmdd_HHMMSS}.txt`                       | One per run: settings, config as written, input checksum, code version, per-contrast histogram range, and whether the run completed | 
| `target_{contrast_name}_gRNA.csv`                     | LFC, Z-corrected LFC, q-value, condition table for target genes                                                                   | 
| `volcano_gRNA_{contrast_name}.png`                    | Volcano plot on gRNA and gene level                                                                                               | 
| `volcano_plot_interactive_{contrast_name}.html`       | Interactive volcano plot on gRNA and gene  CSV                                                                                    | 

Each plot generated by this pipeline is saved alongside a `*_data.csv` file containing its underlying data.

## Pipeline Steps
1. Normalise raw count data as proportion per million <br />
2. Guide-level fold changes <br />
3. Control-based calibration <br />
4. Guide-agreement adjustment <br />
5. Multiple-test correction <br />
6. Hit reporting and contextual validation <br />
7. normalisation of LFC <br />
