# CRISPR Analysis via Locally calIBrated Reference (CALIBR)

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
git clone https://github.com/X/X_code.git
cd X/
```

Create and activate a virtual environment:
```
python -m venv venv
source venv/bin/activate      # On Windows: venv\Scripts\activate
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
3. Run: `python3 calibr.py --config my_config.yaml`

### Configuration parameters

| Parameter         | Description                                               | Example                               |
|-------------------|-----------------------------------------------------------|---------------------------------------|
| `input_counts`    | Path to the counts CSV                                    | `demo_data/WT_....csv`                |
| `output_dir`      | Directory where results are written                       | `output`                              |
| `rep`             | Number of replicates per condition                        | `4`                                   |
| `keep_counts`     | Number of gRNAs per gene                                  | `[1, 4, 10]`                          |
| `alf`             | Significance level (tail of the distribution)             | `0.06`                                |
| `target_type`     | Choice of target genes                                    | `GOI`                                 |
| `control_type`    | Choice of control genes to base Z correction on           | `Zero-expressed genes`                |
| `baseline`        | Baseline condition                                        | `T0`                                  |
| `cond`            | Treatment condition                                       | `5FU`                                 |
| `rep_pairs`       | Name of replicates (leave '' if no replicates present)    | `['60159','60160','60161','60162']`   |
| `st`              | Start of x axis                                           | `-10`                                 |
| `end`             | End of x axis                                             | `10`                                  |
| `step`            | Bin sizes for histograms                                  | `0.05`                                |
| `thr_lfchz`       | LFC threshold for volcano plot                            | `1.7`                                 |

### Outputs
| File Name                                                     | Description                                                                                                                       |
|---------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------|
| `distri_gene_types_{condition_vs_baseline}.png`               | Individual LFC distributions for all the different gene types (eg. genes of interests, zero expressed genes, non-targeting etc)   | 
| `distri_gene_types_overlaid_{condition_vs_baseline}.png`      | Overlaid LFC distributions for all the different gene types (eg. genes of interests, zero expressed genes, non-targeting etc)     |    
| `distri_{control_type}_{condition_vs_baseline}.png`           | Distribution of LFC, Z-corrected LFC and median Z-corrected LFC of chosen control genes                                           |    
| `gene_distribution_gRNA.png`                                  | Number of gRNAs attributed to each gene                                                                                           | 
| `gene_variability_boxplots.html`                              | Interactive box plots of variability between gRNAs for each gene                                                                  | 
| `gRNA_counts_normalisation.png`                               | Total raw (top) and normalised (bottom) read count summed over all gRNAs in each condition                                        | 
| `normalised_counts_mean_med.png`                              | Mean versus median of the normalised gRNA counts for each of condition                                                            | 
| `normalised_counts.csv`                                       | Normalised counts in csv format                                                                                                   | 
| `p_controls_{control_type}_{condition_vs_baseline}.png`       | P-value curve for control gRNAs                                                                                                   | 
| `p_distri_targ_cont_{condition_vs_baseline}.png`              | Control-calibrated p-value curves of target vs control gRNAs                                                                      | 
| `target_{condition_vs_baseline}_gRNA.csv`                     | LFC, Z-corrected LFC, q-value, condition table for target genes                                                                   | 
| `target_distri_LFC_Z_corrected_{condition_vs_baseline}.png`   | LFC and Z-corrected LFC gRNA distribution for target genes                                                                        | 
| `volcano_gRNA_{condition_vs_baseline}.png`                    | Volcano plot on gRNA and gene level                                                                                               | 
| `volcano_plot_interactive_{condition_vs_baseline}.html`       | Interactive volcano plot on gRNA and gene  CSV                                                                                    | 


## Pipeline Steps
### Module 1
Normalise raw count data as proportion per million <br />
Counts number of gRNA per gene/ intergenic <br />
Compute mean and median of normalised counts <br />

### Module 2
Compute and analyse LFC distributions of targets and controls (target, intergenic and zero expressed genes) <br />
Compute q-values, critical values, mean and standard deviation for zero expressed genes controls. <br />
Compute z for targets and all controls, based on mean and standard deviation of zero expressed genes and controls. <br />
Determine significant gRNA enrichment/depletion accounting for both FDR and FNR rates <br />
