# Scripts to support large-scale Boltz-1 protein complex predictions

A suite of tools for high-throughput protein structure prediction and analysis using Boltz-1, the open Alphafold3 alternative (https://github.com/jwohlwend/boltz). 

This repository contains utilities for generating Boltz-1 prediction instructions (MSAs and yaml prediction instruction files) of empiric n-mer complexes, starting from fasta files.

This software enables one-time MSA generation for each unique protein per complex (instead of repeatedly for each complex prediction), saving time and resources.


## Tools Included

1. **YAML Protein Pipeline** (`yaml_protein_pipeline.py`): parses FASTA files, generates MSAs, and generates all possible combinations for a specified number of elements (any n-mer - dimers, trimers, tetramers). 
Output are MSA and yaml-files with prediction instructions per complex, and a bash script to run Boltz-1.

2. **Complex Analysis Utility** (`complex_analysis.py`): analyzes protein complex prediction results, extracting metrics and generating visualizations.

Between the YAML generation and complex analysis, the generated run_predictions.sh script should be executed in a computing environment with Boltz-1 installed and the output file hierarchy available.

Note: upcoming improvements:
- Avoid null bytes in writing yaml files so repair function can be removed
- Find a better way to track protein identity of chains within complexes so no secondary mapping step is needed

## Installation

### Prerequisites

- Python 3.8 or higher
- Boltz can be installed in a Python 3.8 environment:
  ```bash
  pip install boltz -U
  ```
- Required Python libraries: pandas, numpy, biopython, plotly, pyyaml

Install the required dependencies:

```bash
pip install pandas numpy biopython plotly pyyaml
```

## YAML Protein Pipeline

### Features

- Parses FASTA files and extracts unique protein sequences
- Generates MSAs for each protein using YAML format
- Creates YAML files for all protein combinations
- Generates a batch script to run Boltz predictions

### Usage

```bash
python yaml_protein_pipeline.py --fasta_files file1.fasta file2.fasta [OPTIONS]
```

### Command Line Arguments

| Argument | Description |
|----------|-------------|
| `--fasta_files` | Path to one or more FASTA files containing protein sequences |
| `--r` | Number of elements in each combination (default: 2) |
| `--with_replacement` | Use combinations with replacement (default: without replacement) |
| `--skip_msa_generation` | Skip MSA generation (useful if MSAs already exist) |
| `--out_dir` | Output directory for all generated files (default: ./yaml_pipeline_output) |
| `--recycling_steps` | Number of recycling steps for Boltz predictions (default: 3) |
| `--diffusion_samples` | Number of diffusion samples for Boltz predictions (default: 1) |
| `--gpu_count` | Number of GPUs to use for Boltz predictions (default: 1) |

### Examples

#### Basic Pipeline with Defaults

```bash
python yaml_protein_pipeline.py --fasta_files proteins.fasta
```

#### Custom Pipeline for Trimers with 4 Recycling Steps

```bash
python yaml_protein_pipeline.py --fasta_files proteins.fasta --r 3 --recycling_steps 4 --out_dir ./trimer_predictions
```

#### Skip MSA Generation (If Already Generated)

```bash
python yaml_protein_pipeline.py --fasta_files proteins.fasta --skip_msa_generation --out_dir ./reuse_msas
```

## Complex Analysis Utility

### Features

- Works with protein complexes of any size (n-mers)
- Maps protein chains to gene names using FASTA sequences
- Extracts confidence scores, pTMs, and interface metrics
- Generates interactive scatter plots of complex quality metrics
- Handles diverse file structures and naming conventions

### Usage

```bash
python complex_analysis.py --fasta_files file1.fasta file2.fasta --predictions_dir /path/to/predictions [OPTIONS]
```

### Command Line Arguments

| Argument | Description |
|----------|-------------|
| `--fasta_files` | Path to one or more FASTA files containing protein sequences |
| `--predictions_dir` | Directory containing the complex predictions |
| `--complex_yamls_dir` | (Optional) Directory containing the complex YAML files |
| `--output_file` | (Optional) Output CSV file for results (default: complex_analysis_results.csv) |
| `--plot_results` | (Optional) Generate an interactive visualization of the results |

### Examples

#### Basic Analysis

```bash
python complex_analysis.py --fasta_files proteins.fasta --predictions_dir ./predictions
```

#### Analysis with Custom YAML Directory and Plotting

```bash
python complex_analysis.py --fasta_files proteins.fasta --predictions_dir ./predictions --complex_yamls_dir ./yaml_files --plot_results
```

## Complete Workflow Example

This example demonstrates a complete workflow from protein sequence to analysis:

```bash
# Step 1: Generate prediction pipeline for protein dimers
python yaml_protein_pipeline.py --fasta_files proteins.fasta --r 2 --out_dir ./dimer_pipeline

# Step 2: Run the generated batch script to perform predictions
bash ./dimer_pipeline/run_predictions.sh

# Step 3: Analyze the prediction results
python complex_analysis.py --fasta_files proteins.fasta --predictions_dir ./dimer_pipeline/predictions --complex_yamls_dir ./dimer_pipeline/yaml_files --plot_results
```

## Output

### Pipeline Output

- MSAs for each protein
- YAML files for protein combinations
- Batch script for running predictions

### Analysis Output

- CSV file with comprehensive analysis results
- Interactive HTML visualization (when `--plot_results` is used)
- Summary statistics in the terminal output
- Detailed log file

## License

[MIT License](LICENSE)