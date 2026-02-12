# genome_annotation_quality_nextflow_pipeline
Nextflow pipeline for large-scale protein structure prediction using Protenix, with automated chunking, model processing, and pLDDT distribution analysis against reference organism datasets. Also Metapredict and PSAURON analsis
# Nextflow Protenix pLDDT Pipeline

## Overview

This repository contains a **Nextflow DSL2 pipeline** for running
**Protenix structure prediction** from FASTA input, processing resulting
models, and generating **pLDDT distribution plots** combined with
reference CSV organism data.

The pipeline is designed for HPC environments (SLURM supported) and
supports chunked inference for large FASTA datasets.

------------------------------------------------------------------------

## Features

-   FASTA → Protenix JSON conversion\
-   JSON chunk splitting for parallel inference\
-   Protenix prediction using ESM embeddings\
-   Automatic model collection and aggregation\
-   pLDDT extraction and aggregation into PKL\
-   Plot generation:
    -   Histogram distributions
    -   KDE density plots (StatsModels + SciPy)
-   Overlay of reference organism CSV dataset

------------------------------------------------------------------------

## Workflow Diagram

    FASTA
      |
      v
    FASTA_TO_JSON
      |
      v
    SPLIT_JSON
      |
      v
    PROTENIX_PREDICT (per chunk, parallel)
      |
      v
    COLLECT_CHUNKS
      |
      v
    PROCESS_MODELS
      |
      v
    PLOT_PLDDT
      |
      v
    Plots (PNG)

------------------------------------------------------------------------

## Requirements

### Software

-   Nextflow ≥ 23
-   Python 3.10
-   CUDA-enabled GPU (recommended)
-   SLURM (optional, for HPC use)

### Python Environment (Protenix)

-   protenix
-   torch
-   esm
-   pandas
-   numpy
-   scipy
-   statsmodels
-   matplotlib

------------------------------------------------------------------------

## Environment Setup (Example)

``` bash
module load python/3.10
source ~/miniconda3/etc/profile.d/conda.sh
conda activate protenix_env

export PROTENIX_CACHE=/home/<USER>/protenix_cache
```

------------------------------------------------------------------------

## Running the Pipeline

### Basic Run

``` bash
nextflow run main.nf -profile slurm   --fasta input.fasta   --chunk_size 100
```

### Parameters

  Parameter        Description
  ---------------- ----------------------------------------
  `--fasta`        Input FASTA file
  `--chunk_size`   Number of sequences per Protenix chunk

------------------------------------------------------------------------

## Output Structure

    results/
     └── plots/
          ├── plddt_hist_<dataset>.png
          ├── plddt_density_statsmodels_<dataset>.png
          └── plddt_density_scipy_<dataset>.png

Intermediate outputs: - `protenix_out/` - `*_all_predictions/` -
`plddt_all_values_<dataset>.pkl`

------------------------------------------------------------------------

## Example Plots

### Histogram Example

(Add screenshot here once available)

    plddt_density_scipy_after_461.png

### KDE Density Example

(Add screenshot here once available)

    plddt_density_scipy_after_461.png

------------------------------------------------------------------------

## CSV Reference Dataset

The plotting step automatically loads:

    plddt_model_organisms.csv

Filtered species: - Mus musculus\
- Drosophila melanogaster\
- Arabidopsis thaliana\
- Saccharomyces cerevisiae\
- Rattus norvegicus\
- Homo sapiens\
- Candida auris (CaurisB8441)\
- Pan troglodytes

------------------------------------------------------------------------

## Pipeline Processes

  Process            Description
  ------------------ ---------------------------------
  FASTA_TO_JSON      Converts FASTA → Protenix JSON
  SPLIT_JSON         Splits JSON into chunk files
  PROTENIX_PREDICT   Runs Protenix inference
  COLLECT_CHUNKS     Merges chunk prediction outputs
  PROCESS_MODELS     Extracts pLDDT values into PKL
  PLOT_PLDDT         Generates plots from PKL + CSV

------------------------------------------------------------------------

## Troubleshooting

### Pipeline finishes too quickly

Check: - Conda environment activated - PROTENIX_CACHE set - CUDA modules
loaded

### Missing PKL files

Check PROCESS_MODELS logs in Nextflow work directory.

### Missing plots

Ensure CSV file is present in project root.

------------------------------------------------------------------------

## Citation

If you use this pipeline, please cite: - Protenix - ESM models -
Nextflow

------------------------------------------------------------------------

## Author

`<Your Name>`{=html}

------------------------------------------------------------------------
