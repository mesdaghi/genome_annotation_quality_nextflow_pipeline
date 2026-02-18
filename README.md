# genome_annotation_quality_nextflow_pipeline

Nextflow pipeline for large-scale protein structure prediction and
disorder analysis using **Protenix**, **Metapredict**, and optional
**PSAURON** downstream analysis.

------------------------------------------------------------------------

## Overview

This repository contains a **Nextflow DSL2 pipeline** for:

-   Large-scale **protein structure prediction** using Protenix
-   **Intrinsic disorder prediction** using Metapredict
-   Automated chunking and HPC parallelisation
-   Model processing and metric extraction
-   Comparative plotting against **reference proteome datasets**

The pipeline is designed for **HPC environments (SLURM supported)** and
supports chunked inference for large FASTA datasets.

------------------------------------------------------------------------

## Features

### Structure Prediction (Protenix)

-   FASTA → Protenix JSON conversion
-   JSON chunk splitting for parallel inference
-   Protenix prediction using ESM embeddings
-   Automatic model collection and aggregation
-   pLDDT extraction into PKL datasets
-   Distribution plotting vs reference organisms

### Disorder Prediction (Metapredict)

-   GPU-accelerated Metapredict disorder prediction
-   Automatic overlay vs reference proteome disorder database
-   Histogram and KDE density plots

### Downstream Analysis

-   Compatible with PSAURON workflows
-   Designed for multi-metric genome annotation QC workflows

------------------------------------------------------------------------

## Workflow Diagram

                    FASTA
                   /     \
                  /       \
             PROTENIX     METAPREDICT
                |              |
             pLDDT PKL        CSV
                |              |
            PLOT_PLDDT     PLOT_METAPREDICT
                |              |
               PNG            PNG

------------------------------------------------------------------------

## Requirements

### Core Software

-   Nextflow ≥ 23
-   Python 3.10
-   CUDA-enabled GPU (recommended for Metapredict)
-   SLURM (optional, HPC environments)

------------------------------------------------------------------------

## Python Environments

### Protenix Environment

-   protenix
-   torch
-   esm
-   pandas
-   numpy
-   scipy
-   statsmodels
-   matplotlib

### Metapredict Environment

-   metapredict
-   numpy
-   pandas
-   matplotlib
-   scipy
-   statsmodels

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
nextflow run main.nf -profile slurm --chunk_size 100
```

### Parameters

  Parameter      Description
  -------------- ----------------------------------------
  --chunk_size   Number of sequences per Protenix chunk

------------------------------------------------------------------------

## Output Structure

    results/
     ├── plots/
     │    ├── plddt_hist_<dataset>.png
     │    ├── plddt_density_statsmodels_<dataset>.png
     │    └── plddt_density_scipy_<dataset>.png
     │
     └── metapredict_plots/
          ├── <dataset>_mean_disorder_hist.png
          ├── <dataset>_mean_disorder_density_statsmodels.png
          └── <dataset>_mean_disorder_density_scipy.png

Intermediate outputs:

-   protenix_out/
-   \*\_all_predictions/
-   plddt_all_values\_`<dataset>`{=html}.pkl
-   `<dataset>`{=html}\_metapredict.csv

------------------------------------------------------------------------

## Example Output Images

### pLDDT Distribution Example



![pLDDT KDE Example](plddt_density_scipy_after_461.png)

### Metapredict Disorder Example

![Metapredict KDE Example](mean_disorder_density_scipy.png)

------------------------------------------------------------------------

## Reference Datasets

### pLDDT Reference

Used for structure confidence comparison.

### Disorder Reference

    bin/combined_proteome_disorder.pkl

Used for overlay comparison of predicted disorder vs reference
proteomes.

------------------------------------------------------------------------

## Pipeline Processes

  Process                Description
  ---------------------- --------------------------------------
  FASTA_TO_JSON          Converts FASTA → Protenix JSON
  SPLIT_JSON             Splits JSON into chunk files
  PROTENIX_PREDICT       Runs Protenix inference
  COLLECT_CHUNKS         Merges chunk prediction outputs
  PROCESS_MODELS         Extracts pLDDT values into PKL
  PLOT_PLDDT             Generates structure confidence plots
  METAPREDICT_DISORDER   Runs disorder prediction
  PLOT_METAPREDICT       Generates disorder comparison plots

------------------------------------------------------------------------

## Troubleshooting

### Pipeline runs slowly

Check: - SLURM queue load - GPU availability - Chunk size vs cluster
capacity

### Missing PKL files

Check PROCESS_MODELS logs in Nextflow work directory.

### Missing plots

Ensure reference PKL exists:

    bin/combined_proteome_disorder.pkl

------------------------------------------------------------------------

## Citation

If you use this pipeline, please cite:

-   Protenix
-   ESM protein language models
-   Metapredict
-   Nextflow

------------------------------------------------------------------------

## Author

**Shahram Mesdaghi**
