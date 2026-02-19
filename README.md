
# genome_annotation_quality_nextflow_pipeline

Nextflow pipeline for large-scale protein structure prediction,
intrinsic disorder analysis, and genome annotation quality assessment
using **Protenix**, **Metapredict**, and **PSAURON**.

------------------------------------------------------------------------

## Overview

This repository contains a **Nextflow DSL2 pipeline** for:

- Large-scale protein structure prediction using Protenix
- Intrinsic disorder prediction using Metapredict
- Protein-coding sequence quality assessment using PSAURON
- Automated chunking and HPC parallelisation
- Model processing and metric extraction
- Comparative plotting against reference proteome datasets

The pipeline is designed for HPC environments (SLURM supported)
and supports GPU acceleration where appropriate.

------------------------------------------------------------------------

## PSAURON Overview

Evaluating the accuracy of protein-coding sequences in genome annotations
is a challenging problem with no broadly applicable universal solution.

PSAURON (Protein Sequence Assessment Using a Reference ORF Network)
is a machine learning–based tool trained on a diverse dataset of over
1000 plant and animal genomes. It assigns a score to coding DNA or protein
sequences reflecting the likelihood that the sequence represents a genuine
protein-coding region.

PSAURON scores enable:

- Genome-wide protein annotation quality assessment
- Rapid identification of potentially spurious protein annotations
- Comparative overlay against large multi-species reference datasets

This pipeline integrates PSAURON scoring directly into structure and
disorder prediction workflows for multi-metric genome annotation QC.

------------------------------------------------------------------------

## Features

### Structure Prediction (Protenix)

- FASTA → Protenix JSON conversion
- JSON chunk splitting for parallel inference
- GPU-accelerated Protenix prediction
- Automatic model collection and aggregation
- pLDDT extraction into PKL datasets
- Distribution plotting vs reference organisms

### Disorder Prediction (Metapredict)

- GPU-accelerated disorder prediction
- Automatic overlay vs reference proteome disorder database
- Histogram and KDE density plots

### Protein-Coding Quality Assessment (PSAURON)

- GPU-accelerated PSAURON scoring
- Overlay of new dataset vs precomputed multi-species reference
- Histogram and KDE density plots
- Automated integration into Nextflow workflow

------------------------------------------------------------------------

## Workflow Diagram

                          FASTA
                   /       |        \
                  /        |         \
             PROTENIX  METAPREDICT      PSAURON
                |           |              |
             pLDDT PKL     CSV            CSV
                |           |              |
            PLOT_PLDDT  PLOT_METAPREDICT  PLOT_PSAURON
                |           |              |
               PNG         PNG            PNG

------------------------------------------------------------------------

## Requirements

### Core Software

- Nextflow ≥ 23
- Python 3.10
- CUDA-enabled GPU (recommended for Protenix, Metapredict, PSAURON)
- SLURM (optional, HPC environments)

------------------------------------------------------------------------

## Running the Pipeline

### Basic Run

```bash
nextflow run main.nf -profile slurm --fasta after_461.fasta --chunk_size 100
```

### Parameters

| Parameter     | Description                                |
|--------------|--------------------------------------------|
| --chunk_size | Number of sequences per Protenix chunk     |
| --fasta      | Path to input FASTA file                   |

------------------------------------------------------------------------

## Output Structure

    results/
     ├── plots/
     │    ├── plddt_hist_<dataset>.png
     │    ├── plddt_density_statsmodels_<dataset>.png
     │    └── plddt_density_scipy_<dataset>.png
     │
     ├── metapredict_plots/
     │    ├── <dataset>_mean_disorder_hist.png
     │    ├── <dataset>_mean_disorder_density_statsmodels.png
     │    └── <dataset>_mean_disorder_density_scipy.png
     │
     └── psauron_plots/
          ├── psauron_hist.png
          ├── psauron_density_statsmodels.png
          └── psauron_density_scipy.png

Intermediate outputs:

- protenix_out/
- *_all_predictions/
- plddt_all_values_<dataset>_all_one.pkl
- <dataset>_metapredict.csv
- <dataset>_psauron.csv

------------------------------------------------------------------------
## Example Output Images

### pLDDT Distribution Example

Place example images inside:

    docs/images/

Then reference them in README like:

![pLDDT KDE Example](plddt_density_scipy_after_461.png)

### Metapredict Disorder Example

![Metapredict KDE Example](mean_disorder_density_scipy.png)

### PSAURON Example

![Metapredict KDE Example](psauron_density_scipy.png)

------------------------------------------------------------------------
## Reference Datasets

### Disorder Reference

    bin/combined_proteome_disorder.pkl

### PSAURON Reference

    bin/combined_psauron_results.csv

Precomputed multi-species PSAURON score distribution used as background
for comparative overlay plots.

------------------------------------------------------------------------

## Pipeline Processes

| Process               | Description                                 |
|-----------------------|---------------------------------------------|
| FASTA_TO_JSON         | Converts FASTA → Protenix JSON              |
| SPLIT_JSON            | Splits JSON into chunk files                |
| PROTENIX_PREDICT      | Runs Protenix inference (GPU)               |
| COLLECT_CHUNKS        | Merges chunk prediction outputs             |
| PROCESS_MODELS        | Extracts pLDDT values into PKL              |
| PLOT_PLDDT            | Generates structure confidence plots        |
| METAPREDICT_DISORDER  | Runs disorder prediction (GPU)              |
| PLOT_METAPREDICT      | Generates disorder comparison plots         |
| PSAURON_RUN           | Runs PSAURON scoring (GPU)                  |
| PLOT_PSAURON          | Generates PSAURON overlay distribution plots|

------------------------------------------------------------------------

## Author

Shahram Mesdaghi
