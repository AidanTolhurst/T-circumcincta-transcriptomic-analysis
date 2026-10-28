# T. circumcincta Transcriptomic Analysis

Repository for R scripts and reproducible analysis for T. circumcincta transcriptomic heatmaps.

## Overview

This repository contains reproducible R scripts for generating publication-ready heatmaps from transcriptomic data for T. circumcincta across life stages and AMP-predicted peptide gene sets. The analyses are designed for consistent figure formatting and export.

## Repository contents

- `PairwiseComparison.R`  
  Pairwise log2FC heatmap comparing life-stage contrasts across AMP predictor scores.

- `TrajectoryComparison.R`  
  Relative-expression trajectory heatmap for AMP-predicted genes across the developmental series.

## Requirements

- R (>= 4.0)
- R packages:
  - `readxl`
  - `dplyr`
  - `ComplexHeatmap`
  - `circlize`
  - `grid`

## Data
Archived on Zenodo: https://doi.org/10.5281/zenodo.23746703
- File: `T-Circumcincta-Transcriptomic-Dataset.xlsx`
- MD5: 3154a00cebc4ccfda21a06fc3a3489db
- Download it and place it in the folder you run the R script from.

## Usage
From the repository root, open the script in RStudio and run it. Each script will:

1. Load the input Excel file.
2. Filter genes based on AMP predictor scores.
3. Build the heatmap matrix.
4. Export a publication-ready PDF.
5. Generate a PNG preview when enabled.


## Output
The scripts generate:

- PDF output for publication figures
- optional PNG preview files
- manuscript-ready figure dimensions and typography settings

## Citation
If you use these scripts in a publication, please cite the repository and the corresponding manuscript when available.

Suggested citation:

Blank for now

## License
This repository is provided for reproducible scientific analysis and may be used and adapted for academic and research purposes. Please contact the repository owner for any reuse beyond standard academic use.
