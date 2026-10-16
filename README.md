# T. circumcincta Transcriptomic Analysis

Repository for R scripts and reproducible analysis for T. circumcincta transcriptomic heatmaps.

## Overview

This repository contains reproducible R scripts for generating publication-ready heatmaps from transcriptomic data for T. circumcincta across life stages and AMP-predicted peptide gene sets. The analyses are designed for consistent figure formatting and export.

## Repository contents

- `R/Heatmap_Protocol.R`  
  Pairwise log2FC heatmap comparing life-stage contrasts across AMP predictor scores.

- `R/Heatmap_Trajectory_Protocol.R`  
  Relative-expression trajectory heatmap for AMP-predicted genes across the developmental series.

## Requirements

- R (>= 4.0)
- R packages:
  - `readxl`
  - `dplyr`
  - `ComplexHeatmap`
  - `circlize`
  - `grid`

## Input data

The scripts expect a merged Excel file containing transcriptomic and AMP predictor columns, for example:

- `merged_data.xlsx`
- sheet: `Sheet1`
- gene identifier column: `gene_id`
- log2FC columns such as `log2FoldChange_L3_v_xL3`, `log2FoldChange_xL3_v_L4`, etc.

The input path is currently defined as a local Windows path in the scripts and should be adjusted to match your environment before running. This data can be found at 10.5281/zenodo.23415056

## Usage

From the repository root, open the script in RStudio and run it. Each script will:

1. Load the input Excel file.
2. Filter genes based on AMP predictor scores.
3. Build the heatmap matrix.
4. Export a publication-ready PDF.
5. Generate a PNG preview when enabled.

Example:

```r
source("R/Heatmap_Protocol.R")
source("R/Heatmap_Trajectory_Protocol.R")
```

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
