# Human skin single-cell RNA-seq analysis

This repository provides a reproducible workflow for reanalysing the public
human skin single-cell RNA-seq dataset GSE153760. The analysis focuses on
immune cells and CLEC10A-positive conventional dendritic cell 2 (cDC2)
populations in healthy controls (HC) and patients with atopic dermatitis (AD).

## Analysis overview

The workflow loads per-sample 10x Genomics matrices, performs quality control
and dimensionality reduction, identifies immune and cDC2 analysis sets, and
generates summaries for cDC2 abundance, gene expression, and diffusion
pseudotime. Analysis parameters are specified in
[`config/analysis_parameters.yaml`](config/analysis_parameters.yaml).

The analysis includes:

- quality control, normalization, highly variable gene selection, principal
  component analysis, UMAP, and Leiden clustering;
- immune-cell identification by PTPRC expression and cDC2 selection using
  canonical marker genes and the configured Leiden clusters;
- HC-versus-AD summaries of CLEC10A-positive cDC2 abundance and expression of
  CCR2, CCR7, and RAMP1; and
- exploratory diffusion pseudotime analysis of the selected cDC2 population.

## Data

GSE153760 is available from [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153760).
The expression matrices are not redistributed in this repository. Prepare one
set of 10x-format `barcodes.tsv.gz`, `features.tsv.gz`, and `matrix.mtx.gz`
files per sample. Sample names containing `_AD` are assigned to the AD group;
all other sample names are assigned to HC.

## Installation

Create the analysis environment with Conda:

```bash
conda env create -f environment.yml
conda activate chs-cgrp-figure8
```

## Run the analysis

Provide absolute paths to the input and output directories:

```bash
python scripts/run_figure8.py \
  --data-dir /path/to/GSE153760_10x_matrices \
  --output-dir /path/to/figure8-results
```

The output directory contains a run manifest and log, tables of cell-level and
sample-level summaries, and generated figures. Input matrices and generated
outputs are excluded from version control.

## Statistical analysis and interpretation

The analysis reports both cell-level expression distributions and
sample-level summaries. Statistical tests and their units are documented in
[`docs/STATISTICAL_DECISIONS.md`](docs/STATISTICAL_DECISIONS.md). The workflow
is observational: it characterizes associations in GSE153760 and does not
establish causal mechanisms in human atopic dermatitis. Diffusion pseudotime
provides a computational ordering of cells and should not be interpreted as a
direct measurement of cell lineage or temporal progression.

Further information on data provenance and reproducibility is available in
[`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md). Citation metadata are
provided in [`CITATION.cff`](CITATION.cff).
