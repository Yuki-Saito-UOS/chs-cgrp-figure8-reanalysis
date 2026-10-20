# Figure 8 public-data reanalysis

This repository contains the code and analysis settings for the human skin
single-cell RNA-seq reanalysis associated with Figure 8 of a CHS-CGRP study.
It is restricted to the public GEO dataset GSE153760. It does not contain
experimental source data, manuscript files, or generated figures and tables.

## Scope

- Load the per-sample 10x matrices from GSE153760.
- Apply QC, normalization, HVG selection, PCA, UMAP, and Leiden clustering.
- Select immune cells and a marker-defined CLEC10A-positive cDC2 analysis set.
- Summarize sample-level cDC2 proportions and expression of CCR2, CCR7, and
  RAMP1; cell-level distributions are descriptive.
- Calculate and plot the exploratory diffusion pseudotime analysis.

The parameter file is config/analysis_parameters.yaml. The reasoning and known
limits are described in docs/REPRODUCIBILITY.md and
docs/STATISTICAL_DECISIONS.md.

## Data

Download the public dataset from
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153760 and prepare one
matching set of *_barcodes.tsv.gz, *_features.tsv.gz, and *_matrix.mtx.gz files
per sample. Data are not redistributed here. Samples whose names contain _AD
are assigned to the AD group; other samples are assigned to HC, as configured
in the runner.

## Run

Create the environment:

    conda env create -f environment.yml
    conda activate chs-cgrp-figure8

Run the analysis with absolute input and output paths:

    python scripts/run_figure8.py --data-dir /path/to/GSE153760_10x_matrices --output-dir /path/to/figure8-results

The runner writes a manifest, log, tables, and PDF figures to the output
directory. These generated files and input data are excluded from Git.

## Interpretation and release status

The public-data analysis is associative and cannot establish causality in
human atopic dermatitis. The cluster selection and pseudotime root are
data-dependent analytical choices and should be rechecked if inputs or upstream
settings change.

Before citing a deposited version, reconcile the statistical methods in the
manuscript with the primary sample-level tests described in the statistical
decisions document and verify the complete run against the final input files.
