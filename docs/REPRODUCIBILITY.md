# Reproducibility and scope

This repository contains only the public human skin scRNA-seq reanalysis
associated with Figure 8. Its input is GEO GSE153760; no unpublished
experimental data, manuscript drafts, microscopy, flow-cytometry files, or
locally generated figure assets are included.

## Analysis stages

1. Load the per-sample 10x matrices and assign HC or AD group labels from sample
   names.
2. Apply the QC and normalization settings in
   config/analysis_parameters.yaml.
3. Generate a global embedding and Leiden clusters, then select immune cells
   and the marker-defined cDC2 subset.
4. Export sample-level summaries and descriptive cell-level distributions for
   the cDC2 fraction and CCR2, CCR7, and RAMP1 expression.
5. Compute the exploratory diffusion pseudotime and export its cell metadata
   and gene-trend table.

The runner records input filenames and sizes, parameters, software versions,
and execution logs in the chosen output directory. These files may include
local paths and are not intended for committing.

## Data provenance

GSE153760 is a public GEO dataset. This repository does not redistribute its
expression matrices. Users must download and prepare the input files from GEO.

## Manuscript alignment and verification boundary

The October 2, 2026 manuscript describes Mann–Whitney U tests with
Benjamini–Hochberg correction for Figure 8B and 8C. The current code's B-panel
significance labels use cell-level tests with adjustment across three genes;
sample-level expression tests are exported separately. The C-panel
significance label uses a sample-level test without adjustment. Consequently,
the code and manuscript statistical descriptions are not yet fully aligned.

The input matrices were not available in the environment used for this
documentation update. No complete run was performed, and the supplied final
panels were not numerically matched to regenerated outputs. Before citing a
release as the exact source of the published Figure 8, specify the experimental
unit and multiple-testing family in the manuscript, align the plotting code
and legend, then rerun the workflow and compare the outputs.
