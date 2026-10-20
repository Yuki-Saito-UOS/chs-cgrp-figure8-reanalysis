# Reproducibility and data provenance

## Data source

The analysis uses the public human skin single-cell RNA-seq dataset GSE153760,
available through [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE153760).
Expression matrices are not included in this repository. Users should obtain
the data from GEO and prepare one 10x-format matrix set per sample before
running the workflow. Sample identifiers containing `_AD` are assigned to
the AD group; remaining samples are assigned to HC.

## Analysis configuration

Versioned parameters are maintained in
[`config/analysis_parameters.yaml`](../config/analysis_parameters.yaml). The
configuration records the dataset accession, quality-control thresholds,
embedding parameters, marker-based cell selection, pseudotime settings, and
statistical procedures. Review the configuration and confirm that the sample
identifiers match the intended group assignments before analysis.

## Execution and outputs

Run the complete workflow using the command in the repository README. The
runner writes the following to the selected output directory:

- `run_manifest.json`, containing the configuration, input file inventory,
  Python and package versions, and execution environment;
- `run.log`, containing the workflow log;
- `tables/`, containing cell-level expression data, sample-level summaries,
  test results, and pseudotime metadata and trends; and
- `figures/`, containing the generated analysis figures.

The manifest may include machine-specific paths. Keep run outputs and input
matrices outside version control unless they have been reviewed and approved
for distribution.

## Reproduction checks

For each run, retain the manifest and log with the output files. Confirm that
all expected samples were loaded, that group assignments are correct, and
that the generated tables and figures correspond to the selected configuration.
The analysis parameters and statistical units are described in
[`STATISTICAL_DECISIONS.md`](STATISTICAL_DECISIONS.md).

## Scope

This workflow reanalyses public observational data. Results describe patterns
in the source dataset and do not establish causal relationships. Cell
selection, clustering, and pseudotime ordering depend on the configured
markers and analysis parameters; alternative settings may produce different
results.
