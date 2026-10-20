# Statistical analysis

## Analysis units

The workflow produces cell-level expression distributions and summaries at
the sample level. Samples are the biological units for HC-versus-AD
comparisons of cDC2 abundance and sample-mean gene expression. Cell-level
expression tests are retained as a separate analysis and do not account for
dependence among cells from the same sample.

## cDC2 abundance

For each sample, the workflow calculates the fraction of PTPRC-positive
leukocytes assigned to the selected CLEC10A-positive cDC2 analysis set. HC and
AD groups are compared using a two-sided Mann–Whitney U test on sample-level
fractions. The resulting p-value is reported without multiple-testing
adjustment.

## Gene expression

CCR2, CCR7, and RAMP1 expression are summarized within the selected cDC2
population.

- Cell-level expression distributions are visualized by group. To limit the
  number of displayed points, up to 2,500 cells per gene and group are
  selected for plotting with a fixed random seed. The violin distributions
  and Mann–Whitney U tests use all eligible cells. Two-sided p-values are
  adjusted across the three genes using the Benjamini–Hochberg procedure;
  adjusted values determine the significance labels on the plots.
- For each sample, mean expression is calculated separately for each gene.
  HC and AD sample means are compared using two-sided Mann–Whitney U tests.
  The resulting p-values are adjusted across the three genes using the
  Benjamini–Hochberg procedure and exported with the sample-level summaries.

Cell-level tests do not account for within-sample dependence and should be
interpreted as exploratory. Sample-level comparisons use the sample as the
unit of analysis.

## Diffusion pseudotime

Diffusion pseudotime is used to describe transcriptional patterns across the
selected cDC2 cell population. The root is selected using the configured
CCR2-high, CCR7-low, and LAMP3-low score. Gene-expression trends and positive
cell fractions are summarized across quantile-based pseudotime bins and
smoothed for visualization. These analyses are descriptive and do not test a
temporal or lineage relationship.

## Multiple testing and reporting

Benjamini–Hochberg adjustment is applied separately to the three gene-level
tests for cell-level expression and to the three gene-level tests for
sample-level mean expression. The cDC2 abundance comparison is a single
unadjusted test. The run output records sample counts, test results, and
analysis parameters for reporting.
