# Statistical decisions

## Manuscript description (October 2, 2026)

The Figure 8 legend specifies a Mann–Whitney U test with Benjamini–Hochberg
correction for panels B and C. The main Statistical analysis section gives
Student's two-tailed t-test or one-way ANOVA followed by Tukey–Kramer unless
otherwise stated; the Figure 8 legend is the stated exception for these
panels. The legend does not specify the experimental unit or the family of
tests included in the correction.

## Experimental unit and current implementation

Cells from the same sample or donor are not independent biological
replicates. Sample-level summaries are therefore the appropriate unit for
HC-versus-AD inference unless a model accounts for donor clustering.

The current code has two different statistical paths for Figure 8B:

- The violin plot and its significance labels use cell-level two-sided
  Mann–Whitney U tests, followed by Benjamini–Hochberg adjustment across
  CCR2, CCR7, and RAMP1. The displayed cells are sampled for plotting, but the
  tests use all cells. These cell-level p-values do not account for donor
  clustering and should not be interpreted as donor-level inference.
- Per-sample mean expression and sample-level Mann–Whitney U p-values are also
  exported, but they are not used for the violin-plot significance labels.
- The cDC2 proportion is summarized and tested per sample using a two-sided
  Mann–Whitney U test. The current code does not adjust this p-value for
  multiple testing.

## Difference from the manuscript

The October 2, 2026 Figure 8 legend says panels B and C use Mann–Whitney U
tests with Benjamini–Hochberg correction. The current implementation instead
uses cell-level tests for the B-panel labels, sample-level unadjusted tests for
the C-panel label, and exports but does not plot sample-level expression
tests. Thus the code, plotted significance labels, and manuscript description
are not yet fully aligned. The correction family and experimental unit in the
manuscript also need to be stated explicitly to make the analysis auditable.

The input matrices were unavailable in the environment used for this update,
so group sample counts, p-values, and regenerated panels could not be compared
with the supplied manuscript figures. Do not treat the current repository
outputs as verified numerical reproductions of Figure 8.

## Interpretation

The human public-data reanalysis is observational. It can describe
associations in GSE153760 but cannot establish that the TRPV1–CGRP–RAMP1
mechanism causes atopic dermatitis in humans.
