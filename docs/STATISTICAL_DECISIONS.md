# Statistical decisions

## Experimental unit

Cells from the same sample or donor are not independent biological replicates.
Primary HC-versus-AD inference therefore uses one summary per sample (or a
model that accounts for donor clustering). Cell-level tests and distributions
are descriptive only.

## Current implementation

- cDC2 fraction: sample-level comparison using a two-sided Mann–Whitney U test.
- Expression of CCR2, CCR7, and RAMP1: sample-level mean expression is exported;
  the expression analysis also reports cell-level descriptive distributions.
- Multiple testing: Benjamini–Hochberg adjustment across the three expression
  genes is specified by the analysis plan.

## Manuscript reconciliation required

The manuscript text and Figure 8 legend have previously described different
statistical methods. Before citing this code as the final analysis, harmonize
the Methods, figure legend, exported statistics, and this document. Report the
sample counts per group, the primary test, the multiple-testing family, and
whether plotted points represent cells or sample-level summaries.

## Interpretation

The human public-data reanalysis is observational. It can describe associations
in GSE153760 but cannot establish that the TRPV1-CGRP-RAMP1 mechanism causes
atopic dermatitis in humans.
