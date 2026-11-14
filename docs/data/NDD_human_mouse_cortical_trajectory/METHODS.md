# Human-mouse cortical developmental trajectory analysis

## Scope

Velmeshev et al. human cortex is compared with wild-type Di Bella et al.
mouse neocortex across early, mid, and late prenatal development and the
neonatal period.

## Replication and inference

Human pseudobulks use individual donors. Mouse pseudobulks use pooled
libraries; most stages have one pool. Cross-species trajectory metrics are
therefore descriptive. Species-by-phase p-values and FDR values are not
reported, and cells are never treated as biological replicates.

For each gene and matched cell class, the standardized root-mean-square
trajectory error was calculated from the four species-specific z-standardized
log1p(CPM) values. The unscaled value is stored as `z_rmse_raw` (bounded by
0 and 2); `z_rmse` is `z_rmse_raw / 2`, bounded by 0 and 1. The composite
descriptive divergence score is

`(1 - Spearman rho)/2 + z_rmse + (1 - direction agreement) + peak-phase distance/3`

and ranges from 0 to 4. Inferential statistics are not implied by this score.

## Future brain-region integration

All primary tables include brain_region, analysis_scope, phase_set, and
integration_schema_version. region_integration_trajectory_metrics.csv is
the stable result-level input for a future cerebrum-cerebellum comparison.
