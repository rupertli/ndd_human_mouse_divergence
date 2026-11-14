# Human-mouse cerebellar four-phase descriptive trajectory analysis

## Purpose

This is a new, non-overwriting recalculation of the existing Cerebellum
analysis using the same four-phase descriptive framework as the Cerebrum
analysis.

## Shared phases

- early_prenatal
- mid_prenatal
- late_prenatal
- neonatal

Juvenile and adult observations remain available in the original analysis
but are intentionally excluded here.

## Metrics

Spearman rho, z-RMSE, adjacent-phase direction agreement, peak-phase
distance, and the identical Cerebrum descriptive divergence formula are
reported. Inferential p-values, FDR values, interaction statistics, and
bootstrap confidence intervals are intentionally blank so the two brain
regions can be compared under one descriptive contract.

The raw standardized trajectory RMSE is stored as `z_rmse_raw` (0–2).
`z_rmse` is normalized as `z_rmse_raw / 2` (0–1) and used in the composite
descriptive divergence score:

`(1 - Spearman rho)/2 + z_rmse + (1 - direction agreement) + peak-phase distance/3`

The resulting score ranges from 0 to 4 and is descriptive, not inferential.

## Replicates

Human pseudobulks use donors; mouse pseudobulks use individual mice.

## Integration

region_integration_trajectory_metrics.csv uses integration schema 1.0 and
is intended for result-level Cerebrum-Cerebellum integration.
