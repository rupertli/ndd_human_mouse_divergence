# Enriched human-mouse cortical developmental trajectories

Four dataset adapters were combined at donor-level pseudobulk, not cell-level integration.
Inputs: velmeshev_2023_adapter.h5ad, di_bella_2021_adapter.h5ad, wang_2025_adapter.h5ad, gao_2025_adapter.h5ad
Phase schema: event; phase order: event_bin_1, event_bin_2, event_bin_3, event_bin_4, event_bin_5
Minimum shared phases: 4
Minimum cells per donor/pool-stage-class pseudobulk: 20

The event schema is an approved exploratory Workman 2013 baseline model.
It does not apply cortical-neurogenesis interaction corrections; bins are analyst choices.
Out-of-range ages are retained in expression tables but excluded from trajectory comparisons.
Original postconception ages are used, not expression-derived synchronized ages.
Original legacy phases are retained for a separate sensitivity run.

Wang donors matching NIH or exact donor IDs in the existing human atlas are excluded
from the combined input; remaining unidentified overlap cannot be ruled out.
Ambiguous multipotent progenitors in the new studies are excluded.
Counts are summed within dataset/donor/age/class; regions within a donor are pooled
and contributing region/study labels remain in the pseudobulk table.
No batch-corrected values are used. Metrics are descriptive, not disease-effect tests.
Regional mixtures, study effects, neuronal enrichment and cell/nucleus assay differences
can influence trajectories; these outputs are not a region-matched confirmatory analysis.
Di Bella pools are not independent animals. No inferential p-values are reported.
Genes missing in any input are excluded globally under the skip policy.
Mouse symbol matching retains the existing capitalization convention, not verified orthology.

## Website score normalization

These website copies retain raw standardized trajectory RMSE in `z_rmse_raw`
and store `z_rmse = z_rmse_raw / 2`. Composite scores and rankings are recomputed
as `(1 - rho)/2 + z_rmse + (1 - direction_agreement) + peak_phase_distance/4`.
Peak distance uses positions in the full five-bin order, even for four-bin
shared trajectories. Undefined correlations retain undefined composite scores.
Website rankings differ from the raw analysis outputs; original analysis files
are not changed. The four-phase sensitivity copies instead divide peak distance by 3.
