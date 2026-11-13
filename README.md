# ndd_human_mouse_divergence
Human and mouse NDD genes divergence analysis

## Static explorer

The entry point is [sites/index.html](sites/index.html). It uses plain
HTML, CSS, JavaScript modules, and local data; there are no CDN dependencies
or application build steps. Serve it over HTTP rather than opening a file URL.
For example, from the repository root:

```powershell
python -m http.server 8000 --bind 127.0.0.1 --directory sites
```

Open `http://127.0.0.1:8000`. With GitHub Pages publishing the repository root,
the explorer is available at `<Pages URL>/sites/`. Alternatively, deploy the
contents of `sites` as the Pages artifact to make it the site root.

All application and data-preparation scripts live in `sites/scripts`.
All required datasets, including both phase-summary source CSVs and generated
trajectory JSON files, live in `sites/data`. Browser data requests resolve
relative to that data directory, including when hosted under a Pages subpath.
The interface uses an 18px base font and high-contrast text in both themes.
Charts scroll horizontally on narrow screens to keep their labels readable.
Blue accents and local inline SVG DNA, brain, and analysis icons distinguish
the report sections. Human trajectories use solid blue lines; mouse
trajectories use dashed neutral lines, preserving a non-color distinction.

### Ranking and details

- Both regions use the full NDD four-phase descriptive datasets, not the
  CHD6/CHD9-only or six-phase cerebellar analyses.
- Choose top N and most/least divergent globally, with independent cell-type
  filters for cortex and cerebellum. Highlighting indicates membership in both
  displayed lists, not statistical significance.
- Maximum mode ranks by the largest cell-type divergence score. Least divergent
  sorts that maximum ascending, not the minimum score across cell types.
  Exact ties are ordered alphabetically.
- Undefined correlations and divergence scores (for example, constant
  trajectories) remain missing, never zero. Summary statistics omit missing
  values as in the source analysis. Genes without a finite score in the
  selected scope are excluded from both ranking directions.
- Detail medians and means cover eligible cell types in maximum mode, or just
  the selected cell type in cell-type mode. The trajectory selector shows one
  cell type at a time and defaults to the maximum-scoring cell type.
- Charts use `log1p(median_cpm)` for each phase, or its species-specific
  z-standardization. Replicate counts are provided; no confidence intervals
  are inferred. Cortical mouse replicates are pools, not individual mice.
- Unequal cell-type coverage limits cross-region comparisons. These rankings
  are descriptive, not differential-expression significance tests.

### Data preparation and validation

Rankings are computed in the browser from each region's
`trajectory_similarity.csv`. Compact `trajectories.json` assets contain only
eligible gene/cell-type combinations, preserving the four phase-level median
CPM values and biological replicate counts from `aligned_phase_summary.csv`.
They are fetched on first opening a region's detail panel and cached for
subsequent selections. To regenerate them (Node.js 18+):

```powershell
node sites\scripts\build-trajectories.mjs
node --test tests\explorer.test.mjs
```

The generator reads both regions' phase-summary CSVs and similarity CSVs
exclusively from `sites/data`, then writes the generated JSON files beside
them. No external source folders are needed to regenerate or serve the site.

Tests compare every gene summary with published rankings and verify that
every chart's phase data reproduce the reported z-RMSE, direction agreement,
peak-phase distance, and composite divergence score.
