# ndd_human_mouse_divergence
Human and mouse NDD genes divergence analysis

## Static explorer

The entry point is [docs/index.html](docs/index.html). It uses plain
HTML, CSS, JavaScript modules, and local data; there are no CDN dependencies
or application build steps. Serve it over HTTP rather than opening a file URL.
For example, from the repository root:

```powershell
python -m http.server 8000 --bind 127.0.0.1 --directory docs
```

Open `http://127.0.0.1:8000`. To publish with GitHub Pages from a branch,
select the `docs` folder as the publishing source; the explorer is then
available at `<Pages URL>/`. If publishing the repository root instead,
use `<Pages URL>/docs/`. Alternatively, deploy the contents of `docs` as the
Pages artifact to make it the site root.

All application and data-preparation scripts live in `docs/scripts`.
All required datasets, including phase-summary source CSVs and generated
trajectory JSON files, live in `docs/data`. Browser data requests resolve
relative to that data directory, including when hosted under a Pages subpath.
Data requests revalidate cached responses with the server so updated CSV schemas
are loaded after data regeneration. If an already-open preview reports a missing
`z_rmse_raw` column, hard-refresh it (Ctrl+F5) to reload the updated application.
The interface uses an 18px base font and high-contrast text in both themes.
Charts scroll horizontally on narrow screens to keep their labels readable.
Blue accents and local inline SVG DNA, brain, and analysis icons distinguish
the report sections. Human trajectories use solid blue lines; mouse
trajectories use dashed neutral lines, preserving a non-color distinction.

### Ranking and details

- Cortex defaults to `NDD_human_mouse_cortical_enriched`: 19,170 gene/class
  comparisons, 3,303 genes, and six lineage-based classes from Velmeshev,
  Di Bella, Wang, and Gao. Cerebellum retains its four-phase full-NDD dataset.
  Original cortical data and enriched legacy-alignment sensitivity outputs
  are preserved alongside the primary results.
- Cortex uses five analyst-selected Workman 2013 baseline event-score bins
  (edges 0.20, 0.35, 0.50, 0.65, 0.80, 1.00), not published biological phases
  or a cortical-neurogenesis-corrected alignment. Astroglial and oligodendrocyte
  lineage comparisons use four shared bins; other classes use five.
  The cerebellar age framework remains different.
- The interface calls harmonized annotation groups "cell classes"; these may
  contain multiple cell types and developmental states. Existing data fields
  and internal identifiers remain unchanged for compatibility.
- Choose top N and most/least divergent globally, with independent cell-class
  filters for cortex and cerebellum. Highlighting indicates membership in both
  displayed lists, not statistical significance.
- Maximum mode ranks by the largest cell-class divergence score. Least divergent
  sorts that maximum ascending, not the minimum score across cell classes.
  Exact ties are ordered alphabetically.
- Undefined correlations and divergence scores (for example, constant
  trajectories) remain missing, never zero. Summary statistics omit missing
  values as in the source analysis. Genes without a finite score in the
  selected scope are excluded from both ranking directions.
- The raw standardized trajectory RMSE is retained as `z_rmse_raw` (0–2);
  `z_rmse` is normalized as `z_rmse_raw / 2` (0–1) before computing the
  composite divergence score, whose range is 0–4. Peak distance is divided
  by the full phase-order span: 4 for enriched cortex and 3 for four-phase
  data, even when a cortical trajectory uses only four shared bins.
  Website scores and rankings therefore differ from raw scientific outputs;
  source analysis files outside this repository are not modified.
- Detail medians and means cover eligible cell classes in maximum mode, or just
  the selected cell class in cell-class mode. The trajectory selector shows one
  cell class at a time and defaults to the maximum-scoring cell class.
- Clicking a gene in either table selects that same gene for both detail panels,
  including when it is outside the other region's displayed top N. Selection
  persists when rankings or filters change. Each panel respects its region's
  cell-class filter and explicitly reports unavailable data rather than showing
  a different gene. If all scores are undefined, the chart defaults to the first
  eligible cell class and summary scores remain not estimable.
- The shared gene-entry field sits between the ranking tables and detail panels.
  Enter a gene symbol from either dataset to inspect it outside the displayed
  top N, or click a table gene to fill the same field. Lookup trims whitespace
  and ignores case while preserving the dataset's canonical spelling. Unknown
  symbols are reported explicitly without replacing the current details.
- Charts use `log1p(median_cpm)` for each phase, or its species-specific
  z-standardization. Replicate counts are provided; no confidence intervals
  are inferred. Each chart uses only the scored shared phases; absent phases
  are not imputed. Cortical mouse replicates mix Gao donors and Di Bella pools.
- Mixed cortical regions, studies, assays, potential unidentified donor overlap,
  and symbol-based rather than validated one-to-one orthology limit inference.
  Known overlapping Wang donors and ambiguous new progenitors were excluded.
  Unequal cell-class coverage and different age frameworks limit cross-region
  comparisons. These rankings are descriptive, not significance or disease-effect tests.

### NDD bulk-analysis data

`docs/data/ASD_Bulk_analysis` (legacy folder name) contains the NDD expression/K-means v2 results
copied from `BulkRNASeq`, preserving `cluster_result` and
`enrichment/GO`, `enrichment/HPO`, and `enrichment/KEGG`.
The updated source files include per-cluster gene lists, timepoint trends,
and enrichment tables for 2,861 genes in 13 clusters.
GO and KEGG CSVs cover all 13 clusters. Only clusters 4 and 11 have supplied
NDD HPO CSVs; other HPO sources are explicitly marked unavailable.
Header-only exports are displayed as empty results, not loading errors.
The separate [bulk-analysis page](docs/bulk-analysis.html), linked from the
developmental explorer, identifies the selected gene's cluster, lists all members
(including the highlighted selected gene), draws the cluster-average expression
trend, then shows HPO, GO, and KEGG results for the whole cluster.
Enrichment is no longer filtered to rows containing the selected gene.
Member buttons select another gene; `?gene=ACTB` supports direct links.
The chart uses the exported `Mean_Expression`, `SE`, and `N`, with mean +/- SE
shading and a values table. The six developmental phases are ordered categories,
not proportional elapsed time or the five cortical event bins. Source values,
units, and p-values are not transformed or recalculated.
The page shares the developmental explorer's typography, blue accents, cards,
and theme variables, with local SVG DNA/analysis icons and section navigation.
No p-values, enrichment thresholds, or cluster assignments are recalculated.
Selecting a gene in either developmental ranking also shows a gene-specific
bulk-analysis link above "How to interpret these rankings" when it has a
copied bulk cluster assignment. The link preserves the selected gene and theme.
Genes outside the bulk cluster dataset do not receive this link; lookup loading
errors are reported separately without disabling developmental rankings.
The same selected-gene navigation card, above "How to interpret these rankings",
is labeled "BulkRNASeq analysis for selected gene" and displays links vertically
with extra spacing. It
also links to BrainSpan profiles and networks with the selected gene, theme, and
batched-page version preserved. This link remains available for genes without
bulk cluster data and when bulk lookup fails. BrainSpan reports unavailable
genes or profiles explicitly.

`bulk-data.json` indexes the supplied NDD gene lists, trend CSV, and enrichment
paths, recording missing HPO sources as null. HPO `_communities.csv` exports
are preferred when supplied; otherwise the original enrichment CSV is used.
Regenerate the index after replacing or adding source files.
GO/KEGG overlap display uses the bundled `gene-entrez-map.json`,
generated locally from `org.Hs.eg.db` without uploading gene lists.
Unmapped symbols remain selectable; unmapped Entrez IDs are explicitly labeled
as IDs. Serving the page needs no R runtime or external annotation service.
To regenerate the index and mapping, install `AnnotationDbi`, `org.Hs.eg.db`,
and `jsonlite`, then run:

```powershell
node docs\scripts\build-bulk-data.mjs
Rscript docs\scripts\build-bulk-gene-map.R
node --test tests\bulk.test.mjs
```

#### Orphanet HPO onset context

The HPO table retains associated-disease onset counts and possible five-bin
cerebral interval overlaps for 31 terms bundled from the earlier ASD export.
NDD terms outside that evidence set are explicitly labeled as lacking bundled
context; no onset is inferred. The evidence generator now targets the supplied
NDD HPO CSVs, so both evidence and summary can be regenerated for the new terms.
This is **indirect disease-onset context, not phenotype onset or a most-likely
phase**. Exact HPO associations are joined to natural history by ORPHAcode;
zero-frequency exclusions are omitted. Every unique clinical entity contributes
once per reported category, so category totals can exceed the number of
entities. Related subtypes/groups are retained and are not independent patient
observations. No restriction to ACTB-related diseases is applied.

The matching Orphadata XML exports are dated 23 June 2026. Source URLs,
retrieval date, MD5 checksums, per-disease onset categories and evidence links
are bundled in `orphanet-hpo-evidence.json`; `orphanet-hpo-context.json`
contains the summary and bin alignment. Attribution: Orphanet / INSERM,
Orphadata Science, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
Definitions follow [Orphadata product description, Table 8, page 32](https://sciences.orphadata.com/docs/OrphadataFreeAccessProductsDescription.pdf).

For overlap calculations, birth is assumed at 266 postconception days and a
year is 365.25 days. Completed-year ranges end at the next birthday; infancy
ends at 24 months. "All ages" is conservatively open-ended, with no peak.
Half-open intervals overlap the exact cortical boundaries in
`alignment_boundaries.csv`: antenatal spans bins 1-3 and earlier ages; neonatal
spans bins 3-4 (only 0.40 day in bin 3 under this assumption); infancy spans
bins 4-5 and later ages; childhood and older categories are beyond bin 5.
Unknown onset stays unknown. This is not a validated clinical/cortical
alignment; gestational age variation can change boundary overlap.

Regenerate with R (`xml2`, `jsonlite`) and Node.js; only public Orphadata
files are downloaded, and no local gene lists are uploaded:

```powershell
Rscript docs\scripts\build-orphanet-evidence.R
node docs\scripts\build-orphanet-context.mjs
node --test tests\bulk.test.mjs
```

An optional argument to the R script selects an existing XML cache directory.
Without one, raw downloads are temporary and removed after generation.

### BrainSpan NDD page

[docs/brainspan.html](docs/brainspan.html) migrates the BrainSpan Shiny interface
to a separate static page, linked from the home page. It shares `styles.css`,
theme detection, typography, cards, and accents. No R process, CDN, build step,
or dynamic API is needed when serving it on GitHub Pages.

`docs/data/BrainSpan_NDD` contains a 3,921-gene SysNDD/SFARI/DBD catalog and
provenance, 524 samples, driver mappings, and 72 batched profile/network CSVs
under `genes`, with up to 100 available genes per batch. The initial
export includes 3,597 expression profiles, 5,328 activity profiles, and 56,740
network edges. Expression and activity values are copied without normalization
and joined to metadata by sample ID. Missing values are blank, not zero.
Both network endpoint symbols must be in the NDD/ASD union; non-NDD targets
from the original Shiny app are intentionally excluded.

The browser loads the catalog and metadata once, then fetches the selected
gene's profile/network batch and caches it for the session. Selecting another
gene in the same batch requires no additional fetch. Only the selected gene's
rows are rendered. Batching replaces thousands of per-gene files without
removing brain-structure filters or aggregating any sample values. Gene lookup ignores
case; `brainspan.html?gene=GRIN2A&scoutTheme=dark` supports direct links.
After replacing per-gene exports with batches, hard-refresh an already-open page
(Ctrl+F5) to load the new HTML and versioned JavaScript modules.
If the preview retains an old page and reports `Invalid BrainSpan gene catalog`,
open `brainspan.html?dataVersion=2`. The home-page link uses this versioned URL
to avoid reusing the cached pre-batching HTML.
Press Draw after changing gene, driver, or profile. Region, age, smoothing,
confidence-band, and developmental-label controls redraw the loaded data.
Networks support node dragging, background panning, zoom, node selection,
keyboard focus, and an accessible edge table. Downloads provide SVG and
filtered-sample/network CSV, rather than server-generated PDF/HTML. Use browser
print to save PDF. Empty profiles and NDD-only networks are reported explicitly.
The original experiment/control comparison widgets had no implemented analysis
and are intentionally omitted.

Browser smoothing uses local quadratic tricube regression on log10(age), with
lower-degree fits for singular neighborhoods and approximate normal 95% bands.
These are not exact ggplot2 LOESS intervals or donor-adjusted statistical tests.
Sparse groups have no smooth; non-estimable confidence bands are omitted.

To regenerate, install `DBI` and `RSQLite` in R, then run from the repository root:

```powershell
Rscript docs\scripts\build-brainspan-data.R
node --test tests\brainspan.test.mjs
```

Default inputs are the sibling workspace's
`brainspan-dev\brainspan-dev\extdata\brainspan.db` and
`scRNASeq\NDD_human_mouse_trajectory_cerebellum_6_Phases\ndd_gene_provenance.csv`.
For another checkout, pass explicit source paths:

```powershell
Rscript docs\scripts\build-brainspan-data.R 'C:\data\brainspan.db' 'C:\data\ndd_gene_provenance.csv'
```

The exporter removes obsolete generated per-gene/batch CSVs after writing the
new batches and catalog. The database is opened read-only. The source files are only needed for regeneration,
not deployment; commit the generated CSVs along with the page and scripts.
Tests check every shard's NDD membership, driver mapping, sample IDs, availability
counts, missing-value handling, filtering, smoothing, and download CSV escaping.

### Data preparation and validation

Rankings are computed in the browser from each region's
`trajectory_similarity.csv`. Compact `trajectories.json` assets contain only
eligible gene/cell-type combinations, preserving the shared phase-level median
CPM values and biological replicate counts from `aligned_phase_summary.csv`.
They are fetched on first opening a region's detail panel and cached for
subsequent selections. To regenerate them (Node.js 18+):

```powershell
node docs\scripts\build-trajectories.mjs
node docs\scripts\build-trajectories.mjs --enriched
node --test tests\explorer.test.mjs
```

The generator reads the original cortex/cerebellum datasets by default;
`--enriched` generates the five-bin cortex and four-phase sensitivity assets.
It reads phase-summary CSVs and similarity CSVs
exclusively from `docs/data`, then writes the generated JSON files beside
them. No external source folders are needed to regenerate or serve the site.

Tests compare every gene summary with published rankings and verify that
every chart's phase data reproduce raw and normalized z-RMSE, direction
agreement, peak-phase distance, and composite divergence score.

To rescale the stored descriptive metrics and regenerate dependent rankings
and regional comparisons (Node.js 18+):

```powershell
node docs\scripts\rescale-divergence.mjs
```

For newly copied enriched results, normalize only those copies and regenerate
their rankings without touching original datasets or old regional comparisons:

```powershell
node docs\scripts\rescale-divergence.mjs --enriched
node docs\scripts\build-trajectories.mjs --enriched
node --test tests\explorer.test.mjs
```

Only web-sized summaries, metrics, rankings, coverage, QC, methods, and generated
trajectories are copied. H5AD files, caches, donor-level pseudobulk matrices, and
large exact-stage tables are intentionally excluded.
