import test from "node:test";
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { PHASES, parseCSV, prepareMetrics, rankGenes, commonGenes, transformTrajectory } from "../sites/scripts/explorer-model.mjs";

const dataRoot = new URL("../sites/data/", import.meta.url);
const folders = ["NDD_human_mouse_cortical_trajectory", "NDD_cerebellum_four_phase_descriptive"];
const read = (folder, file) => readFile(new URL(`${folder}/${file}`, dataRoot), "utf8");
const close = (actual, expected) => assert.ok(Math.abs(actual - expected) < 1e-10, `${actual} != ${expected}`);

test("CSV handles quoted commas, escaped quotes, CRLF, and malformed input", () => {
  assert.deepEqual(parseCSV('gene,label\r\nA,"a, b"\r\nB,"a ""quote"""\r\n'), [
    { gene: "A", label: "a, b" }, { gene: "B", label: 'a "quote"' },
  ]);
  assert.throws(() => parseCSV('a,b\n"unfinished'), /Unterminated/);
  assert.throws(() => parseCSV("a,b\none"), /fields/);
});

test("max ranking, least-divergent reversal, per-cell filtering, summaries, ties, and overlap", () => {
  const make = (gene, cell, score, rho) => ({
    gene, cell_label: cell, descriptive_divergence_score: score,
    spearman_rho: rho, z_rmse: score / 3, direction_agreement: .5,
  });
  const rows = [make("A", "one", 4, -.5), make("A", "two", 0, 1), make("B", "one", 2, .5), make("C", "one", 2, .5)];
  const maximum = rankGenes(rows);
  assert.deepEqual(maximum.map(row => row.gene), ["A", "B", "C"]);
  assert.equal(maximum[0].median_spearman_rho, .25);
  assert.equal(maximum[0].mean_descriptive_divergence_score, 2);
  assert.equal(maximum[0].maxCell, "one");
  assert.deepEqual(rankGenes(rows, "all", "similar").map(row => row.gene), ["B", "C", "A"]);
  assert.equal(rankGenes(rows, "two")[0].score, 0);
  assert.equal(rankGenes(rows, "two")[0].n_tested_cell_labels, 1);
  assert.deepEqual([...commonGenes(maximum.slice(0, 2), maximum.slice(1))], ["B"]);
  assert.deepEqual(transformTrajectory([4, 4, 4, 4], "z"), [0, 0, 0, 0]);
});

test("undefined scores are not misrepresented as similarity and missing correlations are omitted from medians", () => {
  const rows = [
    { gene: "A", cell_label: "one", descriptive_divergence_score: null, spearman_rho: null, z_rmse: 1, direction_agreement: 0 },
    { gene: "A", cell_label: "two", descriptive_divergence_score: 3, spearman_rho: -.5, z_rmse: 1.5, direction_agreement: 0 },
    { gene: "B", cell_label: "one", descriptive_divergence_score: null, spearman_rho: null, z_rmse: 1, direction_agreement: 0 },
  ];
  const ranked = rankGenes(rows, "all", "similar");
  assert.equal(ranked.length, 1);
  assert.equal(ranked[0].gene, "A");
  assert.equal(ranked[0].n_tested_cell_labels, 2);
  assert.equal(ranked[0].n_scored_cell_labels, 1);
  assert.equal(ranked[0].median_spearman_rho, -.5);
  assert.equal(ranked[0].median_z_rmse, 1.25);
  assert.equal(ranked[0].mean_descriptive_divergence_score, 3);
  assert.equal(rankGenes(rows, "one", "similar").length, 0);
});

const metrics = await Promise.all(folders.map(async folder => prepareMetrics(parseCSV(await read(folder, "trajectory_similarity.csv")))));

for (const [i, folder] of folders.entries()) {
  test(`${folder}: all generated rankings reproduce published gene summaries`, async () => {
    const published = parseCSV(await read(folder, "gene_difference_ranking.csv")).filter(row => row.max_descriptive_divergence_score !== "");
    const ranked = rankGenes(metrics[i]);
    assert.equal(ranked.length, published.length);
    const byGene = new Map(published.map(row => [row.gene, row]));
    for (const item of ranked) {
      const source = byGene.get(item.gene);
      for (const key of ["n_tested_cell_labels", "median_spearman_rho", "median_z_rmse", "median_direction_agreement", "mean_descriptive_divergence_score"]) {
        close(item[key], Number(source[key]));
      }
      close(item.score, Number(source.max_descriptive_divergence_score));
    }
  });

  test(`${folder}: all chart trajectories reproduce source z-RMSE, directions, peak distances, and scores`, async () => {
    const data = JSON.parse(await read(folder, "trajectories.json"));
    assert.deepEqual(data.phases, PHASES);
    for (const row of metrics[i]) {
      const series = data.trajectories[row.cell_label][row.gene];
      const h = series.human.median_cpm;
      const m = series.mouse.median_cpm;
      for (const species of ["human", "mouse"]) {
        assert.equal(series[species].median_cpm.length, 4);
        assert.ok(series[species].n_replicates.every(n => Number.isInteger(n) && n > 0));
      }
      const hz = transformTrajectory(h, "z");
      const mz = transformTrajectory(m, "z");
      const rmse = Math.sqrt(hz.reduce((sum, value, j) => sum + (value - mz[j]) ** 2, 0) / 4);
      const agreement = h.slice(1).filter((value, j) => Math.sign(value - h[j]) === Math.sign(m[j + 1] - m[j])).length / 3;
      const peak = Math.abs(h.indexOf(Math.max(...h)) - m.indexOf(Math.max(...m)));
      close(rmse, row.z_rmse_raw);
      close(rmse / 2, row.z_rmse);
      assert.ok(row.z_rmse >= 0 && row.z_rmse <= 1);
      close(agreement, row.direction_agreement);
      close(peak, row.peak_phase_distance);
      if (row.spearman_rho === null) {
        assert.equal(row.descriptive_divergence_score, null);
        assert.ok(new Set(h).size === 1 || new Set(m).size === 1);
      } else {
        close((1 - row.spearman_rho) / 2 + rmse / 2 + (1 - agreement) + peak / 3, row.descriptive_divergence_score);
        assert.ok(row.descriptive_divergence_score >= 0 && row.descriptive_divergence_score <= 4);
      }
    }
  });
}

test("four-phase default and top-50 overlap reproduce the regional comparison", async () => {
  const ranked = metrics.map(rows => rankGenes(rows));
  assert.equal(commonGenes(...ranked.map(rows => rows.slice(0, 20))).size, 0);
  assert.deepEqual([...commonGenes(...ranked.map(rows => rows.slice(0, 50)))].sort(), ["PIGH", "PLCH1", "TMTC3"]);

  const comparisonFolder = "NDD_cerebrum_cerebellum_four_phase_comparison";
  const comparison = parseCSV(await read(comparisonFolder, "gene_region_rank_comparison.csv"))
    .filter(row => row.available_in_both === "true");
  const summary = JSON.parse(await read(comparisonFolder, "rank_overlap_summary.json"));
  const rankWithinShared = field => new Map(
    [...comparison]
      .sort((left, right) => Number(left[field]) - Number(right[field]))
      .map((row, index) => [row.gene, index + 1]),
  );
  const sharedRanks = {
    cerebrum: rankWithinShared("cerebrum_rank"),
    cerebellum: rankWithinShared("cerebellum_rank"),
  };
  const sharedCerebrumRanks = comparison.map(row => sharedRanks.cerebrum.get(row.gene));
  const sharedCerebellumRanks = comparison.map(row => sharedRanks.cerebellum.get(row.gene));
  const mean = values => values.reduce((sum, value) => sum + value, 0) / values.length;
  const meanCerebrumRank = mean(sharedCerebrumRanks);
  const meanCerebellumRank = mean(sharedCerebellumRanks);
  const covariance = sharedCerebrumRanks.reduce((sum, value, index) =>
    sum + (value - meanCerebrumRank) * (sharedCerebellumRanks[index] - meanCerebellumRank), 0);
  const varianceCerebrum = sharedCerebrumRanks.reduce((sum, value) =>
    sum + (value - meanCerebrumRank) ** 2, 0);
  const varianceCerebellum = sharedCerebellumRanks.reduce((sum, value) =>
    sum + (value - meanCerebellumRank) ** 2, 0);
  close(summary.rank_spearman_rho, covariance / Math.sqrt(varianceCerebrum * varianceCerebellum));
});
