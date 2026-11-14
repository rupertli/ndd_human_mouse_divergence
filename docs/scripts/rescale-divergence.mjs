import { readFile, writeFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import path from "node:path";
import { parseCSV } from "./explorer-model.mjs";

const root = fileURLToPath(new URL("../data/", import.meta.url));
const regions = {
  cerebrum: "NDD_human_mouse_cortical_trajectory",
  cerebellum: "NDD_cerebellum_four_phase_descriptive",
};
const comparisonFolder = "NDD_cerebrum_cerebellum_four_phase_comparison";

function number(value, field) {
  const parsed = Number(value);
  if (value === "" || !Number.isFinite(parsed)) {
    throw new Error(`Missing or invalid ${field}: ${value}`);
  }
  return parsed;
}

function format(value) {
  if (value == null || !Number.isFinite(value)) return "";
  return String(value);
}

function csvEscape(value) {
  const text = value == null ? "" : String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function serializeCSV(rows, columns) {
  return [
    columns.map(csvEscape).join(","),
    ...rows.map(row => columns.map(column => csvEscape(row[column])).join(",")),
  ].join("\n") + "\n";
}

async function readCSV(file) {
  return parseCSV(await readFile(file, "utf8"));
}

async function writeCSV(file, rows, columns = Object.keys(rows[0] ?? {})) {
  await writeFile(file, serializeCSV(rows, columns), "utf8");
}

function rescaleRows(rows) {
  return rows.map(row => {
    const raw = row.z_rmse_raw === undefined || row.z_rmse_raw === ""
      ? number(row.z_rmse, "z_rmse")
      : number(row.z_rmse_raw, "z_rmse_raw");
    if (raw < -1e-12 || raw > 2 + 1e-9) {
      throw new Error(`Raw z-RMSE is outside [0, 2] for ${row.gene}/${row.cell_label}: ${raw}`);
    }
    const normalized = raw / 2;
    const rho = row.spearman_rho === "" ? null : number(row.spearman_rho, "spearman_rho");
    const score = rho == null
      ? null
      : (1 - Math.max(-1, Math.min(1, rho))) / 2
        + normalized
        + (1 - number(row.direction_agreement, "direction_agreement"))
        + number(row.peak_phase_distance, "peak_phase_distance") / 3;
    row.z_rmse_raw = format(raw);
    row.z_rmse = format(normalized);
    row.descriptive_divergence_score = format(score);
    if ("difference_priority_score" in row) {
      row.difference_priority_score = format(score);
    }
    return row;
  });
}

function median(values) {
  const sorted = values.filter(Number.isFinite).sort((a, b) => a - b);
  if (!sorted.length) return null;
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2
    ? sorted[middle]
    : (sorted[middle - 1] + sorted[middle]) / 2;
}

function summarizeRanking(rows) {
  const groups = new Map();
  for (const row of rows) {
    const key = `${row.resolution}\u0000${row.gene}`;
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(row);
  }
  const result = [];
  for (const group of groups.values()) {
    const first = group[0];
    const finite = group
      .map(row => row.descriptive_divergence_score === ""
        ? null : number(row.descriptive_divergence_score, "descriptive_divergence_score"))
      .filter(value => value != null);
    const rhos = group.map(row => row.spearman_rho === "" ? null : Number(row.spearman_rho));
    const zRmses = group.map(row => row.z_rmse === "" ? null : Number(row.z_rmse));
    const directions = group.map(row => row.direction_agreement === "" ? null : Number(row.direction_agreement));
    const peaks = group.map(row => row.peak_phase_distance === "" ? null : Number(row.peak_phase_distance));
    const summary = {
      brain_region: first.brain_region,
      resolution: first.resolution,
      gene: first.gene,
      n_tested_cell_labels: group.length,
      median_spearman_rho: format(median(rhos)),
      median_z_rmse: format(median(zRmses)),
      median_direction_agreement: format(median(directions)),
      max_peak_phase_distance: format(peaks.length ? Math.max(...peaks.filter(Number.isFinite)) : null),
      max_descriptive_divergence_score: format(finite.length ? Math.max(...finite) : null),
      mean_descriptive_divergence_score: format(
        finite.length ? finite.reduce((sum, value) => sum + value, 0) / finite.length : null,
      ),
      inference_status: first.inference_status,
      analysis_scope: first.analysis_scope,
      phase_set: first.phase_set,
      integration_schema_version: first.integration_schema_version,
    };
    result.push(summary);
  }
  return result.sort((a, b) => {
    const left = a.max_descriptive_divergence_score === ""
      ? null : Number(a.max_descriptive_divergence_score);
    const right = b.max_descriptive_divergence_score === ""
      ? null : Number(b.max_descriptive_divergence_score);
    return a.resolution.localeCompare(b.resolution)
      || (left == null && right != null ? 1 : 0)
      || (left != null && right == null ? -1 : 0)
      || (left != null && right != null ? right - left : 0)
      || a.gene.localeCompare(b.gene);
  });
}

function rankForComparison(rows) {
  const sorted = [...rows].sort((a, b) => {
    const left = a.max_descriptive_divergence_score === ""
      ? null : Number(a.max_descriptive_divergence_score);
    const right = b.max_descriptive_divergence_score === ""
      ? null : Number(b.max_descriptive_divergence_score);
    if (left == null && right != null) return 1;
    if (left != null && right == null) return -1;
    if (left != null && right != null && left !== right) return right - left;
    return a.gene.localeCompare(b.gene);
  });
  const finiteCount = sorted.filter(row => row.max_descriptive_divergence_score !== "").length;
  return new Map(sorted.map((row, index) => {
    const rank = index + 1;
    const finite = row.max_descriptive_divergence_score !== "";
    return [row.gene, {
      row,
      rank,
      percentile: finite && finiteCount > 1
        ? 1 - (rank - 1) / finiteCount
        : finite ? 1 : null,
      finite,
    }];
  }));
}

function pearson(left, right) {
  const meanLeft = left.reduce((sum, value) => sum + value, 0) / left.length;
  const meanRight = right.reduce((sum, value) => sum + value, 0) / right.length;
  let covariance = 0;
  let varianceLeft = 0;
  let varianceRight = 0;
  for (let i = 0; i < left.length; i++) {
    const dx = left[i] - meanLeft;
    const dy = right[i] - meanRight;
    covariance += dx * dy;
    varianceLeft += dx * dx;
    varianceRight += dy * dy;
  }
  return covariance / Math.sqrt(varianceLeft * varianceRight);
}

function logGamma(z) {
  const coefficients = [
    676.5203681218851, -1259.1392167224028, 771.3234287776531,
    -176.6150291621406, 12.507343278686905, -0.13857109526572012,
    9.984369578019572e-6, 1.5056327351493116e-7,
  ];
  if (z < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - logGamma(1 - z);
  const x = z - 1;
  let sum = 0.9999999999998099;
  for (let i = 0; i < coefficients.length; i++) sum += coefficients[i] / (x + i + 1);
  const t = x + coefficients.length - 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(sum);
}

function betaContinuedFraction(a, b, x) {
  const maxIterations = 200;
  const epsilon = 3e-14;
  const tiny = 1e-300;
  const qab = a + b;
  const qap = a + 1;
  const qam = a - 1;
  let c = 1;
  let d = 1 - qab * x / qap;
  if (Math.abs(d) < tiny) d = tiny;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= maxIterations; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d;
    if (Math.abs(d) < tiny) d = tiny;
    c = 1 + aa / c;
    if (Math.abs(c) < tiny) c = tiny;
    d = 1 / d;
    const delta = d * c;
    h *= delta;
    if (Math.abs(delta - 1) < epsilon) return h;
  }
  throw new Error("Incomplete beta calculation did not converge.");
}

function regularizedBeta(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const factor = Math.exp(
    logGamma(a + b) - logGamma(a) - logGamma(b)
    + a * Math.log(x) + b * Math.log1p(-x),
  );
  if (x < (a + 1) / (a + b + 2)) {
    return factor * betaContinuedFraction(a, b, x) / a;
  }
  return 1 - factor * betaContinuedFraction(b, a, 1 - x) / b;
}

function spearmanPValue(rho, n) {
  if (n < 3 || Math.abs(rho) >= 1) return Math.abs(rho) >= 1 ? 0 : null;
  const degrees = n - 2;
  const tSquared = rho * rho * degrees / (1 - rho * rho);
  return regularizedBeta(degrees / (degrees + tSquared), degrees / 2, 0.5);
}

function buildRegionalComparison(rankings) {
  const ranked = {
    cerebrum: rankForComparison(rankings.cerebrum.filter(row => row.resolution === "cell_type")),
    cerebellum: rankForComparison(rankings.cerebellum.filter(row => row.resolution === "cell_type")),
  };
  const genes = [...new Set([...ranked.cerebrum.keys(), ...ranked.cerebellum.keys()])].sort();
  const rows = genes.map(gene => {
    const output = { gene };
    for (const region of ["cerebrum", "cerebellum"]) {
      const item = ranked[region].get(gene);
      const row = item?.row;
      for (const [outputField, sourceField] of [
        ["max_descriptive_divergence_score", "max_descriptive_divergence_score"],
        ["median_spearman_rho", "median_spearman_rho"],
        ["median_z_rmse", "median_z_rmse"],
        ["median_direction_agreement", "median_direction_agreement"],
        ["max_peak_phase_distance", "max_peak_phase_distance"],
      ]) output[`${region}_${outputField}`] = row?.[sourceField] ?? "";
      output[`${region}_rank`] = item ? String(item.rank) : "";
      output[`${region}_percentile`] = item?.finite ? format(item.percentile) : "";
    }
    const left = ranked.cerebrum.get(gene);
    const right = ranked.cerebellum.get(gene);
    const both = Boolean(left && right);
    const bothScored = Boolean(left?.finite && right?.finite);
    output.regional_percentile_difference = bothScored
      ? format(left.percentile - right.percentile)
      : "";
    output.available_in_both = String(both);
    return output;
  });
  rows.sort((a, b) =>
    Math.abs(Number(b.regional_percentile_difference || 0))
      - Math.abs(Number(a.regional_percentile_difference || 0))
    || a.gene.localeCompare(b.gene));

  const common = genes.filter(gene =>
    ranked.cerebrum.has(gene) && ranked.cerebellum.has(gene));
  const commonRank = region => new Map(
    [...common]
      .sort((left, right) =>
        ranked[region].get(left).rank - ranked[region].get(right).rank)
      .map((gene, index) => [gene, index + 1]),
  );
  const commonRanks = {
    cerebrum: commonRank("cerebrum"),
    cerebellum: commonRank("cerebellum"),
  };
  const rho = pearson(
    common.map(gene => commonRanks.cerebrum.get(gene)),
    common.map(gene => commonRanks.cerebellum.get(gene)),
  );
  const topOverlap = {};
  for (const size of [5, 10, 20, 50, 100]) {
    const left = [...ranked.cerebrum.entries()]
      .filter(([, value]) => value.finite)
      .sort((a, b) => a[1].rank - b[1].rank)
      .slice(0, size)
      .map(([gene]) => gene);
    const right = [...ranked.cerebellum.entries()]
      .filter(([, value]) => value.finite)
      .sort((a, b) => a[1].rank - b[1].rank)
      .slice(0, size)
      .map(([gene]) => gene);
    const intersection = left.filter(gene => right.includes(gene)).sort();
    const unionSize = new Set([...left, ...right]).size;
    topOverlap[size] = {
      count: intersection.length,
      jaccard: unionSize ? intersection.length / unionSize : 0,
      genes: intersection,
    };
  }
  const cerebrumTop = [...ranked.cerebrum.entries()]
    .filter(([, value]) => value.finite)
    .sort((a, b) => a[1].rank - b[1].rank)
    .slice(0, 5);
  const cerebellumTop = [...ranked.cerebellum.entries()]
    .filter(([, value]) => value.finite)
    .sort((a, b) => a[1].rank - b[1].rank)
    .slice(0, 5);
  const summary = {
    cerebrum_ranked_genes: ranked.cerebrum.size,
    cerebellum_ranked_genes: ranked.cerebellum.size,
    common_ranked_genes: common.length,
    rank_spearman_rho: rho,
    rank_spearman_p_descriptive_only: spearmanPValue(rho, common.length),
    top_overlap: topOverlap,
    cerebrum_top5_in_cerebellum: Object.fromEntries(cerebrumTop.map(([gene]) =>
      [gene, ranked.cerebellum.get(gene)?.finite ? ranked.cerebellum.get(gene).rank : null])),
    cerebellum_top5_in_cerebrum: Object.fromEntries(cerebellumTop.map(([gene]) =>
      [gene, ranked.cerebrum.get(gene)?.finite ? ranked.cerebrum.get(gene).rank : null])),
  };
  return { rows, summary };
}

const similarityByRegion = {};
const rankings = {};
for (const [region, folder] of Object.entries(regions)) {
  const directory = path.join(root, folder);
  const similarityPath = path.join(directory, "trajectory_similarity.csv");
  const similarity = rescaleRows(await readCSV(similarityPath));
  const similarityColumns = Object.keys(similarity[0]);
  if (!similarityColumns.includes("z_rmse_raw")) {
    similarityColumns.splice(similarityColumns.indexOf("z_rmse") + 1, 0, "z_rmse_raw");
  }
  await writeCSV(similarityPath, similarity, similarityColumns);
  similarityByRegion[region] = similarity;

  const integrationPath = path.join(directory, "region_integration_trajectory_metrics.csv");
  const integration = rescaleRows(await readCSV(integrationPath));
  const integrationColumns = Object.keys(integration[0]);
  if (!integrationColumns.includes("z_rmse_raw")) {
    integrationColumns.splice(integrationColumns.indexOf("z_rmse") + 1, 0, "z_rmse_raw");
  }
  await writeCSV(integrationPath, integration, integrationColumns);

  const rankingPath = path.join(directory, "gene_difference_ranking.csv");
  const oldRanking = await readCSV(rankingPath);
  const ranking = summarizeRanking(similarity);
  rankings[region] = ranking;
  await writeCSV(rankingPath, ranking, Object.keys(oldRanking[0]));
}

const comparison = buildRegionalComparison(rankings);
const comparisonDirectory = path.join(root, comparisonFolder);
const comparisonPath = path.join(comparisonDirectory, "gene_region_rank_comparison.csv");
const previousComparison = await readCSV(comparisonPath);
await writeCSV(comparisonPath, comparison.rows, Object.keys(previousComparison[0]));
await writeFile(
  path.join(comparisonDirectory, "rank_overlap_summary.json"),
  `${JSON.stringify(comparison.summary, null, 2)}\n`,
  "utf8",
);

const corticalDirectory = path.join(root, regions.cerebrum);
const topGenesPath = path.join(corticalDirectory, "top5_most_different_ndd_genes.csv");
const topGenes = await readCSV(topGenesPath);
const calculatedTop = rankings.cerebrum
  .filter(row => row.resolution === "cell_type" && row.max_descriptive_divergence_score !== "")
  .slice(0, 5);
const existingTop = new Set(topGenes.map(row => row.gene));
if (calculatedTop.some(row => !existingTop.has(row.gene))) {
  throw new Error("The new top five differs from the curated table; refresh its provenance annotations before writing.");
}
const topMap = new Map(calculatedTop.map((row, index) => [row.gene, { ...row, rank: index + 1 }]));
const updatedTopGenes = topGenes.map(row => ({ ...row, ...topMap.get(row.gene) }));
await writeCSV(topGenesPath, updatedTopGenes, Object.keys(topGenes[0]));

const evidencePath = path.join(corticalDirectory, "top5_lineage_trajectory_evidence.csv");
const evidence = rescaleRows(await readCSV(evidencePath));
const topRanks = new Map(calculatedTop.map((row, index) => [row.gene, String(index + 1)]));
for (const row of evidence) row.rank = topRanks.get(row.gene) ?? "";
evidence.sort((a, b) =>
  Number(a.rank || Infinity) - Number(b.rank || Infinity)
  || Number(b.descriptive_divergence_score || -Infinity) - Number(a.descriptive_divergence_score || -Infinity)
  || a.cell_label.localeCompare(b.cell_label));
const evidenceColumns = Object.keys(evidence[0]);
if (!evidenceColumns.includes("z_rmse_raw")) {
  evidenceColumns.splice(evidenceColumns.indexOf("z_rmse") + 1, 0, "z_rmse_raw");
}
await writeCSV(evidencePath, evidence, evidenceColumns);

console.log("Rescaled z-RMSE to [0, 1] as raw z-RMSE / 2.");
for (const [region, rows] of Object.entries(rankings)) {
  console.log(`${region}: ${rows.filter(row => row.max_descriptive_divergence_score !== "").length} scored genes; top 5 ${rows.slice(0, 5).map(row => row.gene).join(", ")}`);
}
