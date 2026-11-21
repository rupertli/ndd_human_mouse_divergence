import { numeric } from "./explorer-model.mjs";

export const STAGES = {
  Classic: [[35, 63, "Embryonic"], [63, 266, "Fetal development"], [266, 650, "Infancy"], [650, 4300, "Childhood"], [4300, 7220, "Adolescence"], [7220, 14900, "Adulthood"]],
  Trimester: [[35, 63, "Embryonic"], [63, 105, "Prenatal 1st trimester"], [105, 175, "Prenatal 2nd trimester"], [175, 266, "Prenatal 3rd trimester"], [266, 7220, "Childhood"], [7220, 14900, "Adulthood"]],
};
export const REGION_FIELDS = ["structure_acronym", "structure_group1", "structure_group2"];

export function selectGeneData(batch, gene) {
  return {
    profiles: batch.profiles.filter(row => row.gene === gene),
    network: batch.network.filter(row => row.source_gene === gene),
  };
}

export function prepareSamples(rows) {
  const ids = new Set();
  return rows.map(row => {
    if (!row.sample_id || ids.has(row.sample_id)) throw new Error("Missing or duplicate BrainSpan sample ID.");
    ids.add(row.sample_id);
    const Days = numeric(row.Days, "Days");
    if (Days <= 0 || REGION_FIELDS.some(field => !row[field])) throw new Error("Invalid BrainSpan sample metadata.");
    return { ...row, Days };
  });
}

export function profilePoints(rows, samples, profile, id) {
  const matches = rows.filter(row => row.profile === profile && row.id === id);
  if (matches.length > 1) throw new Error(`Duplicate ${profile} profile for ${id}.`);
  if (!matches.length) return [];
  return samples.flatMap(sample => {
    const value = matches[0][sample.sample_id];
    if (!(sample.sample_id in matches[0])) throw new Error(`Missing profile sample ${sample.sample_id}.`);
    if (value === "") return [];
    return [{ ...sample, value: numeric(value, sample.sample_id) }];
  });
}

export function filterPoints(points, options) {
  const { start, end, regions } = options;
  if (!Number.isFinite(start) || !Number.isFinite(end) || start <= 0 || end < start) throw new Error("Enter a valid positive age range (start ≤ end).");
  return points.filter(point => point.Days >= start && point.Days <= end &&
    REGION_FIELDS.every(field => regions[field].includes("All") || regions[field].includes(point[field])));
}

function solve(matrix, rhs) {
  const rows = matrix.map((row, i) => [...row, rhs[i]]);
  for (let i = 0; i < rhs.length; i++) {
    const pivot = rows.slice(i).reduce((best, row, offset) => Math.abs(row[i]) > Math.abs(rows[best][i]) ? i + offset : best, i);
    if (Math.abs(rows[pivot][i]) < 1e-12) return null;
    [rows[i], rows[pivot]] = [rows[pivot], rows[i]];
    const divisor = rows[i][i];
    rows[i] = rows[i].map(value => value / divisor);
    for (let j = 0; j < rhs.length; j++) {
      if (j === i) continue;
      const factor = rows[j][i];
      rows[j] = rows[j].map((value, k) => value - factor * rows[i][k]);
    }
  }
  return rows.map(row => row.at(-1));
}

// Local quadratic tricube regression on log10(age). Bands are approximate,
// not ggplot2's predict.loess confidence intervals.
export function smoothProfile(points, span, steps = 60) {
  if (!(span >= .05 && span <= 1)) throw new Error("Smooth span must be between 0.05 and 1.");
  if (points.length < 4 || new Set(points.map(point => point.Days)).size < 3) return [];
  const data = points.map(point => ({ x: Math.log10(point.Days), y: point.value }));
  const min = Math.min(...data.map(p => p.x)), max = Math.max(...data.map(p => p.x));
  const count = Math.min(data.length, Math.max(4, Math.ceil(span * data.length)));
  return Array.from({ length: steps }, (_, i) => {
    const x = min + (max - min) * i / (steps - 1);
    const neighbors = [...data].sort((a, b) => Math.abs(a.x - x) - Math.abs(b.x - x));
    const radius = Math.abs(neighbors[count - 1].x - x) * 1.000001 || (max - min) / 100;
    const local = data.map(p => ({ ...p, z: (p.x - x) / radius, w: (1 - Math.min(1, Math.abs(p.x - x) / radius) ** 3) ** 3 })).filter(p => p.w > 0);
    let coefficients, matrix, degree;
    for (degree = 2; degree >= 0; degree--) {
      matrix = Array.from({ length: degree + 1 }, (_, a) => Array.from({ length: degree + 1 }, (_, b) => local.reduce((sum, p) => sum + p.w * p.z ** (a + b), 0)));
      coefficients = solve(matrix, Array.from({ length: degree + 1 }, (_, a) => local.reduce((sum, p) => sum + p.w * p.z ** a * p.y, 0)));
      if (coefficients) break;
    }
    if (!coefficients) throw new Error("Unable to fit local regression.");
    const influence = solve(matrix, coefficients.map((_, i) => i === 0 ? 1 : 0));
    const weightSum = local.reduce((sum, p) => sum + p.w, 0);
    const dof = weightSum - degree - 1;
    const variance = dof > 0 ? local.reduce((sum, p) => sum + p.w * (p.y - coefficients.reduce((fit, c, j) => fit + c * p.z ** j, 0)) ** 2, 0) / dof : null;
    const se = variance === null ? null : Math.sqrt(variance * local.reduce((sum, p) => sum + (p.w * influence.reduce((v, c, j) => v + c * p.z ** j, 0)) ** 2, 0));
    return { Days: 10 ** x, value: coefficients[0], low: se === null ? null : coefficients[0] - 1.96 * se, high: se === null ? null : coefficients[0] + 1.96 * se };
  });
}

export function driverNetwork(rows, driver, genes) {
  return rows.filter(row => row.source === driver).map(row => {
    if (!genes.has(row.source_gene) || !genes.has(row.target_gene)) throw new Error("Network contains a non-NDD gene.");
    const MI = numeric(row.MI, "MI");
    if (MI < 0) throw new Error("Network MI must not be negative.");
    return { ...row, MI, pearson: row.pearson === "" ? null : numeric(row.pearson, "pearson"), spearman: row.spearman === "" ? null : numeric(row.spearman, "spearman") };
  });
}

export function toCSV(rows, columns) {
  const quote = value => `"${String(value ?? "").replaceAll('"', '""')}"`;
  return [columns.map(quote).join(","), ...rows.map(row => columns.map(column => quote(row[column])).join(","))].join("\r\n") + "\r\n";
}
