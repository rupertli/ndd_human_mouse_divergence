export const PHASES = ["early_prenatal", "mid_prenatal", "late_prenatal", "neonatal"];

export function parseCSV(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const char = text[i];
    if (char === '"') {
      if (quoted && text[i + 1] === '"') {
        field += '"';
        i++;
      } else {
        quoted = !quoted;
      }
    } else if (char === "," && !quoted) {
      row.push(field);
      field = "";
    } else if ((char === "\n" || char === "\r") && !quoted) {
      row.push(field);
      if (row.some(value => value !== "")) rows.push(row);
      row = [];
      field = "";
      if (char === "\r" && text[i + 1] === "\n") i++;
    } else {
      field += char;
    }
  }
  if (quoted) throw new Error("Unterminated quoted CSV field.");
  row.push(field);
  if (row.some(value => value !== "")) rows.push(row);
  if (!rows.length) throw new Error("CSV is empty.");
  const headers = rows.shift();
  headers[0] = headers[0].replace(/^\uFEFF/, "");
  return rows.map((values, index) => {
    if (values.length !== headers.length) throw new Error(`CSV row ${index + 2} has an unexpected number of fields.`);
    return Object.fromEntries(headers.map((header, i) => [header, values[i]]));
  });
}

export function numeric(value, name) {
  if (value === "" || value == null || !Number.isFinite(Number(value))) {
    throw new Error(`Missing or invalid numeric value: ${name}.`);
  }
  return Number(value);
}

export function prepareMetrics(records) {
  const seen = new Set();
  const columns = ["spearman_rho", "z_rmse", "z_rmse_raw", "direction_agreement", "peak_phase_distance", "descriptive_divergence_score"];
  return records.filter(row => row.resolution === "cell_type").map(row => {
    const key = JSON.stringify([row.cell_label, row.gene]);
    if (!row.cell_label || !row.gene || seen.has(key)) throw new Error(`Missing or duplicate gene/cell-type entry: ${key}.`);
    seen.add(key);
    if (row.common_phases !== PHASES.join("|")) throw new Error(`Expected four shared phases for ${row.gene}.`);
    return {
      ...row,
      ...Object.fromEntries(columns.map(column => [
        column,
        row[column] === "" && ["spearman_rho", "descriptive_divergence_score"].includes(column)
          ? null : numeric(row[column], column),
      ])),
    };
  });
}

export function median(values) {
  const ordered = values.filter(Number.isFinite).sort((a, b) => a - b);
  if (!ordered.length) return null;
  const middle = Math.floor(ordered.length / 2);
  return ordered.length % 2 ? ordered[middle] : (ordered[middle - 1] + ordered[middle]) / 2;
}

export function summarizeGene(rows) {
  if (!rows.length) throw new Error("Cannot summarize an empty gene group.");
  const scored = rows.filter(row => Number.isFinite(row.descriptive_divergence_score));
  const scores = scored.map(row => row.descriptive_divergence_score);
  const maxRow = [...scored].sort((a, b) =>
    b.descriptive_divergence_score - a.descriptive_divergence_score || a.cell_label.localeCompare(b.cell_label))[0];
  return {
    gene: rows[0].gene,
    n_tested_cell_labels: rows.length,
    n_scored_cell_labels: scored.length,
    median_spearman_rho: median(rows.map(row => row.spearman_rho)),
    median_z_rmse: median(rows.map(row => row.z_rmse)),
    median_direction_agreement: median(rows.map(row => row.direction_agreement)),
    mean_descriptive_divergence_score: scores.length ? scores.reduce((sum, score) => sum + score, 0) / scores.length : null,
    score: scores.length ? Math.max(...scores) : null,
    maxCell: maxRow?.cell_label ?? null,
    rows,
  };
}

export function rankGenes(metrics, cellType = "all", direction = "divergent") {
  const groups = new Map();
  for (const row of metrics) {
    if (cellType !== "all" && row.cell_label !== cellType) continue;
    if (!groups.has(row.gene)) groups.set(row.gene, []);
    groups.get(row.gene).push(row);
  }
  const sign = direction === "similar" ? 1 : -1;
  return [...groups.values()].map(summarizeGene).filter(item => item.score !== null).sort((a, b) =>
    sign * (a.score - b.score) || a.gene.localeCompare(b.gene));
}

export function commonGenes(left, right) {
  const rightGenes = new Set(right.map(row => row.gene));
  return new Set(left.filter(row => rightGenes.has(row.gene)).map(row => row.gene));
}

export function transformTrajectory(cpm, scale) {
  const log = cpm.map(value => Math.log1p(value));
  if (scale !== "z") return log;
  const mean = log.reduce((sum, value) => sum + value, 0) / log.length;
  const sd = Math.sqrt(log.reduce((sum, value) => sum + (value - mean) ** 2, 0) / log.length);
  return log.map(value => sd ? (value - mean) / sd : 0);
}
