import { readFile, writeFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import path from "node:path";
import { PHASES, parseCSV, prepareMetrics, numeric } from "./explorer-model.mjs";

const root = fileURLToPath(new URL("../", import.meta.url));
const inputs = [
  "NDD_human_mouse_cortical_trajectory",
  "NDD_cerebellum_four_phase_descriptive",
];
for (const folder of inputs) {
  const target = path.join(root, "data", folder);
  const sourcePath = path.join(target, "aligned_phase_summary.csv");
  const metrics = prepareMetrics(parseCSV(await readFile(path.join(target, "trajectory_similarity.csv"), "utf8")));
  const trajectories = Object.create(null);
  for (const row of metrics) {
    trajectories[row.cell_label] ??= Object.create(null);
    trajectories[row.cell_label][row.gene] = {};
  }
  const phases = parseCSV(await readFile(sourcePath, "utf8"));
  for (const row of phases) {
    if (row.resolution !== "cell_type") continue;
    const gene = trajectories[row.cell_label]?.[row.gene];
    const phaseIndex = PHASES.indexOf(row.phase);
    if (!gene || phaseIndex < 0 || !["human", "mouse"].includes(row.species)) continue;
    const species = gene[row.species] ??= { median_cpm: Array(4).fill(null), n_replicates: Array(4).fill(null) };
    if (species.median_cpm[phaseIndex] !== null) throw new Error(`Duplicate phase: ${row.gene}/${row.cell_label}/${row.species}/${row.phase}`);
    const cpm = numeric(row.median_cpm, "median_cpm");
    const replicates = numeric(row.n_replicates, "n_replicates");
    if (cpm < 0 || !Number.isInteger(replicates) || replicates < 1) throw new Error("Invalid phase CPM or replicate count.");
    species.median_cpm[phaseIndex] = cpm;
    species.n_replicates[phaseIndex] = replicates;
  }
  for (const row of metrics) {
    for (const species of ["human", "mouse"]) {
      const series = trajectories[row.cell_label][row.gene][species];
      if (!series || series.median_cpm.includes(null)) throw new Error(`Incomplete trajectory: ${row.gene}/${row.cell_label}/${species}`);
    }
  }
  const output = {
    phases: PHASES,
    source: "aligned_phase_summary.csv",
    expression: "Median CPM across pseudobulks with at least 20 cells; log1p is applied in the browser.",
    trajectories,
  };
  await writeFile(path.join(target, "trajectories.json"), JSON.stringify(output));
  console.log(`${folder}: ${metrics.length} gene/cell-type trajectories exported.`);
}
