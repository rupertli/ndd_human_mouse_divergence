import { PHASES, parseCSV, prepareMetrics, rankGenes, commonGenes, transformTrajectory } from "./explorer-model.mjs";

const dataRoot = new URL("../data/", import.meta.url);
const configs = [
  { id: "cerebrum", title: "Cerebrum · cerebral cortex", folder: "NDD_human_mouse_cortical_trajectory", replicate: "mouse pools" },
  { id: "cerebellum", title: "Cerebellum", folder: "NDD_cerebellum_four_phase_descriptive", replicate: "individual mice" },
];
const phaseLabels = ["Early prenatal", "Mid prenatal", "Late prenatal", "Neonatal"];
const controls = {
  n: document.querySelector("#top-n"),
  direction: document.querySelector("#direction"),
  overlap: document.querySelector("#overlap"),
  error: document.querySelector("#error"),
};
const regions = [];
const format = value => Number.isFinite(value) ? value.toFixed(3) : "Not estimable";

function el(tag, text, className) {
  const element = document.createElement(tag);
  if (text !== undefined) element.textContent = text;
  if (className) element.className = className;
  return element;
}

function option(value, label) {
  const element = el("option", label);
  element.value = value;
  return element;
}

function icon(name) {
  const svg = svgElement("svg", { class: "icon", "aria-hidden": "true", focusable: "false" });
  svg.append(svgElement("use", { href: `#icon-${name}` }));
  return svg;
}

async function load(url, json = false) {
  const response = await fetch(new URL(url, dataRoot));
  if (!response.ok) throw new Error(`Could not load ${url} (HTTP ${response.status}).`);
  return json ? response.json() : response.text();
}

function showError(error) {
  controls.error.hidden = false;
  controls.error.textContent = `${error.message} Serve the sites folder over HTTP or GitHub Pages; opening index.html directly cannot load data.`;
}

function createRegion(config, metrics) {
  const state = { ...config, metrics, cellType: "all", selected: null, chartCell: null, scale: "log", trajectoriesPromise: null };
  const section = el("section", undefined, "region");
  section.id = config.id;
  const card = el("div", undefined, "card");
  const header = el("div", undefined, "region-header");
  const heading = el("h2", undefined, "icon-heading");
  heading.append(icon("brain"), el("span", config.title));
  header.append(heading, el("span", `${rankGenes(metrics).length.toLocaleString()} rankable genes`, "gene-count"));
  const label = el("label", "Rank genes within");
  const select = el("select", undefined, "cell-select");
  select.id = `${config.id}-cell-type`;
  select.append(option("all", "Maximum across cell types"));
  const cells = [...new Set(metrics.map(row => row.cell_label))].sort();
  for (const cell of cells) select.append(option(cell, cell));
  label.append(select);
  state.scope = el("p", undefined, "scope");
  const wrapper = el("div", undefined, "list-wrap");
  wrapper.tabIndex = 0;
  wrapper.setAttribute("role", "region");
  wrapper.setAttribute("aria-label", `${config.title} gene rankings`);
  const table = el("table", undefined, "rankings");
  const head = el("thead");
  const headings = el("tr");
  for (const text of ["Rank", "Gene", "Divergence"]) {
    const th = el("th", text);
    th.scope = "col";
    headings.append(th);
  }
  head.append(headings);
  state.tbody = el("tbody");
  table.append(head, state.tbody);
  wrapper.append(table);
  state.detail = el("section", undefined, "card detail");
  state.detail.setAttribute("aria-label", `${config.title} selected gene details`);
  card.append(header, label, state.scope, wrapper);
  section.append(card, state.detail);
  document.querySelector("#regions").append(section);
  select.addEventListener("change", () => {
    state.cellType = select.value;
    state.chartCell = null;
    refresh();
  });
  return state;
}

function refresh() {
  if (!controls.n.reportValidity()) return;
  const n = Number(controls.n.value);
  for (const state of regions) {
    state.ranking = rankGenes(state.metrics, state.cellType, controls.direction.value);
    state.visible = state.ranking.slice(0, n);
  }
  const shared = commonGenes(regions[0].visible, regions[1].visible);
  controls.overlap.textContent = `${shared.size} shared ${shared.size === 1 ? "gene" : "genes"} · ${regions[0].visible.length} cortex / ${regions[1].visible.length} cerebellum shown`;
  for (const state of regions) {
    state.scope.textContent = state.cellType === "all"
      ? "Ranking by maximum divergence across eligible cell types. Least divergent = lowest maximum."
      : `Ranking within ${state.cellType}. Summary values below refer only to this cell type.`;
    if (!state.visible.some(row => row.gene === state.selected)) {
      state.selected = state.visible[0]?.gene ?? null;
      state.chartCell = null;
    }
    state.tbody.replaceChildren();
    for (const [index, item] of state.visible.entries()) {
      const tr = el("tr");
      tr.classList.toggle("common", shared.has(item.gene));
      tr.classList.toggle("selected", item.gene === state.selected);
      const gene = el("td");
      const button = el("button", item.gene, "gene-button");
      button.type = "button";
      button.setAttribute("aria-pressed", String(item.gene === state.selected));
      if (shared.has(item.gene)) button.append(el("span", "Shared", "badge"));
      button.addEventListener("click", () => {
        state.selected = item.gene;
        state.chartCell = null;
        for (const row of state.tbody.rows) {
          const selected = row === tr;
          row.classList.toggle("selected", selected);
          row.querySelector("button").setAttribute("aria-pressed", String(selected));
        }
        renderDetail(state);
      });
      gene.append(button);
      const score = el("td", item.score.toFixed(6));
      score.title = String(item.score);
      tr.append(el("td", String(index + 1)), gene, score);
      state.tbody.append(tr);
    }
    renderDetail(state);
  }
}

function renderDetail(state) {
  const item = state.visible.find(row => row.gene === state.selected);
  state.detail.replaceChildren();
  if (!item) {
    state.detail.append(el("p", "No eligible genes for this selection."));
    return;
  }
  const title = el("div", undefined, "detail-title");
  const geneHeading = el("h3", undefined, "gene-name icon-heading");
  geneHeading.append(icon("dna"), el("span", item.gene));
  title.append(geneHeading, el("span", `${item.n_scored_cell_labels}/${item.n_tested_cell_labels} cell types scored`, "muted"));
  const context = state.cellType === "all"
    ? "Summary across all eligible cell types; the chart shows one cell type at a time."
    : "Single-cell-type summary: medians and mean equal this cell type's values.";
  const metrics = el("dl", undefined, "metrics");
  for (const key of ["median_spearman_rho", "median_z_rmse", "median_direction_agreement", "mean_descriptive_divergence_score"]) {
    const metric = el("div", undefined, "metric");
    metric.append(el("dt", key), el("dd", format(item[key])));
    metrics.append(metric);
  }
  const chartControls = el("div", undefined, "chart-controls");
  const cellLabel = el("label", "Trajectory cell type");
  const cellSelect = el("select");
  cellSelect.id = `${state.id}-chart-cell`;
  for (const row of [...item.rows].sort((a, b) => a.cell_label.localeCompare(b.cell_label))) {
    cellSelect.append(option(row.cell_label, row.cell_label));
  }
  if (!item.rows.some(row => row.cell_label === state.chartCell)) state.chartCell = item.maxCell;
  cellSelect.value = state.chartCell;
  cellLabel.append(cellSelect);
  const scaleLabel = el("label", "Expression scale");
  const scaleSelect = el("select");
  scaleSelect.id = `${state.id}-scale`;
  scaleSelect.append(option("log", "log1p(median CPM)"), option("z", "Z-standardized log1p"));
  scaleSelect.value = state.scale;
  scaleLabel.append(scaleSelect);
  chartControls.append(cellLabel, scaleLabel);
  const chart = el("div");
  const chartHeading = el("h4", undefined, "icon-heading chart-heading");
  chartHeading.append(icon("analysis"), el("span", "Developmental expression"));
  state.detail.append(title, el("p", context, "scope"), metrics, chartHeading, chartControls, chart);
  cellSelect.addEventListener("change", () => {
    state.chartCell = cellSelect.value;
    renderChart(state, item, chart);
  });
  scaleSelect.addEventListener("change", () => {
    state.scale = scaleSelect.value;
    renderChart(state, item, chart);
  });
  renderChart(state, item, chart);
}

function svgElement(tag, attributes = {}, text) {
  const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, String(value));
  if (text !== undefined) node.textContent = text;
  return node;
}

function drawChart(series, state, gene) {
  const svg = svgElement("svg", { viewBox: "0 0 600 320", class: "chart", role: "img", "aria-label": `${gene}: human and mouse ${state.chartCell} trajectories, ${state.scale === "z" ? "z-standardized log expression" : "log1p median CPM"}` });
  svg.append(svgElement("title", {}, `${gene} human–mouse developmental trajectories`));
  const transformed = Object.fromEntries(["human", "mouse"].map(species => [species, transformTrajectory(series[species].median_cpm, state.scale)]));
  const values = [...transformed.human, ...transformed.mouse];
  let low = state.scale === "log" ? 0 : Math.min(...values);
  let high = Math.max(...values);
  const padding = (high - low || 1) * .12;
  high += padding;
  if (state.scale === "z") low -= padding;
  const x = i => 64 + i * 160;
  const y = value => 256 - (value - low) / (high - low) * 208;
  for (let i = 0; i <= 4; i++) {
    const value = low + (high - low) * i / 4;
    svg.append(svgElement("line", { x1: 64, x2: 544, y1: y(value), y2: y(value), class: "grid" }));
    svg.append(svgElement("text", { x: 54, y: y(value) + 4, "text-anchor": "end" }, value.toFixed(1)));
  }
  svg.append(svgElement("text", { x: 64, y: 22 }, state.scale === "z" ? "Z-standardized log1p(median CPM)" : "log1p(median CPM)"));
  phaseLabels.forEach((label, index) => {
    const text = svgElement("text", { x: x(index), y: 281, "text-anchor": "middle" });
    const words = label.split(" ");
    text.append(svgElement("tspan", { x: x(index) }, words[0]));
    if (words[1]) text.append(svgElement("tspan", { x: x(index), dy: 16 }, words[1]));
    svg.append(text);
  });
  for (const species of ["human", "mouse"]) {
    svg.append(svgElement("polyline", { points: transformed[species].map((value, i) => `${x(i)},${y(value)}`).join(" "), class: species }));
    transformed[species].forEach((value, i) => {
      const point = svgElement("circle", { cx: x(i), cy: y(value), r: 4, class: species });
      point.append(svgElement("title", {}, `${species}, ${phaseLabels[i]}: ${value.toFixed(4)}; median CPM ${series[species].median_cpm[i].toFixed(4)}; ${series[species].n_replicates[i]} replicates`));
      svg.append(point);
    });
  }
  return svg;
}

async function renderChart(state, item, container) {
  const cell = state.chartCell;
  const token = Symbol();
  container.renderToken = token;
  container.replaceChildren(el("p", "Loading phase-level expression…", "muted"));
  try {
    state.trajectoriesPromise ??= load(`${state.folder}/trajectories.json`, true).catch(error => {
      state.trajectoriesPromise = null;
      throw error;
    });
    const data = await state.trajectoriesPromise;
    if (!container.isConnected || container.renderToken !== token) return;
    if (data.phases.join("|") !== PHASES.join("|")) throw new Error("Trajectory data use unexpected developmental phases.");
    const series = data.trajectories[cell]?.[item.gene];
    for (const species of ["human", "mouse"]) {
      const values = series?.[species];
      if (!values || values.median_cpm.length !== 4 || values.n_replicates.length !== 4 ||
          !values.median_cpm.every(value => Number.isFinite(value) && value >= 0) ||
          !values.n_replicates.every(value => Number.isInteger(value) && value >= 1)) {
        throw new Error(`Incomplete phase-level data for ${item.gene} / ${cell} / ${species}.`);
      }
    }
    const row = item.rows.find(entry => entry.cell_label === cell);
    const legend = el("div", undefined, "legend");
    legend.append(el("span", "Human — solid", "human"), el("span", "Mouse – – dashed", "mouse"));
    const note = el("p", `Plotted cell type: ${cell}. Spearman rho ${format(row.spearman_rho)} · z-RMSE ${format(row.z_rmse)} · direction agreement ${format(row.direction_agreement)} · divergence ${format(row.descriptive_divergence_score)}.`, "chart-note");
    const wrapper = el("div", undefined, "replicates");
    const table = el("table");
    const caption = el("caption", `Biological replicates per phase: human donors / ${state.replicate}`);
    const head = el("thead");
    const labels = el("tr");
    for (const label of ["Species", ...phaseLabels]) {
      const th = el("th", label);
      th.scope = "col";
      labels.append(th);
    }
    head.append(labels);
    const body = el("tbody");
    for (const species of ["human", "mouse"]) {
      const tr = el("tr");
      const th = el("th", species === "human" ? "Human" : "Mouse");
      th.scope = "row";
      tr.append(th, ...series[species].n_replicates.map(value => el("td", String(value))));
      body.append(tr);
    }
    table.append(caption, head, body);
    wrapper.append(table);
    const chartScroll = el("div", undefined, "chart-scroll");
    chartScroll.tabIndex = 0;
    chartScroll.setAttribute("role", "region");
    chartScroll.setAttribute("aria-label", "Developmental trajectory chart; scroll horizontally on narrow screens");
    chartScroll.append(drawChart(series, state, item.gene));
    container.replaceChildren(chartScroll, legend, note, wrapper);
  } catch (error) {
    if (!container.isConnected || container.renderToken !== token) return;
    const retry = el("button", "Retry chart");
    retry.type = "button";
    retry.addEventListener("click", () => renderChart(state, item, container));
    container.replaceChildren(el("p", error.message, "error"), retry);
  }
}

try {
  const loaded = await Promise.all(configs.map(async config => {
    const text = await load(`${config.folder}/trajectory_similarity.csv`);
    const metrics = prepareMetrics(parseCSV(text));
    if (!metrics.length) throw new Error(`No eligible cell-type comparisons in ${config.title}.`);
    return { config, metrics };
  }));
  for (const { config, metrics } of loaded) regions.push(createRegion(config, metrics));
  controls.n.max = String(Math.max(...regions.map(state => new Set(state.metrics.map(row => row.gene)).size)));
  controls.n.disabled = false;
  controls.direction.disabled = false;
  controls.n.addEventListener("change", refresh);
  controls.direction.addEventListener("change", refresh);
  refresh();
} catch (error) {
  controls.overlap.textContent = "Data could not be loaded.";
  showError(error);
}
