import { parseCSV } from "./explorer-model.mjs";
import { prepareClusters, clusterTerms, clusterTrend } from "./bulk-model.mjs";

const dataRoot = new URL("../data/ASD_Bulk_analysis/", import.meta.url);
const input = document.querySelector("#gene-search");
const status = document.querySelector("#bulk-status");
const errorBox = document.querySelector("#bulk-error");
const containers = Object.fromEntries(["HPO", "KEGG", "GO"].map(kind =>
  [kind, document.querySelector(`#${kind.toLowerCase()}-results`)]));
const cache = new Map();
let clusters;
let bulkData;
let symbolsById;
let selectionToken = 0;
let onsetPromise = null;

function el(tag, text, className) {
  const element = document.createElement(tag);
  if (text !== undefined) element.textContent = text;
  if (className) element.className = className;
  return element;
}

async function load(file, json = false) {
  const response = await fetch(new URL(file, dataRoot), { cache: "no-cache" });
  if (!response.ok) throw new Error(`Could not load ${file} (HTTP ${response.status}).`);
  return json ? response.json() : parseCSV(await response.text());
}

function enrichment(cluster, kind) {
  const key = `${cluster}-${kind}`;
  if (!cache.has(key)) {
    cache.set(key, load(bulkData.enrichment[cluster][kind]).catch(error => {
      cache.delete(key);
      throw error;
    }));
  }
  return cache.get(key);
}

function showError(error) {
  errorBox.textContent = error.message;
  if (location.protocol === "file:") errorBox.textContent += " Serve the docs folder over HTTP, not a file URL.";
  errorBox.hidden = false;
}

function termLink(row, kind) {
  const id = kind === "HPO" ? row.Term.match(/HP:\d{7}/)?.[0] : row.ID;
  const valid = kind === "HPO" ? /^HP:\d{7}$/.test(id) :
    kind === "KEGG" ? /^hsa\d{5}$/.test(id) : /^GO:\d{7}$/.test(id);
  if (!valid) throw new Error(`Invalid ${kind} term identifier.`);
  const link = el("a", id);
  link.href = kind === "HPO" ? `https://hpo.jax.org/browse/term/${id}` :
    kind === "KEGG" ? `https://www.kegg.jp/pathway/${id}` : `https://amigo.geneontology.org/amigo/term/${id}`;
  link.target = "_blank";
  link.rel = "noopener noreferrer";
  return link;
}

function onsetCells(id, context) {
  const term = context.terms?.[id];
  if (!term) return [
    el("td", "No bundled Orphanet onset context for this HPO term.", "onset-cell"),
    el("td", "Unknown: no bundled onset evidence; bins are not inferred.", "onset-cell"),
  ];
  const onset = el("td", undefined, "onset-cell");
  onset.append(el("p", `Most reported disease-onset categories: ${term.most_reported_categories.join(", ") || "Unknown"}. Not the phenotype's onset.`, "scope"));
  onset.append(el("p", `${term.entities_with_onset}/${term.associated_entities} associated clinical entities have onset data; ${term.entities_without_onset} have none.`));
  const distribution = el("ul");
  for (const item of term.distribution) distribution.append(el("li", `${item.category}: ${item.count} entities`));
  onset.append(distribution);
  const sources = el("details");
  sources.append(el("summary", "Example associated diseases and sources"));
  for (const example of term.examples) {
    const paragraph = el("p");
    const link = el("a", `ORPHA:${example.orpha_code} - ${example.name}`);
    link.href = example.source_url;
    link.target = "_blank";
    link.rel = "noopener noreferrer";
    paragraph.append(link);
    sources.append(paragraph);
  }
  const download = el("a", "All per-term disease/onset evidence (JSON)");
  download.href = new URL("orphanet-hpo-evidence.json", dataRoot).href;
  sources.append(download);
  onset.append(sources);
  const alignment = el("td", undefined, "onset-cell");
  const possible = new Set(term.distribution.flatMap(item => item.bins));
  const strip = el("div", undefined, "event-bin-strip");
  strip.setAttribute("aria-label", "Possible interval overlaps, not a peak or probability");
  for (let i = 1; i <= 5; i++) {
    const span = el("span", `Bin ${i}`, possible.has(`event_bin_${i}`) ? "bin-overlap" : "");
    span.title = `${context.boundaries_pc_days[i - 1].toFixed(1)}-${context.boundaries_pc_days[i].toFixed(1)} postconception days${possible.has(`event_bin_${i}`) ? "; possible interval overlap" : "; no overlap from reported categories"}`;
    strip.append(span);
  }
  alignment.append(strip);
  const intervals = el("ul");
  for (const item of term.distribution) {
    const flags = [
      item.before_window ? "also before bin 1" : "",
      item.after_window ? "also after bin 5" : "",
    ].filter(Boolean);
    intervals.append(el("li", `${item.category}: ${item.bins.map(bin => bin.replace("event_bin_", "bin ")).join(", ") || "no covered bins"}${flags.length ? `; ${flags.join("; ")}` : ""}.`));
  }
  if (!term.distribution.length) alignment.append(el("p", "Unknown: no disease-onset data available."));
  alignment.append(intervals);
  return [onset, alignment];
}

function renderTerms(container, rows, kind, cluster, context = null) {
  if (!rows.length) {
    container.replaceChildren(el("p", `No ${kind} enrichment rows were exported for cluster ${cluster}. This does not mean its genes lack ${kind} annotations.`, "bulk-empty"));
    return;
  }
  const wrapper = el("div", undefined, "bulk-table-wrap");
  wrapper.tabIndex = 0;
  wrapper.setAttribute("role", "region");
  wrapper.setAttribute("aria-label", `${kind} enrichment for cluster ${cluster}; scroll to view all rows and columns`);
  const table = el("table");
  table.append(el("caption", `${rows.length} ${kind} ${rows.length === 1 ? "result" : "results"} for cluster ${cluster}, ordered by adjusted p-value. Not filtered by selected-gene membership.`));
  const head = el("thead");
  const headings = el("tr");
  const extra = kind === "HPO" ? "HPO community" : kind === "GO" ? "Ontology" : "Category";
  const labels = ["ID", kind === "KEGG" ? "Pathway" : "Term", extra, "Adjusted p-value", "Overlap genes"];
  if (context) labels.push("Associated-disease onset (indirect)", "Possible cerebral bins (indirect)");
  for (const label of labels) {
    const th = el("th", label);
    th.scope = "col";
    headings.append(th);
  }
  head.append(headings);
  const body = el("tbody");
  for (const row of rows) {
    const tr = el("tr");
    const id = el("td");
    id.append(termLink(row, kind));
    const description = kind === "HPO" ? row.Term : row.Description;
    const category = kind === "HPO" ? row.community || "Not provided" :
      kind === "GO" ? row.ONTOLOGY || "Not provided" : row.category || "Not provided";
    const p = Number(row[kind === "HPO" ? "Adjusted.P.value" : "p.adjust"]);
    const genes = kind === "HPO" ? row.Genes.split(";").map(g => g.trim()) :
      row.geneID.split("/").map(id => symbolsById.get(id)?.length ? symbolsById.get(id).join(" / ") : `Entrez ${id}`);
    tr.append(id, el("td", description), el("td", category), el("td", p.toExponential(3)), el("td", genes.join("; ")));
    if (context) tr.append(...onsetCells(row.Term.match(/HP:\d{7}/)?.[0], context));
    body.append(tr);
  }
  table.append(head, body);
  wrapper.append(table);
  container.replaceChildren(wrapper);
  if (context) container.prepend(el("p", "Orphanet disease-onset context, not phenotype onset. Counts are unweighted clinical entities, not patients; one entity may contribute to several categories. Shaded bins indicate any interval overlap, not the most likely phase.", "scope"));
}

function svgElement(tag, attributes = {}, text) {
  const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, value);
  if (text !== undefined) node.textContent = text;
  return node;
}

function renderTrend(container, rows, cluster) {
  const x = i => 100 + i * 144;
  const lower = Math.min(0, ...rows.map(row => row.mean - row.se));
  const upper = Math.max(0, ...rows.map(row => row.mean + row.se));
  const span = upper - lower || 1;
  const min = lower - span * .05;
  const max = upper + span * .05;
  const y = value => 280 - (value - min) / (max - min) * 230;
  const svg = svgElement("svg", {
    viewBox: "0 0 900 390", class: "trend-chart", role: "img",
    "aria-label": `Cluster ${cluster}: mean expression across developmental phases, with mean plus or minus one standard error`,
  });
  svg.append(svgElement("title", {}, `Cluster ${cluster} expression trend`));
  svg.append(svgElement("text", { x: 100, y: 25 }, "Mean expression (exported units)"));
  for (let i = 0; i <= 4; i++) {
    const value = min + (max - min) * i / 4;
    svg.append(svgElement("line", { x1: 100, x2: 820, y1: y(value), y2: y(value), class: "trend-grid" }));
    svg.append(svgElement("text", { x: 88, y: y(value) + 4, "text-anchor": "end" }, Number(value.toPrecision(3)).toString()));
  }
  const band = [
    ...rows.map((row, i) => `${x(i)},${y(row.mean + row.se)}`),
    ...rows.map((row, i) => `${x(i)},${y(row.mean - row.se)}`).reverse(),
  ];
  svg.append(svgElement("polygon", { points: band.join(" "), class: "trend-band" }));
  svg.append(svgElement("polyline", { points: rows.map((row, i) => `${x(i)},${y(row.mean)}`).join(" "), class: "trend-line" }));
  rows.forEach((row, i) => {
    const point = svgElement("circle", { cx: x(i), cy: y(row.mean), r: 5, class: "trend-point" });
    point.append(svgElement("title", {}, `${row.phase}: mean ${row.mean}; SE ${row.se}; N ${row.n}`));
    svg.append(point);
    const label = svgElement("text", { x: x(i), y: 310, "text-anchor": "middle" });
    row.phase.replaceAll("_", " ").split(" ").forEach((word, j) => {
      label.append(svgElement("tspan", { x: x(i), dy: j ? 18 : 0 }, word));
    });
    svg.append(label);
  });
  svg.append(svgElement("text", { x: 460, y: 375, "text-anchor": "middle" }, "Developmental phase (ordered categories, not proportional elapsed time)"));
  const chart = el("div", undefined, "trend-chart-wrap");
  chart.append(svg);
  const details = el("details");
  details.append(el("summary", "View expression trend values"));
  const table = el("table");
  table.append(el("caption", `Exported cluster ${cluster} expression summaries`));
  const head = el("thead");
  const headings = el("tr");
  for (const label of ["Developmental phase", "Mean expression", "SE", "N"]) {
    const th = el("th", label);
    th.scope = "col";
    headings.append(th);
  }
  head.append(headings);
  const body = el("tbody");
  for (const row of rows) {
    const tr = el("tr");
    for (const value of [row.phase, row.mean, row.se, row.n]) tr.append(el("td", value));
    body.append(tr);
  }
  table.append(head, body);
  const wrapper = el("div", undefined, "bulk-table-wrap");
  wrapper.append(table);
  details.append(wrapper);
  const download = el("a", "Download cluster expression trends (CSV)");
  download.href = new URL(bulkData.trends, dataRoot).href;
  container.replaceChildren(
    el("p", `Cluster ${cluster}: exported mean expression; shaded band = mean ± SE. These are cluster averages, not the selected gene's expression.`, "hint"),
    chart, details, download,
  );
}

async function showTrend(container, cluster, token) {
  container.replaceChildren(el("p", "Loading cluster expression trend...", "bulk-empty"));
  try {
    if (!cache.has("trends")) {
      cache.set("trends", load(bulkData.trends).catch(error => {
        cache.delete("trends");
        throw error;
      }));
    }
    const records = await cache.get("trends");
    if (token !== selectionToken) return 0;
    renderTrend(container, clusterTrend(records, cluster), cluster);
    return 0;
  } catch (error) {
    if (token !== selectionToken) return 0;
    const retry = el("button", "Retry expression trend");
    retry.type = "button";
    retry.addEventListener("click", selectGene);
    container.replaceChildren(el("p", error.message, "error"), retry);
    return 1;
  }
}

async function selectGene() {
  const token = ++selectionToken;
  const requested = input.value.trim();
  const gene = [...clusters.byGene.keys()].find(symbol => symbol.toUpperCase() === requested.toUpperCase()) ?? requested;
  const cluster = clusters.byGene.get(gene);
  errorBox.hidden = true;
  const peers = document.querySelector("#cluster-peers");
  const trend = document.querySelector("#cluster-trend");
  peers.replaceChildren();
  trend.replaceChildren();
  for (const container of Object.values(containers)) container.replaceChildren();
  if (!cluster) {
    status.textContent = "Gene not found in this dataset.";
    showError(new Error(`No expression/K-means v2 cluster assignment for "${gene}". Choose a gene from the suggestions.`));
    return;
  }
  input.value = gene;
  const url = new URL(location.href);
  url.searchParams.set("gene", gene);
  history.replaceState(null, "", url);
  const members = clusters.byCluster.get(cluster);
  const summary = `${gene} | Cluster ${cluster} | ${members.length} ${members.length === 1 ? "gene" : "genes"} in cluster.`;
  status.textContent = `${summary} Loading cluster data...`;
  peers.append(el("p", `${gene} belongs to cluster ${cluster}. All ${members.length} cluster members are listed below, including the selected gene.`));
  const list = el("div", undefined, "peer-list");
  list.tabIndex = 0;
  list.setAttribute("role", "region");
  list.setAttribute("aria-label", `All ${members.length} genes in cluster ${cluster}; scroll to view all members`);
  for (const peer of members) {
    const button = el("button", peer, "peer-button");
    button.type = "button";
    button.setAttribute("aria-pressed", String(peer === gene));
    button.addEventListener("click", () => {
      input.value = peer;
      selectGene();
    });
    list.append(button);
  }
  peers.append(list);
  let failures = 0;
  const results = await Promise.all([showTrend(trend, cluster, token), ...Object.entries(containers).map(async ([kind, container]) => {
    container.append(el("p", `Loading ${kind} enrichment...`, "bulk-empty"));
    if (bulkData.enrichment[cluster][kind] === null) {
      container.replaceChildren(el("p", `${kind} source data are unavailable for cluster ${cluster}: no ${kind} result file was supplied. No terms can be displayed.`, "bulk-empty"));
      return;
    }
    try {
      const records = await enrichment(cluster, kind);
      if (token !== selectionToken) return;
      const rows = clusterTerms(records, kind);
      renderTerms(container, rows, kind, cluster);
      if (kind === "HPO" && rows.length) {
        try {
          onsetPromise ??= load("orphanet-hpo-context.json", true).catch(error => {
            onsetPromise = null;
            throw error;
          });
          const context = await onsetPromise;
          if (token !== selectionToken) return;
          renderTerms(container, rows, kind, cluster, context);
        } catch (error) {
          if (token !== selectionToken) return;
          failures++;
          const retry = el("button", "Retry Orphanet onset context");
          retry.type = "button";
          retry.addEventListener("click", selectGene);
          container.append(el("p", `Orphanet onset context could not be loaded: ${error.message}`, "error"), retry);
        }
      }
    } catch (error) {
      if (token !== selectionToken) return;
      failures++;
      const retry = el("button", `Retry ${kind} loading`);
      retry.type = "button";
      retry.addEventListener("click", selectGene);
      container.replaceChildren(el("p", error.message, "error"), retry);
    }
  })]);
  failures += results[0];
  if (token === selectionToken) status.textContent = `${summary}${failures ? ` ${failures} data sources failed; see errors below.` : " Cluster data loaded."}`;
}

document.querySelector("#gene-form").addEventListener("submit", event => {
  event.preventDefault();
  selectGene();
});

try {
  const [data, mapping] = await Promise.all([
    load("bulk-data.json", true), load("gene-entrez-map.json", true),
  ]);
  bulkData = data;
  clusters = prepareClusters(data.assignments);
  if (!mapping.symbols_by_id || typeof mapping.symbols_by_id !== "object") throw new Error("Missing Entrez-to-symbol mapping. Regenerate the bulk gene map.");
  symbolsById = new Map(Object.entries(mapping.symbols_by_id));
  const options = document.querySelector("#gene-options");
  for (const gene of [...clusters.byGene.keys()].sort()) {
    const option = el("option");
    option.value = gene;
    options.append(option);
  }
  input.disabled = false;
  document.querySelector("#show-gene").disabled = false;
  input.value = new URLSearchParams(location.search).get("gene") ?? "ACTB";
  await selectGene();
} catch (error) {
  status.textContent = "Bulk data could not be loaded.";
  showError(error);
}
