import { parseCSV } from "./explorer-model.mjs";
import { prepareClusters, clusterPeers, geneTerms } from "./bulk-model.mjs";

const dataRoot = new URL("../data/ASD_Bulk_analysis/", import.meta.url);
const input = document.querySelector("#gene-search");
const status = document.querySelector("#bulk-status");
const errorBox = document.querySelector("#bulk-error");
const containers = Object.fromEntries(["HPO", "KEGG", "GO"].map(kind =>
  [kind, document.querySelector(`#${kind.toLowerCase()}-results`)]));
const cache = new Map();
let clusters;
let geneMap;
let symbolsById;
let selectionToken = 0;

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
    const names = {
      HPO: `expression_kmeans_HPO_enrichment_resultv2_Cluster${cluster}_communities.csv`,
      KEGG: `expression_kmeans_KEGG_enrichv2_Cluster${cluster}.csv`,
      GO: `expression_kmeans_GO_enrichment_resultv2_Cluster${cluster}.csv`,
    };
    cache.set(key, load(`enrichment/${kind}/${names[kind]}`).catch(error => {
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

function renderTerms(container, rows, kind, gene, cluster) {
  container.replaceChildren();
  if (!rows.length) {
    container.append(el("p", `No exported ${kind} enrichment rows contain ${gene} in cluster ${cluster}. This does not mean the gene lacks ${kind} annotations.`, "bulk-empty"));
    return;
  }
  const wrapper = el("div", undefined, "bulk-table-wrap");
  wrapper.tabIndex = 0;
  wrapper.setAttribute("role", "region");
  wrapper.setAttribute("aria-label", `${kind} enrichment for ${gene}; scroll to view all rows and columns`);
  const table = el("table");
  table.append(el("caption", `${rows.length} ${kind} ${rows.length === 1 ? "result" : "results"} involving ${gene} in cluster ${cluster}, ordered by adjusted p-value.`));
  const head = el("thead");
  const headings = el("tr");
  const extra = kind === "HPO" ? "HPO community" : kind === "GO" ? "Ontology" : "Category";
  for (const label of ["ID", kind === "KEGG" ? "Pathway" : "Term", extra, "Adjusted p-value", "Overlap genes"]) {
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
      row.geneID.split("/").map(id => symbolsById.get(id)?.join(" / ") ?? `Entrez ${id}`);
    tr.append(id, el("td", description), el("td", category), el("td", p.toExponential(3)), el("td", genes.join("; ")));
    body.append(tr);
  }
  table.append(head, body);
  wrapper.append(table);
  container.append(wrapper);
}

async function selectGene() {
  const token = ++selectionToken;
  const gene = input.value.trim().toUpperCase();
  const cluster = clusters.byGene.get(gene);
  errorBox.hidden = true;
  const peers = document.querySelector("#cluster-peers");
  peers.replaceChildren();
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
  const others = clusterPeers(clusters, gene);
  const summary = `${gene} | Cluster ${cluster} | ${others.length + 1} ${others.length ? "genes" : "gene"} in cluster.`;
  status.textContent = `${summary} Loading enrichment...`;
  peers.append(el("p", `${gene} is in cluster ${cluster}. ${others.length} other ${others.length === 1 ? "gene shares" : "genes share"} this cluster.`));
  if (!others.length) peers.append(el("p", "This is a singleton cluster; there are no other genes.", "bulk-empty"));
  const list = el("div", undefined, "peer-list");
  for (const peer of others) {
    const button = el("button", peer, "peer-button");
    button.type = "button";
    button.addEventListener("click", () => {
      input.value = peer;
      selectGene();
    });
    list.append(button);
  }
  peers.append(list);
  let failures = 0;
  await Promise.all(Object.entries(containers).map(async ([kind, container]) => {
    container.append(el("p", `Loading ${kind} enrichment...`, "bulk-empty"));
    if (kind === "HPO" && cluster === 10) {
      container.replaceChildren(el("p", "HPO source data are unavailable for cluster 10: no HPO result file was exported. No terms can be displayed.", "bulk-empty"));
      return;
    }
    try {
      const records = await enrichment(cluster, kind);
      if (token !== selectionToken) return;
      renderTerms(container, geneTerms(records, kind, gene, geneMap), kind, gene, cluster);
    } catch (error) {
      if (token !== selectionToken) return;
      failures++;
      const retry = el("button", `Retry ${kind} loading`);
      retry.type = "button";
      retry.addEventListener("click", selectGene);
      container.replaceChildren(el("p", error.message, "error"), retry);
    }
  }));
  if (token === selectionToken) status.textContent = `${summary}${failures ? ` ${failures} enrichment sources failed; see errors below.` : " Enrichment loaded."}`;
}

document.querySelector("#gene-form").addEventListener("submit", event => {
  event.preventDefault();
  selectGene();
});

try {
  const [assignments, mapping] = await Promise.all([
    load("cluster_result/expression_kmeans_all_clusters_genesv2.csv"), load("gene-entrez-map.json", true),
  ]);
  clusters = prepareClusters(assignments);
  geneMap = mapping.genes;
  symbolsById = new Map();
  const options = document.querySelector("#gene-options");
  for (const gene of [...clusters.byGene.keys()].sort()) {
    const ids = geneMap?.[gene];
    if (!Array.isArray(ids) || !ids.length || !ids.every(id => /^\d+$/.test(id))) throw new Error(`Missing or invalid Entrez mapping for ${gene}.`);
    for (const id of ids) {
      if (!symbolsById.has(id)) symbolsById.set(id, []);
      symbolsById.get(id).push(gene);
    }
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
