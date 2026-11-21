import { parseCSV, resolveGene } from "./explorer-model.mjs";
import { STAGES, REGION_FIELDS, selectGeneData, prepareSamples, profilePoints, filterPoints, smoothProfile, driverNetwork, toCSV } from "./brainspan-model.mjs?v=2";

const $ = id => document.getElementById(id);
const root = new URL("../data/BrainSpan_NDD/", import.meta.url);
const ns = "http://www.w3.org/2000/svg";
const cache = new Map();
let catalog, samples, drivers, applied = null, plotted = [], edges = [], request = 0;

function element(tag, attributes = {}, text) {
  const node = document.createElementNS(ns, tag);
  for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, value);
  if (text !== undefined) node.textContent = text;
  return node;
}
function option(value, label = value) {
  const node = document.createElement("option");
  node.value = value;
  node.textContent = label;
  return node;
}
function showError(error) {
  $("error").textContent = error.message;
  $("error").hidden = false;
}
function clearError() { $("error").hidden = true; }
async function csv(file) {
  const response = await fetch(new URL(file, root), { cache: "no-cache" });
  if (!response.ok) throw new Error(`Could not load ${file}: HTTP ${response.status}.`);
  return parseCSV(await response.text());
}
function currentGene() { return resolveGene(new Set(catalog.keys()), $("gene").value); }
function updateDrivers() {
  const gene = currentGene();
  const previous = $("driver").value;
  const choices = drivers.filter(row => row.gene === gene);
  $("driver").replaceChildren(...choices.map(row => option(row.driver_id, `${row.driver_id} (${row.function_type})`)));
  if (choices.some(row => row.driver_id === previous)) $("driver").value = previous;
  else if (choices.some(row => row.driver_id === `${gene}_TF`)) $("driver").value = `${gene}_TF`;
  $("driver").disabled = !choices.length;
}
function setBusy(busy) {
  for (const id of ["gene", "profile", "draw"]) $(id).disabled = busy;
  $("driver").disabled = busy || !$("driver").options.length;
  $("filters").disabled = busy || !applied;
}
function clearResults() {
  applied = null;
  plotted = [];
  edges = [];
  for (const id of ["pattern", "legend", "network", "edges"]) $(id).replaceChildren();
  for (const id of ["pattern-svg", "pattern-csv", "network-svg", "network-csv", "zoom-in", "zoom-out", "reset-network", "node"]) $(id).disabled = true;
  $("node").replaceChildren();
}

async function draw(event) {
  event?.preventDefault();
  clearError();
  const gene = currentGene();
  if (!gene) return showError(new Error("Unknown NDD gene. Choose a symbol from the catalog."));
  $("gene").value = gene;
  const choice = { gene, driver: $("driver").value, profile: $("profile").value };
  const manifest = catalog.get(gene);
  const token = ++request;
  clearResults();
  setBusy(true);
  $("status").textContent = `Loading ${gene} CSVs…`;
  try {
    let data = cache.get(manifest.file);
    if (!data) {
      if (manifest.file) {
        const [profiles, network] = await Promise.all([csv(`genes/${manifest.file}-profiles.csv`), csv(`genes/${manifest.file}-network.csv`)]);
        data = { profiles, network };
      } else data = { profiles: [], network: [] };
      if (manifest.file) cache.set(manifest.file, data);
    }
    if (token !== request) return;
    applied = { ...choice, ...selectGeneData(data, gene) };
    renderPattern();
    renderNetwork();
    $("status").textContent = `${gene} · ${manifest.sources} · ${manifest.expression === "TRUE" ? "expression available" : "no expression"} · ${manifest.activity_drivers} activity drivers · ${manifest.network_edges} NDD-only edges across all drivers.`;
  } catch (error) {
    clearResults();
    $("status").textContent = `Could not display ${gene}.`;
    showError(error);
  } finally {
    if (token === request) setBusy(false);
  }
}

function renderPattern() {
  $("pattern").replaceChildren();
  $("legend").replaceChildren();
  $("pattern-note").textContent = "";
  plotted = [];
  $("pattern-svg").disabled = true;
  $("pattern-csv").disabled = true;
  if (!applied) return;
  const id = applied.profile === "Expression" ? applied.gene : applied.driver;
  $("pattern-title").textContent = `${applied.profile} for ${id || applied.gene}`;
  const points = profilePoints(applied.profiles, samples, applied.profile, id);
  if (!points.length) {
    $("pattern-note").textContent = `No ${applied.profile.toLowerCase()} data available for ${id || applied.gene}.`;
    return;
  }
  const start = Number($("start").value), end = Number($("end").value);
  const regions = Object.fromEntries(REGION_FIELDS.map(field => [field, [...$(field).selectedOptions].map(o => o.value)]));
  plotted = filterPoints(points, { start, end, regions });
  if (!plotted.length) {
    $("pattern-note").textContent = "No samples match the selected regions and age range.";
    return;
  }
  const field = $("group").value;
  const groups = new Map();
  for (const point of plotted) {
    const key = field === "none" ? "All samples" : point[field];
    if (!groups.has(key)) groups.set(key, []);
    groups.get(key).push(point);
  }
  const series = [...groups].sort(([a], [b]) => a.localeCompare(b)).map(([name, points]) => ({ name, points, smooth: smoothProfile(points, Number($("span").value)) }));
  const confidence = $("confidence").checked;
  const values = plotted.map(p => p.value);
  for (const item of series) for (const p of item.smooth) {
    values.push(p.value);
    if (confidence && p.low !== null) values.push(p.low, p.high);
  }
  let low = Math.min(...values), high = Math.max(...values);
  const padding = (high - low || Math.abs(high) || 1) * .1;
  low -= padding; high += padding;
  const width = 960, height = 580, left = 90, right = 28, top = 125, bottom = 95;
  const xmin = Math.log10(Math.min(...plotted.map(p => p.Days))), xmax = Math.log10(Math.max(...plotted.map(p => p.Days)));
  const x0 = xmin === xmax ? xmin - .05 : xmin, x1 = xmin === xmax ? xmax + .05 : xmax;
  const x = days => left + (Math.log10(days) - x0) / (x1 - x0) * (width - left - right);
  const y = value => top + (high - value) / (high - low) * (height - top - bottom);
  const svg = element("svg", { viewBox: `0 0 ${width} ${height}`, class: "chart brainspan-chart", role: "img", "aria-label": `${applied.profile} for ${id}: ${plotted.length} observations versus age` });
  svg.append(element("title", {}, `${applied.profile} for ${id}`));
  for (let i = 0; i <= 5; i++) {
    const value = low + (high - low) * i / 5;
    svg.append(element("line", { x1: left, x2: width - right, y1: y(value), y2: y(value), class: "grid" }),
      element("text", { x: left - 10, y: y(value) + 6, "text-anchor": "end" }, Number(value.toPrecision(3)).toString()));
  }
  const ticks = [[56, "8pcw"], [63, "9pcw"], [84, "12pcw"], [91, "13pcw"], [105, "15pcw"], [175, "25pcw"], [266, "38pcw"], [645, "1yr"], [1375, "3yr"], [4300, "11yr"], [7220, "19yr"], [14900, "40yr"]];
  // Dense source age ticks are thinned only for readability, not the data.
  let previousX = -Infinity;
  for (const [days, label] of ticks.filter(([days]) => Math.log10(days) >= x0 && Math.log10(days) <= x1)) {
    if (x(days) - previousX < 40) continue;
    previousX = x(days);
    svg.append(element("text", { transform: `translate(${x(days)},${height - bottom + 22}) rotate(45)`, "text-anchor": "start" }, label));
  }
  for (const [a, b, label] of STAGES[$("style").value]) {
    if (a < start || b > end || Math.log10(b) < x0 || Math.log10(a) > x1) continue;
    if (Math.log10(b) <= x1) svg.append(element("line", { x1: x(b), x2: x(b), y1: top, y2: height - bottom, class: "boundary" }));
    const middle = Math.sqrt(a * b);
    if (Math.log10(middle) >= x0 && Math.log10(middle) <= x1) svg.append(element("text", { transform: `translate(${x(middle)},${top - 8}) rotate(-45)`, class: "stage-label" }, label));
  }
  if (Math.log10(266) >= x0 && Math.log10(266) <= x1) svg.append(element("line", { x1: x(266), x2: x(266), y1: top, y2: height - bottom, class: "birth" }));
  svg.append(element("line", { x1: left, x2: width - right, y1: height - bottom, y2: height - bottom, class: "axis" }));
  const path = points => points.map((p, i) => `${i ? "L" : "M"}${x(p.Days)},${y(p.value)}`).join(" ");
  series.forEach((item, index) => {
    const group = element("g", { class: `series-${index % 6}` });
    const band = item.smooth.filter(p => p.low !== null);
    if (confidence && band.length > 1) group.append(element("path", { class: "band", d: path(band.map(p => ({ ...p, value: p.high }))) + " " + path([...band].reverse().map(p => ({ ...p, value: p.low }))).replace(/^M/, "L") + " Z" }));
    for (const point of item.points) {
      const mark = field === "none" ? element("circle", { cx: x(point.Days), cy: y(point.value), r: 3.5, class: "point" }) :
        element("text", { x: x(point.Days), y: y(point.value) + 4, "text-anchor": "middle", style: "fill:currentColor;font-size:12px", class: "point" }, String(index + 1));
      mark.append(element("title", {}, `${point.structure_acronym} · ${point.age} · ${point.gender} · donor ${point.donor_id} · sample ${point.sample_id} · ${point.value}`));
      group.append(mark);
    }
    if (item.smooth.length) group.append(element("path", { class: "smooth", d: path(item.smooth), "stroke-dasharray": index < 6 ? "none" : `${3 + Math.floor(index / 6) * 2} 4` }));
    svg.append(group);
    const label = document.createElement("span");
    label.textContent = `${field === "none" ? "" : `${index + 1} · `}${item.name} (${item.points.length})`;
    $("legend").append(label);
  });
  svg.append(element("text", { x: (left + width - right) / 2, y: height - 8, "text-anchor": "middle" }, "Age (logarithmic; pcw = postconception weeks)"),
    element("text", { transform: `translate(24,${(top + height - bottom) / 2}) rotate(-90)`, "text-anchor": "middle" }, `${applied.profile} level`));
  $("pattern").append(svg);
  const unsmoothed = series.filter(item => !item.smooth.length).length;
  const missingBands = series.some(item => item.smooth.some(point => point.low === null));
  $("pattern-note").textContent = `${plotted.length} samples from ${new Set(plotted.map(p => p.donor_id)).size} donors. ${unsmoothed ? `${unsmoothed} region group(s) have insufficient ages/samples for smoothing. ` : ""}${confidence && missingBands ? "Some confidence bands are not estimable. " : ""}Hover over an observation for its sample details.`;
  $("pattern-svg").disabled = false;
  $("pattern-csv").disabled = false;
}

function renderNetwork() {
  edges = applied.driver ? driverNetwork(applied.network, applied.driver, new Set(catalog.keys())) : [];
  $("network-title").textContent = `Network for ${applied.driver || applied.gene}`;
  if (!edges.length) {
    $("network-note").textContent = applied.driver ? "No NDD-to-NDD network edges are available for this driver." : "No driver is available for this gene.";
    return;
  }
  const ids = [...new Set(edges.flatMap(row => [row.source, row.target]))];
  $("network-note").textContent = `${ids.length} nodes · ${edges.length} NDD-only directed edges. Non-NDD targets are excluded.`;
  $("node").replaceChildren(option("", "No highlight"), ...ids.map(id => option(id)));
  const svg = element("svg", { viewBox: "0 0 960 600", class: "chart brainspan-chart network-chart", role: "group", "aria-label": `Interactive NDD network for ${applied.driver}` });
  const defs = element("defs"), marker = element("marker", { id: "network-arrow", viewBox: "0 0 10 10", refX: 9, refY: 5, markerWidth: 5, markerHeight: 5, orient: "auto" });
  marker.append(element("path", { d: "M0 0 L10 5 L0 10 Z" }));
  defs.append(marker); svg.append(defs);
  const viewport = element("g");
  svg.append(viewport);
  const targets = ids.filter(id => id !== applied.driver);
  const positions = new Map(ids.map(id => {
    const index = targets.indexOf(id);
    const angle = index * 2 * Math.PI / (targets.length || 1);
    return [id, id === applied.driver ? { x: 480, y: 300 } : { x: 480 + 340 * Math.cos(angle), y: 300 + 220 * Math.sin(angle) }];
  }));
  const maxMI = Math.max(...edges.map(row => row.MI)) || 1;
  const edgeNodes = edges.map(row => {
    const path = element("path", { class: "edge", "stroke-width": 1 + 4 * row.MI / maxMI, "marker-end": "url(#network-arrow)" });
    path.append(element("title", {}, `${row.source} → ${row.target} · MI: ${row.MI} · Pearson: ${row.pearson ?? "NA"} · Spearman: ${row.spearman ?? "NA"}`));
    viewport.append(path);
    return { row, path };
  });
  const nodeNodes = ids.map(id => {
    const source = id === applied.driver;
    const node = element("g", { class: `node${source ? " source" : ""}`, tabindex: 0, role: "button", "aria-label": `Highlight ${id}` });
    node.append(element("circle", { r: source ? 22 : 12 }),
      element("text", { y: source ? -32 : -20, "text-anchor": "middle" }, id),
      element("title", {}, id));
    viewport.append(node);
    return { id, node };
  });
  function highlight(id) {
    const adjacent = new Set([id]);
    for (const { row, path } of edgeNodes) {
      const selected = row.source === id || row.target === id;
      path.classList.toggle("highlight", selected);
      if (selected) { adjacent.add(row.source); adjacent.add(row.target); }
    }
    for (const item of nodeNodes) item.node.classList.toggle("highlight", Boolean(id) && adjacent.has(item.id));
  }
  let scale = 1, tx = 0, ty = 0;
  const view = () => viewport.setAttribute("transform", `translate(${tx},${ty}) scale(${scale})`);
  function layout() {
    for (const { id, node } of nodeNodes) {
      const p = positions.get(id);
      node.setAttribute("transform", `translate(${p.x},${p.y})`);
    }
    for (const { row, path } of edgeNodes) {
      const a = positions.get(row.source), b = positions.get(row.target);
      const dx = b.x - a.x, dy = b.y - a.y, distance = Math.hypot(dx, dy);
      path.setAttribute("d", distance ? `M${a.x + dx * 24 / distance},${a.y + dy * 24 / distance} L${b.x - dx * 16 / distance},${b.y - dy * 16 / distance}` : `M${a.x + 15},${a.y - 15} C${a.x + 80},${a.y - 90} ${a.x - 80},${a.y - 90} ${a.x - 15},${a.y - 15}`);
    }
  }
  function location(event) {
    const point = svg.createSVGPoint();
    point.x = event.clientX; point.y = event.clientY;
    return point.matrixTransform(svg.getScreenCTM().inverse());
  }
  let drag = null;
  svg.addEventListener("pointerdown", event => {
    const item = nodeNodes.find(item => item.node.contains(event.target));
    drag = { id: item?.id, origin: location(event), tx, ty };
    svg.setPointerCapture(event.pointerId);
    if (item) highlight(item.id);
  });
  svg.addEventListener("pointermove", event => {
    if (!drag) return;
    const p = location(event);
    if (drag.id) {
      positions.set(drag.id, { x: (p.x - tx) / scale, y: (p.y - ty) / scale });
      layout();
    } else { tx = drag.tx + p.x - drag.origin.x; ty = drag.ty + p.y - drag.origin.y; view(); }
  });
  for (const type of ["pointerup", "pointercancel", "lostpointercapture"]) svg.addEventListener(type, () => { drag = null; });
  for (const { id, node } of nodeNodes) {
    node.addEventListener("mouseenter", () => highlight(id));
    node.addEventListener("mouseleave", () => { if (!drag) highlight($("node").value); });
    node.addEventListener("focus", () => highlight(id));
    node.addEventListener("blur", () => highlight($("node").value));
    node.addEventListener("keydown", event => {
      if (event.key === "Enter" || event.key === " ") { event.preventDefault(); $("node").value = id; highlight(id); }
    });
  }
  $("node").onchange = () => highlight($("node").value);
  function zoom(factor) {
    const next = Math.max(.25, Math.min(4, scale * factor));
    tx = 480 - (480 - tx) * next / scale;
    ty = 300 - (300 - ty) * next / scale;
    scale = next; view();
  }
  $("zoom-in").onclick = () => zoom(1.25);
  $("zoom-out").onclick = () => zoom(.8);
  $("reset-network").onclick = () => { scale = 1; tx = 0; ty = 0; view(); };
  layout(); view();
  $("network").append(svg);
  const table = document.createElement("table"), head = document.createElement("thead"), body = document.createElement("tbody");
  const columns = ["source", "target", "MI", "pearson", "spearman"];
  const heading = document.createElement("tr");
  for (const column of columns) { const cell = document.createElement("th"); cell.scope = "col"; cell.textContent = column; heading.append(cell); }
  head.append(heading);
  for (const row of edges) {
    const tr = document.createElement("tr");
    for (const column of columns) { const cell = document.createElement("td"); cell.textContent = row[column] ?? "NA"; tr.append(cell); }
    body.append(tr);
  }
  table.append(head, body); $("edges").append(table);
  for (const id of ["network-svg", "network-csv", "zoom-in", "zoom-out", "reset-network", "node"]) $(id).disabled = false;
}

function download(text, type, name) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const link = document.createElement("a");
  link.href = url; link.download = name;
  document.body.append(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function saveSVG(container, name) {
  const original = $(container).querySelector("svg");
  if (!original) return;
  const clone = original.cloneNode(true);
  const originals = [original, ...original.querySelectorAll("*")], copies = [clone, ...clone.querySelectorAll("*")];
  originals.forEach((node, i) => {
    const style = getComputedStyle(node);
    for (const key of ["color", "fill", "stroke", "stroke-width", "stroke-dasharray", "opacity", "font-family", "font-size", "font-weight"]) copies[i].style.setProperty(key, style.getPropertyValue(key));
  });
  clone.setAttribute("xmlns", ns);
  download(new XMLSerializer().serializeToString(clone), "image/svg+xml", name);
}
function guarded(action) {
  clearError();
  try { action(); } catch (error) { showError(error); }
}

$("selection").addEventListener("submit", draw);
$("gene").addEventListener("input", () => { if (catalog) updateDrivers(); });
$("filters").addEventListener("input", () => {
  $("span-value").textContent = Number($("span").value).toFixed(2);
  guarded(renderPattern);
});
for (const name of ["pattern", "network"]) {
  $(name + "-tab").addEventListener("click", () => {
    for (const other of ["pattern", "network"]) {
      $(other + "-tab").setAttribute("aria-selected", String(name === other));
      $(other + "-tab").tabIndex = name === other ? 0 : -1;
      $(other + "-panel").hidden = name !== other;
    }
  });
  $(name + "-tab").addEventListener("keydown", event => {
    if (["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) {
      event.preventDefault();
      const target = event.key === "Home" ? "pattern" : event.key === "End" ? "network" : name === "pattern" ? "network" : "pattern";
      $(target + "-tab").click(); $(target + "-tab").focus();
    }
  });
}
function profileId() { return applied.profile === "Expression" ? applied.gene : applied.driver; }
$("pattern-svg").onclick = () => guarded(() => saveSVG("pattern", `BrainSpan_${profileId()}_${applied.profile}.svg`));
$("network-svg").onclick = () => guarded(() => saveSVG("network", `Network_${applied.driver}.svg`));
$("pattern-csv").onclick = () => {
  const rows = plotted.map(point => ({ gene: applied.gene, profile: applied.profile, profile_id: profileId(), ...point }));
  download(toCSV(rows, ["gene", "profile", "profile_id", "sample_id", "donor_id", "age", "gender", "Days", "age_group", ...REGION_FIELDS, "value"]), "text/csv", `BrainSpan_${profileId()}_${applied.profile}_samples.csv`);
};
$("network-csv").onclick = () => download(toCSV(edges, ["source", "target", "source_gene", "target_gene", "MI", "pearson", "spearman"]), "text/csv", `Network_${applied.driver}.csv`);
for (const link of document.querySelectorAll("nav a")) {
  const url = new URL(link.href);
  url.searchParams.set("scoutTheme", document.documentElement.dataset.theme);
  link.href = url.href;
}

async function init() {
  try {
    const [manifest, sampleRows, driverRows] = await Promise.all([csv("manifest.csv"), csv("samples.csv"), csv("drivers.csv")]);
    catalog = new Map();
    for (const row of manifest) {
      if (!row.gene || catalog.has(row.gene) || (row.file && !/^batch-\d{5}$/.test(row.file))) throw new Error("Invalid BrainSpan gene catalog.");
      catalog.set(row.gene, row);
    }
    samples = prepareSamples(sampleRows);
    drivers = driverRows;
    if (!catalog.size || !samples.length || drivers.some(row => !catalog.has(row.gene))) throw new Error("Invalid BrainSpan catalog or driver mappings.");
    $("genes").replaceChildren(...manifest.map(row => option(row.gene)));
    for (const field of REGION_FIELDS) {
      $(field).replaceChildren(option("All"), ...[...new Set(samples.map(row => row[field]))].sort().map(value => option(value)));
      $(field).value = "All";
    }
    $("start").value = Math.min(...samples.map(row => row.Days));
    $("end").value = Math.max(...samples.map(row => row.Days));
    const query = new URLSearchParams(location.search).get("gene");
    if (query) $("gene").value = query;
    updateDrivers(); setBusy(false);
    if (currentGene()) await draw();
    else { $("status").textContent = `${catalog.size} NDD genes loaded. Choose a gene and press Draw.`; showError(new Error("The requested gene is not in the NDD catalog.")); }
  } catch (error) {
    $("status").textContent = "BrainSpan data could not be loaded. Check the data files and serve this page over HTTP.";
    showError(error);
  }
}
init();
