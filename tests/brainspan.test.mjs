import test from "node:test";
import assert from "node:assert/strict";
import { readFile, readdir } from "node:fs/promises";
import { parseCSV } from "../docs/scripts/explorer-model.mjs";
import { selectGeneData, prepareSamples, profilePoints, filterPoints, smoothProfile, driverNetwork, toCSV } from "../docs/scripts/brainspan-model.mjs";

const root = new URL("../docs/data/BrainSpan_NDD/", import.meta.url);
const csv = async file => parseCSV(await readFile(new URL(file, root), "utf8"));
const sample = (sample_id, Days, region = "A") => ({
  sample_id, Days: String(Days), donor_id: "donor", structure_acronym: region, structure_group1: "cortex", structure_group2: "frontal",
});

test("profiles join by sample ID, retain missingness and reject invalid or duplicated values", () => {
  const samples = prepareSamples([sample("second", 20), sample("first", 10), sample("missing", 30)]);
  const rows = [{ profile: "Expression", id: "GENE", first: "5", second: "-2", missing: "" }];
  assert.deepEqual(profilePoints(rows, samples, "Expression", "GENE").map(p => [p.sample_id, p.value]), [["second", -2], ["first", 5]]);
  assert.deepEqual(profilePoints(rows, samples, "Activity", "GENE_TF"), []);
  assert.throws(() => profilePoints([...rows, ...rows], samples, "Expression", "GENE"), /Duplicate/);
  assert.throws(() => profilePoints([{ ...rows[0], second: "NA" }], samples, "Expression", "GENE"), /invalid/);
  assert.throws(() => profilePoints(rows, prepareSamples([sample("unknown", 10)]), "Expression", "GENE"), /Missing profile sample/);
  assert.throws(() => prepareSamples([sample("same", 10), sample("same", 20)]), /duplicate/);
  assert.throws(() => prepareSamples([sample("bad", 0)]), /Invalid/);
});

test("region filters intersect, All overrides, empty selections and age boundaries behave explicitly", () => {
  const points = prepareSamples([sample("1", 10), sample("2", 20, "B"), sample("3", 30)]);
  const all = { structure_acronym: ["All"], structure_group1: ["All"], structure_group2: ["All"] };
  assert.deepEqual(filterPoints(points, { start: 10, end: 20, regions: all }).map(p => p.sample_id), ["1", "2"]);
  assert.deepEqual(filterPoints(points, { start: 10, end: 30, regions: { ...all, structure_acronym: ["A"] } }).map(p => p.sample_id), ["1", "3"]);
  assert.equal(filterPoints(points, { start: 10, end: 30, regions: { ...all, structure_acronym: ["A", "All"] } }).length, 3);
  assert.equal(filterPoints(points, { start: 10, end: 30, regions: { ...all, structure_group1: [] } }).length, 0);
  assert.equal(filterPoints(points, { start: 40, end: 50, regions: all }).length, 0);
  assert.throws(() => filterPoints(points, { start: 30, end: 10, regions: all }), /valid positive age/);
});

test("local regression recovers a quadratic in log age and handles sparse or repeated-age samples", () => {
  const points = Array.from({ length: 24 }, (_, i) => {
    const x = 1 + i / 24;
    return { Days: 10 ** x, value: 2 + 3 * x - x ** 2 };
  });
  for (const fit of smoothProfile(points, .6)) {
    const x = Math.log10(fit.Days);
    assert.ok(Math.abs(fit.value - (2 + 3 * x - x ** 2)) < 1e-8);
    assert.ok(fit.low <= fit.value && fit.high >= fit.value);
  }
  assert.deepEqual(smoothProfile(points.slice(0, 3), .6), []);
  assert.deepEqual(smoothProfile(Array.from({ length: 10 }, () => ({ Days: 10, value: 2 })), .6), []);
  const repeated = [10, 10, 10, 20, 20, 30, 30, 40].map((Days, i) => ({ Days, value: i }));
  assert.ok(smoothProfile(repeated, .05).every(row => Number.isFinite(row.value)));
  assert.throws(() => smoothProfile(points, 0), /span/);
});

test("networks select the exact driver and reject non-NDD endpoints; download CSV quotes safely", () => {
  const rows = [
    { source: "A_TF", target: "B", source_gene: "A", target_gene: "B", MI: ".4", pearson: "", spearman: "-.3" },
    { source: "A_SIG", target: "B", source_gene: "A", target_gene: "B", MI: ".5", pearson: ".1", spearman: ".2" },
  ];
  const selected = driverNetwork(rows, "A_TF", new Set(["A", "B"]));
  assert.equal(selected.length, 1);
  assert.equal(selected[0].pearson, null);
  assert.equal(selected[0].spearman, -.3);
  assert.throws(() => driverNetwork(rows, "A_TF", new Set(["A"])), /non-NDD/);
  assert.throws(() => driverNetwork([{ ...rows[0], MI: "-1" }], "A_TF", new Set(["A", "B"])), /negative/);
  const fixture = [{ id: 'A,"B"\nC', value: null }];
  assert.deepEqual(parseCSV(toCSV(fixture, ["id", "value"])), [{ id: fixture[0].id, value: "" }]);
});

test("shared batches isolate a selected gene's profiles and source edges, not edges targeting that gene", () => {
  const batch = {
    profiles: [{ gene: "A", profile: "Expression", id: "A" }, { gene: "B", profile: "Expression", id: "B" }],
    network: [{ source_gene: "A", target_gene: "B" }, { source_gene: "B", target_gene: "A" }],
  };
  assert.deepEqual(selectGeneData(batch, "A"), { profiles: [batch.profiles[0]], network: [batch.network[0]] });
  assert.deepEqual(selectGeneData(batch, "MISSING"), { profiles: [], network: [] });
});

test("exported batches retain all sample values, isolate catalog genes and use at most 100 genes per pair of files", async () => {
  const manifest = await csv("manifest.csv"), provenance = await csv("ndd_gene_provenance.csv");
  const samples = prepareSamples(await csv("samples.csv")), mappings = await csv("drivers.csv");
  const genes = new Set(provenance.map(row => row.gene));
  assert.deepEqual(manifest.map(row => row.gene), [...genes].sort());
  const expectedFiles = [];
  const batches = new Map();
  for (const file of new Set(manifest.map(row => row.file).filter(Boolean))) {
    assert.match(file, /^batch-\d{5}$/);
    const members = new Set(manifest.filter(row => row.file === file).map(row => row.gene));
    assert.ok(members.size <= 100);
    expectedFiles.push(`${file}-profiles.csv`, `${file}-network.csv`);
    const data = { profiles: await csv(`genes/${file}-profiles.csv`), network: await csv(`genes/${file}-network.csv`) };
    assert.ok(data.profiles.every(row => members.has(row.gene)));
    assert.ok(data.network.every(row => members.has(row.source_gene)));
    batches.set(file, data);
  }
  assert.equal(expectedFiles.length, Math.ceil(manifest.filter(row => row.file).length / 100) * 2);
  let expressionCount = 0, activityCount = 0, edgeCount = 0;
  for (const row of manifest) {
    const drivers = mappings.filter(driver => driver.gene === row.gene).map(driver => driver.driver_id);
    if (!row.file) {
      assert.equal(row.expression, "FALSE");
      assert.equal(Number(row.activity_drivers), 0);
      assert.equal(Number(row.network_edges), 0);
      continue;
    }
    const { profiles, network } = selectGeneData(batches.get(row.file), row.gene);
    const seen = new Set();
    for (const profile of profiles) {
      assert.ok(["Expression", "Activity"].includes(profile.profile));
      assert.equal(seen.has(`${profile.profile}:${profile.id}`), false);
      seen.add(`${profile.profile}:${profile.id}`);
      assert.deepEqual(Object.keys(profile).filter(key => !["gene", "profile", "id"].includes(key)).sort(), samples.map(p => p.sample_id).sort());
      if (profile.profile === "Expression") { assert.equal(profile.id, row.gene); expressionCount++; }
      else { assert.ok(drivers.includes(profile.id)); activityCount++; }
      assert.ok(profilePoints([profile], samples, profile.profile, profile.id).every(p => Number.isFinite(p.value)));
    }
    assert.equal(profiles.filter(p => p.profile === "Expression").length, row.expression === "TRUE" ? 1 : 0);
    assert.equal(profiles.filter(p => p.profile === "Activity").length, Number(row.activity_drivers));
    assert.equal(network.length, Number(row.network_edges));
    for (const edge of network) {
      assert.ok(drivers.includes(edge.source));
      assert.equal(edge.source_gene, row.gene);
      assert.ok(genes.has(edge.target_gene));
    }
    for (const driver of drivers) driverNetwork(network, driver, genes);
    edgeCount += network.length;
  }
  assert.deepEqual((await readdir(new URL("genes/", root))).sort(), expectedFiles.sort(), "No stale or unlisted shards");
  assert.ok(expressionCount > 0 && activityCount > 0 && edgeCount > 0, "Exports include all three data types");
});
