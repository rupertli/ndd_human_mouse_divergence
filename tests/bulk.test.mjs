import test from "node:test";
import assert from "node:assert/strict";
import { readFile, readdir } from "node:fs/promises";
import { parseCSV } from "../docs/scripts/explorer-model.mjs";
import { prepareClusters, clusterPeers, geneTerms, clusterTerms, clusterTrend, BULK_PHASES, bulkGeneUrl } from "../docs/scripts/bulk-model.mjs";
import { alignOnset, summarizeOnset } from "../docs/scripts/orphanet-model.mjs";

const root = new URL("../docs/data/ASD_Bulk_analysis/", import.meta.url);
const csv = async file => parseCSV(await readFile(new URL(file, root), "utf8"));

const boundaries = [49.94887588317583, 87.26940995423321, 152.4749011684023,
  266.4002827394653, 465.4478219027315, 979.4586293554893];

test("Orphanet intervals retain all overlapping bins and flag outside-window ages", () => {
  assert.deepEqual(alignOnset("Antenatal", boundaries), {
    bins: ["event_bin_1", "event_bin_2", "event_bin_3"], before_window: true, after_window: false,
  });
  assert.deepEqual(alignOnset("Neonatal", boundaries), {
    bins: ["event_bin_3", "event_bin_4"], before_window: false, after_window: false,
  });
  assert.deepEqual(alignOnset("Infancy", boundaries), {
    bins: ["event_bin_4", "event_bin_5"], before_window: false, after_window: true,
  });
  assert.deepEqual(alignOnset("Childhood", boundaries), {
    bins: [], before_window: false, after_window: true,
  });
  assert.deepEqual(alignOnset("All ages", boundaries).bins, ["event_bin_3", "event_bin_4", "event_bin_5"]);
  assert.throws(() => alignOnset("Unknown", boundaries), /Unsupported/);
  assert.throws(() => alignOnset("Infancy", [1, 2]), /six increasing/);
  const exact = [0, 266, 294, 400, 500, 1000];
  assert.deepEqual(alignOnset("Neonatal", exact).bins, ["event_bin_2"]);
});

test("disease onset distributions retain ties, multi-category entities, and unknowns without inventing a phenotype peak", () => {
  const diseases = [
    { orpha_code: "1", onset_categories: ["Infancy", "Childhood", "Infancy"] },
    { orpha_code: "2", onset_categories: ["Infancy"] },
    { orpha_code: "3", onset_categories: ["Childhood"] },
    { orpha_code: "4", onset_categories: ["No data available"] },
    { orpha_code: "5", onset_categories: [] },
  ];
  const result = summarizeOnset(diseases, boundaries);
  assert.equal(result.associated_entities, 5);
  assert.equal(result.entities_with_onset, 3);
  assert.equal(result.entities_without_onset, 2);
  assert.deepEqual(result.most_reported_categories, ["Infancy", "Childhood"]);
  assert.deepEqual(result.distribution.map(item => item.count), [2, 2]);
  assert.deepEqual(summarizeOnset([], boundaries).distribution, []);
  assert.throws(() => summarizeOnset([diseases[0], diseases[0]], boundaries), /duplicate/);
});

test("bundled HPO onset context reproduces its evidence and exact cortical boundaries", async () => {
  const context = JSON.parse(await readFile(new URL("orphanet-hpo-context.json", root), "utf8"));
  const evidence = JSON.parse(await readFile(new URL("orphanet-hpo-evidence.json", root), "utf8"));
  assert.deepEqual(context.boundaries_pc_days, boundaries);
  assert.deepEqual(Object.keys(context.terms).sort(), Object.keys(evidence.terms).sort());
  assert.ok(Object.keys(context.terms).length > 0);
  assert.equal(context.sources.phenotypes.date.slice(0, 10), context.sources.ages.date.slice(0, 10));
  for (const id of Object.keys(evidence.terms)) {
    const raw = evidence.terms[id];
    const computed = summarizeOnset(raw.diseases, boundaries);
    for (const [key, value] of Object.entries(computed)) assert.deepEqual(context.terms[id][key], value);
    for (const disease of raw.diseases) {
      assert.notEqual(disease.phenotype_frequency, "Excluded (0%)");
      assert.equal(disease.source_url, `https://www.orpha.net/en/disease/detail/${disease.orpha_code}`);
    }
  }
});

test("selected-gene bulk links require exact membership and preserve the Pages subpath and theme", () => {
  const clusters = prepareClusters([{ Gene: "KDM5C", Cluster: "12" }, { Gene: "GENE+1", Cluster: "1" }]);
  const base = "https://example.org/report/docs/index.html?scoutTheme=dark";
  assert.equal(bulkGeneUrl(clusters, "KDM5C", base),
    "https://example.org/report/docs/bulk-analysis.html?gene=KDM5C&scoutTheme=dark");
  assert.equal(bulkGeneUrl(clusters, "KDM5", base), null);
  assert.equal(bulkGeneUrl(clusters, null, base), null);
  const encoded = new URL(bulkGeneUrl(clusters, "GENE+1", base));
  assert.equal(encoded.searchParams.get("gene"), "GENE+1");
});

test("cluster peers exclude the selected gene and other clusters, including singletons", () => {
  const clusters = prepareClusters([
    { Gene: "ABC", Cluster: "1" }, { Gene: "ABC2", Cluster: "1" }, { Gene: "OTHER", Cluster: "2" },
  ]);
  assert.deepEqual(clusterPeers(clusters, "ABC"), ["ABC2"]);
  assert.deepEqual(clusterPeers(clusters, "OTHER"), []);
  assert.throws(() => clusterPeers(clusters, "UNKNOWN"), /No cluster/);
  assert.throws(() => prepareClusters([{ Gene: "A", Cluster: "x" }]), /invalid/);
  assert.throws(() => prepareClusters([{ Gene: "A", Cluster: "1" }, { Gene: "A", Cluster: "2" }]), /duplicate/);
});

test("HPO membership matches whole symbols; GO/KEGG match whole Entrez IDs and sort adjusted p-values", () => {
  const mapping = { ABC: ["60", "600"] };
  const hpo = [
    { Term: "Unrelated", Genes: "ABC2;OTHER", "Adjusted.P.value": ".01" },
    { Term: "Match", Genes: "OTHER; ABC", "Adjusted.P.value": ".04" },
  ];
  assert.deepEqual(geneTerms(hpo, "HPO", "ABC", mapping).map(row => row.Term), ["Match"]);
  for (const kind of ["GO", "KEGG"]) {
    const rows = [
      { ID: "a", Description: "No prefix match", geneID: "6000/160", "p.adjust": ".001" },
      { ID: "b", Description: "Match", geneID: "60/12", "p.adjust": ".04" },
      { ID: "c", Description: "Second mapped ID", geneID: "600", "p.adjust": ".02" },
    ];
    assert.deepEqual(geneTerms(rows, kind, "ABC", mapping).map(row => row.ID), ["c", "b"]);
    assert.throws(() => geneTerms(rows, kind, "ABC", {}), /mapping/);
    assert.throws(() => geneTerms([{ ...rows[0], "p.adjust": "" }], kind, "ABC", mapping), /Invalid/);
    assert.deepEqual(geneTerms([], kind, "ABC", mapping), []);
  }
});

test("cluster enrichment retains terms that do not contain the selected gene and does not mutate exports", () => {
  const rows = [
    { Term: "Contains ABC", Genes: "ABC", "Adjusted.P.value": ".04" },
    { Term: "Does not contain ABC", Genes: "OTHER", "Adjusted.P.value": ".01" },
  ];
  assert.deepEqual(clusterTerms(rows, "HPO").map(row => row.Term), ["Does not contain ABC", "Contains ABC"]);
  assert.equal(rows[0].Term, "Contains ABC");
  assert.deepEqual(clusterTerms([], "GO"), []);
  assert.throws(() => clusterTerms([{ ...rows[0], "Adjusted.P.value": "" }], "HPO"), /Invalid/);
});

test("cluster trends preserve exported means and SE, order phases, and reject incomplete or invalid data", () => {
  const rows = BULK_PHASES.map((TimePoint, i) => ({
    Cluster: "Cluster2", TimePoint, Mean_Expression: String(i + 1), SE: ".5", N: "16",
  }));
  const result = clusterTrend([...rows].reverse(), 2);
  assert.deepEqual(result.map(row => row.phase), BULK_PHASES);
  assert.deepEqual(result[0], { phase: "Embryo", mean: 1, se: .5, n: 16 });
  assert.throws(() => clusterTrend(rows, 1), /No expression trend/);
  assert.throws(() => clusterTrend(rows.slice(1), 2), /Incomplete/);
  assert.throws(() => clusterTrend([...rows, rows[0]], 2), /duplicate/);
  for (const patch of [{ SE: "-1" }, { Mean_Expression: "" }, { N: "0" }, { TimePoint: "unknown" }]) {
    assert.throws(() => clusterTrend([{ ...rows[0], ...patch }, ...rows.slice(1)], 2), /Invalid/);
  }
  const flat = rows.map(row => ({ ...row, Mean_Expression: "0", SE: "0" }));
  assert.ok(clusterTrend(flat, 2).every(row => row.mean === 0 && row.se === 0));
});

test("NDD index preserves all supplied cluster members, trends, and cluster-wide enrichment", async () => {
  const index = JSON.parse(await readFile(new URL("bulk-data.json", root), "utf8"));
  const clusters = prepareClusters(index.assignments);
  const mapping = JSON.parse(await readFile(new URL("gene-entrez-map.json", root), "utf8"));
  const trends = await csv(index.trends);
  const files = await readdir(new URL("cluster_result/", root));
  const geneFiles = files.filter(file => /^expression_kmeans_ndd_genesv2_Cluster\d+\.txt$/.test(file));
  assert.equal(clusters.byCluster.size, geneFiles.length);
  assert.ok(clusters.byGene.size > 2800);
  assert.deepEqual(mapping.genes.ACTB, ["60"]);
  assert.ok(clusters.byGene.has("CCDC115"), "Unmapped symbols remain selectable");
  for (const file of geneFiles) {
    const cluster = Number(file.match(/Cluster(\d+)/)[1]);
    const members = (await readFile(new URL(`cluster_result/${file}`, root), "utf8")).trim().split(/\r?\n/).sort();
    assert.deepEqual([...clusters.byCluster.get(cluster)].sort(), members);
  }
  for (const [cluster, genes] of clusters.byCluster) {
    const trend = clusterTrend(trends, cluster);
    assert.equal(trend.length, 6);
    for (const point of trend) {
      const source = trends.find(row => row.Cluster === `Cluster${cluster}` && row.TimePoint === point.phase);
      assert.equal(point.mean, Number(source.Mean_Expression));
      assert.equal(point.se, Number(source.SE));
    }
    assert.equal(clusterPeers(clusters, genes[0]).length + 1, genes.length);
    for (const kind of ["HPO", "KEGG", "GO"]) {
      const file = index.enrichment[cluster][kind];
      if (file === null) {
        assert.equal(kind, "HPO");
        assert.ok(![4, 11].includes(cluster));
        continue;
      }
      const rows = await csv(file);
      assert.equal(clusterTerms(rows, kind).length, rows.length);
      if (kind === "HPO") {
        assert.ok(file.endsWith("_communities.csv"), `Community export not selected for cluster ${cluster}`);
        const original = await csv(file.replace("_communities.csv", ".csv"));
        assert.equal(rows.length, original.length);
        const byTerm = new Map(original.map(row => [row.Term, row]));
        for (const row of rows) {
          assert.ok(row.community?.trim(), `Missing community for ${row.Term}`);
          const source = byTerm.get(row.Term);
          assert.ok(source, `Unexpected HPO term ${row.Term}`);
          assert.equal(row.Genes, source.Genes);
          assert.equal(Number(row["Adjusted.P.value"]), Number(source["Adjusted.P.value"]));
        }
      } else {
        for (const row of rows) {
          for (const id of row.geneID.split("/")) {
            assert.ok(Array.isArray(mapping.symbols_by_id[id]), `Missing Entrez entry: ${id}`);
          }
        }
      }
    }
  }
});
