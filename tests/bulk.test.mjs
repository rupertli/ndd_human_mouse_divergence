import test from "node:test";
import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { parseCSV } from "../docs/scripts/explorer-model.mjs";
import { prepareClusters, clusterPeers, geneTerms, bulkGeneUrl } from "../docs/scripts/bulk-model.mjs";

const root = new URL("../docs/data/ASD_Bulk_analysis/", import.meta.url);
const csv = async file => parseCSV(await readFile(new URL(file, root), "utf8"));

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

test("bundled data support all 318 cluster genes with correct ACTB and empty HPO results", async () => {
  const clusters = prepareClusters(await csv("cluster_result/expression_kmeans_all_clusters_genesv2.csv"));
  const mapping = JSON.parse(await readFile(new URL("gene-entrez-map.json", root), "utf8")).genes;
  assert.equal(clusters.byGene.size, 318);
  assert.equal(clusters.byCluster.size, 12);
  assert.equal(clusters.byGene.get("ACTB"), 7);
  assert.deepEqual(clusterPeers(clusters, "ACTB"), []);
  assert.deepEqual(mapping.ACTB, ["60"]);
  for (const gene of clusters.byGene.keys()) {
    assert.ok(mapping[gene]?.length, `Missing mapping for ${gene}`);
    assert.ok(mapping[gene].every(id => /^\d+$/.test(id)));
  }
  for (const [cluster, genes] of clusters.byCluster) {
    for (const kind of ["HPO", "KEGG", "GO"]) {
      if (kind === "HPO" && cluster === 10) continue;
      const file = kind === "HPO" ? `HPO_enrichment_resultv2_Cluster${cluster}_communities.csv` :
        kind === "GO" ? `GO_enrichment_resultv2_Cluster${cluster}.csv` : `KEGG_enrichv2_Cluster${cluster}.csv`;
      const rows = await csv(`enrichment/${kind}/expression_kmeans_${file}`);
      if (kind === "HPO" && cluster !== 7) assert.equal(rows.length, 0);
      if (cluster === 7) {
        assert.ok(rows.length > 0);
        assert.equal(geneTerms(rows, kind, "ACTB", mapping).length, rows.length);
      }
      for (const gene of genes) {
        const found = geneTerms(rows, kind, gene, mapping);
        for (const row of found) {
          assert.ok(kind === "HPO" ? row.Genes.split(";").map(g => g.trim()).includes(gene) :
            mapping[gene].some(id => row.geneID.split("/").includes(id)));
        }
      }
    }
  }
});
