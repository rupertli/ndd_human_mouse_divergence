import { readdir, readFile, access, writeFile } from "node:fs/promises";
import { prepareClusters, clusterTerms, clusterTrend } from "./bulk-model.mjs";
import { parseCSV } from "./explorer-model.mjs";

const root = new URL("../data/ASD_Bulk_analysis/", import.meta.url);
const files = await readdir(new URL("cluster_result/", root));
const assignments = [];
for (const file of files.sort()) {
  const match = file.match(/^expression_kmeans_ndd_genesv2_Cluster(\d+)\.txt$/);
  if (!match) continue;
  const genes = (await readFile(new URL(`cluster_result/${file}`, root), "utf8"))
    .split(/\r?\n/).map(gene => gene.trim()).filter(Boolean);
  if (!genes.length) throw new Error(`Empty cluster gene list: ${file}`);
  assignments.push(...genes.map(Gene => ({ Gene, Cluster: Number(match[1]) })));
}
const clusters = prepareClusters(assignments);
const trends = "cluster_result/cluster_expression_ndd_kmeans_trends_by_timepointv2.csv";
const trendRows = parseCSV(await readFile(new URL(trends, root), "utf8"));
const hpoFiles = new Set(await readdir(new URL("enrichment/HPO/", root)));
const enrichment = {};
for (const cluster of [...clusters.byCluster.keys()].sort((a, b) => a - b)) {
  clusterTrend(trendRows, cluster);
  const hpo = `expression_kmeans_ndd_HPO_enrichment_resultv2_Cluster${cluster}.csv`;
  const communities = hpo.replace(/\.csv$/, "_communities.csv");
  const hpoSource = hpoFiles.has(communities) ? communities : hpoFiles.has(hpo) ? hpo : null;
  enrichment[cluster] = {
    HPO: hpoSource ? `enrichment/HPO/${hpoSource}` : null,
    GO: `enrichment/GO/expression_kmeans_ndd_GO_enrichment_resultv2_Cluster${cluster}.csv`,
    KEGG: `enrichment/KEGG/expression_kmeans_ndd_KEGG_enrichv2_Cluster${cluster}.csv`,
  };
  for (const [kind, file] of Object.entries(enrichment[cluster])) {
    if (!file) continue;
    await access(new URL(file, root));
    clusterTerms(parseCSV(await readFile(new URL(file, root), "utf8")), kind);
  }
}
await writeFile(new URL("bulk-data.json", root), JSON.stringify({
  assignments, trends, enrichment,
}, null, 2) + "\n");
console.log(`Built bulk data index: ${clusters.byGene.size} genes, ${clusters.byCluster.size} clusters.`);
