export function prepareClusters(records) {
  const byGene = new Map();
  const byCluster = new Map();
  for (const row of records) {
    const gene = row.Gene?.trim();
    const cluster = Number(row.Cluster);
    if (!gene || !Number.isInteger(cluster) || cluster < 1 || byGene.has(gene)) {
      throw new Error("Missing, invalid, or duplicate cluster assignment.");
    }
    byGene.set(gene, cluster);
    if (!byCluster.has(cluster)) byCluster.set(cluster, []);
    byCluster.get(cluster).push(gene);
  }
  if (!byGene.size) throw new Error("Cluster assignments are empty.");
  for (const genes of byCluster.values()) genes.sort((a, b) => a.localeCompare(b));
  return { byGene, byCluster };
}

export function geneTerms(records, kind, gene, geneMap) {
  const ids = geneMap[gene];
  if (!Array.isArray(ids) || !ids.length || !ids.every(id => /^\d+$/.test(id))) {
    throw new Error(`Missing or invalid Entrez mapping for ${gene}.`);
  }
  const column = kind === "HPO" ? "Genes" : "geneID";
  const pColumn = kind === "HPO" ? "Adjusted.P.value" : "p.adjust";
  for (const row of records) {
    if (!(column in row) || !(pColumn in row) ||
        !(kind === "HPO" ? row.Term : row.ID && row.Description) ||
        row[pColumn] === "" || !Number.isFinite(Number(row[pColumn])) ||
        Number(row[pColumn]) < 0 || Number(row[pColumn]) > 1) {
      throw new Error(`Invalid ${kind} enrichment row.`);
    }
  }
  return records.filter(row => {
    const members = row[column].split(kind === "HPO" ? ";" : "/").map(value => value.trim());
    return kind === "HPO" ? members.includes(gene) : ids.some(id => members.includes(id));
  }).sort((a, b) => Number(a[pColumn]) - Number(b[pColumn]) ||
    (a.Term ?? a.ID).localeCompare(b.Term ?? b.ID));
}

export function clusterPeers(clusters, gene) {
  const cluster = clusters.byGene.get(gene);
  if (!cluster) throw new Error(`No cluster assignment for ${gene}.`);
  return clusters.byCluster.get(cluster).filter(peer => peer !== gene);
}

export function bulkGeneUrl(clusters, gene, pageUrl) {
  if (!clusters.byGene.has(gene)) return null;
  const current = new URL(pageUrl);
  const url = new URL("bulk-analysis.html", current);
  url.searchParams.set("gene", gene);
  const theme = current.searchParams.get("scoutTheme");
  if (theme) url.searchParams.set("scoutTheme", theme);
  return url.href;
}
