import { readFile, writeFile } from "node:fs/promises";
import { parseCSV } from "./explorer-model.mjs";
import { BIRTH_PC_DAYS, summarizeOnset } from "./orphanet-model.mjs";

const root = new URL("../data/ASD_Bulk_analysis/", import.meta.url);
const [raw, text] = await Promise.all([
  readFile(new URL("orphanet-hpo-evidence.json", root), "utf8"),
  readFile(new URL("../NDD_human_mouse_cortical_enriched/alignment_boundaries.csv", root), "utf8"),
]);
const evidence = JSON.parse(raw);
const boundaries = parseCSV(text).map(row => Number(row.human_pc_days));
const output = {
  sources: evidence.sources,
  attribution: evidence.attribution,
  license_url: evidence.license_url,
  birth_pc_days: BIRTH_PC_DAYS,
  boundaries_pc_days: boundaries,
  method: "Exact HPO-to-ORPHA associations excluding frequency 0%; one count per unique clinical entity per listed disease-onset category. Categories are not mutually exclusive. No phenotype-onset or patient-incidence inference.",
  terms: Object.fromEntries(Object.entries(evidence.terms).map(([id, term]) =>
    [id, { label: term.label, ...summarizeOnset(term.diseases, boundaries),
      examples: term.diseases.slice(0, 3).map(disease => ({
        orpha_code: disease.orpha_code, name: disease.name, source_url: disease.source_url,
      })) }])),
};
await writeFile(new URL("orphanet-hpo-context.json", root), `${JSON.stringify(output, null, 2)}\n`);
console.log(`Built onset context and interval overlaps for ${Object.keys(output.terms).length} HPO terms.`);
