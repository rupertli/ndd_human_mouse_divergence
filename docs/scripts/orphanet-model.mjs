export const ONSET_CATEGORIES = ["Antenatal", "Neonatal", "Infancy", "Childhood", "Adolescent", "Adult", "Elderly", "All ages"];
export const BIRTH_PC_DAYS = 266;
const YEAR_DAYS = 365.25;
export const ONSET_INTERVALS = {
  Antenatal: [0, BIRTH_PC_DAYS],
  Neonatal: [BIRTH_PC_DAYS, BIRTH_PC_DAYS + 28],
  Infancy: [BIRTH_PC_DAYS + 28, BIRTH_PC_DAYS + 2 * YEAR_DAYS],
  Childhood: [BIRTH_PC_DAYS + 2 * YEAR_DAYS, BIRTH_PC_DAYS + 12 * YEAR_DAYS],
  Adolescent: [BIRTH_PC_DAYS + 12 * YEAR_DAYS, BIRTH_PC_DAYS + 19 * YEAR_DAYS],
  Adult: [BIRTH_PC_DAYS + 19 * YEAR_DAYS, BIRTH_PC_DAYS + 66 * YEAR_DAYS],
  Elderly: [BIRTH_PC_DAYS + 66 * YEAR_DAYS, Infinity],
  "All ages": [BIRTH_PC_DAYS, Infinity],
};

export function alignOnset(category, boundaries) {
  const interval = ONSET_INTERVALS[category];
  if (!interval) throw new Error(`Unsupported Orphanet onset category: ${category}.`);
  if (boundaries.length !== 6 || !boundaries.every((value, i) =>
    Number.isFinite(value) && (i === 0 || value > boundaries[i - 1]))) {
    throw new Error("Expected six increasing cerebral event-bin boundaries.");
  }
  const [start, end] = interval;
  return {
    bins: boundaries.slice(0, -1).flatMap((lower, i) =>
      start < boundaries[i + 1] && end > lower ? [`event_bin_${i + 1}`] : []),
    before_window: start < boundaries[0],
    after_window: end > boundaries.at(-1),
  };
}

export function summarizeOnset(diseases, boundaries) {
  const counts = Object.fromEntries(ONSET_CATEGORIES.map(category => [category, 0]));
  const seen = new Set();
  let unknown = 0;
  for (const disease of diseases) {
    if (!/^\d+$/.test(disease.orpha_code) || seen.has(disease.orpha_code) || !Array.isArray(disease.onset_categories)) {
      throw new Error("Missing or duplicate Orphanet clinical entity.");
    }
    seen.add(disease.orpha_code);
    const categories = [...new Set(disease.onset_categories)];
    if (categories.some(category => !ONSET_CATEGORIES.includes(category) && category !== "No data available")) {
      throw new Error("Unexpected Orphanet onset category.");
    }
    const known = categories.filter(category => category !== "No data available");
    if (!known.length) unknown++;
    for (const category of known) counts[category]++;
  }
  const distribution = ONSET_CATEGORIES.filter(category => counts[category]).map(category => ({
    category, count: counts[category], ...alignOnset(category, boundaries),
  }));
  const max = Math.max(0, ...distribution.map(item => item.count));
  return {
    associated_entities: diseases.length,
    entities_with_onset: diseases.length - unknown,
    entities_without_onset: unknown,
    most_reported_categories: distribution.filter(item => item.count === max).map(item => item.category),
    distribution,
  };
}
