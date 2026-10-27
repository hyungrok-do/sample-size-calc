export const DEFAULT_POWER = { sampleSizes: [100, 200, 400], repetitions: 500, alpha: .05,
  side: 'two-sided', seed: 2026, populationSize: 20000, calibrationEdge: '',
  calibrationLower: -3, calibrationUpper: 3, checkNull: true, targetPower: .8 };
const fail = (message) => { throw new Error(message); };
export function parsePowerSettings(raw = {}) {
  const x = { ...DEFAULT_POWER, ...raw };
  const sizes = Array.isArray(x.sampleSizes) ? x.sampleSizes.map(Number) : String(x.sampleSizes).trim().split(/[\s,]+/).filter(Boolean).map(Number);
  x.sampleSizes = [...new Set(sizes)].sort((a, b) => a - b);
  if (!x.sampleSizes.length || x.sampleSizes.length > 10 || x.sampleSizes.some((n) => !Number.isInteger(n) || n < 4 || n > 50000)) fail('Choose 1–10 total sample sizes, each a whole number from 4 to 50,000.');
  for (const [key, min, max] of [['repetitions', 10, 10000], ['populationSize', 1000, 50000], ['seed', 0, 4294967295]]) {
    x[key] = Number(x[key]); if (!Number.isInteger(x[key]) || x[key] < min || x[key] > max) fail(`${key} must be an integer from ${min} to ${max}.`);
  }
  for (const key of ['alpha', 'targetPower']) { x[key] = Number(x[key]); if (!(x[key] > 0 && x[key] < 1)) fail(`${key} must be strictly between 0 and 1.`); }
  for (const key of ['calibrationLower', 'calibrationUpper']) { x[key] = Number(x[key]); if (!Number.isFinite(x[key])) fail('Calibration bounds must be finite.'); }
  if (x.calibrationLower >= x.calibrationUpper) fail('The lower calibration bound must be smaller than the upper bound.');
  if (!['two-sided', 'greater', 'less'].includes(x.side)) fail('Choose a two-sided, greater-than, or less-than test.');
  if (typeof x.checkNull !== 'boolean') fail('The null-check setting must be true or false.');
  if (typeof x.calibrationEdge !== 'string') fail('Choose a calibration arrow.');
  return Object.fromEntries(Object.keys(DEFAULT_POWER).map((key) => [key, x[key]]));
}
