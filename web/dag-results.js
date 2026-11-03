// Reading a power answer for people: where the target power is crossed, which sizes to simulate next, whether
// each tested size's test can be trusted, and a paragraph for a protocol. Pure functions of an engine answer.
// The estimated size is an interpolation between simulated sizes and is always labelled as one; only a run
// at that size measures its power.
const MIN_SIZE = 4, MAX_SIZE = 50000;

// Inverse standard normal CDF (Acklam's rational approximation; relative error below 1.2e-9).
export function probit(p) {
  const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239];
  const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572];
  const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783];
  const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
  const tail = (q) => (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  if (p < 0.02425) return tail(Math.sqrt(-2 * Math.log(p)));
  if (p > 1 - 0.02425) return -tail(Math.sqrt(-2 * Math.log(1 - p)));
  const q = p - 0.5, r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
}

const tested = (answer) => answer.results.filter((row) => row.alternative.valid > 0).sort((a, b) => a.n - b.n);
const clamp = (n) => Math.min(MAX_SIZE, Math.max(MIN_SIZE, Math.round(n)));

// Where power reaches the target: 'between' two tested sizes (with an interpolated size, linear in probit power
// against log size), 'below' every tested size (all already reach it), 'above' them (none does) or 'none'.
export function estimateSize(answer) {
  const rows = tested(answer), target = answer.power.targetPower;
  if (!rows.length) return { status: 'none' };
  const i = rows.findIndex((row) => row.alternative.probability >= target);
  if (i === 0) return { status: 'below', upper: rows[0] };
  if (i < 0) return { status: 'above', lower: rows.at(-1) };
  const lower = rows[i - 1], upper = rows[i];
  const z = (p) => probit(Math.min(0.995, Math.max(0.005, p)));
  const zl = z(lower.alternative.probability), zu = z(upper.alternative.probability), zt = z(target);
  const share = zu > zl ? (zt - zl) / (zu - zl) : 1;
  const n = Math.ceil(Math.exp(Math.log(lower.n) + Math.min(1, Math.max(0, share)) * (Math.log(upper.n) - Math.log(lower.n))));
  return { status: 'between', n: Math.min(upper.n, Math.max(lower.n + 1, n)), lower, upper };
}

// The sizes to simulate next: around the estimate, or further out when the target lies beyond the tested range.
export function nextSizes(estimate) {
  const unique = (list) => [...new Set(list.map(clamp))].sort((a, b) => a - b);
  if (estimate.status === 'between') return unique([estimate.n * 0.85, estimate.n, estimate.n * 1.15]);
  if (estimate.status === 'above') return estimate.lower.n >= MAX_SIZE ? [] : unique([2, 4, 8].map((k) => estimate.lower.n * k));
  if (estimate.status === 'below') return estimate.upper.n <= MIN_SIZE ? [] : unique([1 / 8, 1 / 4, 1 / 2].map((k) => estimate.upper.n * k));
  return [];
}

// Reasons a tested size's power may not mean what it says. Each check fires only when the evidence is clear
// against the simulation's own Monte Carlo error, so a correct test is not flagged by chance.
export function rowChecks(answer, row) {
  const alpha = answer.power.alpha, checks = [], alt = row.alternative;
  const positive = ['risk_ratio', 'odds_ratio', 'exponentiated_coefficient'].includes(answer.estimand.measure);
  if (row.null?.valid > 0 && row.null.interval[0] > alpha)
    checks.push({ kind: 'size', short: `Type I error ${pct(row.null.probability)}`,
      text: `the test rejects ${pct(row.null.probability)} of studies with no effect (α = ${alpha})` });
  if (answer.estimand.kind !== 'joint' && alt.coverage != null && alt.coverageValid >= 50) {
    const nominal = 1 - alpha;
    if (alt.coverage < nominal - 3 * Math.sqrt(nominal * alpha / alt.coverageValid))
      checks.push({ kind: 'coverage', short: `Coverage ${pct(alt.coverage)}`,
        text: `${pct(nominal)} intervals contain the true effect in only ${pct(alt.coverage)} of studies` });
  }
  const truth = answer.alternative.truth?.estimate;
  if (!positive && answer.estimand.kind !== 'joint' && Number.isFinite(alt.bias) && Number.isFinite(truth) && alt.valid >= 30) {
    const effect = Math.abs(truth - answer.estimand.nullValue), noise = 3 * (alt.estimateSD ?? 0) / Math.sqrt(alt.valid);
    if (effect > 0 && Math.abs(alt.bias) > 0.1 * effect && Math.abs(alt.bias) > noise)
      checks.push({ kind: 'bias', short: `Bias ${signed(alt.bias)}`,
        text: `estimates are off by ${signed(alt.bias)} on average (${Math.round(100 * Math.abs(alt.bias) / effect)}% of the effect)` });
  }
  const failed = alt.failures + (row.null?.failures ?? 0), total = alt.repetitions + (row.null?.repetitions ?? 0);
  if (failed > 0.05 * total) checks.push({ kind: 'failures', short: `${failed} failed`, text: `${failed} of ${total} simulated studies could not be analyzed` });
  return checks;
}

const pct = (p) => `${(100 * p).toFixed(1)}%`;
const signed = (x) => `${x > 0 ? '+' : '−'}${Number(Math.abs(x).toPrecision(2))}`;
const number = (x) => Number(Number(x).toPrecision(4)).toString();

// The effect the run used, in words a reader can check against the graph.
export function effectText(answer) {
  const truth = answer.alternative.truth;
  if (answer.estimand.kind === 'joint') return 'the predictor block specified in the model';
  if (!truth || !Number.isFinite(truth.estimate)) return 'the effect implied by the model';
  const exact = truth.exact || !(truth.mcSE > 1e-9 * Math.max(1, Math.abs(truth.estimate)));
  return exact ? number(truth.estimate) : `${number(truth.estimate)} (simulated true value, Monte Carlo SE ${Number(truth.mcSE.toPrecision(2))})`;
}

// A paragraph for a protocol or grant. Every number in it comes from the answer.
export function summaryText(answer, question) {
  const p = answer.power, unit = answer.graph.study ? answer.sizeUnit : 'participants in total';
  const sides = p.side === 'two-sided' ? 'two-sided' : 'one-sided';
  const lines = [`Power was estimated by simulation from a data generating model specified as a directed acyclic graph of ${answer.graph.nodes.length} variables.`,
    `Question: ${question}`,
    `Effect assumed: ${effectText(answer)}${answer.alternative.target != null ? ` (the model was calibrated to a planning effect of ${number(answer.alternative.target)})` : ''}; null value ${answer.estimand.nullValue}.`,
    `Analysis: ${answer.analysis.label}, ${sides} test at α = ${p.alpha}.`,
    `Each sample size was evaluated with ${p.repetitions} simulated studies${p.seed != null ? ` (seed ${p.seed})` : ''}; failed analyses count as non-rejections.`];
  for (const row of tested(answer)) {
    const a = row.alternative;
    lines.push(`With ${row.n} ${unit}, estimated power is ${pct(a.probability)} (95% Monte Carlo interval ${pct(a.interval[0])}–${pct(a.interval[1])})`
      + (row.null?.valid ? `; the rejection rate with no effect is ${pct(row.null.probability)}.` : '.'));
  }
  const first = answer.results.find((row) => row.n === answer.firstTestedN);
  lines.push(first ? `The smallest tested size reaching ${pct(p.targetPower)} power is ${first.n} ${unit}.` : `No tested size reached ${pct(p.targetPower)} power.`);
  lines.push('These results are conditional on the assumed model, effect sizes and analysis; the intervals reflect simulation error only.');
  return lines.join(' ');
}
