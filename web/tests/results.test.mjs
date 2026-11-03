import assert from 'node:assert/strict';
import test from 'node:test';
import { probit, estimateSize, nextSizes, rowChecks, effectText, summaryText } from '../dag-results.js';

// a minimal engine answer: rows of (n, power) with optional null rejection, coverage and bias
const scenario = (rows, extra = {}) => ({
  power: { targetPower: .8, alpha: .05, side: 'two-sided', repetitions: 500, seed: 7 },
  estimand: { kind: 'marginal', measure: 'mean_difference', nullValue: 0 },
  alternative: { truth: { estimate: .5, mcSE: 0, exact: true }, target: null },
  analysis: { label: 'Welch t test of the marginal mean difference' },
  graph: { nodes: [{}, {}] },
  firstTestedN: rows.find((r) => r.power >= .8)?.n ?? null,
  results: rows.map((r) => ({ n: r.n,
    alternative: { valid: 500, repetitions: 500, failures: r.failures ?? 0, probability: r.power, interval: [r.power - .03, r.power + .03],
      coverage: r.coverage ?? .95, coverageValid: 500, bias: r.bias ?? 0, estimateSD: .2 },
    ...(r.null == null ? {} : { null: { valid: 500, repetitions: 500, failures: 0, probability: r.null, interval: [r.null - .02, r.null + .02] } }) })),
  ...extra,
});

test('probit inverts the normal distribution', () => {
  for (const [p, z] of [[.5, 0], [.8, 0.8416212335729143], [.975, 1.959963984540054], [.01, -2.326347874040841], [.999, 3.090232306167813]])
    assert.ok(Math.abs(probit(p) - z) < 1e-8, `${p}`);
});

test('the estimated size lies between the tested sizes that bracket the target', () => {
  const answer = scenario([{ n: 100, power: .698 }, { n: 200, power: .928 }, { n: 400, power: .996 }]);
  const e = estimateSize(answer);
  assert.equal(e.status, 'between');
  assert.equal(e.lower.n, 100); assert.equal(e.upper.n, 200);
  assert.ok(e.n > 100 && e.n < 200, String(e.n));
  // two-sample normal theory for these powers puts 80% near 130 in total
  assert.ok(Math.abs(e.n - 130) <= 8, String(e.n));
  const next = nextSizes(e);
  assert.equal(next.length, 3); assert.ok(next[0] < e.n && next[2] > e.n && next.includes(e.n));
});

test('targets outside the tested range suggest the next sizes to simulate', () => {
  const low = estimateSize(scenario([{ n: 100, power: .3 }, { n: 200, power: .5 }]));
  assert.equal(low.status, 'above');
  assert.deepEqual(nextSizes(low), [400, 800, 1600]);
  const high = estimateSize(scenario([{ n: 100, power: .9 }, { n: 200, power: .99 }]));
  assert.equal(high.status, 'below');
  assert.deepEqual(nextSizes(high), [13, 25, 50]);
  assert.deepEqual(nextSizes(estimateSize(scenario([{ n: 50000, power: .4 }]))), [], 'nothing beyond the largest allowed size');
  assert.equal(estimateSize({ ...scenario([]), results: [] }).status, 'none');
});

test('checks flag a test that over-rejects, poor coverage, bias and failures, and stay quiet otherwise', () => {
  const clean = scenario([{ n: 200, power: .9, null: .052 }]);
  assert.deepEqual(rowChecks(clean, clean.results[0]), []);
  const bad = scenario([{ n: 200, power: .9, null: .246, coverage: .46, bias: .3, failures: 80 }]);
  const kinds = rowChecks(bad, bad.results[0]).map((c) => c.kind);
  assert.deepEqual(kinds, ['size', 'coverage', 'bias', 'failures']);
  assert.match(rowChecks(bad, bad.results[0])[0].short, /Type I error 24\.6%/);
  // ratio measures report bias on the ratio scale, so bias is not judged there
  const ratio = scenario([{ n: 200, power: .9, bias: .3 }], { estimand: { kind: 'marginal', measure: 'risk_ratio', nullValue: 1 } });
  assert.deepEqual(rowChecks(ratio, ratio.results[0]).map((c) => c.kind), []);
});

test('the protocol text states the question, effect, analysis and every tested size', () => {
  const answer = scenario([{ n: 100, power: .698, null: .04 }, { n: 200, power: .928, null: .05 }]);
  const text = summaryText(answer, 'Compare Outcome if Treatment were set to 1 versus 0, using the mean difference.');
  for (const part of ['Compare Outcome', 'Effect assumed: 0.5', 'Welch t test', 'α = 0.05', '500 simulated studies (seed 7)',
    'With 100 participants in total, estimated power is 69.8%', 'rejection rate with no effect is 5.0%', 'smallest tested size reaching 80.0% power is 200'])
    assert.ok(text.includes(part), part);
  assert.equal(effectText({ ...answer, alternative: { truth: { estimate: .4982, mcSE: .004, exact: false } } }), '0.4982 (simulated true value, Monte Carlo SE 0.004)');
  assert.equal(effectText({ ...answer, alternative: { truth: { estimate: .5, mcSE: 3e-19, exact: false } } }), '0.5', 'no Monte Carlo noise is shown for an exact effect');
});
