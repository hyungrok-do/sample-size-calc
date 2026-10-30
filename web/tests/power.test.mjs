import assert from 'node:assert/strict';
import test from 'node:test';
import { readFile } from 'node:fs/promises';
import { exampleGraph, makeNode, DISTRIBUTIONS } from '../dag-model.js';
import { exportDesign, parseDesign } from '../estimand-model.js';
import { compilePowerPlan, powerCompatibility, parsePowerSettings, calibrationEdges } from '../power-model.js';
import { native, simulateGraph, analyzeReplicate, populationEffect, calibrateGraph, noEffectGraph, runPower } from './dag-engine.mjs';

const near = (actual, expected, tolerance = 1e-9) => assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} differs from ${expected} by more than ${tolerance}`);
const target = (extra = {}) => ({ version: 1, kind: 'marginal', population: 'Graph population', definition: '',
  exposure: 'A', outcome: 'Y', measure: 'mean_difference', reference: 0, active: 1, horizon: null,
  nullValue: 0, alternative: null, ...extra });
const settings = (extra = {}) => parsePowerSettings({ sampleSizes: [100, 300], repetitions: 100, populationSize: 2000, ...extra });


test('paired interventions follow mediated paths and are invariant to skipped exposure draws', () => {
  const graph = exampleGraph();
  const lo = simulateGraph(graph, { n: 200, seed: 14, paired: true, interventions: { A: 0 }, observe: false });
  const hi = simulateGraph(graph, { n: 200, seed: 14, paired: true, interventions: { A: 1 }, observe: false });
  assert.deepEqual(lo.columns.X, hi.columns.X);
  for (let i = 0; i < 200; i++) { near(hi.columns.M[i] - lo.columns.M[i], .8); near(hi.columns.Y[i] - lo.columns.Y[i], .78); }
  assert.throws(() => simulateGraph(graph, { interventions: { A: .5 } }), /binary intervention/);
  assert.throws(() => simulateGraph(graph, { interventions: { Z: 0 } }), /existing variable/);
});

test('the saved estimand chooses the estimator, test scale, and survival horizon', () => {
  const graph = exampleGraph();
  assert.equal(compilePowerPlan(graph, target()).method, 'welch');
  assert.equal(compilePowerPlan(exampleGraph('binary'), target({ measure: 'risk_difference' })).method, 'risk_score');
  assert.equal(compilePowerPlan(exampleGraph('binary'), target({ measure: 'risk_ratio', nullValue: 1 })).testScale, 'log');
  const survival = exampleGraph('survival');
  assert.equal(compilePowerPlan(survival, target({ outcome: 'T', measure: 'rmst_difference', horizon: 14 })).method, 'km_rmst');
  assert.match(compilePowerPlan(survival, target({ outcome: 'T', measure: 'survival_difference', horizon: 14 })).label, /14/);
  assert.equal(compilePowerPlan(survival, target({ outcome: 'T', kind: 'coefficient', measure: 'coefficient' })).method, 'cox');
  assert.equal(powerCompatibility(graph, null).supported, false);
  assert.equal(powerCompatibility(graph, target({ kind: 'custom', measure: 'custom', definition: 'A custom functional' })).supported, false);
});

test('a confounded comparison uses AIPW, adjusting for its parents and not for mediators', () => {
  const graph = exampleGraph();
  graph.edges.push({ source: 'X', target: 'A', coefficient: .6 });
  const plan = compilePowerPlan(graph, target());
  assert.equal(plan.method, 'aipw');
  assert.deepEqual(plan.terms, [{ left: 'X' }]);
  assert.deepEqual(plan.propensityTerms, [{ left: 'X' }]);
  assert.match(plan.label, /AIPW/);
  assert.match(plan.formula, /propensity A ~ X \(logistic\); outcome Y ~ X within each arm/);
  const engine = native({ op: 'plan', graph, estimand: target() });
  for (const key of ['method', 'label', 'formula', 'testScale', 'outcomeModel']) assert.equal(engine[key], plan[key]);
  assert.deepEqual(engine.terms, plan.terms);
  assert.deepEqual(engine.propensityTerms, plan.propensityTerms);
  assert.doesNotThrow(() => compilePowerPlan(graph, target({ kind: 'coefficient', measure: 'coefficient' })));
  assert.throws(() => compilePowerPlan(exampleGraph(), target({ exposure: 'X' })), /binary comparison/);
  const survival = exampleGraph('survival');
  survival.edges.push({ source: 'X', target: 'A', coefficient: .6 });
  assert.throws(() => compilePowerPlan(survival, target({ outcome: 'T', measure: 'rmst_difference', horizon: 14 })), /censoring-adjusted/);
});

test('AIPW power recovers the total effect of a confounded exposure in WebAssembly', () => {
  const graph = exampleGraph();
  graph.edges.push({ source: 'X', target: 'A', coefficient: .8 });
  const truth = populationEffect(graph, target(), { n: 2000 });
  near(truth.estimate, .78, 1e-9);
  const answer = runPower(graph, target(), settings({ sampleSizes: [200], repetitions: 200 }));
  assert.equal(answer.analysis.method, 'aipw');
  const alt = answer.results[0].alternative;
  assert.equal(alt.failures, 0);
  near(alt.meanEstimate, .78, 4 * alt.estimateSD / Math.sqrt(200));
  assert.ok(answer.results[0].null.probability < .12);
});

test('KM power rejects dependent censoring and horizons beyond administrative support', () => {
  const graph = exampleGraph('survival');
  const estimand = target({ outcome: 'T', measure: 'rmst_difference', horizon: 14 });
  graph.edges.push({ source: 'X', target: 'C', coefficient: .3 });
  assert.throws(() => compilePowerPlan(graph, estimand), /Unadjusted KM/);
  assert.doesNotThrow(() => compilePowerPlan(graph, { ...estimand, kind: 'coefficient', measure: 'coefficient', horizon: null }));
  graph.edges.push({ source: 'T', target: 'C', coefficient: .01 });
  assert.throws(() => compilePowerPlan(graph, { ...estimand, kind: 'coefficient', measure: 'coefficient', horizon: null }), /descendant/);
  const administrative = exampleGraph('survival');
  administrative.nodes.find((n) => n.id === 'T').censoring = { mode: 'administrative', time: 10, variable: '' };
  assert.throws(() => compilePowerPlan(administrative, estimand), /horizon exceeds/);
});


test('survival analyses use observed time and event status, never latent outcome values', () => {
  const graph = exampleGraph('survival'), estimand = target({ outcome: 'T', measure: 'rmst_difference', horizon: 10 });
  const data = simulateGraph(graph, { n: 1000, seed: 34 });
  const first = analyzeReplicate(graph, estimand, data);
  data.columns.T.fill(1000000);
  assert.deepEqual(analyzeReplicate(graph, estimand, data), first);
  const coefficient = { ...estimand, kind: 'coefficient', measure: 'coefficient', horizon: null };
  const a = analyzeReplicate(graph, coefficient, data);
  data.columns.T.fill(.01);
  assert.deepEqual(analyzeReplicate(graph, coefficient, data), a);
});

test('total-effect calibration includes mediation and preserves the drawn graph', () => {
  const graph = exampleGraph(), original = JSON.stringify(graph), estimand = target({ alternative: .5 });
  near(populationEffect(graph, estimand, { n: 1000, seed: 20 }).estimate, .78);
  const calibrated = calibrateGraph(graph, estimand, .5, settings(), 34);
  near(calibrated.truth.estimate, .5, .005);
  near(calibrated.parameter.used, .5 - .6 * .8, .005);
  assert.equal(JSON.stringify(graph), original);
  assert.equal(calibrationEdges(graph, estimand)[0].source, 'A');
  assert.equal(calibrationEdges(graph, estimand)[0].target, 'Y');
  assert.throws(() => calibrateGraph(graph, estimand, 20, settings({ calibrationLower: -1, calibrationUpper: 1 }), 34), /cannot be reached/);
});

test('the no-effect null removes all treatment paths and interactions, not just the direct arrow', () => {
  const graph = exampleGraph(), estimand = target();
  graph.nodes.at(-1).interactions.push({ left: 'A', right: 'X', coefficient: .2 });
  const nullGraph = noEffectGraph(graph, estimand);
  assert.equal(nullGraph.edges.find((e) => e.source === 'A' && e.target === 'M').coefficient, 0);
  assert.equal(nullGraph.edges.find((e) => e.source === 'A' && e.target === 'Y').coefficient, 0);
  assert.equal(nullGraph.nodes.at(-1).interactions[0].coefficient, 0);
  assert.equal(nullGraph.edges.find((e) => e.source === 'M' && e.target === 'Y').coefficient, .6);
  near(populationEffect(nullGraph, estimand, { n: 1000, seed: 22 }).estimate, 0);
});

test('conditional ratio planning calibrates and tests the correct coefficient scale', () => {
  const graph = exampleGraph('survival'), estimand = target({ kind: 'coefficient', outcome: 'T', measure: 'exponentiated_coefficient', nullValue: 1, alternative: .65 });
  const calibrated = calibrateGraph(graph, estimand, .65, settings(), 3);
  near(calibrated.parameter.used, Math.log(.65)); near(calibrated.truth.estimate, .65);
  assert.equal(calibrated.truth.exact, true);
});

test('power increases with sample size, null rejection is near alpha, and runs are reproducible', () => {
  const graph = exampleGraph(), estimand = target({ alternative: .5 });
  const options = settings({ sampleSizes: [80, 400], repetitions: 500, populationSize: 1000 });
  const first = runPower(graph, estimand, options);
  assert.ok(first.results[1].alternative.probability > first.results[0].alternative.probability + .3);
  for (const row of first.results) {
    assert.ok(row.null.probability > .015 && row.null.probability < .09);
    assert.equal(row.alternative.failures, 0);
    assert.equal(row.alternative.repetitions, options.repetitions);
  }
  assert.equal(first.firstTestedN, 400);
  const restored = parseDesign(JSON.stringify(first));
  assert.deepEqual(restored.graph, graph); assert.deepEqual(restored.estimand, estimand);
  assert.deepEqual(restored.power, options);
  assert.deepEqual(runPower(graph, estimand, options), first);
});

test('analysis failures stay in the denominator rather than inflating power', () => {
  const graph = exampleGraph('binary'); graph.nodes.find((n) => n.id === 'A').intercept = Math.log(.01 / .99);
  const result = runPower(graph, target({ measure: 'risk_difference' }), settings({ sampleSizes: [10], repetitions: 100, checkNull: false }));
  const row = result.results[0].alternative;
  assert.ok(row.failures > 80); assert.equal(row.valid + row.failures, 100);
  assert.equal(row.probability, row.rejected / 100);
  assert.ok(row.failureReasons.empty_arm > 0);
});

test('editing the estimand, horizon or rejection direction changes the actual analysis', () => {
  const graph = exampleGraph('survival'), base = target({ outcome: 'T', measure: 'rmst_difference', horizon: 10 });
  const data = simulateGraph(graph, { n: 1000, seed: 36 });
  const at10 = analyzeReplicate(graph, base, data), at14 = analyzeReplicate(graph, { ...base, horizon: 14 }, data);
  assert.notEqual(at10.estimate, at14.estimate); assert.notEqual(at10.statistic, at14.statistic);
  const survival = analyzeReplicate(graph, { ...base, measure: 'survival_difference' }, data);
  assert.notEqual(survival.estimate, at10.estimate);
  assert.notEqual(analyzeReplicate(graph, base, data, 'greater').pValue, analyzeReplicate(graph, base, data, 'less').pValue);
});

test('saved designs preserve numerical power settings alongside the graph and estimand', () => {
  const graph = exampleGraph(), estimand = target(), configuration = settings({ sampleSizes: '80, 200', repetitions: 200 });
  const design = exportDesign(graph, estimand, configuration), parsed = parseDesign(JSON.stringify(design));
  assert.deepEqual(parsed.power, configuration); assert.deepEqual(parsed.estimand, estimand);
  assert.deepEqual(parsed.graph, graph);
  for (const invalid of [{ repetitions: 0 }, { sampleSizes: '100, NaN' }, { alpha: 0 }, { seed: -.1 }, { calibrationLower: 4, calibrationUpper: 3 }])
    assert.throws(() => parsePowerSettings(invalid));
});


test('WASM regression estimates agree with independently fitted fixed datasets', async () => {
  const reference = JSON.parse(await readFile(new URL('../../core/tests/fixtures/dag_reference.json', import.meta.url), 'utf8'));
  for (const item of reference.cases) {
    const fit = analyzeReplicate(item.graph, item.estimand, { n: item.n, columns: item.columns });
    near(fit.estimate, item.reference.coefficient, 2e-6); near(fit.se, item.reference.se, 2e-6);
    near(fit.pValue, item.reference.p, 1e-8);
  }
});

test('UI analysis hints agree with the executable C++ plans', () => {
  for (const [graph, estimand] of [
    [exampleGraph(), target()],
    [exampleGraph('binary'), target({ measure: 'risk_difference' })],
    [exampleGraph('survival'), target({ outcome: 'T', measure: 'rmst_difference', horizon: 14 })],
    [exampleGraph('survival'), target({ outcome: 'T', kind: 'coefficient', measure: 'exponentiated_coefficient', nullValue: 1 })],
  ]) assert.deepEqual(native({ op: 'plan', graph, estimand }), compilePowerPlan(graph, estimand));
});
