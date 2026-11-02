import test from 'node:test';
import assert from 'node:assert/strict';
import { exampleGraph, parseGraph, renameNode, removeNode, validateGraph, dataCsv } from '../dag-model.js';
import { exportDesign, parseDesign } from '../estimand-model.js';
import { defaultTarget } from '../dag-power.js';
import { contrastWeights } from '../dag-study.js';
import { compilePowerPlan } from '../power-model.js';
import { native } from './dag-engine.mjs';
const close = (a, b, tol = 1e-8) => assert.ok(Math.abs(a - b) < tol, `${a} versus ${b}`);
const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;
const variance = (a) => { const m = mean(a); return a.reduce((x, y) => x + (y - m) ** 2, 0) / (a.length - 1); };
const simulate = (graph, n = 300, seed = 42) => native({ op: 'simulate', graph, options: { n, seed } });
const analyze = (graph, estimand, data, analysis = {}) => native({ op: 'analyze', graph, estimand, data, analysis });
const population = (graph, estimand, n = 1000) => native({ op: 'population', graph, estimand, options: { n, seed: 17 } });
const target = (g, contrast = 'final') => ({ ...defaultTarget(g), repeated: { contrast } });

test('study structure, observation rules, scopes and contrasts survive saved designs and renaming', () => {
  let g = exampleGraph('repeated');
  g.nodes.at(-1).lags = [{ source: 'Y', lag: 1, coefficient: .2, initial: 0 }];
  g.study.observation = { ...g.study.observation, mode: 'dropout', historyVariable: 'Y', historyCoefficient: .2 };
  const e = target(g, 'change');
  const saved = parseDesign(exportDesign(g, e));
  assert.deepEqual(saved.graph, parseGraph(g)); assert.deepEqual(saved.estimand.repeated, e.repeated);
  g = renameNode(saved.graph, 'Y', 'Response');
  assert.equal(g.study.observation.historyVariable, 'Response'); assert.equal(g.nodes.at(-1).lags[0].source, 'Response');
  const without = removeNode(g, 'Response'); assert.equal(without.study.observation.historyVariable, ''); validateGraph(without);
});

test('structured simulation preserves subject/cluster scope, identifiers and seeded prefixes', () => {
  const g = exampleGraph('clustered-repeated'), d = simulate(g, 6), larger = simulate(g, 9);
  assert.equal(d.units, 6); assert.equal(d.subjects, 72); assert.equal(d.n, 288); assert.equal(d.sizeUnit, 'clusters');
  for (const key of Object.keys(d.columns)) assert.deepEqual(d.columns[key], larger.columns[key].slice(0, d.n));
  for (let i = 1; i < d.n; i++) {
    if (d.columns.subject_id[i] === d.columns.subject_id[i - 1]) assert.equal(d.columns.X[i], d.columns.X[i - 1]);
    if (d.columns.cluster_id[i] === d.columns.cluster_id[i - 1]) assert.equal(d.columns.A[i], d.columns.A[i - 1]);
    assert.equal(d.columns.Time[i], g.study.times[i % 4]);
  }
  const reorder = { ...g, nodes: [...g.nodes].reverse(), edges: [...g.edges].reverse() };
  assert.deepEqual(simulate(reorder, 6).columns, d.columns);
  assert.match(dataCsv(g, d), /subject_id/); assert.match(dataCsv(g, d), /Y_recorded/);
});

test('irregular-time covariance matches random intercept plus continuous-time residual correlation', () => {
  const g = exampleGraph('repeated'); g.study.times = [0, 1, 3];
  for (const edge of g.edges) edge.coefficient = 0;
  const y = g.nodes.at(-1); y.sd = 1.2; y.interactions = []; y.repeated = { subjectSD: .7, residualCorrelation: .4 };
  const d = simulate(g, 8000), v = [0, 1, 2].map((j) => d.columns.Y.filter((_, i) => i % 3 === j));
  close(variance(v[0]), .49 + 1.44, .08);
  for (const j of [1, 2]) {
    const m0 = mean(v[0]), mj = mean(v[j]);
    const cov = mean(v[0].map((x, i) => (x - m0) * (v[j][i] - mj)));
    close(cov, .49 + 1.44 * .4 ** g.study.times[j], .06);
  }
});

test('lagged outcome uses the previous visit and an explicit initial value', () => {
  const g = exampleGraph('repeated'); const base = simulate(g, 20);
  g.nodes.at(-1).lags = [{ source: 'Y', lag: 1, coefficient: .25, initial: 2 }];
  const lagged = simulate(g, 20);
  for (let i = 0; i < lagged.n; i++) close(lagged.columns.Y[i], base.columns.Y[i] + .25 * (i % 4 ? lagged.columns.Y[i - 1] : 2));
});

test('MMRM estimates and variance agree with independent complete-data contrast calculations', () => {
  const g = exampleGraph('repeated'), d = simulate(g, 400);
  for (const contrast of ['final', 'average', 'change', 'slope']) {
    const e = target(g, contrast), w = contrastWeights(g, e), arms = [[], []];
    for (let i = 0; i < d.n; i += 4) arms[d.columns.A[i]].push(w.reduce((x, v, j) => x + v * d.columns.Y[i + j], 0));
    const f = analyze(g, e, d, { method: 'mmrm' });
    const pooled = ((arms[0].length - 1) * variance(arms[0]) + (arms[1].length - 1) * variance(arms[1])) / (d.subjects - 2);
    close(f.estimate, mean(arms[1]) - mean(arms[0]));
    close(f.se, Math.sqrt(pooled * (1 / arms[0].length + 1 / arms[1].length)));
  }
});

test('linear mixed-model trend agrees with independently computed subject slopes', () => {
  const g = exampleGraph('repeated'), e = target(g, 'slope'), d = simulate(g, 140), w = contrastWeights(g, e), arms = [[], []];
  for (let i = 0; i < d.n; i += 4) arms[d.columns.A[i]].push(w.reduce((x, v, j) => x + v * d.columns.Y[i + j], 0));
  const f = analyze(g, e, d, { method: 'lmm', timeBasis: 'linear', randomSlope: true });
  close(f.estimate, mean(arms[1]) - mean(arms[0]), 1e-6);
  assert.ok(f.se > 0 && f.df > 0);
});

test('clustered trajectory variance agrees with an independent cluster influence calculation', () => {
  const g = exampleGraph('clustered-repeated'), e = target(g), d = simulate(g, 45), groups = [[], []];
  for (let r = 3; r < d.n; r += 4) groups[d.columns.A[r]].push(r);
  const means = groups.map((rows) => mean(rows.map((r) => d.columns.Y[r]))), influence = Array(d.units).fill(0);
  for (const arm of [0, 1]) for (const r of groups[arm]) influence[d.columns.cluster_id[r]] += (arm ? 1 : -1) * (d.columns.Y[r] - means[arm]) / groups[arm].length;
  const f = analyze(g, e, d);
  close(f.estimate, means[1] - means[0]);
  close(f.se, Math.sqrt(influence.reduce((s, x) => s + x * x, 0) * d.units / (d.units - 1)));
  assert.equal(f.df, 43);
});

test('dropout affects observation masks and never exposes latent outcomes to MMRM', () => {
  const g = exampleGraph('repeated');
  g.study.observation = { mode: 'dropout', intercept: -2, historyVariable: 'Y', historyCoefficient: .2 };
  const d = simulate(g, 500), e = target(g), f = analyze(g, e, d);
  let missing = 0;
  for (let i = 0; i < d.n; i++) {
    if (!d.columns.Y_observed[i]) { assert.equal(d.columns.Y_recorded[i], null); d.columns.Y[i] = 1e10; missing++; }
    if (i % 4 && !d.columns.Y_observed[i - 1]) assert.equal(d.columns.Y_observed[i], 0);
  }
  assert.ok(missing > 0); assert.deepEqual(analyze(g, e, d), f);
});

test('intermittent missing visits work with LMM and are rejected by monotone MMRM', () => {
  const g = exampleGraph('repeated'); g.study.observation.mode = 'intermittent'; g.study.observation.intercept = -1;
  const d = simulate(g, 160), e = target(g);
  assert.ok(analyze(g, e, d, { method: 'lmm', timeBasis: 'linear' }).se > 0);
  assert.throws(() => analyze(g, e, d, { method: 'mmrm' }), /monotone/);
});

test('population targets, effect calibration and null construction use the selected trajectory contrast', () => {
  const g = exampleGraph('repeated');
  for (const [contrast, truth] of [['final', .72], ['average', .315], ['change', .72], ['slope', .06]]) close(population(g, target(g, contrast)).estimate, truth);
  const e = target(g);
  const calibrated = native({ op: 'calibrate', graph: g, estimand: e, target: .4, settings: { populationSize: 1000 } });
  close(population(calibrated.graph, e).estimate, .4, .005);
  const nullGraph = native({ op: 'nullGraph', graph: g, estimand: e }); close(population(nullGraph, e).estimate, 0);
});

test('binary repeated outcomes retain a marginal risk-difference target', () => {
  const g = exampleGraph('repeated-binary'), e = target(g), plan = native({ op: 'plan', graph: g, estimand: e });
  assert.equal(plan.method, 'trajectory'); assert.equal(e.measure, 'risk_difference');
  const truth = population(g, e, 10000); assert.ok(truth.estimate > .1 && truth.estimate < .2);
  assert.throws(() => native({ op: 'plan', graph: g, estimand: e, analysis: { method: 'lmm' } }), /continuous/);
});

test('study validation rejects invalid levels, schedules, identifiers and incompatible estimators', () => {
  const g = exampleGraph('repeated');
  for (const edit of [(x) => { x.study.times = [0, 3, 2]; }, (x) => { x.nodes[0].scope = 'cluster'; },
    (x) => { x.edges.push({ source: 'Y', target: 'X', coefficient: .2 }); },
    (x) => { x.nodes.at(-1).repeated.residualCorrelation = 1; }, (x) => { x.nodes.at(-1).lags = [{ source: 'Y', lag: 0, coefficient: 1 }]; }]) {
    const bad = structuredClone(g); edit(bad); assert.throws(() => validateGraph(bad)); assert.throws(() => simulate(bad, 10));
  }
  const clustered = exampleGraph('clustered-repeated'); assert.throws(() => native({ op: 'plan', graph: clustered, estimand: target(clustered), analysis: { method: 'mmrm' } }), /independent subjects/);
  const huge = exampleGraph('clustered-repeated'); assert.throws(() => simulate(huge, 50000), /measurement rows/);
  const d = simulate(g, 30); d.columns.A[1] = 1 - d.columns.A[0]; assert.throws(() => analyze(g, target(g), d), /changed within/);
});

test('browser and native structured planning agree on the estimator', () => {
  for (const example of ['repeated', 'clustered-repeated', 'repeated-binary']) {
    const g = exampleGraph(example), e = target(g);
    for (const analysis of [{ method: 'auto' }, { method: 'trajectory' }]) assert.equal(compilePowerPlan(g, e, analysis).method, native({ op: 'plan', graph: g, estimand: e, analysis }).method);
  }
});

test('repeated power preserves the null rate, coverage, independent-unit axis and thread reproducibility', () => {
  const g = exampleGraph('repeated'), e = target(g);
  const request = { op: 'power', graph: g, estimand: e, settings: { sampleSizes: [50, 150], repetitions: 250, populationSize: 1000, seed: 709 } };
  const r = native({ ...request, threads: 1 }), mt = native({ ...request, threads: 4 });
  assert.deepEqual(r.results, mt.results); assert.equal(r.sizeUnit, 'subjects'); assert.match(r.assumptions, /conditional/);
  assert.ok(r.results[1].alternative.probability > r.results[0].alternative.probability + .15);
  for (const row of r.results) {
    assert.equal(row.measurementRows, row.n * 4); assert.equal(row.subjects, row.n);
    assert.equal(row.alternative.failures, 0); assert.equal(row.null.failures, 0);
    assert.ok(row.null.probability < .10 && row.null.probability > .005);
    assert.ok(row.alternative.coverage > .88);
  }
  const grouped = exampleGraph('clustered-repeated');
  const result = native({ op: 'power', graph: grouped, estimand: target(grouped), settings: { sampleSizes: [30], repetitions: 50, populationSize: 1000, seed: 17 } });
  assert.equal(result.sizeUnit, 'clusters'); assert.equal(result.results[0].subjects, 360); assert.equal(result.results[0].measurementRows, 1440);
});
