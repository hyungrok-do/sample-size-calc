import assert from 'node:assert/strict';
import test from 'node:test';
import { makeNode, parseGraph, renameNode, removeEdge, validateGraph } from '../dag-model.js';
import { parseAnalysis, transformValue, parseFunction } from '../dag-functions.js';
import { compilePowerPlan } from '../power-model.js';
import { exportDesign, parseDesign } from '../estimand-model.js';
import { native, simulateGraph, runPower } from './dag-engine.mjs';

const near = (x, y, tol = 1e-8) => assert.ok(Math.abs(x - y) <= tol, `${x} != ${y} within ${tol}`);
const target = (extra = {}) => ({ version: 1, kind: 'marginal', population: 'Simulated population', definition: '',
  exposure: 'A', outcome: 'Y', measure: 'mean_difference', reference: 0, active: 1, horizon: null, nullValue: 0, alternative: null, ...extra });
const settings = (extra = {}) => ({ sampleSizes: [400], repetitions: 300, populationSize: 20000, checkNull: true, seed: 941, ...extra });
const quadratic = () => ({ version: 1, nodes: [makeNode('A'), { ...makeNode('Y'), functions: [{ left: 'A', transform: { kind: 'square' }, coefficient: .8 }] }],
  edges: [{ source: 'A', target: 'Y', coefficient: .1 }] });
const confounded = (binary = false) => ({ version: 1,
  nodes: [{ ...makeNode('X'), family: 'uniform', width: 4 }, { ...makeNode('A', 'binary'), intercept: -.5, functions: [{ left: 'X', transform: { kind: 'square' }, coefficient: .3 }] },
    { ...makeNode('Y', binary ? 'binary' : 'continuous'), intercept: binary ? -1 : 0, functions: [{ left: 'X', transform: { kind: 'square' }, coefficient: .7 }] }],
  edges: [{ source: 'X', target: 'A', coefficient: .2 }, { source: 'X', target: 'Y', coefficient: .2 }, { source: 'A', target: 'Y', coefficient: .4 }] });

test('nonlinear generation agrees with explicit functions and transformed interactions', () => {
  const graph = { version: 1, nodes: [makeNode('X'), makeNode('Z'), makeNode('Y')],
    edges: [{ source: 'X', target: 'Y', coefficient: .2 }, { source: 'Z', target: 'Y', coefficient: -.1 }] };
  const functions = [{ kind: 'square', center: .4, scale: 2 }, { kind: 'threshold', knot: .2 }, { kind: 'hinge', knot: -.5 },
    { kind: 'saturation', center: -.2, scale: .8 }, { kind: 'spline', knots: [-2, -1, .5, 2], index: 0 }, { kind: 'spline', knots: [-2, -1, .5, 2], index: 1 }];
  const base = simulateGraph(graph, { n: 300, seed: 4 });
  for (const transform of functions) {
    graph.nodes[2].functions = [{ left: 'X', transform, right: 'Z', rightTransform: { kind: 'saturation', center: 0, scale: 1 }, coefficient: .6 }];
    const data = simulateGraph(graph, { n: 300, seed: 4 });
    for (let i = 0; i < 300; i++) near(data.columns.Y[i] - base.columns.Y[i], .6 * transformValue(data.columns.X[i], transform) * Math.tanh(data.columns.Z[i]));
  }
  for (const transform of functions.filter((t) => t.kind === 'spline')) {
    near(transformValue(20, transform) - transformValue(10, transform), transformValue(30, transform) - transformValue(20, transform));
    near(transformValue(-10, transform), 0);
  }
});

test('response terms validate, rename, delete and survive design export', () => {
  const graph = parseGraph(quadratic());
  assert.deepEqual(parseGraph(graph), graph);
  const renamed = renameNode(graph, 'A', 'Dose');
  assert.equal(renamed.nodes[1].functions[0].left, 'Dose');
  assert.equal(removeEdge(renamed, 'Dose', 'Y').nodes[1].functions.length, 0);
  const config = settings({ analysis: { method: 'gformula', basis: 'linear', models: { Y: { basis: 'spline', knots: [-2, 0, 2] } } } });
  const design = exportDesign(graph, target(), config);
  assert.deepEqual({ version: 1, ...parseDesign(JSON.stringify(design)) }, design);
  for (const bad of [{ kind: 'expression', code: '1' }, { kind: 'spline', knots: [0, 0, 1] }, { kind: 'saturation', scale: 0 }]) {
    const g = structuredClone(graph); g.nodes[1].functions[0].transform = bad;
    assert.throws(() => validateGraph(g)); assert.throws(() => native({ op: 'simulate', graph: g }));
  }
  const bad = structuredClone(graph); bad.nodes[1].functions[0].left = 'Y';
  assert.throws(() => validateGraph(bad)); assert.throws(() => native({ op: 'simulate', graph: bad }));
});

test('the UI and C++ agree on independent nonlinear analysis specifications', () => {
  for (const graph of [quadratic(), confounded(), confounded(true)]) {
    const e = target(graph.nodes.at(-1).kind === 'binary' ? { measure: 'risk_difference' } : {});
    for (const basis of ['match', 'linear', 'quadratic', 'spline']) {
      const analysis = parseAnalysis({ method: graph.nodes[1].id === 'A' ? 'aipw' : 'gformula', basis, knots: [-2, 0, 2] });
      assert.deepEqual(native({ op: 'plan', graph, estimand: e, analysis }), compilePowerPlan(graph, e, analysis));
    }
  }
  for (const method of ['aipw', 'tmle']) {
    const graph = confounded(true), e = target({ measure: 'risk_difference' }), analysis = parseAnalysis({ method, folds: 3, basis: 'linear', propensity: { basis: 'quadratic' } });
    assert.deepEqual(native({ op: 'plan', graph, estimand: e, analysis }), compilePowerPlan(graph, e, analysis));
  }
  for (const measure of ['overall_effect', 'nonlinearity']) {
    const graph = quadratic(), e = target({ kind: 'joint', measure }), analysis = parseAnalysis({ basis: 'quadratic' });
    assert.deepEqual(native({ op: 'plan', graph, estimand: e, analysis }), compilePowerPlan(graph, e, analysis));
  }
});

test('quadratic truth is recovered by g-computation and exposes a misspecified linear analysis', () => {
  const graph = quadratic(), e = target();
  const matched = runPower(graph, e, settings({ analysis: { method: 'gformula', basis: 'quadratic' } }));
  const linear = runPower(graph, e, settings({ analysis: { method: 'gformula', basis: 'linear' } }));
  const m = matched.results[0].alternative, l = linear.results[0].alternative;
  near(matched.alternative.truth.estimate, .9);
  near(m.meanEstimate, .9, .025); assert.ok(m.coverage > .90 && m.coverage < .99);
  assert.ok(Math.abs(l.bias) > .65 && l.coverage < .1);
  assert.ok(matched.results[0].null.probability < .10);
  assert.deepEqual(graph, quadratic());
});

test('joint tests detect a U-shaped response while a linear block has little power', () => {
  const graph = quadratic(); graph.edges[0].coefficient = 0;
  const e = target({ kind: 'joint', measure: 'overall_effect' });
  const curved = runPower(graph, e, settings({ analysis: { basis: 'quadratic' } }));
  const linear = runPower(graph, e, settings({ analysis: { basis: 'linear' } }));
  assert.ok(curved.results[0].alternative.probability > .95);
  assert.ok(linear.results[0].alternative.probability < .12);
  assert.ok(curved.results[0].null.probability < .10);
  assert.equal(curved.results[0].alternative.coverage, null);
  const nonlinearity = runPower(graph, target({ kind: 'joint', measure: 'nonlinearity' }), settings({ analysis: { basis: 'quadratic' } }));
  assert.ok(nonlinearity.results[0].alternative.probability > .95);
  assert.throws(() => compilePowerPlan(graph, target({ kind: 'joint', measure: 'nonlinearity' }), { basis: 'linear' }), /no terms/);
});

test('cross-fitted AIPW exposes double misspecification and recovers with either nuisance model', () => {
  const graph = confounded(), e = target();
  for (const [outcome, propensity] of [['quadratic', 'quadratic'], ['linear', 'quadratic'], ['quadratic', 'linear'], ['linear', 'linear']]) {
    const answer = runPower(graph, e, settings({ sampleSizes: [800], analysis: { method: 'aipw', basis: outcome, propensity: { basis: propensity }, folds: 2 } }));
    const a = answer.results[0].alternative;
    if (outcome === 'linear' && propensity === 'linear') assert.ok(Math.abs(a.bias) > .15);
    else near(a.meanEstimate, .4, .045);
    assert.ok(a.valid > 290);
  }
});

test('cross-fitted binary TMLE gives bounded risk contrasts and calibrated inference', () => {
  const graph = confounded(true), e = target({ measure: 'risk_difference' });
  const answer = runPower(graph, e, settings({ sampleSizes: [800], analysis: { method: 'tmle', basis: 'quadratic', folds: 2 } }));
  const a = answer.results[0].alternative;
  assert.ok(a.valid > 290); assert.ok(Math.abs(a.bias) < .015);
  assert.ok(a.coverage > .90 && a.coverage < .99);
  assert.ok(answer.results[0].null.probability < .10);
  assert.equal(answer.analysis.method, 'tmle');
  for (const seed of [3, 29]) {
    const data = simulateGraph(graph, { n: 500, seed });
    const result = native({ op: 'analyze', graph, estimand: e, data, analysis: { method: 'tmle', basis: 'quadratic', folds: 2 } });
    assert.ok(result.estimate >= -1 && result.estimate <= 1);
  }
});

test('the average-effect reference retains nonlinear subgroup effects', () => {
  const graph = { version: 1, nodes: [makeNode('X'), makeNode('A', 'binary'), { ...makeNode('Y'),
    functions: [{ left: 'X', transform: { kind: 'square' }, right: 'A', coefficient: .5 }] }],
    edges: [{ source: 'A', target: 'Y', coefficient: .3 }, { source: 'X', target: 'Y', coefficient: .1 }] };
  const answer = runPower(graph, target(), settings({ nullMode: 'average_effect', populationSize: 50000 }));
  const reference = answer.nullScenario;
  assert.equal(reference.graph.nodes[2].functions[0].coefficient, .5);
  near(reference.graph.edges[0].coefficient, -.5, .015);
  assert.ok(Math.abs(reference.truth.estimate) < .02);
  assert.equal(reference.exactNull, false);
  assert.ok(answer.results[0].null.probability < .10);
  const sharp = native({ op: 'nullGraph', graph, estimand: target() });
  assert.equal(sharp.nodes[2].functions[0].coefficient, 0);
});
