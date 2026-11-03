import assert from 'node:assert/strict';
import test from 'node:test';
import { native, simulateGraph } from './dag-engine.mjs';
import { SCENARIOS, scenarioDesign } from '../dag-scenarios.js';
import { compilePowerPlan, measuredTerms } from '../power-model.js';
import { parseGraph, validateGraph, equation, makeNode } from '../dag-model.js';
import { nodeTerms } from '../dag-functions.js';
import { validateEstimand } from '../estimand-model.js';

test('every scenario loads as a valid design, and its analysis plan matches the engine', () => {
  for (const s of SCENARIOS) {
    const d = scenarioDesign(s.id);
    if (s.id === 'empty') continue;
    assert.ok(d.estimand, s.id);
    assert.doesNotThrow(() => validateEstimand(d.graph, d.estimand), s.id);
    const plan = native({ op: 'plan', graph: d.graph, estimand: d.estimand, analysis: d.power.analysis });
    // structured designs are planned by the engine's study planner, mirrored only in outline
    if (!d.graph.study) assert.deepEqual(plan, compilePowerPlan(d.graph, d.estimand, d.power.analysis), s.id);
    else assert.equal(plan.method, compilePowerPlan(d.graph, d.estimand, d.power.analysis).method, s.id);
  }
});

test('Student t, Beta and negative binomial variables are kept with their own parameters', () => {
  const graph = parseGraph({ version: 1, nodes: [
    { ...makeNode('T', 'continuous', 10, 10), family: 'student', df: 4 },
    { ...makeNode('B', 'continuous', 10, 90), family: 'beta', precision: 6 },
    { ...makeNode('C', 'count', 10, 170), family: 'negbin', overdispersion: .4 },
    { ...makeNode('N', 'continuous', 10, 250), df: 9 }], edges: [] });
  const [t, b, c, n] = graph.nodes;
  assert.equal(t.df, 4); assert.equal(b.precision, 6); assert.equal(c.overdispersion, .4);
  assert.equal(n.df, undefined, 'a parameter its distribution does not use is not kept');
  assert.match(equation(graph, t), /Student t\(4 df\)/);
  assert.match(equation(graph, b), /Beta\(μφ, \(1 − μ\)φ\), φ = 6/);
  assert.match(equation(graph, c), /NegBin\(μ, Var = μ \+ 0\.4 μ²\)/);
  assert.throws(() => validateGraph({ ...graph, nodes: [{ ...t, df: 1 }] }), /more than 1 degree of freedom/);
  // simulated values respect each distribution's support
  const data = simulateGraph(graph, { n: 2000, seed: 3 });
  assert.ok(data.columns.B.every((x) => x > 0 && x < 1));
  assert.ok(data.columns.C.every((x) => Number.isInteger(x) && x >= 0));
});

test('an unmeasured variable is drawn dashed in equations, simulated, and replaced in analyses by its measured causes', () => {
  const d = scenarioDesign('unmeasured-confounding');
  assert.equal(d.graph.nodes.find((n) => n.id === 'U').measured, false);
  assert.equal(d.graph.nodes.find((n) => n.id === 'X').measured, undefined, 'measured variables carry no flag');
  assert.match(equation(d.graph, d.graph.nodes.find((n) => n.id === 'U')), /\(unmeasured\)$/);
  const plan = compilePowerPlan(d.graph, d.estimand, d.power.analysis);
  assert.deepEqual(plan.propensityTerms, [{ left: 'X' }]);
  assert.match(plan.notes[0], /U is unmeasured/);
  assert.ok('U' in simulateGraph(d.graph, { n: 10, seed: 1 }).columns, 'the data still carry it');
  // through a chain of unmeasured variables to the measured ones behind it
  const chain = parseGraph({ version: 1, nodes: [makeNode('Z', 'continuous', 0, 0), { ...makeNode('U', 'continuous', 0, 100), measured: false },
    { ...makeNode('V', 'continuous', 0, 200), measured: false }, makeNode('Y', 'continuous', 0, 300)],
  edges: [{ source: 'Z', target: 'U', coefficient: 1 }, { source: 'U', target: 'V', coefficient: 1 }, { source: 'V', target: 'Y', coefficient: 1 }] });
  assert.deepEqual(measuredTerms(chain, nodeTerms(chain, chain.nodes[3])), [{ left: 'Z' }]);
  assert.throws(() => validateEstimand(d.graph, { ...d.estimand, exposure: 'U' }), /must be measured/);
  assert.throws(() => compilePowerPlan(d.graph, d.estimand, { method: 'aipw', basis: 'custom', terms: [{ left: 'U' }] }), /cannot use U: it is unmeasured/);
});

test('the interaction test covers only the terms that multiply the compared variable', () => {
  const d = scenarioDesign('effect-modification');
  const plan = compilePowerPlan(d.graph, d.estimand, d.power.analysis);
  assert.equal(plan.label, 'Wald test of the interaction (1 df)');
  assert.deepEqual(plan.testedIndices.map((i) => plan.terms[i - 1]), [{ left: 'A', right: 'X' }]);
});
