import assert from 'node:assert/strict';
import test from 'node:test';
import { exampleGraph, makeNode, removeNode } from '../dag-model.js';
import { availableMeasures, validateEstimand, describeEstimand, renameEstimand, exportDesign, parseDesign } from '../estimand-model.js';
import { questionSummary, studySummary } from '../dag-workflow.js';

const target = (extra = {}) => ({ version: 1, kind: 'marginal', population: 'Eligible participants represented by the DGP',
  definition: '', outcome: 'Y', exposure: 'A', measure: 'mean_difference', reference: 0, active: 1,
  horizon: null, nullValue: 0, alternative: null, ...extra });

test('plain-language questions preserve comparison direction and distinguish intervention from coefficient targets', () => {
  const graph = exampleGraph('trial');
  assert.match(questionSummary(graph, target({ active: 0, reference: 1 })), /Treatment were set to 0 versus 1.*mean difference/);
  const coefficient = questionSummary(graph, target({ kind: 'coefficient', measure: 'coefficient' }));
  assert.match(coefficient, /coefficient of Treatment in the model for Outcome/);
  assert.doesNotMatch(coefficient, /were set|mean difference/);
  assert.doesNotThrow(() => validateEstimand(graph, target()));
});

test('plain-language questions preserve scheduled visit, change, trend and custom-weight targets', () => {
  const graph = exampleGraph('repeated');
  const describe = (repeated) => questionSummary(graph, target({ repeated }));
  assert.match(describe({ contrast: 'final', visit: 1 }), /Outcome at time 3/);
  assert.match(describe({ contrast: 'change' }), /change in Outcome from time 0 to 12/);
  assert.match(describe({ contrast: 'slope' }), /linear trend in Outcome over time/);
  assert.match(describe({ contrast: 'average' }), /averaged across scheduled visits/);
  assert.match(describe({ contrast: 'custom', weights: [-1, 0, 0, 1] }), /custom visit weights/);
  assert.match(studySummary(exampleGraph('clustered-repeated')), /size counts clusters/);
});

test('plain-language survival questions distinguish a time point from restricted follow-up', () => {
  const graph = exampleGraph('survival');
  assert.match(questionSummary(graph, target({ outcome: 'T', horizon: 28, measure: 'survival_difference' })), /at time 28.*survival probability difference/);
  assert.match(questionSummary(graph, target({ outcome: 'T', horizon: 28, measure: 'rmst_difference' })), /through time 28.*restricted mean survival time difference/);
});

test('a marginal total effect is defined independently of the graph’s mediator labels', () => {
  const graph = exampleGraph(), estimand = target();
  assert.doesNotThrow(() => validateEstimand(graph, estimand));
  const description = describeEstimand(graph, estimand);
  assert.match(description.formula, /do\(A = 1\)/);
  assert.match(description.formula, /do\(A = 0\)/);
  assert.match(description.note, /Mediators follow their generating models/);
  assert.match(description.hypothesis, /θ = 0/);
  assert.equal(graph.nodes.at(-1).interactions.length, 0, 'specifying a target does not modify the DGP');
});

test('conditional coefficients require an existing arrow and distinguish coefficient and ratio scales', () => {
  const graph = exampleGraph();
  const coefficient = target({ kind: 'coefficient', measure: 'coefficient', alternative: .3 });
  const description = describeEstimand(graph, coefficient);
  assert.match(description.formula, /current DGP value = 0.3/);
  assert.match(description.note, /conditional structural coefficient/);
  assert.match(description.note, /rather than a marginal total effect/);
  assert.throws(() => validateEstimand(graph, { ...coefficient, exposure: 'Y', outcome: 'A' }), /requires an arrow/);
  assert.throws(() => validateEstimand(graph, { ...coefficient, measure: 'exponentiated_coefficient', nullValue: 1 }), /compatible/);
  const survival = exampleGraph('survival');
  const hazardRatio = target({ kind: 'coefficient', outcome: 'T', measure: 'exponentiated_coefficient', nullValue: 1 });
  assert.match(describeEstimand(survival, hazardRatio).formula, /current DGP value = 0.7/i);
});

test('binary contrasts validate ratio and difference scales', () => {
  const graph = exampleGraph('binary');
  assert.ok(availableMeasures(graph, 'marginal', 'Y').includes('risk_ratio'));
  assert.throws(() => validateEstimand(graph, target({ measure: 'risk_difference', alternative: 1.2 })), /between/);
  assert.throws(() => validateEstimand(graph, target({ measure: 'risk_ratio', nullValue: 0 })), /positive/);
  assert.doesNotThrow(() => validateEstimand(graph, target({ measure: 'risk_ratio', nullValue: 1, alternative: .75 })));
  assert.throws(() => validateEstimand(graph, target({ active: 2 })), /Binary comparison/);
  assert.throws(() => validateEstimand(graph, target({ active: 0 })), /different/);
  assert.throws(() => validateEstimand(graph, target({ alternative: 0 })), /differ from/);
});

test('survival probability and RMST contrasts require a horizon and target latent event time', () => {
  const graph = exampleGraph('survival');
  const rmst = target({ outcome: 'T', measure: 'rmst_difference', horizon: 28, alternative: 2 });
  const description = describeEstimand(graph, rmst);
  assert.match(description.formula, /min\(T, 28\)/);
  assert.match(description.note, /latent event times/);
  assert.throws(() => validateEstimand(graph, { ...rmst, horizon: null }), /positive time horizon/);
  assert.throws(() => validateEstimand(graph, { ...rmst, alternative: 30 }), /RMST horizon/);
  assert.match(describeEstimand(graph, { ...rmst, measure: 'survival_difference', alternative: .1 }).formula, /T > 28/);
});

test('generic intervention variables can be continuous or count variables with valid levels', () => {
  const graph = exampleGraph();
  assert.doesNotThrow(() => validateEstimand(graph, target({ exposure: 'X', reference: -1, active: 2 })));
  const count = { ...graph, nodes: [...graph.nodes, makeNode('Z', 'count')] };
  assert.doesNotThrow(() => validateEstimand(count, target({ exposure: 'Z', active: 3 })));
  assert.throws(() => validateEstimand(count, target({ exposure: 'Z', active: .5 })), /nonnegative integers/);
});

test('a custom definition can be recorded without assigning special roles to variables', () => {
  const custom = target({ kind: 'custom', measure: 'custom', outcome: '', exposure: '', definition: 'Difference in the 75th percentile of event time between specified interventions.' });
  assert.doesNotThrow(() => validateEstimand(exampleGraph(), custom));
  assert.match(describeEstimand(exampleGraph(), custom).note, /estimator.*specified explicitly/);
  assert.throws(() => validateEstimand(exampleGraph(), { ...custom, definition: '' }), /custom estimand definition/);
});

test('design export stores the estimand separately from the DGP and accepts standalone DAG files', () => {
  const graph = exampleGraph(), estimand = target({ alternative: .78 });
  const design = exportDesign(graph, estimand);
  assert.deepEqual(parseDesign(JSON.stringify(design)), { graph, estimand });
  assert.deepEqual(parseDesign(JSON.stringify(graph)), { graph, estimand: null });
  const renamed = renameEstimand(estimand, 'Y', 'Outcome');
  assert.equal(renamed.outcome, 'Outcome'); assert.equal(renamed.exposure, 'A');
  assert.equal(renameEstimand(null, 'Y', 'Outcome'), null);
  assert.throws(() => validateEstimand(removeNode(graph, 'Y'), estimand), /Select an outcome/);
  assert.throws(() => parseDesign({ version: 2, graph }), /Unsupported/);
  const extra = { ...design, estimand: { ...estimand, arbitraryCode: 'not executable' } };
  assert.deepEqual(parseDesign(extra).estimand, estimand);
});
