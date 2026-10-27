import assert from 'node:assert/strict';
import test from 'node:test';
import { readFile } from 'node:fs/promises';
import { addEdge, removeEdge, removeNode, renameNode, emptyGraph, exampleGraph, makeNode, parseGraph, validateGraph,
  topologicalOrder, simulateGraph, equation, dataCsv, DISTRIBUTIONS } from '../dag-model.js';
const near = (actual, expected, tolerance, label = '') => assert.ok(Math.abs(actual - expected) < tolerance,
  `${label}: expected ${expected} ± ${tolerance}, got ${actual}`);
const simple = (family = 'gaussian', intercept = 0, coefficient = 0) => ({ version: 1,
  nodes: [makeNode('A', 'binary'), { ...makeNode('Y', DISTRIBUTIONS[family].kind), family, intercept }],
  edges: [{ source: 'A', target: 'Y', coefficient }] });
function groupMean(result, id, arm) {
  const values = result.columns[id], a = result.columns.A;
  let sum = 0, n = 0;
  for (let i = 0; i < values.length; i++) if (a[i] === arm) { sum += values[i]; n++; }
  return sum / n;
}

test('any variable DAG is allowed: types and labels do not prescribe trial roles', () => {
  assert.doesNotThrow(() => validateGraph(emptyGraph()));
  assert.throws(() => simulateGraph(emptyGraph()), /at least one variable/);
  const graph = { version: 1, nodes: [makeNode('X'), { ...makeNode('A', 'binary'), label: 'Treatment' }], edges: [] };
  assert.doesNotThrow(() => simulateGraph(addEdge(graph, 'X', 'A'), { n: 1 }));
  assert.doesNotThrow(() => simulateGraph({ version: 1, nodes: [makeNode('X')], edges: [] }));
  assert.equal(removeNode({ version: 1, nodes: [makeNode('X')], edges: [] }, 'X').nodes.length, 0);
});
test('topological generation rejects cycles, duplicates and missing endpoints', () => {
  const graph = exampleGraph();
  for (const nodes of [graph.nodes, [...graph.nodes].reverse()]) {
    const order = topologicalOrder({ ...graph, nodes });
    for (const edge of graph.edges) assert.ok(order.indexOf(edge.source) < order.indexOf(edge.target));
  }
  assert.throws(() => addEdge(graph, 'Y', 'M'), /cycle/);
  assert.throws(() => addEdge(graph, 'Y', 'Y'), /itself/);
  assert.throws(() => addEdge(graph, 'X', 'M'), /already connected/);
  assert.throws(() => addEdge(graph, 'Z', 'Y'), /missing variable/);
});
test('model JSON round-trip and variable renaming preserve all dependencies', () => {
  const graph = exampleGraph();
  graph.nodes.at(-1).interactions.push({ left: 'A', right: 'X', coefficient: 0.2 });
  assert.deepEqual(parseGraph(JSON.stringify(graph)), graph);
  assert.deepEqual(parseGraph({ ...graph, expression: 'arbitrary', nodes: graph.nodes.map((n) => ({ ...n, code: 'arbitrary' })) }), graph);
  const renamed = renameNode(graph, 'X', 'Age');
  assert.equal(renamed.nodes[0].id, 'Age');
  assert.ok(renamed.edges.some((e) => e.source === 'Age'));
  assert.equal(renamed.nodes.at(-1).interactions[0].right, 'Age');
  const survival = renameNode(exampleGraph('survival'), 'C', 'Censor');
  assert.equal(survival.nodes.find((n) => n.id === 'T').censoring.variable, 'Censor');
  assert.deepEqual(parseGraph(JSON.stringify(survival)), survival);
  assert.throws(() => renameNode(graph, 'X', 'Y'), /unique/);
  assert.throws(() => parseGraph({ ...graph, version: 2 }), /version/);
  assert.throws(() => parseGraph('not json'), SyntaxError);
  for (const id of ['x-y', '1x', '__proto__', '', null, undefined, true]) {
    const invalid = simple(); invalid.nodes[1].id = id;
    assert.throws(() => parseGraph(invalid), /IDs/);
  }
});
test('invalid distributions, parameters and numerical requests fail without clipping', () => {
  for (const sd of [0, -1, Infinity]) {
    const g = simple(); g.nodes[1].sd = sd;
    assert.throws(() => simulateGraph(g), /sd must be positive/);
  }
  const mismatch = simple(); mismatch.nodes[1].kind = 'survival';
  assert.throws(() => validateGraph(mismatch), /compatible/);
  for (const n of [0, 2.5, 50001, NaN]) assert.throws(() => simulateGraph(simple(), { n }), /Sample size/);
  for (const seed of [-1, 0.5, 4294967296, NaN]) assert.throws(() => simulateGraph(simple(), { seed }), /Seed/);
  for (const intercept of [Math.log(1001), -1000, 1000]) assert.throws(() => simulateGraph(simple('poisson', intercept)), /Poisson rate/);
  const huge = simple(); huge.nodes[1].intercept = Infinity;
  assert.throws(() => validateGraph(huge), /finite intercept/);
});
test('seeded data and CSV exports are reproducible for both seed endpoints', () => {
  const graph = exampleGraph();
  const a = simulateGraph(graph, { n: 250, seed: 0 });
  assert.deepEqual(a, simulateGraph(parseGraph(JSON.stringify(graph)), { n: 250, seed: 0 }));
  assert.notDeepEqual(a.columns.Y, simulateGraph(graph, { n: 250, seed: 4294967295 }).columns.Y);
  const csv = dataCsv(graph, a).trim().split('\n');
  assert.equal(csv.length, 251); assert.equal(csv[0], 'X,A,M,Y');
  assert.deepEqual(csv[1].split(',').map(Number), graph.nodes.map((node) => a.columns[node.id][0]));
});
test('Gaussian mediation agrees with analytic marginal mean, variance and total effect', () => {
  const result = simulateGraph(exampleGraph(), { n: 50000, seed: 371 });
  const y = result.summaries.find((s) => s.id === 'Y');
  near(groupMean(result, 'Y', 1) - groupMean(result, 'Y', 0), .78, .04, 'total effect');
  near(y.mean, .39, .025, 'mean'); near(y.sd ** 2, 1.9217, .05, 'variance');
  near(groupMean(result, 'X', 1) - groupMean(result, 'X', 0), 0, .04, 'independent root variables');
});
test('Bernoulli conditional probabilities follow the logistic link', () => {
  const result = simulateGraph(simple('bernoulli', -1, Math.log(2)), { n: 50000, seed: 19 });
  near(groupMean(result, 'Y', 0), 1 / (1 + Math.exp(1)), .015);
  near(groupMean(result, 'Y', 1), 1 / (1 + Math.exp(1 - Math.log(2))), .015);
  assert.ok(result.columns.Y.every((v) => v === 0 || v === 1));
  assert.ok(simulateGraph(simple('bernoulli', 1000), { n: 100 }).columns.Y.every((v) => v === 1));
  assert.ok(simulateGraph(simple('bernoulli', -1000), { n: 100 }).columns.Y.every((v) => v === 0));
});
test('Poisson means and variances match log-link specifications', () => {
  const result = simulateGraph(simple('poisson', Math.log(3), Math.log(2)), { n: 50000, seed: 7 });
  near(groupMean(result, 'Y', 0), 3, .06); near(groupMean(result, 'Y', 1), 6, .08);
  assert.ok(result.columns.Y.every((v) => Number.isInteger(v) && v >= 0));
  const high = simulateGraph(simple('poisson', Math.log(100)), { n: 50000, seed: 8 }).summaries[1];
  near(high.mean, 100, .3); near(high.sd ** 2, 100, 2);
});
test('Uniform, Log-normal and Gamma distributions match their parameterizations', () => {
  const uniform = simple('uniform', 3); uniform.nodes[1].width = 4;
  const u = simulateGraph(uniform, { n: 50000, seed: 14 });
  near(u.summaries[1].mean, 3, .02); near(u.summaries[1].sd ** 2, 16 / 12, .025);
  assert.ok(u.columns.Y.every((v) => v > 1 && v < 5));
  const lognormal = simple('lognormal', .2); lognormal.nodes[1].sd = .6;
  const l = simulateGraph(lognormal, { n: 50000, seed: 15 });
  near(l.summaries[1].mean, Math.exp(.2 + .36 / 2), .025);
  for (const shape of [.5, 4]) {
    const g = simple('gamma', .4); g.nodes[1].shape = shape;
    const result = simulateGraph(g, { n: 50000, seed: 16 }).summaries[1];
    near(result.mean, Math.exp(.4), .035); near(result.sd ** 2, Math.exp(.8) / shape, .12);
  }
});
test('Weibull and exponential survival obey the specified PH survival functions', () => {
  for (const family of ['weibull', 'exponential']) {
    const graph = simple(family, 0, Math.log(2));
    graph.nodes[1].shape = 1.7; graph.nodes[1].scale = 2.5;
    const result = simulateGraph(graph, { n: 50000, seed: 28 });
    const shape = family === 'exponential' ? 1 : 1.7;
    for (const arm of [0, 1]) {
      let count = 0, alive = 0;
      for (let i = 0; i < result.n; i++) if (result.columns.A[i] === arm) { count++; alive += result.columns.Y[i] > 2; }
      near(alive / count, Math.exp(-((2 / 2.5) ** shape) * (arm ? 2 : 1)), .015, family);
    }
    assert.ok(result.columns.Y.every((v) => v > 0));
    assert.deepEqual(result.columns.Y_time, result.columns.Y);
    assert.ok(result.columns.Y_event.every((v) => v === 1));
  }
});
test('administrative and variable censoring produce consistent observed times and indicators', () => {
  const graph = simple('exponential');
  graph.nodes[1].scale = 5; graph.nodes[1].censoring = { mode: 'administrative', time: 3, variable: '' };
  const result = simulateGraph(graph, { n: 50000, seed: 35 });
  near(result.survival[0].events / result.n, 1 - Math.exp(-3 / 5), .015);
  for (let i = 0; i < result.n; i++) {
    assert.equal(result.columns.Y_time[i], Math.min(result.columns.Y[i], 3));
    assert.equal(result.columns.Y_event[i], Number(result.columns.Y[i] <= 3));
  }
  const variable = exampleGraph('survival');
  const v = simulateGraph(variable, { n: 1000, seed: 36 });
  for (let i = 0; i < v.n; i++) {
    assert.equal(v.columns.T_time[i], Math.min(v.columns.T[i], v.columns.C[i]));
    assert.equal(v.columns.T_event[i], Number(v.columns.T[i] <= v.columns.C[i]));
  }
  assert.ok(v.survival[0].events > 0 && v.survival[0].censored > 0);
  assert.equal(dataCsv(variable, v).split('\n')[0], 'X,A,T,C,T_time,T_event,C_time,C_event');
  const removed = removeNode(variable, 'C');
  assert.equal(removed.nodes.find((n) => n.id === 'T').censoring.mode, 'none');
  assert.doesNotThrow(() => validateGraph(removed));
  const collision = { ...variable, nodes: [...variable.nodes, makeNode('T_time')] };
  assert.throws(() => validateGraph(collision), /conflicts/);
  const bad = simple('exponential'); bad.nodes.push({ ...makeNode('C'), intercept: -100 });
  bad.nodes[1].censoring = { mode: 'variable', variable: 'C', time: 20 };
  assert.throws(() => simulateGraph(bad), /censoring times must be nonnegative/);
});
test('interactions use both parents, and deleting dependencies removes invalid terms', () => {
  const graph = exampleGraph();
  const y = graph.nodes.at(-1); y.sd = 1e-10;
  y.interactions = [{ left: 'A', right: 'X', coefficient: -.7 }];
  const { columns: c } = simulateGraph(graph, { n: 100, seed: 3 });
  for (let i = 0; i < c.Y.length; i++) near(c.Y[i], .4 * c.X[i] + .6 * c.M[i] + .3 * c.A[i] - .7 * c.X[i] * c.A[i], 1e-8);
  assert.match(equation(graph, y), /− 0.7 × A × X/);
  const invalid = structuredClone(graph); invalid.nodes.at(-1).interactions[0].left = 'Y';
  assert.throws(() => validateGraph(invalid), /two different parents/);
  const duplicate = structuredClone(graph); duplicate.nodes.at(-1).interactions.push({ left: 'X', right: 'A', coefficient: .5 });
  assert.throws(() => validateGraph(duplicate), /interaction already/);
  assert.equal(removeNode(graph, 'X').nodes.at(-1).interactions.length, 0);
  assert.equal(removeEdge(graph, 'X', 'Y').nodes.at(-1).interactions.length, 0);
  assert.equal(graph.nodes.at(-1).interactions.length, 1, 'original graph is unchanged');
});
test('the simulator page is independent and its offline header matches its policy', async () => {
  const base = new URL('../', import.meta.url);
  const [html, sw, app] = await Promise.all(['dag.html', 'sw.js', 'app.js'].map((name) => readFile(new URL(name, base), 'utf8')));
  const meta = /http-equiv="Content-Security-Policy" content="([^"]+)"/.exec(html)[1];
  const expression = /const DAG_POLICY = ((?:"[^"]*"\s*\+?\s*)+);/.exec(sw)[1];
  assert.equal([...expression.matchAll(/"([^"]*)"/g)].map((m) => m[1]).join(''), meta + "; frame-ancestors 'none'");
  assert.doesNotMatch(html, /src="app\.js|unsafe-eval|wasm-unsafe-eval|id="dag-add-role"/);
  assert.doesNotMatch(app, /dag|mountDag/);
  for (const name of ['dag.html', 'dag.css', 'dag-entry.js', 'dag.js', 'dag-model.js', 'dag-worker.js']) assert.ok(sw.includes(`${name}?v=dev`));
});
