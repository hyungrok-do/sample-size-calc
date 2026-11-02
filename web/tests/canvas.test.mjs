// The DAG canvas's geometry and reading of the graph (web/dag-canvas.js): plates for the study's levels, random
// effects drawn on them, moving variables between levels, switching levels on and off, layout, which arrows may
// be drawn, inferred roles and the view.
import assert from 'node:assert/strict';
import test from 'node:test';
import { exampleGraph, makeNode, validateGraph, parseGraph, addEdge, removeEdge } from '../dag-model.js';
import { NODE_W, NODE_H, CANVAS_LIMIT, ASPECT, plateLayout, plateAt, changeScope, setStructure, edgePath, loopPath, layout,
  freeSpot, connectTargets, inferRoles, contentBounds, fitView, zoomView, viewBox } from '../dag-canvas.js';

const inside = (inner, outer) => inner.x >= outer.x && inner.y >= outer.y && inner.x + inner.w <= outer.x + outer.w && inner.y + inner.h <= outer.y + outer.h;
const box = (n) => ({ x: n.x, y: n.y, w: NODE_W, h: NODE_H });

test('independent designs have no plates; structured designs nest clusters ⊃ subjects ⊃ visits', () => {
  assert.deepEqual(plateLayout(exampleGraph('mediation')), { plates: [], effects: [] });
  const repeated = exampleGraph('repeated');
  const { plates } = plateLayout(repeated);
  assert.deepEqual(plates.map((p) => p.level), ['subject', 'visit']);
  const [subject, visit] = plates;
  assert.ok(inside(visit, subject), 'the visit plate sits inside the subject plate');
  assert.deepEqual(visit.members.sort(), ['Time', 'Y']);
  for (const id of visit.members) assert.ok(inside(box(repeated.nodes.find((n) => n.id === id)), visit));
  assert.match(visit.label, /t = 0, 3, 6, 12/);
  const clustered = exampleGraph('clustered-repeated');
  const nested = plateLayout(clustered).plates;
  assert.deepEqual(nested.map((p) => p.level), ['cluster', 'subject', 'visit']);
  assert.ok(inside(nested[1], nested[0]) && inside(nested[2], nested[1]));
  assert.deepEqual(nested[0].members, ['A'], 'the randomized treatment is a cluster variable');
  assert.match(nested[0].label, /8, 12, 16/);
});

test('an empty visit plate is offered as a drop target; one visit with no visit variables has none', () => {
  const graph = setStructure(exampleGraph('trial'), { clusters: false, visits: true });
  const { plates } = plateLayout(graph);
  const visit = plates.find((p) => p.level === 'visit');
  assert.ok(visit.empty);
  for (const n of graph.nodes) assert.ok(!inside(box(n), visit), 'the empty plate does not cover a box');
  const clusteredOnly = setStructure(exampleGraph('trial'), { clusters: true, visits: false });
  assert.deepEqual(plateLayout(clusteredOnly).plates.map((p) => p.level), ['cluster', 'subject']);
});

test('random effects are drawn on the plate they vary over, above their variable', () => {
  const graph = exampleGraph('clustered-repeated');
  const y = graph.nodes.find((n) => n.id === 'Y');
  y.repeated = { ...y.repeated, subjectSD: 0.5, slopeSD: 0.1, clusterSD: 0.3 };
  const { plates, effects } = plateLayout(graph);
  const [cluster, subject, visit] = plates;
  assert.deepEqual(effects.map((e) => e.key).sort(), ['clusterSD', 'slopeSD', 'subjectSD']);
  const at = (e) => ({ x: e.cx - 14, y: e.cy - 14, w: 28, h: 28 });
  for (const e of effects.filter((e) => e.level === 'subject')) {
    assert.ok(inside(at(e), subject) && !inside(at(e), visit), `${e.key} lies between the subject and visit plates`);
    assert.ok(e.cy < y.y);
  }
  const c = effects.find((e) => e.key === 'clusterSD');
  assert.ok(inside(at(c), cluster) && !inside(at(c), subject));
  assert.match(c.text, /Cluster random intercept on Y, SD 0.3/);
});

test('a point is in the innermost plate around it', () => {
  const { plates } = plateLayout(exampleGraph('clustered-repeated'));
  const [cluster, subject, visit] = plates;
  assert.equal(plateAt(plates, { x: visit.x + 5, y: visit.y + 5 }), 'visit');
  assert.equal(plateAt(plates, { x: subject.x + 5, y: subject.y + 5 }), 'subject');
  assert.equal(plateAt(plates, { x: cluster.x + 5, y: cluster.y + 5 }), 'cluster');
  assert.equal(plateAt(plates, { x: cluster.x - 50, y: cluster.y }), null);
});

test('moving a variable between levels keeps the graph valid and clears visit-only settings', () => {
  const graph = exampleGraph('repeated');
  const y = graph.nodes.find((n) => n.id === 'Y');
  y.repeated = { ...y.repeated, subjectSD: 0.4 }; y.lags = [{ source: 'Y', lag: 1, coefficient: 0.2, initial: 0 }];
  // Y depends on visit time, so it cannot leave the visits until that arrow goes
  assert.throws(() => validateGraph(changeScope(graph, 'Y', 'subject')), /finer observation level/);
  const moved = changeScope(removeEdge(graph, 'Time', 'Y'), 'Y', 'subject');
  const after = moved.nodes.find((n) => n.id === 'Y');
  assert.equal(after.scope, 'subject');
  assert.equal(after.repeated.subjectSD, 0);
  assert.deepEqual(after.lags, []);
  assert.doesNotThrow(() => validateGraph(moved));
  const x = changeScope(graph, 'X', 'visit');
  assert.equal(x.nodes.find((n) => n.id === 'X').scope, 'visit');
  assert.doesNotThrow(() => validateGraph(x));
  assert.equal(changeScope(graph, 'Time', 'subject').nodes.find((n) => n.id === 'Time').scope, 'time', 'visit time stays on the schedule');
  // a variable cannot depend on a finer level: X at visits cannot feed a subject-level A
  const fed = { ...graph, edges: [...graph.edges, { source: 'X', target: 'A', coefficient: 0.2 }] };
  assert.throws(() => validateGraph(changeScope(fed, 'X', 'visit')), /finer observation level/);
});

test('levels switch on and off without breaking the graph', () => {
  const trial = exampleGraph('mediation');
  const visits = setStructure(trial, { clusters: false, visits: true });
  assert.deepEqual(visits.study.times, [0, 3, 6, 12]);
  assert.deepEqual(visits.study.clusterSizes, []);
  assert.ok(visits.nodes.every((n) => n.scope === 'subject'));
  assert.doesNotThrow(() => validateGraph(visits));
  const both = setStructure(visits, { clusters: true, visits: true });
  assert.deepEqual(both.study.clusterSizes, [10]);
  assert.doesNotThrow(() => validateGraph(both));
  // switching visits off moves visit variables to the subject and keeps the first time
  const repeated = exampleGraph('clustered-repeated');
  const y = repeated.nodes.find((n) => n.id === 'Y');
  y.repeated = { ...y.repeated, subjectSD: 0.5, clusterSD: 0.2, residualCorrelation: 0.3 };
  const once = setStructure(repeated, { clusters: true, visits: false });
  assert.deepEqual(once.study.times, [0]);
  assert.equal(once.nodes.find((n) => n.id === 'Y').scope, 'subject');
  assert.equal(once.nodes.find((n) => n.id === 'Y').repeated.subjectSD, 0);
  assert.doesNotThrow(() => validateGraph(once));
  // and back on starts the standard schedule; clusters off moves the cluster treatment to the subject
  const noClusters = setStructure(repeated, { clusters: false, visits: true });
  assert.equal(noClusters.nodes.find((n) => n.id === 'A').scope, 'subject');
  assert.equal(noClusters.nodes.find((n) => n.id === 'Y').repeated.clusterSD, 0);
  assert.deepEqual(noClusters.study.times, [0, 3, 6, 12]);
  assert.doesNotThrow(() => validateGraph(noClusters));
  const none = setStructure(repeated, { clusters: false, visits: false });
  assert.equal(none.study, undefined);
  assert.ok(none.nodes.every((n) => n.scope === undefined && n.lags === undefined && n.repeated === undefined));
  assert.doesNotThrow(() => validateGraph(none));
  // a design saved after the change reads back the same
  assert.deepEqual(parseGraph(JSON.stringify(both)).study, both.study);
});

test('arrows leave the side facing the other box; a lag can bow beside a same-visit arrow; a self-lag loops', () => {
  const a = { x: 0, y: 0 }, b = { x: 400, y: 0 };
  const { d, mid } = edgePath(a, b);
  assert.match(d, new RegExp(`^M ${NODE_W} ${NODE_H / 2} `));
  assert.match(d, / 400 52$/);
  assert.equal(mid[0], (NODE_W + 400) / 2);
  assert.notEqual(edgePath(a, b, 40).mid[1], mid[1]);
  const below = edgePath(a, { x: 0, y: 300 });
  assert.match(below.d, new RegExp(`^M ${NODE_W / 2} ${NODE_H} `));
  const loop = loopPath({ x: 100, y: 200 });
  assert.ok(loop.mid[1] < 200);
});

test('tidy layout: columns by generation order, levels left to right, inside the canvas', () => {
  const tidy = layout(exampleGraph('mediation'));
  const at = (id) => tidy.nodes.find((n) => n.id === id);
  assert.equal(at('X').x, at('A').x, 'roots share the first column');
  assert.ok(at('M').x > at('A').x && at('Y').x > at('M').x);
  assert.notEqual(at('X').y, at('A').y);
  assert.doesNotThrow(() => validateGraph(tidy));
  const clustered = layout(exampleGraph('clustered-repeated'));
  const pos = (id) => clustered.nodes.find((n) => n.id === id);
  assert.ok(pos('A').x < pos('X').x, 'cluster variables before subject variables');
  assert.ok(pos('X').x < pos('Time').x && pos('X').x < pos('Y').x, 'subject variables before visit variables');
  for (const n of clustered.nodes) assert.ok(n.x >= 0 && n.y >= 0 && n.x <= CANVAS_LIMIT.x && n.y <= CANVAS_LIMIT.y);
  // the laid-out plates do not cover a box of another level
  const { plates } = plateLayout(clustered);
  const visit = plates.find((p) => p.level === 'visit');
  for (const n of clustered.nodes.filter((n) => !visit.members.includes(n.id))) assert.ok(!inside(box(n), visit), `${n.id} outside the visit plate`);
  assert.doesNotThrow(() => validateGraph(clustered));
});

test('a new box goes to a free spot', () => {
  const graph = exampleGraph('mediation');
  const spot = freeSpot(graph, { x: graph.nodes[0].x, y: graph.nodes[0].y });
  for (const n of graph.nodes) assert.ok(Math.abs(n.x - spot.x) >= NODE_W || Math.abs(n.y - spot.y) >= NODE_H);
  assert.deepEqual(freeSpot({ nodes: [] }, { x: 50, y: 60 }), { x: 50, y: 60 });
});

test('arrow targets: no self-arrows, duplicates, cycles, visit-time parents or finer-to-coarser links', () => {
  const graph = exampleGraph('mediation');  // X → M, X → Y, A → M, M → Y, A → Y
  const fromY = connectTargets(graph, 'Y');
  assert.match(fromY.get('M').reason, /cycle/);
  assert.match(fromY.get('Y').reason, /two different/);
  const fromA = connectTargets(graph, 'A');
  assert.match(fromA.get('M').reason, /already/);
  assert.equal(fromA.get('X').ok, true);
  // the same answer as the model's own check
  for (const [id, { ok }] of fromA) {
    let valid = true;
    try { addEdge(graph, 'A', id); } catch { valid = false; }
    assert.equal(ok, valid, `A → ${id}`);
  }
  const repeated = exampleGraph('clustered-repeated');
  repeated.nodes.push({ ...makeNode('Z', 'continuous', 600, 40), scope: 'subject' });
  assert.match(connectTargets(repeated, 'Y').get('Z').reason, /cannot affect Z/);
  assert.throws(() => addEdge(repeated, 'Y', 'Z'), /finer observation level/);
  assert.match(connectTargets(repeated, 'X').get('Time').reason, /schedule/);
  assert.equal(connectTargets(repeated, 'A').get('X').ok, true, 'a cluster treatment can affect a subject variable');
});

test('inferred roles for the comparison', () => {
  const graph = { version: 1, nodes: ['C', 'A', 'M', 'Y', 'Z', 'P', 'K', 'Q'].map((id) => makeNode(id)), edges: [
    { source: 'C', target: 'A', coefficient: 1 }, { source: 'C', target: 'Y', coefficient: 1 },
    { source: 'A', target: 'M', coefficient: 1 }, { source: 'M', target: 'Y', coefficient: 1 },
    { source: 'Z', target: 'A', coefficient: 1 }, { source: 'P', target: 'Y', coefficient: 1 },
    { source: 'A', target: 'K', coefficient: 1 }, { source: 'Y', target: 'K', coefficient: 1 }] };
  const roles = inferRoles(graph, 'A', 'Y');
  assert.deepEqual(Object.fromEntries(roles), { C: 'confounder', M: 'mediator', Z: 'instrument', P: 'prognostic', K: 'collider' });
  assert.equal(roles.has('Q'), false);
  assert.equal(inferRoles(graph, 'A', 'A').size, 0);
  assert.equal(inferRoles(graph, 'A', 'missing').size, 0);
});

test('the view fits the content, zooms about a point and stays within limits', () => {
  const graph = exampleGraph('clustered-repeated');
  const { plates, effects } = plateLayout(graph);
  const content = contentBounds(graph, plates, effects);
  for (const p of plates) assert.ok(inside(p, content));
  const view = fitView(content);
  assert.ok(view.w >= content.w && view.w * ASPECT >= content.h - 1e-9);
  assert.ok(view.x <= content.x && view.y <= content.y);
  assert.equal(fitView({ x: 0, y: 0, w: 100, h: 50 }).w, 1100, 'a small graph is not blown up');
  const zoomed = zoomView({ x: 0, y: 0, w: 1000 }, 2, { x: 500, y: 100 });
  assert.equal(zoomed.w, 500);
  assert.equal(zoomed.x, 250, 'the point under the pointer stays put');
  assert.equal(zoomed.y, 50);
  assert.equal(zoomView({ x: 0, y: 0, w: 1000 }, 1000).w, 360);
  assert.equal(viewBox({ x: 1, y: 2, w: 1100 }), '1.0 2.0 1100.0 580.0');
});

test('positions beyond the old 920 × 440 area are valid; negative or huge ones are not', () => {
  const graph = { version: 1, nodes: [makeNode('A', 'binary', 2000, 1200)], edges: [] };
  assert.doesNotThrow(() => validateGraph(graph));
  assert.throws(() => validateGraph({ ...graph, nodes: [makeNode('A', 'binary', -1, 0)] }), /canvas position/);
  assert.throws(() => validateGraph({ ...graph, nodes: [makeNode('A', 'binary', CANVAS_LIMIT.x + 1, 0)] }), /canvas position/);
});
