import assert from 'node:assert/strict';
import test from 'node:test';
import createEngine from '../engine/ssc_web.mjs';
import { guideData } from '../guide-data.js';
import { GuideState } from '../guide-state.js';
import { prepareSelection, selectionHash, parseSelectionHash } from '../selection.js';

const M = await createEngine();
function take(ptr) {
  try { return JSON.parse(M.UTF8ToString(ptr)); } finally { M._ssc_free(ptr); }
}
const catalog = take(M._ssc_catalog_json());
const schemas = new Map(catalog.procedures.map(({ id }) => {
  const ptr = M._malloc(M.lengthBytesUTF8(id) + 1);
  try {
    M.stringToUTF8(id, ptr, M.lengthBytesUTF8(id) + 1);
    return [id, take(M._ssc_procedure_json(ptr))];
  } finally { M._free(ptr); }
}));
const routes = [];
function visit(node, path = []) {
  if (node.type === 'question') node.choices.forEach((choice, i) => visit(choice.next, [...path, [node.node_id, String(i)]]));
  else routes.push({ node, path });
}
visit(guideData.root);
function follow(route) {
  const state = new GuideState(guideData);
  route.path.forEach(([id, choice]) => state.choose(id, choice));
  return state;
}

test('every engine procedure is reachable and every route transfers valid schema options', () => {
  const seen = new Set();
  for (const route of routes) {
    const state = follow(route);
    assert.equal(state.snapshot().selection, null, 'additional requirements must be answered');
    state.toggle('none', true);
    const { selection } = state.snapshot();
    if (route.node.type !== 'candidate') { assert.equal(selection, null); continue; }
    seen.add(selection.procedure);
    assert.deepEqual(selection.inputs, route.node.proposed_option_values);
    assert.ok(prepareSelection(schemas.get(selection.procedure), selection));
    assert.deepEqual(parseSelectionHash(selectionHash(selection)), selection);
  }
  assert.deepEqual([...seen].sort(), [...schemas.keys()].sort());
});

test('additional requirements are checked together, with no silent fallback', () => {
  const route = routes.find((r) => r.node.procedure_id === 'means.two.t');
  const state = follow(route);
  state.toggle('cluster', true);
  state.toggle('interim', true);
  assert.deepEqual([...state.extra], ['cluster', 'interim']);
  assert.equal(state.snapshot().selection, null);
  assert.equal(state.snapshot().blockers.length, 2);
  state.toggle('none', true);
  assert.deepEqual([...state.extra], ['none']);
  assert.ok(state.snapshot().selection);
  state.toggle('unsure', true);
  assert.equal(state.snapshot().selection, null);
});

test('changing an ancestor clears dependent choices, retaining additional requirements', () => {
  const route = routes.find((r) => r.node.procedure_id === 'prop.two');
  const state = follow(route);
  state.toggle('cluster', true);
  state.choose(guideData.root.node_id, 'unsure');
  assert.equal(state.answers.size, 1);
  assert.equal(state.snapshot().unsure, true);
  assert.deepEqual([...state.extra], ['cluster']);
  state.choose(guideData.root.node_id, route.path[0][1]);
  assert.equal(state.snapshot().node.type, 'question');
  state.reset();
  assert.equal(state.answers.size, 0);
  assert.equal(state.extra.size, 0);
});

test('unknown, conditional, or numerical presets fail instead of silently applying defaults', () => {
  const proc = schemas.get('prop.two');
  assert.throws(() => prepareSelection(proc, { procedure: proc.id, inputs: { test: 'unknown' } }));
  assert.throws(() => prepareSelection(proc, { procedure: proc.id, inputs: { p1: '0.7' } }));
  assert.throws(() => prepareSelection(proc, { procedure: proc.id, solve_for: 'unknown' }));
  assert.throws(() => prepareSelection(proc, { procedure: 'unknown' }));
  const conditional = { id: 'conditional', solve_for: [{ value: 'n' }], params: [
    { id: 'analysis', type: 'select', default: 'a', options: [{ value: 'a' }, { value: 'b' }] },
    { id: 'test', type: 'select', default: 't', visible: [[{ key: 'analysis', values: ['b'] }]], options: [{ value: 't' }] },
  ] };
  assert.throws(() => prepareSelection(conditional, { procedure: conditional.id, inputs: { test: 't' } }));
  assert.ok(prepareSelection(conditional, { procedure: conditional.id, inputs: { test: 't', analysis: 'b' } }));
  assert.equal(parseSelectionHash('#%').procedure, '');
});

test('Fisher explicitly requests exact enumeration', () => {
  const route = routes.find((r) => r.node.procedure_id === 'prop.two' && r.node.proposed_option_values.test === 'fisher');
  assert.equal(route.node.proposed_option_values.method, 'exact');
});
