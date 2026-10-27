import assert from 'node:assert/strict';
import test from 'node:test';
import createThreadedEngine from '../engine/ssc_web_mt.mjs';
import { exampleGraph } from '../dag-model.js';
import { native } from './dag-engine.mjs';

test('the two shipped WASM builds return identical power and report native progress', async () => {
  const M = await createThreadedEngine();
  const request = { op: 'power', graph: exampleGraph(), estimand: { version: 1, kind: 'marginal',
    population: 'Graph population', definition: '', exposure: 'A', outcome: 'Y', measure: 'mean_difference',
    reference: 0, active: 1, horizon: null, nullValue: 0, alternative: .5 },
  settings: { sampleSizes: [100, 300], repetitions: 80, populationSize: 1000 },
  threads: Math.min(4, M._ssc_hardware_threads()) };
  const expected = native(request), progress = [];
  M.onDagProgress = (event) => progress.push(event);
  const text = JSON.stringify(request), size = M.lengthBytesUTF8(text) + 1, input = M._malloc(size);
  M.stringToUTF8(text, input, size); let ptr;
  try {
    ptr = M._ssc_dag_json(input, 0);
    const actual = JSON.parse(M.UTF8ToString(ptr));
    assert.equal(actual.ok, true, actual.error);
    assert.equal(actual.result.engine.threads, request.threads);
    assert.deepEqual(actual.result.results, expected.results);
    assert.deepEqual(actual.result.alternative, expected.alternative);
    assert.deepEqual(actual.result.nullScenario, expected.nullScenario);
    assert.ok(progress.some((event) => event.stage === 'calibration'));
    assert.equal(progress.at(-1).completed, progress.at(-1).total);
    assert.equal(progress.at(-1).completed, 320);
  } finally { if (ptr) M._ssc_free(ptr); M._free(input); M.onDagProgress = null; }
});
