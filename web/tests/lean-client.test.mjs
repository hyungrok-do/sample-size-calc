import test from 'node:test';
import assert from 'node:assert/strict';
import { BrowserLeanCompiler } from '../lean-wasm-client.js';

test('fatal runtime errors discard the worker and the next request initializes a fresh one', async () => {
  const oldWorker = globalThis.Worker, oldIsolation = globalThis.crossOriginIsolated;
  const oldStorage = globalThis.sessionStorage;
  const workers = [];
  globalThis.crossOriginIsolated = true;
  globalThis.sessionStorage = { removeItem() {} };
  globalThis.Worker = class {
    constructor() { workers.push(this); }
    postMessage(data) {
      if (data.type === 'initialize') queueMicrotask(() => this.onmessage({ data: { type: 'ready' } }));
      else this.request = data;
    }
    terminate() { this.terminated = true; }
  };
  const compiler = new BrowserLeanCompiler(() => {});
  try {
    await compiler.initialize();
    const rejected = assert.rejects(compiler.compile('import SSC'), /stack overflow/);
    await new Promise(resolve => setImmediate(resolve));
    workers[0].onmessage({ data: { type: 'error', id: workers[0].request.id,
      fatal: true, message: 'stack overflow' } });
    await rejected;
    assert.equal(workers[0].terminated, true);
    assert.equal(compiler.ready, false);
    await compiler.initialize();
    assert.equal(workers.length, 2);
    assert.equal(compiler.ready, true);
  } finally {
    compiler.stop();
    globalThis.Worker = oldWorker;
    globalThis.crossOriginIsolated = oldIsolation;
    globalThis.sessionStorage = oldStorage;
  }
});
