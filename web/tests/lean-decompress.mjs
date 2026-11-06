// Exercise the real, pinned Brotli WASM decoder in its temporary worker.
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { brotliCompressSync } from 'node:zlib';
import { chromium } from 'playwright';

const root = resolve(process.argv[2] || new URL('..', import.meta.url).pathname);
// The override permits this small protocol test before packaging is complete.
// CI uses only the decoder copied into the authenticated compiler package.
const decoderRoot = process.env.SSC_BROTLI_TEST_DECODER_ROOT;
const manifest = decoderRoot ? null : JSON.parse(await readFile(resolve(root, 'lean-wasm/manifest.json'), 'utf8'));
const modulePath = decoderRoot ? resolve(decoderRoot, 'brotli_dec_wasm.js')
  : resolve(root, 'lean-wasm', manifest.decoder.module.file);
const wasmPath = decoderRoot ? resolve(decoderRoot, 'brotli_dec_wasm_bg.wasm')
  : resolve(root, 'lean-wasm', manifest.decoder.wasm.file);
const fixture = Buffer.from('Lean αβγ kernel replay and ordinary source.\n'.repeat(800));
const resources = new Map([
  ['/lean-decompress-worker.js', ['text/javascript', await readFile(resolve(root, 'lean-decompress-worker.js'))]],
  ['/brotli_wasm.js', ['text/javascript', await readFile(modulePath)]],
  ['/brotli_wasm_bg.wasm', ['application/wasm', await readFile(wasmPath)]],
  ['/fixture.br', ['application/octet-stream', brotliCompressSync(fixture)]],
  ['/empty.br', ['application/octet-stream', brotliCompressSync(Buffer.alloc(0))]],
]);
const server = createServer((req, res) => {
  res.setHeader('Content-Security-Policy', "default-src 'none'; script-src 'self' 'wasm-unsafe-eval'; worker-src 'self'; connect-src 'self'");
  const path = new URL(req.url, 'http://localhost').pathname;
  if (path === '/') return res.end('<!doctype html><title>Lean asset decompression</title>');
  const resource = resources.get(path);
  if (!resource) return res.writeHead(404).end();
  res.setHeader('Content-Type', resource[0]);
  res.end(resource[1]);
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const browser = await chromium.launch({
  executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined,
  args: ['--disable-dev-shm-usage'],
});
try {
  const page = await browser.newPage();
  const failures = [];
  page.on('pageerror', error => failures.push(error.message));
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  const result = await page.evaluate(async expectedBytes => {
    const worker = new Worker('/lean-decompress-worker.js');
    const pending = new Map();
    let serial = 0;
    worker.onmessage = ({ data }) => {
      const waiting = pending.get(data.id);
      if (!waiting) throw Error('Unexpected decompression response ID.');
      pending.delete(data.id);
      clearTimeout(waiting.timer);
      waiting.resolve(data);
    };
    const rpc = (data, transfer = []) => new Promise((resolve, reject) => {
      const id = ++serial;
      const timer = setTimeout(() => reject(Error('Decompression worker timed out.')), 10000);
      pending.set(id, { resolve, reject, timer });
      worker.postMessage({ ...data, id }, transfer);
    });
    try {
      const early = await rpc({ type: 'decompress', bytes: new Uint8Array(), expectedBytes: 0 });
      const wrongOrigin = await rpc({ type: 'initialize', moduleURL: 'https://example.invalid/decoder.js', wasmBinary: new Uint8Array() });
      const wasmBinary = new Uint8Array(await (await fetch('/brotli_wasm_bg.wasm')).arrayBuffer());
      const encoded = new Uint8Array(await (await fetch('/fixture.br')).arrayBuffer());
      const queuedBytes = encoded.slice();
      // Queue a decode before asynchronous module initialization has completed.
      const initialized = rpc({ type: 'initialize', moduleURL: '/brotli_wasm.js', wasmBinary }, [wasmBinary.buffer]);
      const queued = rpc({ type: 'decompress', bytes: queuedBytes, expectedBytes }, [queuedBytes.buffer]);
      const ready = await initialized;
      const decoded = await queued;
      const wrongSize = await rpc({ type: 'decompress', bytes: encoded.slice(), expectedBytes: expectedBytes + 1 });
      const invalidSize = await rpc({ type: 'decompress', bytes: encoded.slice(), expectedBytes: -1 });
      const invalidBytes = await rpc({ type: 'decompress', bytes: [], expectedBytes });
      const corrupt = await rpc({ type: 'decompress', bytes: new TextEncoder().encode('not brotli'), expectedBytes });
      const emptyBytes = new Uint8Array(await (await fetch('/empty.br')).arrayBuffer());
      const empty = await rpc({ type: 'decompress', bytes: emptyBytes, expectedBytes: 0 });
      const retry = await rpc({ type: 'decompress', bytes: encoded.slice(), expectedBytes });
      return {
        early, wrongOrigin, ready, wrongSize, invalidSize, invalidBytes, corrupt,
        decodedText: decoded.type === 'result' ? new TextDecoder().decode(decoded.bytes) : decoded,
        retryText: retry.type === 'result' ? new TextDecoder().decode(retry.bytes) : retry,
        emptyType: empty.type, emptyLength: empty.bytes?.byteLength,
        transferredInputDetached: queuedBytes.byteLength === 0,
        outputIsBytes: decoded.bytes instanceof Uint8Array,
      };
    } finally {
      for (const waiting of pending.values()) clearTimeout(waiting.timer);
      worker.terminate();
    }
  }, fixture.byteLength);
  assert.equal(result.early.type, 'error');
  assert.match(result.early.message, /not initialized/);
  assert.equal(result.wrongOrigin.type, 'error');
  assert.match(result.wrongOrigin.message, /same origin/);
  assert.equal(result.ready.type, 'ready');
  assert.equal(result.decodedText, fixture.toString());
  assert.equal(result.retryText, fixture.toString());
  assert.equal(result.outputIsBytes, true);
  assert.equal(result.transferredInputDetached, true);
  assert.equal(result.emptyType, 'result');
  assert.equal(result.emptyLength, 0);
  for (const key of ['wrongSize', 'invalidSize', 'invalidBytes', 'corrupt'])
    assert.equal(result[key].type, 'error', key);
  assert.match(result.wrongSize.message, /size/);
  assert.deepEqual(failures, []);
  console.log('Pinned Brotli worker: initialization, transferred bytes, corruption/size rejection and recovery passed.');
} finally {
  await browser.close();
  await new Promise(resolve => server.close(resolve));
}
