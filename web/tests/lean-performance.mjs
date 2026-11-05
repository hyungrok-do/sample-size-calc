// Reproducible browser measurements. The production deadline remains 120 s.
// An explicit larger deadline is for diagnosing baselines that already time out.
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { createHash } from 'node:crypto';
import { readFile, writeFile } from 'node:fs/promises';
import { extname, resolve, sep } from 'node:path';
import { chromium } from 'playwright';

const root = resolve(process.argv[2] || new URL('..', import.meta.url).pathname);
const deadlineMs = Number(process.env.SSC_PERFORMANCE_DEADLINE_MS || 120000);
assert.ok(Number.isSafeInteger(deadlineMs) && deadlineMs > 0 && deadlineMs <= 600000);
const files = process.argv.slice(3);
const cases = files.length ? await Promise.all(files.map(async file => ({ name: file, source: await readFile(file, 'utf8') }))) : [
  ...[100, 500, 1000].map(n => ({ name: `${n} declarations`,
    source: 'import SSC\n' + Array.from({ length: n }, (_, i) => `def value${i} : Nat := ${i}`).join('\n') })),
  ...await Promise.all(['bernoulli-observations', 'paired-categorical', 'dependent-observations'].map(async name => ({
    name, source: await readFile(new URL(`../../examples/lean-derivations/${name}.lean`, import.meta.url), 'utf8') })))
];
const server = createServer(async (req, res) => {
  res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
  res.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
  const pathname = new URL(req.url, 'http://localhost').pathname;
  if (pathname === '/') return res.end('<!doctype html><title>Lean performance</title>');
  const file = resolve(root, '.' + pathname);
  if (!file.startsWith(root + sep)) return res.writeHead(404).end();
  try {
    res.setHeader('Content-Type', ({ '.js': 'text/javascript', '.json': 'application/json' })[extname(file)] || 'application/octet-stream');
    res.end(await readFile(file));
  } catch { res.writeHead(404).end(); }
});
await new Promise(r => server.listen(0, '127.0.0.1', r));
const browser = await chromium.launch({ executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined,
  args: ['--disable-dev-shm-usage'] });
const page = await browser.newPage();
const browserCdp = await browser.newBrowserCDPSession();
const cpuSnapshot = async () => (await browserCdp.send('SystemInfo.getProcessInfo')).processInfo;
const rendererCpu = (before, after) => {
  const pairs = after.filter(p => p.type === 'renderer').map(p => [before.find(old => old.id === p.id), p]);
  // Some sandboxes expose process IDs but return zero for every CPU counter.
  // A missing/replaced process is not a measured zero-cost compilation either.
  if (!pairs.length || pairs.some(([p]) => !p) ||
      pairs.every(([p, q]) => p.cpuTime === 0 && q.cpuTime === 0)) return null;
  return pairs.reduce((total, [p, q]) => total + Math.max(0, q.cpuTime - p.cpuTime), 0);
};
const measurements = { browser: browser.version(), deadlineMs, diagnosticDeadline: deadlineMs > 120000,
  cpuProfiled: Boolean(process.env.SSC_PERFORMANCE_CPU), cases: [] };
page.on('console', message => { if (message.text().startsWith('Performance:')) console.log(message.text()); });
try {
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  const started = performance.now();
  await page.evaluate(() => new Promise((resolve, reject) => {
    const worker = window.performanceWorker = new Worker('/lean-wasm-worker.js');
    const timer = setTimeout(() => reject(Error('Compiler initialization timed out')), 300000);
    worker.onerror = event => reject(Error(event.message));
    worker.onmessage = ({ data }) => {
      if (data.type === 'progress') console.log('Performance: initialization', data.text);
      if (data.type === 'ready') { clearTimeout(timer); resolve(); }
      if (data.type === 'error') { clearTimeout(timer); reject(Error(data.message)); }
    };
    worker.postMessage({ type: 'initialize' });
  }));
  measurements.initializationMs = performance.now() - started;
  console.log('Performance: initialization complete', Math.round(measurements.initializationMs), 'ms');
  let profileCommand;
  if (process.env.SSC_PERFORMANCE_CPU) {
    const cdp = browserCdp;
    const { targetInfos } = await cdp.send('Target.getTargets');
    const target = targetInfos.find(t => t.type === 'worker' && t.url.includes('lean-wasm-worker.js'));
    assert.ok(target, 'Compiler worker target is available for profiling');
    const { sessionId } = await cdp.send('Target.attachToTarget', { targetId: target.targetId, flatten: false });
    let serial = 0;
    const pending = new Map();
    cdp.on('Target.receivedMessageFromTarget', event => {
      if (event.sessionId !== sessionId) return;
      const response = JSON.parse(event.message), waiter = pending.get(response.id);
      if (!waiter) return;
      pending.delete(response.id);
      if (response.error) waiter.reject(Error(response.error.message)); else waiter.resolve(response.result);
    });
    profileCommand = (method, params = {}) => new Promise((resolve, reject) => {
      const id = ++serial; pending.set(id, { resolve, reject });
      cdp.send('Target.sendMessageToTarget', { sessionId, message: JSON.stringify({ id, method, params }) }).catch(reject);
    });
    await profileCommand('Profiler.enable');
    await profileCommand('Profiler.setSamplingInterval', { interval: 1000 });
  }
  const compile = async source => page.evaluate(({ source, deadlineMs }) => new Promise((resolve, reject) => {
    const worker = window.performanceWorker;
    const timer = setTimeout(() => { worker.terminate(); reject(Error(`Benchmark exceeded ${deadlineMs} ms`)); }, deadlineMs);
    worker.onmessage = ({ data }) => {
      if (data.type === 'progress') console.log('Performance:', data.text);
      if (data.type === 'result') { clearTimeout(timer); resolve(data.value); }
      if (data.type === 'error') { clearTimeout(timer); reject(Error(data.message)); }
    };
    worker.postMessage({ type: 'compile', id: 1, source });
  }), { source, deadlineMs });
  const fast = [];
  for (let i = 0; i < 31; i++) {
    const started = performance.now();
    const result = await compile('import SSC\ndef warm : Nat := 7');
    assert.equal(result.ok, true);
    if (!i) {
      measurements.leanVersion = result.leanVersion;
      measurements.packageVersion = result.packageVersion;
    }
    if (i) fast.push(performance.now() - started);
  }
  fast.sort((a, b) => a - b);
  measurements.warm = { runs: fast.length, medianMs: fast[15], p95Ms: fast[28], maxMs: fast[29] };
  console.log('Performance: warm requests', JSON.stringify(measurements.warm));
  for (const { name, source } of cases) {
    console.log('Performance: starting', name);
    const input = { name, sourceBytes: Buffer.byteLength(source),
      sourceSha256: createHash('sha256').update(source).digest('hex') };
    if (profileCommand) await profileCommand('Profiler.start');
    const cpuBefore = await cpuSnapshot();
    const started = performance.now();
    let result;
    try { result = await compile(source); }
    catch (error) {
      measurements.cases.push({ ...input, elapsedMs: performance.now() - started,
        rendererCpuSeconds: rendererCpu(cpuBefore, await cpuSnapshot()), error: error.message });
      throw error;
    }
    const elapsedMs = performance.now() - started;
    const rendererCpuSeconds = rendererCpu(cpuBefore, await cpuSnapshot());
    if (profileCommand) {
      const { profile } = await profileCommand('Profiler.stop');
      await writeFile(`${process.env.SSC_PERFORMANCE_CPU}.${measurements.cases.length}.cpuprofile`, JSON.stringify(profile));
    }
    assert.equal(result.ok, true);
    const finiteProfiles = result.diagnostics.filter(d => d.severity === 'information' && d.data.startsWith('finite_derive:')).map(d => d.data);
    measurements.cases.push({ ...input, elapsedMs, rendererCpuSeconds,
      packageVersion: result.packageVersion,
      declarations: result.results.length, phases: result.performance,
      calculationMessages: result.diagnostics.filter(d => d.severity === 'information' && d.data.startsWith('SSC_CALCULATION ')).map(d => d.data),
      ...(finiteProfiles.length ? { finiteProfiles } : {}) });
    console.log('Performance:', name, Math.round(measurements.cases.at(-1).elapsedMs), 'ms');
  }
} finally {
  if (process.env.SSC_PERFORMANCE_REPORT)
    await writeFile(process.env.SSC_PERFORMANCE_REPORT, JSON.stringify(measurements, null, 2) + '\n');
  console.log(JSON.stringify(measurements, null, 2));
  await browser.close();
  await new Promise(r => server.close(r));
}
