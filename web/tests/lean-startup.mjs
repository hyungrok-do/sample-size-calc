// Real browser startup phases, content-cache restarts, and offline recovery.
// Pass a second root for a paired comparison: warm both packages, then alternate
// cached workers in the same browser/page. Warmups are excluded from summaries.
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { createServer } from 'node:http';
import { readFile, writeFile } from 'node:fs/promises';
import { extname, resolve, sep } from 'node:path';
import { chromium } from 'playwright';

const roots = [resolve(process.argv[2] || new URL('..', import.meta.url).pathname),
  ...(process.argv[3] ? [resolve(process.argv[3])] : [])];
const paired = roots.length === 2;
const variants = roots.map((root, i) => ({ root, name: paired ? ['before', 'after'][i] : 'single',
  prefix: paired ? `/${['before', 'after'][i]}/` : '/' }));
const requests = [];
const server = createServer(async (req, res) => {
  res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
  res.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
  const url = new URL(req.url, 'http://localhost');
  if (url.pathname === '/') return res.end('<!doctype html><title>Lean startup benchmark</title>');
  const variant = variants.find(v => url.pathname.startsWith(v.prefix));
  if (!variant) return res.writeHead(404).end();
  const file = resolve(variant.root, url.pathname.slice(variant.prefix.length));
  if (!file.startsWith(variant.root + sep)) return res.writeHead(404).end();
  try {
    const bytes = await readFile(file);
    requests.push({ variant: variant.name, path: url.pathname, bytes: bytes.length });
    res.setHeader('Content-Type', ({ '.js': 'text/javascript', '.json': 'application/json' })[extname(file)] || 'application/octet-stream');
    // Versioned worker/glue scripts can restart offline as on a static host.
    if (url.searchParams.has('v')) res.setHeader('Cache-Control', 'public, max-age=31536000, immutable');
    res.end(bytes);
  } catch { res.writeHead(404).end(); }
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const browser = await chromium.launch({
  executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined,
  args: ['--disable-dev-shm-usage'],
});
const context = await browser.newContext();
const page = await context.newPage();
page.on('console', message => {
  if (message.text().startsWith('Startup:')) console.log(message.text());
});
const identity = async ({ root, name }) => {
  const bytes = await readFile(resolve(root, 'lean-wasm/manifest.json'));
  const manifest = JSON.parse(bytes);
  return { name, leanVersion: manifest.leanVersion, packageVersion: manifest.version, hasDecoder: Boolean(manifest.decoder),
    manifestSha256: createHash('sha256').update(bytes).digest('hex'),
    workerSha256: createHash('sha256').update(await readFile(resolve(root, 'lean-wasm-worker.js'))).digest('hex') };
};
const identities = await Promise.all(variants.map(identity));
const report = {
  browser: browser.version(), ...(paired ? { variants: identities } : identities[0]),
  measurement: paired
    ? 'Local HTTP server, one browser/page. Each package receives one unmeasured warmup, then three alternating cached restarts and an offline restart. Browser-level WASM compilation caches may be reused. Warmup results are retained but excluded from summaries.'
    : 'Local HTTP server. First empty-cache worker followed by restarted workers in the same browser context; browser-level WASM compilation caches may be reused.',
  settings: { paired, cachedRuns: 3, versionedHttpCache: true, kernelProofAfterEachStartup: true,
    startupTimeoutMs: 300000, sameBrowserContext: true },
  runs: [],
};
// This ordinary dynamic evaluation needs an imported IR body that is not in
// the pinned runtime's native export table. Keeping it catches unsafe package
// trimming even when elementary proofs still happen to pass.
const startupSource = `import SSC
open Lean
unsafe def dynamicCore : CoreM Unit := do
  let env ← getEnv
  let opts ← getOptions
  match env.evalConst (Array Syntax → Syntax) opts \`\`Lean.mkNullNode false with
  | .ok f => logInfo m!"DYNAMIC_CORE_IR {(f #[]).getArgs.size}"
  | .error e => throwError "{e}"
#eval dynamicCore
theorem startup_valid (n : Nat) : n + 0 = n := by simp
`;
const stopAll = () => page.evaluate(() => Object.values(window.compilers || {}).forEach(c => c.stop()));
const runStartup = async (variant, mode, measured = true) => {
  await stopAll();
  const requestStart = requests.length;
  const run = await page.evaluate(async ({ name, mode, measured, source }) => {
    const compiler = window.compilers[name];
    const started = performance.now();
    const ready = await compiler.initialize();
    const roundTripMs = performance.now() - started;
    const proof = await compiler.compile(source);
    return { variant: name, mode, measured, roundTripMs, initialization: ready.initialization,
      proof: proof.results.find(d => d.name === 'startup_valid'), compileMs: proof.performance.totalMs,
      dynamicOutput: [proof.output, ...proof.diagnostics.map(d => d.data)].join('\n') };
  }, { name: variant.name, mode, measured, source: startupSource });
  assert.equal(run.proof?.kind, 'theorem', `${variant.name}/${mode}: startup must leave an operational checked compiler`);
  assert.match(run.dynamicOutput, /DYNAMIC_CORE_IR 0/);
  delete run.dynamicOutput;
  assert.ok(!run.proof.axioms.some(a => /sorry|trustCompiler|ofReduce/.test(a)));
  assert.ok(run.initialization && run.initialization.totalMs > 0);
  for (const [name, value] of Object.entries(run.initialization)) {
    if (name.endsWith('Ms')) assert.ok(Number.isFinite(value) && value >= 0, name);
  }
  const newRequests = requests.slice(requestStart);
  run.network = { requests: newRequests.length,
    responseBytes: newRequests.reduce((n, r) => n + r.bytes, 0),
    packRequests: newRequests.filter(r => /\.pack$/.test(r.path)).length };
  if (mode.startsWith('cached-') || mode === 'offline-cache')
    assert.equal(run.network.packRequests, 0, 'Verified content cache must avoid downloading packs again');
  if (mode === 'offline-cache') assert.equal(run.network.requests, 0);
  delete run.proof;
  report.runs.push(run);
  console.log('Startup result:', JSON.stringify(run));
};
try {
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  await page.evaluate(async variants => {
    window.compilers = {};
    for (const { name, prefix } of variants) {
      const { BrowserLeanCompiler } = await import(prefix + 'lean-wasm-client.js');
      window.compilers[name] = new BrowserLeanCompiler(text => console.log('Startup:', name, text));
    }
  }, variants.map(({ name, prefix }) => ({ name, prefix })));
  report.features = await page.evaluate(() => ({
    crossOriginIsolated,
    brotliDecompression: (() => { try { new DecompressionStream('brotli'); return true; } catch { return false; } })(),
  }));
  if (paired) {
    for (const variant of variants) await runStartup(variant, 'warmup', false);
    for (let round = 1; round <= 3; round++) {
      const order = round % 2 ? variants : [...variants].reverse();
      for (const variant of order) await runStartup(variant, `cached-${round}`);
    }
  } else {
    for (const mode of ['empty-cache', 'cached-1', 'cached-2', 'cached-3'])
      await runStartup(variants[0], mode);
  }
  await context.setOffline(true);
  for (const variant of variants) await runStartup(variant, 'offline-cache');
  report.cached = Object.fromEntries(variants.map(variant => {
    const runs = report.runs.filter(r => r.variant === variant.name && r.mode.startsWith('cached-'));
    const median = values => [...values].sort((a, b) => a - b)[1];
    const phases = Object.keys(runs[0].initialization).filter(name => name.endsWith('Ms'));
    return [variant.name, { medianMs: median(runs.map(r => r.roundTripMs)),
      phaseMedians: Object.fromEntries(phases.map(name => [name, median(runs.map(r => r.initialization[name]))])) }];
  }));
  if (!paired) report.cachedMedianMs = report.cached.single.medianMs;
  // Corruption is tested after all timings, so the benchmark cannot accidentally
  // measure a repaired cache. The worker must reject before running bad bytes.
  await context.setOffline(false);
  await stopAll();
  const checked = variants.at(-1);
  await page.evaluate(async prefix => {
    const cache = await caches.open('ssc-lean-content-v1');
    const key = (await cache.keys()).find(request => {
      const path = new URL(request.url).pathname;
      return path.startsWith(prefix + 'lean-wasm/') && path.includes('lean.wasm') && path.endsWith('.pack');
    });
    if (!key) throw Error('No cached runtime asset to validate');
    await cache.put(key, new Response('intentionally corrupted cache entry'));
  }, checked.prefix);
  await assert.rejects(page.evaluate(name => window.compilers[name].initialize(), checked.name), /damaged|checksum/i);
  report.corruptCacheRejected = checked.name;
  // A stopped compiler must also terminate its temporary decoder worker.
  // Run this only where the selected package/browser needs that helper.
  if (identities.at(-1).hasDecoder && !report.features.brotliDecompression) {
    await stopAll();
    const cdp = await browser.newBrowserCDPSession();
    await page.evaluate(name => {
      window.cancelledInitialization = window.compilers[name].initialize().then(
        () => ({ ready: true }), error => ({ ready: false, message: error.message }));
    }, checked.name);
    let helper;
    const deadline = Date.now() + 15000;
    while (Date.now() < deadline) {
      const { targetInfos } = await cdp.send('Target.getTargets');
      helper = targetInfos.find(t => t.type === 'worker' && t.url.includes(checked.prefix + 'lean-decompress-worker.js'));
      if (helper) break;
      await new Promise(resolve => setTimeout(resolve, 50));
    }
    assert.ok(helper, 'Temporary decoder worker must be active before cancellation');
    await stopAll();
    const cancelled = await page.evaluate(() => window.cancelledInitialization);
    assert.equal(cancelled.ready, false);
    const closedDeadline = Date.now() + 5000;
    let remains = true;
    while (Date.now() < closedDeadline) {
      const { targetInfos } = await cdp.send('Target.getTargets');
      remains = targetInfos.some(t => t.targetId === helper.targetId);
      if (!remains) break;
      await new Promise(resolve => setTimeout(resolve, 50));
    }
    assert.equal(remains, false, 'Stopping the compiler must terminate the decoder child');
    await cdp.detach();
    report.initializationCancellation = { variant: checked.name, childTerminated: true };
    await runStartup(checked, 'cancellation-recovery', false);
  }
  // Emulate an upgrade from the original gzip-only package, reusing exact
  // authenticated bytes already downloaded for the comparison baseline.
  if (paired && identities.at(-1).hasDecoder) {
    await stopAll();
    const copied = await page.evaluate(async ({ before, after }) => {
      const cache = await caches.open('ssc-lean-content-v1');
      const keys = await cache.keys();
      for (const key of keys) {
        const path = new URL(key.url).pathname;
        if (path.startsWith(after + 'lean-wasm/') && path.endsWith('.br.pack')) await cache.delete(key);
      }
      let count = 0;
      for (const key of keys) {
        const url = new URL(key.url);
        if (!url.pathname.startsWith(before + 'lean-wasm/') || !url.pathname.endsWith('.pack') || url.pathname.endsWith('.br.pack')) continue;
        const response = await cache.match(key);
        url.pathname = after + url.pathname.slice(before.length);
        await cache.put(url.href, response);
        count++;
      }
      return count;
    }, { before: variants[0].prefix, after: checked.prefix });
    assert.ok(copied > 0);
    await runStartup(checked, 'gzip-cache-migration', false);
    const migration = report.runs.at(-1);
    assert.equal(migration.initialization.compression, 'gzip');
    assert.equal(migration.initialization.decoder, 'native');
    assert.equal(migration.network.packRequests, 0, 'Existing gzip assets must not trigger replacement downloads');
    report.gzipCacheMigration = { copiedAssets: copied, noPackDownloads: true };
  }
  console.log('Lean startup, cache/offline integrity, imported IR, initialization cancellation and gzip migration checks passed.');
} finally {
  if (process.env.SSC_STARTUP_REPORT)
    await writeFile(process.env.SSC_STARTUP_REPORT, JSON.stringify(report, null, 2) + '\n');
  await stopAll().catch(() => {});
  await context.close();
  await browser.close();
  await new Promise(resolve => server.close(resolve));
}
