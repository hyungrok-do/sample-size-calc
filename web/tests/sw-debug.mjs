// Diagnostic only (debug branch): why WebKit keeps no offline copies.  Serves the assembled site as site.mjs does,
// with the service worker instrumented to report each caching step to the page's console, and checks the Cache
// API from the page itself, in a new context (memory) and in a persistent one (disk).
// node tests/sw-debug.mjs <site dir>; SSC_BROWSER picks the browser.
import { createServer } from 'node:http';
import { mkdtemp, readFile, readdir, rm } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { extname, join, resolve } from 'node:path';
import { chromium, firefox, webkit } from 'playwright';

const site = resolve(process.argv[2] ?? '.');
const files = new Map();
async function walk(dir, rel = '') {
  for (const e of await readdir(join(dir, rel), { withFileTypes: true })) {
    const r = rel ? `${rel}/${e.name}` : e.name;
    if (e.isDirectory()) await walk(dir, r); else files.set(r, await readFile(join(dir, r)));
  }
}
await walk(site);

// the service worker, instrumented
let sw = files.get('sw.js').toString('utf8');
const report = `
const __log = [];
function report(msg) {
  __log.push(msg);
  self.clients.matchAll({ includeUncontrolled: true, type: 'all' })
    .then((cs) => cs.forEach((c) => c.postMessage('[sw] ' + msg))).catch(() => {});
}
self.addEventListener('message', (event) => {
  (async () => {
    const view = {};
    try { for (const name of await caches.keys()) view[name] = (await (await caches.open(name)).keys()).length; }
    catch (e) { view.error = String(e); }
    event.source?.postMessage('[sw] state ' + JSON.stringify({ caches: view, log: __log.slice(-40) }));
  })();
});
self.addEventListener('error', (e) => report('error ' + e.message));
self.addEventListener('unhandledrejection', (e) => report('unhandledrejection ' + e.reason));
report('script evaluated');
`;
function swap(old, now) {
  if (!sw.includes(old)) throw new Error(`sw.js no longer has: ${old}`);
  sw = sw.replace(old, now);
}
swap("event.waitUntil(caches.open(CACHE).then((cache) => cache.addAll(fresh)).catch(() => {}));",
  "report('install'); event.waitUntil(caches.open(CACHE).then((cache) => cache.addAll(fresh))" +
  ".then(() => report('install addAll ok'), (e) => report('install addAll FAILED ' + e)));");
swap("await self.clients.claim();", "await self.clients.claim(); report('activated, caches ' + JSON.stringify(await caches.keys()));");
swap("event.waitUntil(caches.open(CACHE).then((cache) => cache.put(key, copy)));",
  "const put = caches.open(CACHE).then((cache) => cache.put(key, copy)).then(() => report('put ok ' + key.split('/').pop())," +
  " (e) => report('put FAILED ' + key.split('/').pop() + ' ' + e));\n" +
  "        try { event.waitUntil(put); } catch (e) { report('waitUntil threw ' + e); }");
swap("} catch (error) {", "} catch (error) { report('fetch failed ' + request.url.split('/').pop() + ' ' + error);");
files.set('sw.js', Buffer.from(report + sw));

const types = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.wasm': 'application/wasm' };
const prefix = '/sample-size-calc/';
const server = createServer((req, res) => {
  const pathname = new URL(req.url, 'http://localhost').pathname;
  const file = pathname.slice(prefix.length) || 'index.html';
  const body = files.get(file);
  if (!pathname.startsWith(prefix) || !body) { res.writeHead(404).end(); return; }
  res.writeHead(200, { 'Content-Type': types[extname(file)] ?? 'application/octet-stream', 'Cache-Control': 'max-age=600' })
    .end(body);
});
await new Promise((done) => server.listen(0, '127.0.0.1', done));
const url = `http://127.0.0.1:${server.address().port}${prefix}`;
const kind = process.env.SSC_BROWSER ?? 'chromium';
const type = { chromium, firefox, webkit }[kind];
const launch = kind === 'chromium' ? { executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined } : {};

async function session(label, context) {
  console.log(`\n===== ${kind}, ${label}`);
  const page = await context.newPage();
  page.on('console', (m) => console.log(`  console.${m.type()}: ${m.text()}`));
  page.on('pageerror', (e) => console.log(`  pageerror: ${e.message}`));
  await page.addInitScript(() => navigator.serviceWorker?.addEventListener('message', (e) => console.log(String(e.data))));
  const ready = () => page.waitForFunction(() => /ready in|did not load/.test(document.getElementById('engine-status')?.textContent),
    null, { timeout: 60000 }).catch((e) => console.log(`  not ready: ${e.message.split('\n')[0]}`));
  const probe = async (when) => {
    const view = await page.evaluate(async () => {
      const out = { controlled: !!navigator.serviceWorker?.controller, status: document.getElementById('engine-status')?.textContent };
      try {
        out.caches = {};
        for (const name of await caches.keys()) out.caches[name] = (await (await caches.open(name)).keys()).length;
      } catch (e) { out.caches = String(e); }
      try {  // the Cache API from the page itself
        const c = await caches.open('probe');
        await c.put('probe-entry', new Response('x'));
        out.pageCache = (await (await c.match('probe-entry'))?.text()) === 'x' ? 'works' : 'no match';
        await caches.delete('probe');
      } catch (e) { out.pageCache = String(e); }
      try { out.estimate = await navigator.storage?.estimate?.(); } catch (e) { out.estimate = String(e); }
      try { out.persisted = await navigator.storage?.persisted?.(); } catch (e) { out.persisted = String(e); }
      navigator.serviceWorker?.controller?.postMessage('state');
      return out;
    }).catch((e) => ({ evaluate: e.message }));
    await page.waitForTimeout(500);
    console.log(`  ${when}: ${JSON.stringify(view)}`);
  };
  await page.goto(url);
  await ready();
  await page.waitForFunction(() => navigator.serviceWorker?.controller, null, { timeout: 30000 })
    .catch(() => console.log('  never controlled'));
  await page.waitForTimeout(1500);
  await probe('after the first visit');
  await page.goto(url);
  await ready();
  await page.waitForTimeout(1500);
  await probe('after the second visit');
  await page.goto(url);
  await ready();
  await page.waitForTimeout(1500);
  await probe('after the third visit');
  await page.close();
}

try {
  const browser = await type.launch(launch);
  await session('new context (memory)', await browser.newContext({ serviceWorkers: 'allow' }));
  await browser.close();
  const profile = await mkdtemp(join(tmpdir(), 'ssc-swdebug-'));
  const persistent = await type.launchPersistentContext(profile, { ...launch, serviceWorkers: 'allow' });
  await session('persistent context (disk)', persistent);
  await persistent.close();
  await rm(profile, { recursive: true, force: true });
} finally {
  server.close();
}
