// Real static-host entry point: the service worker supplies isolation headers,
// then an offline page reload must start a new Lean runtime and check new source.
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { extname, resolve, sep } from 'node:path';
import { chromium } from 'playwright';

const root = resolve(process.argv[2] || new URL('..', import.meta.url).pathname);
const prefix = '/sample-size-calc/';
const requests = [], errors = [];
const mime = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript',
  '.wasm': 'application/wasm', '.json': 'application/json', '.css': 'text/css' };
const server = createServer(async (req, res) => {
  const path = new URL(req.url, 'http://localhost').pathname;
  requests.push({ path, method: req.method });
  if (path.includes('/api/')) return res.writeHead(503).end('No compiler backend');
  if (!path.startsWith(prefix)) return res.writeHead(404).end();
  const file = resolve(root, path.slice(prefix.length) || 'index.html');
  if (!file.startsWith(root + sep)) return res.writeHead(404).end();
  try {
    // Deliberately no COOP/COEP headers and no invented engine artifacts. Missing
    // manual-calculator files can make install-time cache.addAll fail; normal
    // service-worker fetch caching must still support the actual Lean page.
    res.setHeader('Content-Type', mime[extname(file)] || 'application/octet-stream');
    res.end(await readFile(file));
  } catch { res.writeHead(404).end(); }
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const browser = await chromium.launch({
  executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined,
  args: ['--disable-dev-shm-usage'],
});
const context = await browser.newContext();
const page = await context.newPage();
page.on('pageerror', error => errors.push(error.message));
const base = `http://127.0.0.1:${server.address().port}${prefix}`;
let progress = '';
const timer = setInterval(async () => {
  try {
    const current = await page.locator('#lean-connection-status').textContent();
    if (current !== progress) { console.log('Static startup:', current); progress = current; }
  } catch {}
}, 10000);
const ready = async () => {
  await page.waitForFunction(() => !document.getElementById('lean-compile').disabled
    || !document.getElementById('lean-setup').hidden, null, { timeout: 300000 });
  assert.equal(await page.locator('#lean-compile').isDisabled(), false,
    await page.locator('#lean-connection-status').textContent());
  assert.equal(await page.evaluate(() => crossOriginIsolated), true);
};
const checkSource = async (name, count, sum) => {
  await page.locator('#lean-source').fill(`import SSC
theorem ${name} (n : Nat) : n + 0 = n := by simp
#eval IO.println ("STATIC_PAGE_VALUE " ++ toString (List.range ${count}).sum)
`);
  await page.locator('#lean-compile').click();
  await page.waitForFunction(() => document.getElementById('lean-status').textContent
    .startsWith('Lean compilation complete.') || !document.getElementById('lean-error').hidden,
  null, { timeout: 120000 });
  assert.equal(await page.locator('#lean-error').isVisible(), false,
    await page.locator('#lean-error').textContent());
  assert.match(await page.locator('#lean-query').textContent(), new RegExp(name));
  assert.match(await page.locator('#lean-proof-status').textContent(), /Theorem checked by the Lean kernel/);
  assert.match(await page.locator('#lean-evaluation').textContent(), new RegExp(`STATIC_PAGE_VALUE ${sum}(?:\\s|$)`));
};
try {
  await page.goto(base + 'lean.html');
  await ready();
  assert.ok(await page.evaluate(() => Boolean(navigator.serviceWorker.controller)),
    'The static page must be controlled by the isolation service worker');
  await checkSource('online_startup_valid', 8, 28);
  const helperUsed = requests.some(r => r.path.endsWith('/lean-decompress-worker.js'));
  if (helperUsed) {
    assert.equal(await page.evaluate(async () => Boolean(await caches.match(
      new URL('lean-decompress-worker.js?v=dev', location.href)))), true,
    'The temporary decompression helper must be available offline');
  }

  // Stop an active request through the actual page control. The handler begins
  // reconnecting; immediate reload discards it and exercises a fresh page/worker.
  await page.locator('#lean-source').fill('import SSC\n#eval IO.sleep 60000');
  await page.locator('#lean-compile').click();
  await page.waitForTimeout(500);
  await context.setOffline(true);
  const beforeOffline = requests.length;
  await page.locator('#lean-stop').click();
  await page.waitForFunction(() => document.getElementById('lean-status').textContent === 'Stopped.');
  await page.reload();
  await ready();
  await checkSource('offline_startup_valid', 7, 21);
  assert.doesNotMatch(await page.locator('#lean-query').textContent(), /online_startup_valid/);
  // Chromium's browser-managed service-worker update checks can bypass the
  // page/context offline emulation. Only that exact script path is exempt:
  // the page, helper, decoder and compiler assets must all come from caches.
  assert.deepEqual(requests.slice(beforeOffline).filter(r =>
    r.path !== prefix + 'sw.js' || r.method !== 'GET'), [],
  'Offline application startup must not fetch page, compiler or decoder assets');
  assert.deepEqual(requests.filter(r => r.path.includes('/api/') || /\/study\.(?:mjs|wasm)$/.test(r.path)), []);
  assert.deepEqual(errors, []);
  console.log(`Static Lean page: isolation, checked source, UI cancellation, offline fresh-worker reload and new proof passed${helperUsed ? ' with cached decompression helper' : ' with native decompression'}.`);
} catch (error) {
  console.error('Static page status:', await page.locator('#lean-connection-status').textContent().catch(() => 'page unavailable'));
  console.error('Static page errors:', errors);
  throw error;
} finally {
  clearInterval(timer);
  await context.close();
  await browser.close();
  await new Promise(resolve => server.close(resolve));
}
