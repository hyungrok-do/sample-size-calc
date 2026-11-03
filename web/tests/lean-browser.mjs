// Real browser -> C++/WASM smoke and regression checks for the independent third tab.
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { resolve, extname } from 'node:path';
import { chromium } from 'playwright';
import { examples } from '../lean-examples.js';

const root = resolve(process.argv[2] || new URL('..', import.meta.url).pathname);
const types = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm' };
const server = createServer(async (req, res) => {
  const path = new URL(req.url, 'http://localhost').pathname;
  const file = resolve(root, `.${path === '/' ? '/index.html' : path}`);
  if (!file.startsWith(`${root}/`)) return res.writeHead(403).end();
  try {
    res.writeHead(200, { 'Content-Type': types[extname(file)] || 'application/octet-stream',
      'Cross-Origin-Opener-Policy': 'same-origin', 'Cross-Origin-Embedder-Policy': 'require-corp' });
    res.end(await readFile(file));
  } catch { res.writeHead(404).end(); }
});
await new Promise(r => server.listen(0, '127.0.0.1', r));
const browser = await chromium.launch({ executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined,
  args: ['--disable-dev-shm-usage', ...JSON.parse(process.env.SSC_CHROMIUM_ARGS || '[]')] });
const context = await browser.newContext({ viewport: { width: 1360, height: 1000 }, serviceWorkers: 'allow' });
const page = await context.newPage(), errors = [], violations = [];
page.on('pageerror', e => errors.push(e.message));
await page.exposeFunction('recordViolation', e => violations.push(e));
await page.addInitScript(() => document.addEventListener('securitypolicyviolation', e => window.recordViolation({ blocked: e.blockedURI, directive: e.violatedDirective })));
const base = `http://127.0.0.1:${server.address().port}`;
async function compile() {
  await page.locator('#lean-compile').click();
  await page.waitForFunction(() => !document.getElementById('lean-solve').disabled || !document.getElementById('lean-error').hidden, null, { timeout: 90000 });
  assert.equal(await page.locator('#lean-error').isVisible(), false, await page.locator('#lean-error').textContent());
}
async function solve() {
  await page.locator('#lean-solve').click();
  await page.waitForFunction(() => document.getElementById('lean-status').textContent === 'Calculation complete.');
}
try {
  await page.goto(`${base}/lean.html?engine=single`);
  assert.equal(await page.locator('.workspace-tabs a').count(), 3);
  assert.equal(await page.locator('#lean-solve').isDisabled(), true);
  await compile(); await solve();
  assert.equal(await page.locator('#lean-answer').textContent(), '63');
  assert.match(await page.locator('#lean-result-detail').textContent(), /126 participants total/);
  assert.match(await page.locator('#lean-cpp').textContent(), /int main\(\)/);
  if (process.env.SCREENSHOTS) await page.screenshot({ path: resolve(process.env.SCREENSHOTS, 'lean-desktop.png'), fullPage: true });
  await page.locator('#lean-target').selectOption('power');
  assert.equal(await page.locator('#lean-solve').isDisabled(), true);
  assert.equal(await page.locator('#lean-result').isVisible(), false);
  await compile(); await solve();
  assert.match(await page.locator('#lean-answer').textContent(), /^80\.1302%$/);
  await page.locator('#lean-target').selectOption('effect'); await compile(); await solve();
  assert.ok(Number(await page.locator('#lean-answer').textContent()) < .5);
  await page.locator('#lean-source').fill(examples.normal.replace('sd := 1', 'sd := 0'));
  await page.locator('#lean-compile').click();
  await page.waitForFunction(() => !document.getElementById('lean-error').hidden);
  assert.equal(await page.locator('#lean-solve').isDisabled(), true);
  assert.match(await page.locator('#lean-error').textContent(), /sd must be strictly positive/);
  await page.locator('#lean-source').fill(examples.normal + '\n#eval 1 + 1');
  await page.locator('#lean-compile').click(); await page.waitForFunction(() => !document.getElementById('lean-error').hidden);
  assert.match(await page.locator('#lean-error').textContent(), /unsupported token|unexpected code/);
  await page.locator('#lean-example').selectOption('influence'); await page.locator('#lean-load').click();
  await page.locator('#lean-target').selectOption('sample_size'); await compile(); await solve();
  assert.equal(await page.locator('#lean-answer').textContent(), '126');
  assert.match(await page.locator('#lean-basis').textContent(), /asymptotic/);
  await page.reload(); assert.equal(await page.locator('#lean-source').inputValue(), examples.influence);
  assert.equal(await page.locator('#lean-solve').isDisabled(), true);
  await page.setViewportSize({ width: 390, height: 844 }); await compile(); await solve();
  assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth), false);
  if (process.env.SCREENSHOTS) await page.screenshot({ path: resolve(process.env.SCREENSHOTS, 'lean-mobile.png'), fullPage: true });
  // The new tab's assets must remain navigable and compilable offline.
  await page.evaluate(async () => { await navigator.serviceWorker.ready; });
  await page.waitForFunction(async () => (await caches.keys()).length > 0);
  await page.waitForFunction(async () => !!(await caches.match('lean-examples.js?v=dev', { ignoreSearch: true })), null, { timeout: 60000 });
  await context.setOffline(true); await page.reload(); await compile(); await solve();
  assert.equal(await page.locator('#lean-answer').textContent(), '126');
  assert.deepEqual(errors, []); assert.deepEqual(violations, []);
  console.log('Lean tab: compile, three solve targets, invalidation, errors, storage, responsive layout and offline execution passed.');
} finally { await context.close(); await browser.close(); server.close(); }
