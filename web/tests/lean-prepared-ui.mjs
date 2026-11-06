// UI lifecycle checks use a deterministic compiler double. Actual Lean execution,
// replay/audit and numerical equivalence are covered by lean-prepared.mjs.
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { resolve, extname, sep } from 'node:path';
import { chromium } from 'playwright';

const root = resolve(new URL('..', import.meta.url).pathname);
const server = createServer(async (req, res) => {
  const pathname = new URL(req.url, 'http://localhost').pathname;
  const file = resolve(root, '.' + pathname);
  if (!file.startsWith(root + sep)) return res.writeHead(404).end();
  try {
    res.setHeader('Content-Type', ({ '.html': 'text/html', '.js': 'text/javascript', '.css': 'text/css' })[extname(file)] || 'application/octet-stream');
    res.end(await readFile(file));
  } catch { res.writeHead(404).end(); }
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const browser = await chromium.launch({
  executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined,
  args: ['--disable-dev-shm-usage'],
});
const page = await browser.newPage({ viewport: { width: 1360, height: 1000 }, serviceWorkers: 'block' });
const errors = [];
page.on('pageerror', error => errors.push(error.message));
await page.route('**/lean-wasm-client.js*', route => route.fulfill({
  contentType: 'text/javascript', body: `
const pause = () => new Promise(resolve => setTimeout(resolve, 30));
const report = source => ({results: [], diagnostics: [], output: '', module: source,
  performance: {totalMs: 250, commandMs: 200, kernelReplayMs: 30, axiomAuditMs: 10,
    declarationRenderMs: 4, traceRenderMs: 2, headerMs: 2, parseMs: 2, commands: []}});
export class BrowserLeanCompiler {
  async initialize() { this.ready = true; }
  async compile(source) { await pause(); return report(source); }
  async prepare(source, entry) {
    await pause();
    if (!entry || entry === 'missing') throw Error('Unknown entry.');
    return {...report(source), prepared: {id: 'prepared', entry}};
  }
  async evaluate(id, parameters) {
    await pause();
    if (parameters.effectMilli === 0) throw Error('effectMilli must be positive.');
    return {ok: true, verification: 'compiled-execution', concreteResultKernelChecked: false,
      value: {selected: parameters.effectMilli === 500 ? 63 : 44, method: 'numerical'},
      performance: {evaluationMs: 1.25, totalMs: 1.5}};
  }
  stop() { this.ready = false; }
}`,
}));
const visibleError = () => page.locator('#lean-error').isVisible();
try {
  await page.goto(`http://127.0.0.1:${server.address().port}/lean.html?example=prepared-normal-search`);
  await page.waitForFunction(() => !document.getElementById('lean-prepare').disabled);
  assert.equal(await page.locator('#lean-repeated').getAttribute('open'), '');
  assert.match(await page.locator('#lean-source').inputValue(), /def sscEvaluate/);
  assert.equal(await page.locator('#lean-evaluate-parameters').isDisabled(), true);
  await page.locator('#lean-prepare').click();
  await page.waitForFunction(() => !document.getElementById('lean-evaluate-parameters').disabled);
  assert.equal(await visibleError(), false);
  assert.equal(await page.locator('#lean-compiled').isVisible(), true);
  assert.match(await page.locator('#lean-prepared-status').textContent(), /Prepared sscEvaluate/);
  await page.locator('#lean-evaluate-parameters').click();
  await page.locator('#lean-prepared-result').waitFor({ state: 'visible' });
  assert.match(await page.locator('#lean-prepared-result').textContent(), /Runtime evaluation.*no new numerical theorem/);
  assert.match(await page.locator('#lean-prepared-timing').textContent(), /1\.25 ms/);
  assert.equal(JSON.parse(await page.locator('#lean-prepared-value').textContent()).selected, 63);
  if (process.env.SSC_PREPARED_UI_SCREENSHOT)
    await page.screenshot({ path: process.env.SSC_PREPARED_UI_SCREENSHOT, fullPage: true });
  await page.setViewportSize({ width: 390, height: 844 });
  assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth), false);
  await page.locator('#lean-inputs').fill('{"effectMilli":600}');
  assert.equal(await page.locator('#lean-prepared-result').isVisible(), false);
  assert.equal(await page.locator('#lean-evaluate-parameters').isDisabled(), false);
  await page.locator('#lean-evaluate-parameters').click();
  await page.locator('#lean-prepared-result').waitFor({ state: 'visible' });
  assert.equal(JSON.parse(await page.locator('#lean-prepared-value').textContent()).selected, 44);
  await page.locator('#lean-inputs').fill('{broken JSON');
  await page.locator('#lean-evaluate-parameters').click();
  await page.locator('#lean-error').waitFor({ state: 'visible' });
  assert.equal(await page.locator('#lean-prepared-result').isVisible(), false);
  assert.equal(await page.locator('#lean-evaluate-parameters').isDisabled(), false);
  await page.locator('#lean-inputs').fill('{"effectMilli":0}');
  await page.locator('#lean-evaluate-parameters').click();
  await page.waitForFunction(() => document.getElementById('lean-error').textContent.includes('positive'));
  await page.locator('#lean-entry').fill('missing');
  assert.equal(await page.locator('#lean-evaluate-parameters').isDisabled(), true);
  await page.locator('#lean-prepare').click();
  await page.waitForFunction(() => document.getElementById('lean-error').textContent.includes('Unknown'));
  assert.equal(await page.locator('#lean-evaluate-parameters').isDisabled(), true);
  await page.locator('#lean-entry').fill('sscEvaluate');
  await page.locator('#lean-prepare').click();
  await page.waitForFunction(() => !document.getElementById('lean-evaluate-parameters').disabled);
  await page.locator('#lean-source').fill('import SSC\ndef replacement := 1');
  assert.equal(await page.locator('#lean-evaluate-parameters').isDisabled(), true);
  assert.equal(await page.locator('#lean-compiled').isVisible(), false);
  await page.locator('#lean-prepare').click();
  await page.waitForFunction(() => !document.getElementById('lean-evaluate-parameters').disabled);
  await page.locator('#lean-compile').click();
  await page.waitForFunction(() => !document.getElementById('lean-compile').disabled);
  assert.equal(await page.locator('#lean-evaluate-parameters').isDisabled(), true);
  assert.deepEqual(errors, []);
  console.log('Prepared UI: parameter edits, separate runtime status, errors, source/entry invalidation, ordinary compile and mobile layout passed.');
} finally {
  await browser.close();
  await new Promise(resolve => server.close(resolve));
}
