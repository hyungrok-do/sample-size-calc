// Exercises the real page and WebAssembly engine; no mocked schemas or calculator.
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { readFile, mkdir } from 'node:fs/promises';
import { dirname, resolve, extname, sep } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright';
import { guideData } from '../guide-data.js';

const web = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const types = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.wasm': 'application/wasm' };
const server = createServer(async (req, res) => {
  const path = resolve(web, '.' + (new URL(req.url, 'http://localhost').pathname === '/' ? '/index.html' : new URL(req.url, 'http://localhost').pathname));
  if (!path.startsWith(web + sep)) { res.writeHead(403).end(); return; }
  try {
    const body = await readFile(path);
    res.writeHead(200, { 'Content-Type': types[extname(path)] ?? 'application/octet-stream', 'Cache-Control': 'no-store' });
    res.end(body);
  } catch { res.writeHead(404).end(); }
});
await new Promise((resolve) => server.listen(0, '127.0.0.1', resolve));
const base = `http://127.0.0.1:${server.address().port}/?engine=single`;
const browser = await chromium.launch({
  executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined,
  args: ['--disable-dev-shm-usage'],
});
try {
  const context = await browser.newContext({ serviceWorkers: 'block', viewport: { width: 1280, height: 900 } });
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (error) => errors.push(error.message));
  // Browsing all methods comes first and opens a method at once, without scrolling past the introduction.
  await page.goto(base);
  await page.waitForSelector('#in-test');
  assert.equal(await page.locator('#mode-browse').getAttribute('aria-pressed'), 'true');
  assert.equal(await page.evaluate(() => document.querySelector('.entry-modes button').id), 'mode-browse');
  assert.equal(await page.locator('#guide').isVisible(), false);
  assert.equal(new URL(page.url()).hash, '#means.two.t');
  assert.equal(await page.evaluate(() => scrollY), 0);
  async function openGuide() {
    await page.locator('#mode-guide').click();
    await page.waitForSelector('#guide-questions input');
  }
  await openGuide();
  assert.equal(await page.locator('#method').isVisible(), false);
  assert.match(await page.locator('#guide .guide-experimental').textContent(), /This is an experimental feature\./);
  assert.equal(await page.locator('#guide-recommendation-title').textContent(), 'Recommended analysis');

  const routes = [];
  function visit(node, path = []) {
    if (node.type === 'question') node.choices.forEach((c, i) => visit(c.next, [...path, [node.node_id, String(i)]]));
    else routes.push({ node, path });
  }
  visit(guideData.root);
  async function select(route) {
    await page.locator('#guide-reset').click();
    await page.evaluate((path) => {
      for (const [id, value] of path) {
        const input = document.querySelector(`#guide input[name="${id}"][value="${value}"]`);
        input.checked = true;
        input.dispatchEvent(new Event('change', { bubbles: true }));
      }
    }, route.path);  // no additional requirement needs confirming
  }
  async function apply() {
    await page.locator('#guide-recommendation-content button').click();
    await page.waitForFunction(() => document.querySelector('#guide-recommendation-content button')?.textContent === 'Reapply this analysis');
    assert.equal(await page.locator('#method').isVisible(), true);
  }
  for (const route of routes) {
    await select(route);
    if (route.node.type !== 'candidate') {
      assert.equal(await page.locator('#guide-recommendation-content button').count(), 0);
      continue;
    }
    await apply();
    assert.equal(new URL(page.url()).hash.split('?')[0], `#${route.node.procedure_id}`);
    for (const [key, value] of Object.entries(route.node.proposed_option_values)) {
      assert.equal(await page.locator(`#in-${key}`).inputValue(), value, `${route.node.node_id}: ${key}`);
      assert.equal(await page.locator(`#in-${key}`).isVisible(), true);
    }
  }
  console.log(`Browser: ${routes.filter((r) => r.node.type === 'candidate').length} routes applied; unsupported routes withheld.`);

  // One answer already lists the analyses that fit; a requirement narrows them; any of them opens directly.
  await page.locator('#guide-reset').click();
  await page.locator(`#guide input[name="${guideData.root.node_id}"][value="0"]`).check();
  const listed = () => page.locator('#guide-recommendation-content .guide-candidate').evaluateAll((bs) => bs.map((b) => b.dataset.procedure));
  const afterOne = await listed();
  assert.ok(afterOne.length > 1, 'analyses listed after one answer');
  assert.equal(await page.locator('#guide-requirements').isVisible(), true);
  await page.locator('input[name="guide-extra"][value="cluster"]').check();
  const cluster = await listed();
  const allowed = guideData.features.find((f) => f.id === 'cluster').allow;
  assert.ok(cluster.length > 0 && cluster.length < afterOne.length && cluster.every((id) => allowed.includes(id)));
  await page.locator('#guide-recommendation-content .guide-candidate').first().click();
  await page.waitForFunction((id) => location.hash.startsWith(`#${id}`), cluster[0]);
  assert.equal(await page.locator('#method').isVisible(), true);
  assert.match(await page.locator('#guide-recommendation-content [role="status"]').textContent(), /open below/);

  const fisher = routes.find((r) => r.node.procedure_id === 'prop.two' && r.node.proposed_option_values.test === 'fisher');
  const score = routes.find((r) => r.node.procedure_id === 'prop.two' && r.node.proposed_option_values.test === 'zpool');
  await select(fisher); await apply();
  await page.locator('#in-p1').fill('0.7');
  await page.locator('#calc').click();
  await page.waitForSelector('#result .statement');
  assert.match(await page.locator('#result').textContent(), /Fisher/);
  await page.locator('#in-test').selectOption('zunpool');
  assert.equal(await page.locator('#result').textContent(), '');
  assert.match(await page.locator('#m-source').textContent(), /settings changed/);
  await apply();
  assert.equal(await page.locator('#in-p1').inputValue(), '0.7');
  assert.equal(await page.locator('#in-test').inputValue(), 'fisher');
  await select(score); await apply();
  assert.equal(await page.locator('#in-test').inputValue(), 'zpool', 'same-procedure change applies');
  await page.locator('input[name="guide-extra"][value="cluster"]').check();
  await page.locator('input[name="guide-extra"][value="interim"]').check();
  assert.equal(await page.locator('#method').isVisible(), false);
  assert.equal(await page.locator('#guide-recommendation-content button').count(), 0);
  assert.equal(await page.locator('input[name="guide-extra"]:checked').count(), 2);
  await page.locator('#guide-reset').click();
  const first = page.locator('#guide-questions input').first();
  await first.focus();
  await page.keyboard.press('Space');
  assert.equal(await first.isChecked(), true);
  await page.keyboard.press('ArrowRight');
  assert.equal(await page.locator('#guide-questions input').nth(1).isChecked(), true);
  await page.locator('#guide-questions input[value="unsure"]').first().check();
  assert.equal(await page.locator('#guide-recommendation-content button').count(), 0);

  // A change in design changes the meaning of CV, so it must not inherit the previous numeric input.
  const crossover = routes.find((r) => r.node.procedure_id === 'means.be.ratio' && r.node.proposed_option_values.design === '2x2');
  const parallel = routes.find((r) => r.node.procedure_id === 'means.be.ratio' && r.node.proposed_option_values.design === 'parallel');
  await select(crossover); await apply();
  await page.locator('#in-cv').fill('0.77');
  await select(parallel); await apply();
  assert.notEqual(await page.locator('#in-cv').inputValue(), '0.77');

  // A changed answer wins over a recommendation still opening asynchronously.
  await select(fisher);
  await page.evaluate(() => {
    document.querySelector('#guide-recommendation-content button').click();
    document.querySelector('#guide-reset').click();
  });
  assert.equal(await page.locator('#method').isVisible(), false);

  await page.goto(base + '#prop.two?test=missing');
  await page.waitForSelector('#navigation-error');
  assert.equal(await page.locator('#method').isVisible(), false);

  // Deep links retain their options, and legacy method-only links still open the calculator directly.
  await page.goto(base + '#prop.two?test=fisher&method=exact');
  await page.waitForSelector('#in-test');
  assert.equal(await page.locator('#in-test').inputValue(), 'fisher');
  await page.locator('#search').fill('Welch');
  await page.locator('#method-list button[data-id="means.two.welch"]').click();
  await page.waitForFunction(() => location.hash === '#means.two.welch');
  await page.goBack();
  await page.waitForFunction(() => document.querySelector('#in-test')?.value === 'fisher');
  await page.goto(base + '#means.two.t');
  await page.waitForSelector('#in-test');
  assert.equal(await page.locator('#in-test').inputValue(), 't');

  // A change made while a simulation is in flight must not receive its old result.
  await page.goto(base + '#sim.means.robust');
  await page.waitForSelector('#in-reps');
  await page.locator('#in-reps').fill('20000');
  await page.locator('#calc').click();
  await page.locator('#in-dist').selectOption('gamma');
  await page.waitForFunction(() => document.querySelector('#calc').getAttribute('aria-disabled') === 'false');
  assert.equal(await page.locator('#result').textContent(), '');

  // Stop stays reachable when returning to the guide hides the running calculation's form.
  await page.locator('#calc').click();
  await page.locator('#mode-guide').click();
  assert.equal(await page.locator('#activity-stop').isVisible(), true);
  await page.locator('#activity-stop').click();
  await page.waitForFunction(() => document.querySelector('#activity-stop').hidden);

  await page.goto(base);
  await page.waitForSelector('#in-test');
  await openGuide();
  const welch = routes.find((r) => r.node.procedure_id === 'means.two.welch');
  await select(welch); await apply();
  const order = await page.evaluate(() => {
    const questions = document.querySelector('#guide-requirements').getBoundingClientRect();
    const rec = document.querySelector('#guide-recommendation').getBoundingClientRect();
    const calc = document.querySelector('#method').getBoundingClientRect();
    return rec.top >= questions.bottom && calc.top >= rec.bottom;
  });
  assert.equal(order, true);
  if (process.env.SSC_SCREENSHOTS) await mkdir(process.env.SSC_SCREENSHOTS, { recursive: true });
  for (const [width, theme] of [[1280, 'light'], [390, 'light'], [390, 'dark']]) {
    await page.setViewportSize({ width, height: 900 });
    await page.emulateMedia({ colorScheme: theme });
    assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), true, `No overflow at ${width}px`);
    if (process.env.SSC_SCREENSHOTS) await page.screenshot({ path: `${process.env.SSC_SCREENSHOTS}/guide-${width}-${theme}.png`, fullPage: true });
  }
  assert.deepEqual(errors, []);
  await context.close();

  // The guide must be included in the first-visit cache, not only work while online.
  const offline = await browser.newContext({ serviceWorkers: 'allow' });
  const offlinePage = await offline.newPage();
  await offlinePage.goto(base);
  await offlinePage.waitForSelector('#guide-questions input', { state: 'attached' });
  await offlinePage.evaluate(async () => {
    await navigator.serviceWorker.ready;
    const cache = await caches.open('ssc-dev');
    const expected = ['guide.js', 'guide-state.js', 'guide-data.js', 'selection.js', 'engine/ssc_web.wasm'];
    for (const file of expected) {
      if (!await cache.match(new URL(file + '?v=dev', location.href))) throw new Error(`Missing cached asset: ${file}`);
    }
  });
  await offline.setOffline(true);
  await offlinePage.reload();
  await offlinePage.waitForSelector('#guide-questions input', { state: 'attached' });
  await offline.close();
  console.log('Browser: calculation, same-method options, stale results, keyboard, history, mobile/dark layout, and offline reload passed.');
} finally {
  await browser.close();
  server.close();
}
