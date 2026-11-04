// The site as published (web/assemble.sh), served as GitHub Pages serves it: under a subpath, without the
// cross-origin isolation headers, and with files a browser may keep for 10 minutes.  First, every file the page can
// load must have been published, with the version in its address.  Then, in a browser: a first visit, the guide,
// the multi-core engine from the second visit (the service worker adds the headers), a start without a network,
// and an update that must never pair the new page with older files.
// node tests/site.mjs <site dir>; SSC_BROWSER picks chromium (default), firefox or webkit.
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { createServer } from 'node:http';
import { existsSync } from 'node:fs';
import { readFile, readdir } from 'node:fs/promises';
import { extname, join, resolve, sep } from 'node:path';
import { chromium, firefox, webkit } from 'playwright';
import { policy, watchViolations } from './csp.mjs';

const site = resolve(process.argv[2] ?? '.');
// SCREENSHOTS=<folder>: also save pictures of the DAG canvas and of the whole DAG page there (CI uploads them for review)
const shots = process.env.SCREENSHOTS;
const shoot = async (page, name) => { if (shots) await page.locator('.dag-scroll:has(#dag-canvas)').screenshot({ path: join(shots, `${name}.png`) }); };
const shootPage = async (page, name) => { if (shots) await page.screenshot({ path: join(shots, `${name}.png`), fullPage: true }); };
assert.ok(existsSync(join(site, 'index.html')), 'usage: node tests/site.mjs <site dir made by web/assemble.sh>');
const text = (file) => readFile(join(site, file), 'utf8');

// ----------------------------------------------------------------------------------------- the files published
const index = await text('index.html');
const version = /src="app\.js\?v=([^"&]+)"/.exec(index)?.[1];
assert.ok(version && version !== 'dev', 'index.html loads app.js with a stamped version');
const missing = [], unversioned = [];
const scripts = (await readdir(site)).filter((f) => f.endsWith('.js'));
const seen = new Set();
function reference(raw, from) {
  // an engine named by a variable: both engines
  for (const ref of raw.includes('${') ? ['ssc_web', 'ssc_web_mt'].map((v) => raw.replace(/\$\{[^}]*\}/g, v)) : [raw]) {
    const [path, query = ''] = ref.replace(/^\.\//, '').split('?');
    // These runtime URLs use the compiler manifest's content hashes and version.
    if (from === 'lean-wasm-worker.js' && ['lean-wasm/', 'manifest.json', 'lean.js'].includes(path)) continue;
    if (path === '' || /^[a-z]+:/.test(path)) continue;
    seen.add(path);
    if (!existsSync(join(site, path))) missing.push(`${path} (from ${from})`);
    // the service worker itself is fetched without the browser's copy, and the page's own address has no version
    if (path !== 'sw.js' && !new URLSearchParams(query).getAll('v').includes(version)) unversioned.push(`${ref} (from ${from})`);
  }
}
for (const file of ['index.html', 'dag.html', 'lean.html', 'study.html'])
  for (const m of (await text(file)).matchAll(/(?:src|href)="([^"#]+)"/g)) if (!/^https?:/.test(m[1])) reference(m[1], file);
for (const file of scripts) {
  const code = await text(file);
  for (const m of code.matchAll(/(?:from\s+|import\(\s*|new URL\(\s*|register\(\s*)(['"`])([^'"`]+)\1/g)) reference(m[2], file);
  if (file === 'sw.js') {
    const list = /const FILES = \[([^\]]*)\]/.exec(code)?.[1] ?? '';
    for (const m of list.matchAll(/'([^']*)'/g)) if (m[1] !== './') reference(m[1], 'sw.js (offline copies)');
  }
}
const compilerManifest = JSON.parse(await text('lean-wasm/manifest.json'));
for (const entry of [...Object.entries(compilerManifest.runtime).map(([file,meta])=>({file,...meta})), ...compilerManifest.packs]) {
  const bytes = await readFile(join(site,'lean-wasm',entry.file));
  assert.equal(createHash('sha256').update(bytes).digest('hex'),entry.sha256,`compiler asset ${entry.file}`);
}
assert.deepEqual(missing, [], 'files the page loads but the site does not have');
assert.deepEqual(unversioned, [], 'files loaded without the version in their address');
for (const file of scripts) assert.ok(seen.has(file), `${file} is published but nothing loads it`);
for (const engine of ['ssc_web', 'ssc_web_mt']) for (const ext of ['mjs', 'wasm']) assert.ok(seen.has(`engine/${engine}.${ext}`));
console.log(`Site: version ${version}; all ${seen.size} files the page loads are published and versioned.`);
// stamping must leave the <style> as its hash in the policy has it
policy(index, await text('sw.js'));
console.log('Site: the Content-Security-Policy lets in the page\'s <style>, and the service worker sends the same policy.');
if (process.argv.includes('--files-only')) process.exit(0);

// ------------------------------------------------------------------------------- served as GitHub Pages serves it
// Two versions to publish in turn: this one, and the same files as the next version (another version string, and
// a marker in the page).
const files = new Map();
for (const file of [...seen, 'index.html']) files.set(file, await readFile(join(site, file)));
const next = `${version}n`;
const nextFiles = new Map([...files].map(([file, body]) => {
  if (!/\.(html|js)$/.test(file) || file.startsWith('engine/')) return [file, body];
  let s = body.toString('utf8').replaceAll(`v=${version}`, `v=${next}`).replaceAll(`'ssc-${version}'`, `'ssc-${next}'`);
  if (file === 'index.html') s = s.replace('</head>', '<meta name="ssc-test" content="next">\n</head>');
  return [file, Buffer.from(s)];
}));
let live = files;
let down = false;  // the server drops every connection, as when the network is gone
// what a failed check reports: the requests dropped while the network was gone (and what kind of request each was),
// requests the browser saw fail, and the errors and warnings in the console
const dropped = [], failed = [], logged = [];
const types = { '.html': 'text/html', '.css': 'text/css', '.js': 'text/javascript', '.mjs': 'text/javascript', '.wasm': 'application/wasm' };
const prefix = '/sample-size-calc/';
const server = createServer((req, res) => {
  if (down) {
    dropped.push(`${req.url} (${req.headers['sec-fetch-dest'] ?? 'no destination'})`);
    req.socket.destroy();
    return;
  }
  const pathname = new URL(req.url, 'http://localhost').pathname;
  if (!pathname.startsWith(prefix)) { res.writeHead(404).end(); return; }
  const file = pathname.slice(prefix.length) || 'index.html';
  const body = live.get(file);
  if (!body || file.includes('..') || file.split('/').includes('')) { res.writeHead(404).end(); return; }
  const etag = `"${createHash('sha1').update(body).digest('hex')}"`;
  const headers = { 'Content-Type': types[extname(file)] ?? 'application/octet-stream', 'Cache-Control': 'max-age=600', ETag: etag };
  if (req.headers['if-none-match'] === etag) { res.writeHead(304, headers).end(); return; }
  res.writeHead(200, headers).end(body);
});
await new Promise((done) => server.listen(0, '127.0.0.1', done));
const url = `http://127.0.0.1:${server.address().port}${prefix}`;

const kind = process.env.SSC_BROWSER ?? 'chromium';
const browser = await { chromium, firefox, webkit }[kind].launch(kind === 'chromium'
  ? { executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined, args: ['--disable-dev-shm-usage'] } : {});
const context = await browser.newContext({ serviceWorkers: 'allow' });
const violations = [];  // what the Content-Security-Policy blocked: nothing, in every visit below
await watchViolations(context, violations);
// the offline copies the service worker keeps, by cache
const stored = (p) => p.evaluate(async () => {
  const kept = {};
  for (const name of await caches.keys())
    kept[name] = (await (await caches.open(name)).keys()).map((r) => new URL(r.url).pathname.split('/').pop() + new URL(r.url).search);
  return kept;
});
let keptOnline = null;  // the copies kept just before the network goes down
async function ready(p) {
  try {
    await p.waitForFunction(() => /ready in/.test(document.getElementById('engine-status')?.textContent)
      && !document.getElementById('method').hidden, null, { timeout: 60000 });
  } catch (error) {
    const state = await p.evaluate(() => ({
      status: document.getElementById('engine-status')?.textContent, controlled: !!navigator.serviceWorker?.controller,
      isolated: self.crossOriginIsolated, title: document.title,
      error: document.getElementById('navigation-error')?.textContent,
    })).catch((e) => ({ page: e.message }));
    const kept = await stored(p).catch((e) => e.message);
    const recent = (list) => list.slice(-15);
    throw new Error(`the page did not become ready: ${JSON.stringify({ ...state, kept, keptOnline,
      dropped: recent(dropped), failed: recent(failed), console: recent(logged) }, null, 1)}`);
  }
}
async function calculates(p) {
  await p.click('#calc');
  await p.waitForSelector('#result .statement', { timeout: 60000 });
}
// the versions in the addresses of the scripts this page loaded, the engine's worker included
const loaded = (p) => p.evaluate(() => performance.getEntriesByType('resource').map((e) => new URL(e.name))
  .filter((u) => /\.(m?js|wasm)$/.test(u.pathname) && !u.pathname.endsWith('/sw.js')).map((u) => u.searchParams.get('v')));
try {
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));
  context.on('requestfailed', (r) => failed.push(`${r.url()}: ${r.failure()?.errorText}`));
  page.on('console', (m) => { if (m.type() === 'error' || m.type() === 'warning') logged.push(`${m.type()}: ${m.text()}`); });

  await page.goto(url);
  await ready(page);
  assert.match(await page.textContent('#engine-status'), /1 thread/, 'first visit: the one-thread engine (no isolation yet)');
  await calculates(page);
  await page.click('#mode-guide');
  await page.check(`#guide input[name="q1"][value="0"]`);
  assert.ok(await page.locator('.guide-candidate').count() > 1, 'the guide lists analyses after one answer');
  await page.locator('.guide-candidate').first().click();
  await page.waitForFunction(() => !document.getElementById('method').hidden);
  console.log(`Site: first visit — ${await page.textContent('#engine-status')}; calculation and guide work.`);

  await page.waitForFunction(() => navigator.serviceWorker.controller !== null, null, { timeout: 30000 });
  const served = await page.goto(url);
  // the service worker's headers (as Chromium reports them), the policy with frame-ancestors among them
  if (kind === 'chromium') assert.match(served.headers()['content-security-policy'] ?? '', /; frame-ancestors 'none'$/);
  await ready(page);
  const isolated = await page.evaluate(() => self.crossOriginIsolated);
  const cores = await page.evaluate(() => navigator.hardwareConcurrency);
  const status = await page.textContent('#engine-status');
  if (kind === 'chromium' && cores > 1) assert.match(status, new RegExp(`${cores} threads`), 'second visit: every core');
  await calculates(page);
  console.log(`Site: second visit — ${isolated ? 'cross-origin isolated' : 'not isolated'}; ${status}`);

  // any other address is answered but not kept, so links cannot fill the site's storage with copies
  await page.goto(`${url}?utm_source=link#prop.two`);
  await ready(page);
  await page.evaluate(() => Promise.all(['engine/ssc_web.wasm?x=1', '?junk=1'].map((a) => fetch(a).then((r) => r.arrayBuffer()))));
  await page.waitForTimeout(500);
  const kept = await page.evaluate(async () => {
    const all = [];
    for (const name of await caches.keys()) all.push(...(await (await caches.open(name)).keys()).map((r) => new URL(r.url).search));
    return all;
  });
  assert.ok(kept.length && kept.every((q) => q === '' || /^\?(engine=ssc_web(_mt)?&)?v=[^&]+$/.test(q)), `kept: ${kept.join(' ')}`);
  // back at the page's own address: the browser's cache keeps each address (10 minutes on GitHub Pages), and the
  // update below reloads this one
  await page.goto(url);
  await ready(page);

  keptOnline = await stored(page);
  // WebKit can lose every copy it keeps when the network goes down while the visit's own work is still under way:
  // about one offline reload in four lost them in CI when the network went down the moment the page was ready,
  // and none of 16 after a second's wait (an awaited update check of the service worker did not help).  The
  // network goes down a second after the page is ready, as for someone who goes offline after using it.
  await page.waitForTimeout(1000);
  down = true;
  await page.reload();
  await ready(page);
  await calculates(page);
  down = false;
  console.log('Site: without a network the page, guide and engine start from the offline copies.');

  // An update: the next version is published while this one is open and kept by the browser.
  live = nextFiles;
  await page.reload();
  await ready(page);
  assert.equal(await page.getAttribute('meta[name="ssc-test"]', 'content'), 'next', 'a reload brings the new page');
  assert.deepEqual([...new Set(await loaded(page))], [next], 'every script of the new page is the new version');
  assert.ok(page.workers().length > 0 && page.workers().every((w) => new URL(w.url()).searchParams.get('v') === next));
  await calculates(page);
  await page.goto(url);
  await ready(page);
  assert.deepEqual([...new Set(await loaded(page))], [next], 'and on the next visit');
  console.log('Site: after an update the page and every file it loads are the new version.');

  // An independent tab: no calculator scripts or engine, with the same offline/versioning guarantees.
  await page.getByRole('link', { name: 'DAG simulator', exact: true }).click();
  await page.waitForSelector('#dag-example option', { state: 'attached' });
  const desktopViewport = page.viewportSize();
  await page.setViewportSize({ width: 390, height: 844 });
  const fitsViewport = async () => assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth), true, 'the panel fits a narrow viewport');
  await fitsViewport();
  // one panel, top to bottom: no steps or tabs, and the graph, the question and power are all on the page
  assert.equal(await page.locator('[data-step], [data-dag-step], #dag-show-all, #dag-step-next').count(), 0, 'no steps or tabs');
  assert.equal(await page.isVisible('#dag-canvas'), true);
  // a blank start opens an empty model; a starting point opens its model, question and analysis together
  await page.selectOption('#dag-example', 'empty');
  await page.click('#dag-example-load');
  assert.equal(await page.locator('#dag-canvas .dag-node').count(), 0);
  await page.selectOption('#dag-example', 'trial');
  await page.click('#dag-example-load');
  assert.match(await page.textContent('#dag-node-count'), /2 variables.*sample size counts observations/);
  await fitsViewport();
  assert.equal(await page.isVisible('#power-run'), true, 'power is on the same panel as the graph');
  assert.match(await page.textContent('#power-question-summary'), /Outcome.*Treatment.*1 versus 0.*mean difference/);
  // the estimand card names the target, its formula, hypotheses and what the analysis estimates
  const estimandCard = await page.textContent('#estimand-card');
  for (const part of [/Average treatment effect \(ATE\) on the mean of Outcome/, /if everyone had Treatment \(A\) = 1/, /E\[Y under do\(A = 1\)\]/, /H₀: θ = 0; H₁: θ ≠ 0/, /Welch/])
    assert.match(estimandCard, part);
  await page.waitForFunction(() => /θ = 0\.5 under the graph as drawn/.test(document.querySelector('#estimand-card')?.textContent ?? ''));
  // the effect the model implies is shown before anything runs
  await page.waitForFunction(() => /mean difference 0\.5/.test(document.querySelector('#estimand-implied')?.textContent ?? ''));
  assert.equal(await page.getAttribute('#estimand-alternative', 'placeholder'), 'As drawn (0.5)');
  assert.equal(await page.isVisible('#estimand-active'), true, 'comparison levels are a primary control');
  assert.equal(await page.isVisible('#analysis-method'), true, 'the selected analysis is visible without a disclosure');
  await page.selectOption('#dag-example', 'repeated');
  await page.click('#dag-example-load');
  assert.equal(await page.isVisible('#dag-study'), true, 'a structured study shows its settings');
  await page.selectOption('#estimand-trajectory', 'change');
  assert.match(await page.textContent('#power-question-summary'), /change in Outcome from time 0 to 12/);
  await fitsViewport();
  await shootPage(page, 'dag-page-narrow');
  await page.reload();
  await page.waitForSelector('#power-run');
  assert.equal(await page.inputValue('#estimand-trajectory'), 'change', 'a reload retains the scientific target');
  await page.setViewportSize(desktopViewport);
  assert.equal(await page.isVisible('#dag-canvas'), true);
  assert.equal(await page.isVisible('#estimand-outcome'), true);
  await shootPage(page, 'dag-page-repeated');
  await page.selectOption('#dag-example', 'mediation');
  await page.click('#dag-example-load');
  await page.waitForSelector('#dag-equations li', { state: 'attached' });
  assert.equal(await page.locator('#dag-canvas .dag-node').count(), 4);
  assert.equal(await page.locator('#dag-canvas .dag-edge').count(), 5);
  assert.equal(await page.locator('#engine-status, #method, #guide').count(), 0);
  // each box is two lines (label; ID and distribution), and the selected variable's settings sketch its values
  // from a sample the engine draws in the background; that engine then closes
  assert.equal(await page.locator('#dag-canvas .dag-node[data-node="X"] text').first().textContent(), 'Covariate');
  assert.equal(await page.locator('#dag-canvas .dag-node[data-node="X"] .dag-node-family').textContent(), 'X · N(0, 1²)');
  await page.selectOption('#dag-node-select', 'M');
  await page.waitForSelector('#dag-node-sketch svg');
  assert.match(await page.textContent('#dag-node-sketch'), /Simulated: mean/);
  assert.match(await page.locator('#dag-canvas .dag-node[data-node="M"] > title').textContent(), /Mediator: on a path/);
  for (let i = 0; i < 50 && page.workers().length; i++) await page.waitForTimeout(100);
  assert.equal(page.workers().length, 0, 'no engine keeps running on the DAG page');
  assert.deepEqual([...new Set(await loaded(page))], [next], 'the simulator scripts use the current version');
  // the graph marks the power analysis's variables and the arrows that carry its effect
  assert.match(await page.textContent('#dag-canvas'), /Compared[\s\S]*Outcome|Outcome[\s\S]*Compared/);
  assert.equal(await page.locator('#dag-canvas .dag-edge.effect').count(), 3, 'A → M, M → Y and A → Y carry the total effect');
  await shoot(page, 'dag-canvas-mediation');
  assert.match(await page.textContent('#dag-paths'), /A → M → Y\s*causal · open/);
  assert.match(await page.textContent('#dag-paths'), /A → M ← X → Y\s*non-causal · blocked at collider M/);
  // a coefficient is changed where it is drawn, and roles are set on the graph (keyboard, menu, variable settings)
  await page.locator('[data-edge="A>Y"] .dag-edge-label').click();
  await page.fill('.dag-inline', '0.35');
  await page.keyboard.press('Enter');
  assert.match(await page.textContent('#dag-equations'), /0\.35 × A/);
  await page.click('#dag-undo');
  assert.match(await page.textContent('#dag-equations'), /0\.3 × A/);
  const box = (id) => page.locator(`#dag-canvas .dag-node[data-node="${id}"]`);
  await box('M').focus();
  await page.keyboard.press('o');
  assert.equal(await page.inputValue('#estimand-outcome'), 'M');
  assert.match(await page.textContent('#dag-paths'), /Paths from A to M/);
  await page.click('#dag-undo');
  assert.equal(await page.inputValue('#estimand-outcome'), 'Y');
  await box('X').click({ button: 'right' });
  await page.getByRole('menuitem', { name: 'Compare this variable', exact: true }).click();
  assert.equal(await page.inputValue('#estimand-exposure'), 'X');
  await page.selectOption('#dag-node-select', 'A');
  await page.click('#dag-role-exposure');
  assert.equal(await page.inputValue('#estimand-exposure'), 'A');
  await page.click('#dag-sim > summary');
  await page.click('#dag-simulate');
  await page.waitForSelector('#dag-results table');
  assert.match(await page.textContent('#dag-results'), /N = 500/);
  const firstData = await page.textContent('#dag-results');
  await page.click('#dag-simulate');
  await page.waitForSelector('#dag-results table');
  assert.equal(await page.textContent('#dag-results'), firstData, 'same graph and seed reproduce the trial');
  // a variable is offered only the new parents that keep the graph acyclic
  await page.selectOption('#dag-node-select', 'A');
  assert.deepEqual(await page.$$eval('#dag-source option', (o) => o.map((x) => x.value)), ['X']);
  await page.getByRole('button', { name: 'Add arrow', exact: true }).click();
  assert.equal(await page.locator('#dag-canvas .dag-edge').count(), 6);
  await page.click('#dag-undo');
  assert.equal(await page.locator('#dag-canvas .dag-edge').count(), 5);
  await page.selectOption('#dag-node-select', 'M');
  assert.equal(await page.locator('#dag-source').count(), 0, 'X and A already point at M, and Y depends on it');
  // an arrow drawn on the canvas that would close a cycle is refused
  await page.getByRole('button', { name: 'Start an arrow from Y', exact: true }).click();
  await page.locator('#dag-canvas .dag-node[data-node="M"]').click();
  assert.match(await page.textContent('#dag-status'), /cycle/);
  assert.equal(await page.locator('#dag-canvas .dag-edge').count(), 5);
  await page.keyboard.press('Escape');
  // an arrow dragged from a box's port onto another box; Ctrl+Z takes it back
  // Drawing replaces node elements; scroll the persistent canvas before measuring them.
  await page.locator('#dag-canvas').scrollIntoViewIfNeeded();
  const portBox = await page.getByRole('button', { name: 'Start an arrow from X', exact: true }).boundingBox();
  const targetBox = await page.locator('#dag-canvas .dag-node[data-node="A"]').locator('.dag-box').boundingBox();
  await page.mouse.move(portBox.x + portBox.width / 2, portBox.y + portBox.height / 2);
  await page.mouse.down();
  await page.mouse.move(targetBox.x + targetBox.width / 2, targetBox.y + targetBox.height / 2, { steps: 8 });
  await page.mouse.up();
  assert.equal(await page.locator('#dag-canvas [data-edge="X>A"]').count(), 1, await page.textContent('#dag-status'));
  await page.keyboard.press('Control+z');
  assert.equal(await page.locator('#dag-canvas .dag-edge').count(), 5);
  await page.selectOption('#dag-node-select', 'Y');
  await page.getByRole('button', { name: 'Add interaction', exact: true }).click();
  assert.match(await page.textContent('#dag-equations'), /0.2 × X × M/);
  assert.equal(await page.locator('#dag-results table').count(), 0, 'a structural edit clears old results');
  assert.equal(await page.isDisabled('#dag-download-csv'), true);
  await page.click('#dag-undo');
  assert.doesNotMatch(await page.textContent('#dag-equations'), /0.2 × X × M/);
  // Nonlinear terms are editable, previewed, exported and undoable; analysis
  // changes are separate from the graph and participate in the same history.
  await page.click('#dag-function-add > summary');
  await page.selectOption('#dag-function-left', 'X');
  await page.selectOption('#dag-function-kind', 'square');
  await page.fill('#dag-function-coefficient', '0.3');
  await page.getByRole('button', { name: 'Add response term', exact: true }).click();
  assert.match(await page.textContent('#dag-equations'), /square\(\(X - 0\) \/ 1\)/);
  await page.click('#dag-curve-panel > summary');
  assert.equal(await page.locator('#dag-curve-preview svg path').count(), 1);
  await page.click('#power-analysis-models > summary');
  await page.selectOption('#analysis-method', 'gformula');
  await page.selectOption('#analysis-basis', 'linear');
  const curveEvent = page.waitForEvent('download');
  await page.click('#dag-export');
  const curveFile = await curveEvent;
  const curveDesign = JSON.parse(await readFile(await curveFile.path(), 'utf8'));
  assert.equal(curveDesign.graph.nodes.find((n) => n.id === 'Y').functions[0].transform.kind, 'square');
  assert.equal(curveDesign.power.analysis.basis, 'linear');
  await page.click('#dag-undo');
  assert.equal(await page.inputValue('#analysis-basis'), 'match');
  await page.click('#dag-undo');
  assert.equal(await page.inputValue('#analysis-method'), 'auto');
  await page.click('#dag-undo');
  assert.doesNotMatch(await page.textContent('#dag-equations'), /square/);
  // an unmeasured variable is drawn dashed and simulated, but no analysis adjusts for it; its settings measure it
  await page.selectOption('#dag-example', 'unmeasured-confounding');
  await page.click('#dag-example-load');
  await page.waitForSelector('#dag-canvas .dag-node.unmeasured');
  assert.equal(await page.locator('#dag-canvas .dag-node.unmeasured .dag-node-family').textContent(), 'U · unmeasured');
  assert.match(await page.textContent('#dag-equations'), /\(unmeasured\)/);
  assert.match(await page.textContent('#dag-paths'), /A ← U → Y\s*back-door · open/);
  await page.selectOption('#dag-node-select', 'U');
  await page.check('#dag-node-measured');
  await page.waitForFunction(() => !document.querySelector('#dag-canvas .dag-node.unmeasured'));
  assert.match(await page.textContent('#dag-paths'), /A ← U → Y\s*back-door · blocked/);
  await page.click('#dag-undo');
  await page.waitForSelector('#dag-canvas .dag-node.unmeasured');
  // a relationship of unknown form: "?" on its arrows, a new shape on request, and TMLE with a Super Learner
  await page.selectOption('#dag-example', 'unknown-confounding');
  await page.click('#dag-example-load');
  await page.waitForFunction(() => /\?/.test(document.querySelector('[data-edge="X1>Y"] .dag-edge-label')?.textContent ?? ''));
  assert.match(await page.textContent('#power-method'), /TMLE.*Super Learner/);
  await page.selectOption('#dag-node-select', 'Y');
  assert.match(await page.textContent('#dag-unknown summary'), /Unknown form: f\(X1, X2\)/);
  assert.match(await page.textContent('#dag-unknown-redraw'), /now #5/);
  await page.click('#dag-unknown-redraw');
  await page.waitForFunction(() => /now #6/.test(document.querySelector('#dag-unknown-redraw')?.textContent ?? ''));
  await page.click('#dag-undo');
  // Clear empties the canvas in one step, and one undo brings the model back
  await page.click('#dag-clear');
  assert.equal(await page.locator('#dag-canvas .dag-node').count(), 0);
  assert.equal(await page.isDisabled('#dag-clear'), true);
  assert.match(await page.textContent('#dag-status'), /Cleared 4 variables/);
  await page.click('#dag-undo');
  assert.equal(await page.locator('#dag-canvas .dag-node').count(), 4);
  await page.selectOption('#dag-example', 'binary');
  await page.click('#dag-example-load');
  await page.click('#dag-simulate');
  await page.waitForSelector('#dag-results table');
  assert.match(await page.textContent('#dag-equations'), /logit/);
  await page.click('#dag-tidy');
  const at = (id) => page.locator(`#dag-canvas .dag-node[data-node="${id}"]`).getAttribute('transform');
  assert.equal((await at('X')).split(' ')[0], (await at('A')).split(' ')[0], 'Tidy puts the two roots in the first column');
  assert.notEqual((await at('Y')).split(' ')[0], (await at('X')).split(' ')[0]);
  // the power panel starts from a target it can run: the binary variable that affects others, on what it affects
  assert.equal(await page.inputValue('#estimand-exposure'), 'A');
  assert.equal(await page.inputValue('#estimand-outcome'), 'Y');
  assert.equal(await page.inputValue('#estimand-measure'), 'risk_difference');
  await page.fill('#estimand-alternative', '0.1');
  assert.match(await page.textContent('#power-estimand'), /do\(A = 1\)/);
  assert.match(await page.textContent('#power-method'), /proportion score/);
  await page.fill('#power-sizes', '60, 120');
  await page.fill('#power-repetitions', '40');
  await page.click('#power-more > summary');
  await page.fill('#power-population', '1000');
  await page.click('#power-run');
  await page.waitForSelector('#power-table');
  assert.equal(await page.locator('#power-table tbody tr').count(), 2);
  assert.equal(await page.locator('#power-plot .power-plot-dot').count(), 2);
  // the results lead with an answer, check each size's test, and offer the next sizes and a protocol paragraph
  assert.match(await page.textContent('#power-answer'), /for 80\.0% power|already gives|Not reached/);
  assert.match(await page.textContent('#power-table thead'), /Checks/);
  assert.match(await page.textContent('#power-summary'), /Question: Compare.*With 60 participants in total, estimated power is/);
  const refine = page.locator('#power-answer .dag-power-next');
  if (await refine.count()) assert.match(await refine.getAttribute('data-sizes'), /^\d+(,\d+)+$/);
  const displayedPower = await page.textContent('#power-table');
  await page.selectOption('#dag-node-select', 'X');
  assert.equal(await page.textContent('#power-table'), displayedPower, 'selecting a variable keeps calculated results');
  assert.match(await page.textContent('#power-calibration'), /planning effect 0.1/);
  const powerEvent = page.waitForEvent('download');
  await page.click('#power-json');
  const powerDownload = await powerEvent;
  const powerReport = JSON.parse(await readFile(await powerDownload.path(), 'utf8'));
  assert.equal(powerReport.estimand.measure, 'risk_difference');
  assert.deepEqual(powerReport.power.sampleSizes, [60, 120]);
  assert.equal(powerReport.results.length, 2);
  assert.equal(powerReport.engine.backend, 'C++');
  assert.match(await page.textContent('#power-engine'), /C\+\+\/WebAssembly/);
  assert.ok(['ssc_web', 'ssc_web_mt'].includes(powerReport.engine.variant));
  await page.fill('#estimand-alternative', '2');
  assert.equal(await page.locator('#power-table').count(), 0, 'an edit of the target clears the results');
  assert.match(await page.textContent('#power-method'), /between −1 and 1/);
  assert.equal(await page.isDisabled('#power-run'), true, 'an invalid target cannot run');
  await page.fill('#estimand-alternative', '0.2');
  assert.equal(await page.isDisabled('#power-run'), false, 'a valid edit can run at once, with no separate save');
  const csvEvent = page.waitForEvent('download');
  await page.click('#dag-download-csv');
  const csvDownload = await csvEvent;
  const csv = await readFile(await csvDownload.path(), 'utf8');
  assert.equal(csv.trim().split('\n').length, 501);
  assert.equal(csv.split('\n')[0], 'X,A,Y');
  const modelEvent = page.waitForEvent('download');
  await page.click('#dag-export');
  const modelDownload = await modelEvent;
  const model = JSON.parse(await readFile(await modelDownload.path(), 'utf8'));
  assert.equal(model.version, 1);
  assert.equal(model.graph.nodes.length, 3);
  assert.equal(model.estimand.measure, 'risk_difference');
  assert.equal(model.estimand.alternative, .2);
  assert.equal(model.power.repetitions, 40);
  await page.selectOption('#dag-example', 'mediation');
  await page.click('#dag-example-load');
  await page.setInputFiles('#dag-import', { name: 'trial-dag.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(model)) });
  await page.waitForFunction(() => document.querySelectorAll('#dag-canvas .dag-node').length === 3);
  assert.equal(await page.inputValue('#estimand-measure'), 'risk_difference');
  assert.equal(await page.inputValue('#power-repetitions'), '40');
  await page.fill('#estimand-alternative', '');
  await page.click('#power-more > summary');
  await page.selectOption('#power-null-mode', 'average_effect');
  assert.equal(await page.isVisible('#power-lower'), true, 'the average-effect reference exposes calibration bounds');
  await page.selectOption('#power-null-mode', 'no_effect');
  await page.click('#power-more > summary');
  await page.click('#power-analysis-models > summary');
  await page.selectOption('#analysis-method', 'tmle');
  await page.selectOption('#analysis-basis', 'quadratic');
  await page.fill('#analysis-folds', '2');
  await page.fill('#power-sizes', '400');
  await page.fill('#power-repetitions', '20');
  await page.click('#power-run');
  await page.waitForSelector('#power-table');
  assert.match(await page.textContent('#power-method'), /TMLE/);
  assert.match(await page.textContent('#power-results'), /CI coverage/);
  // A built-in U-shaped example supplies a ready-to-run dose contrast and
  // can instead be tested for its nonlinear regression block.
  await page.selectOption('#analysis-method', 'auto');
  await page.fill('#analysis-folds', '1');
  await page.selectOption('#dag-example', 'nonlinear-dose');
  await page.click('#dag-example-load');
  await page.selectOption('#estimand-kind', 'joint');
  await page.selectOption('#estimand-measure', 'nonlinearity');
  await page.click('#power-run');
  await page.waitForSelector('#power-table');
  assert.match(await page.textContent('#power-method'), /Joint Wald test of nonlinearity/);
  assert.match(await page.textContent('#power-truth'), /predictor-block coefficients/);
  // Repeated data retain independent-unit sizes, visit contrasts and a separately planned analysis.
  await page.selectOption('#dag-example', 'repeated');
  await page.click('#dag-example-load');
  assert.match(await page.textContent('#power-method'), /MMRM/);
  await page.click('#dag-timeline > summary');
  assert.equal(await page.locator('#dag-timeline rect').count(), 10, 'two shared variables and eight visit variables');
  await page.selectOption('#estimand-trajectory', 'slope');
  await page.click('#power-analysis-models > summary');
  await page.selectOption('#analysis-method', 'lmm');
  await page.selectOption('#analysis-time-basis', 'linear');
  await page.check('#analysis-random-slope');
  await page.fill('#power-sizes', '60');
  await page.fill('#power-repetitions', '20');
  await page.reload();
  await page.waitForSelector('#dag-equations li', { state: 'attached' });
  assert.equal(await page.inputValue('#estimand-trajectory'), 'slope');
  assert.equal(await page.inputValue('#analysis-time-basis'), 'linear');
  assert.equal(await page.isChecked('#analysis-random-slope'), true);
  await page.click('#power-run');
  await page.waitForSelector('#power-table');
  assert.match(await page.textContent('#power-method'), /Linear mixed-model/);
  assert.deepEqual(await page.locator('#power-table tbody tr:first-child td').allTextContents().then((v) => v.slice(0, 3)), ['60', '60', '240']);
  assert.match(await page.textContent('#power-results'), /conditional on the generating/);
  await page.selectOption('#dag-example', 'clustered-repeated');
  await page.click('#dag-example-load');
  await page.fill('#power-sizes', '12, 24');
  await page.fill('#power-repetitions', '20');
  await page.click('#power-run');
  await page.waitForSelector('#power-table');
  assert.deepEqual(await page.locator('#power-table tbody tr:first-child td').allTextContents().then((v) => v.slice(0, 3)), ['12', '144', '576']);
  assert.match(await page.textContent('#power-table thead'), /Number of clusters/);
  assert.match(await page.textContent('#power-method'), /cluster-robust/);
  // the study's levels are drawn as nested plates, with the shared effects on the outcome as circles
  assert.deepEqual(await page.locator('#dag-canvas .dag-plate').evaluateAll((g) => g.map((p) => p.dataset.plate)), ['cluster', 'subject', 'visit']);
  assert.equal(await page.getAttribute('#dag-level-clusters', 'aria-pressed'), 'true');
  assert.equal(await page.locator('#dag-canvas .dag-random').count(), 3, 'subject intercept, subject slope and cluster effects');
  await shoot(page, 'dag-canvas-clustered-repeated');
  // Scenarios load the question, estimator and independent-unit sizes together.
  assert.equal(await page.locator('#dag-example optgroup[label="TMLE"] option').count(), 6);
  await page.selectOption('#dag-example', 'tmle-repeated');
  assert.match(await page.textContent('#dag-scenario-description'), /no time-varying treatment or dropout/);
  await page.click('#dag-example-load');
  assert.equal(await page.inputValue('#analysis-method'), 'tmle');
  assert.equal(await page.inputValue('#estimand-trajectory'), 'average');
  await page.click('#power-analysis-models > summary');
  assert.equal(await page.isVisible('#analysis-folds'), true);
  assert.equal(await page.isVisible('#analysis-time-basis'), false);
  assert.equal(await page.inputValue('#analysis-folds'), '3');
  await page.fill('#power-sizes', '200');
  await page.fill('#power-repetitions', '20');
  await page.click('#power-run');
  await page.waitForSelector('#power-table');
  assert.match(await page.textContent('#power-results'), /largest inverse weight/);
  await page.selectOption('#dag-example', 'tmle-bounded');
  await page.click('#dag-example-load');
  assert.equal(await page.inputValue('#analysis-outcome-bounds'), '0, 100');
  await page.selectOption('#dag-example', 'mar-dropout');
  await page.click('#dag-example-load');
  await page.selectOption('#analysis-method', 'tmle');
  assert.match(await page.textContent('#power-method'), /requires complete outcomes/);
  assert.equal(await page.isEnabled('#power-run'), false);
  await page.setInputFiles('#dag-import', { name: 'trial-dag.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify(model)) });
  await page.selectOption('#dag-example', 'survival');
  await page.click('#dag-example-load');
  await page.selectOption('#dag-node-select', 'T');
  assert.equal(await page.inputValue('#dag-node-kind'), 'survival');
  assert.equal(await page.inputValue('#dag-censor-variable'), 'C');
  await page.selectOption('#estimand-exposure', 'A');
  await page.selectOption('#estimand-outcome', 'T');
  await page.selectOption('#estimand-measure', 'rmst_difference');
  assert.equal(await page.inputValue('#power-sizes'), '100, 200, 400', 'the scenario supplies its own size grid');
  await page.fill('#power-sizes', '120, 240');
  await page.fill('#power-repetitions', '20');
  await page.fill('#estimand-horizon', '28');
  assert.match(await page.textContent('#power-estimand'), /min\(T, 28\)/);
  // The repeated-study reload above restores the design with disclosures closed.
  if (!await page.isVisible('#dag-simulate')) await page.click('#dag-sim > summary');
  await page.click('#dag-simulate');
  await page.waitForSelector('#dag-results table');
  assert.match(await page.textContent('#dag-results'), /Censored %/);
  assert.match(await page.textContent('#dag-results'), /T_time/);
  // The graph is kept independently of calculator inputs, including when the page is reopened offline.
  down = true;
  await page.reload();
  await page.waitForSelector('#dag-equations li', { state: 'attached' });
  assert.equal(await page.locator('#dag-canvas .dag-node').count(), 4);
  assert.equal(await page.inputValue('#estimand-measure'), 'rmst_difference');
  assert.equal(await page.inputValue('#estimand-horizon'), '28');
  assert.equal(await page.inputValue('#power-sizes'), '120, 240');
  assert.equal(await page.inputValue('#power-repetitions'), '20');
  await page.click('#power-run');
  await page.waitForSelector('#power-table');
  assert.match(await page.textContent('#power-method'), /RMST/);
  assert.deepEqual(await page.locator('#power-table tbody tr td:first-child').allTextContents(), ['120', '240'], 'the cached worker calculates the saved size grid offline');
  await page.click('#dag-sim > summary');
  await page.click('#dag-simulate');
  await page.waitForSelector('#dag-results table');
  down = false;
  await page.getByRole('link', { name: 'Calculator', exact: true }).click();
  await ready(page);
  await calculates(page);
  console.log('Site: independent DAG tab, seeded datasets, survival and censoring, cycle checks, interactions, undo, CSV/JSON interchange, offline use, and return to calculator work.');
  assert.deepEqual(errors, []);
  assert.deepEqual(violations, [], 'the Content-Security-Policy blocked nothing');
  console.log(`Site: every check passed in ${kind}.`);
} finally {
  await browser.close();
  server.close();
}

