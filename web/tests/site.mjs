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
    if (path === '' || /^[a-z]+:/.test(path)) continue;
    seen.add(path);
    if (!existsSync(join(site, path))) missing.push(`${path} (from ${from})`);
    // the service worker itself is fetched without the browser's copy, and the page's own address has no version
    if (path !== 'sw.js' && !new URLSearchParams(query).getAll('v').includes(version)) unversioned.push(`${ref} (from ${from})`);
  }
}
for (const file of ['index.html', 'dag.html'])
  for (const m of (await text(file)).matchAll(/(?:src|href)="([^"#]+)"/g)) if (!/^https?:/.test(m[1])) reference(m[1], file);
for (const file of scripts) {
  const code = await text(file);
  for (const m of code.matchAll(/(?:from\s+|import\(\s*|new URL\(\s*|register\(\s*)(['"`])([^'"`]+)\1/g)) reference(m[2], file);
  if (file === 'sw.js') {
    const list = /const FILES = \[([^\]]*)\]/.exec(code)?.[1] ?? '';
    for (const m of list.matchAll(/'([^']*)'/g)) if (m[1] !== './') reference(m[1], 'sw.js (offline copies)');
  }
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
  await page.waitForSelector('#dag-add');
  await page.selectOption('#dag-example', 'mediation');
  await page.click('#dag-example-load');
  await page.waitForSelector('#dag-equations li');
  assert.equal(await page.locator('#dag-canvas .dag-node').count(), 4);
  assert.equal(await page.locator('#dag-canvas .dag-edge').count(), 5);
  assert.equal(await page.locator('#engine-status, #method, #guide').count(), 0);
  assert.equal(page.workers().length, 0, 'the DAG page does not start the calculator engine');
  assert.deepEqual([...new Set(await loaded(page))], [next], 'the simulator scripts use the current version');
  await page.click('#dag-simulate');
  await page.waitForSelector('#dag-results table');
  assert.match(await page.textContent('#dag-results'), /N = 500/);
  const firstData = await page.textContent('#dag-results');
  await page.click('#dag-simulate');
  await page.waitForSelector('#dag-results table');
  assert.equal(await page.textContent('#dag-results'), firstData, 'same graph and seed reproduce the trial');
  await page.selectOption('#dag-source', 'Y');
  await page.selectOption('#dag-target', 'M');
  await page.getByRole('button', { name: 'Add arrow', exact: true }).click();
  assert.match(await page.textContent('#dag-status'), /cycle/);
  assert.equal(await page.locator('#dag-canvas .dag-edge').count(), 5);
  await page.getByRole('button', { name: 'Add interaction', exact: true }).click();
  assert.match(await page.textContent('#dag-equations'), /0.2 × X × M/);
  assert.equal(await page.locator('#dag-results table').count(), 0, 'a structural edit clears old results');
  assert.equal(await page.isDisabled('#dag-download-csv'), true);
  await page.click('#dag-undo');
  assert.doesNotMatch(await page.textContent('#dag-equations'), /0.2 × X × M/);
  await page.selectOption('#dag-example', 'binary');
  await page.click('#dag-example-load');
  await page.click('#dag-simulate');
  await page.waitForSelector('#dag-results table');
  assert.match(await page.textContent('#dag-equations'), /logit/);
  await page.selectOption('#estimand-exposure', 'A');
  await page.selectOption('#estimand-outcome', 'Y');
  await page.selectOption('#estimand-measure', 'risk_difference');
  await page.fill('#estimand-alternative', '0.1');
  await page.getByRole('button', { name: 'Save estimand', exact: true }).click();
  assert.match(await page.textContent('#estimand-preview'), /do\(A = 1\)/);
  assert.match(await page.textContent('#power-method'), /proportion score/);
  await page.fill('#power-sizes', '60, 120');
  await page.fill('#power-repetitions', '40');
  await page.getByText('Population calibration precision', { exact: true }).click();
  await page.fill('#power-population', '1000');
  await page.click('#power-run');
  await page.waitForSelector('#power-table');
  assert.equal(await page.locator('#power-table tbody tr').count(), 2);
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
  await page.fill('#estimand-alternative', '0.2');
  assert.equal(await page.locator('#power-table').count(), 0, 'an estimand edit invalidates power');
  assert.equal(await page.isDisabled('#power-run'), true, 'the draft target cannot run power');
  await page.getByRole('button', { name: 'Save estimand', exact: true }).click();
  assert.equal(await page.isDisabled('#power-run'), false);
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
  await page.selectOption('#dag-example', 'survival');
  await page.click('#dag-example-load');
  await page.selectOption('#dag-node-select', 'T');
  assert.equal(await page.inputValue('#dag-node-kind'), 'survival');
  assert.equal(await page.inputValue('#dag-censor-variable'), 'C');
  await page.selectOption('#estimand-exposure', 'A');
  await page.selectOption('#estimand-outcome', 'T');
  await page.selectOption('#estimand-measure', 'rmst_difference');
  await page.fill('#estimand-horizon', '28');
  await page.getByRole('button', { name: 'Save estimand', exact: true }).click();
  assert.match(await page.textContent('#estimand-preview'), /min\(T, 28\)/);
  await page.click('#dag-simulate');
  await page.waitForSelector('#dag-results table');
  assert.match(await page.textContent('#dag-results'), /Survival observations/);
  assert.match(await page.textContent('#dag-results'), /T_time/);
  // The graph is kept independently of calculator inputs, including when the page is reopened offline.
  down = true;
  await page.reload();
  await page.waitForSelector('#dag-equations li');
  assert.equal(await page.locator('#dag-canvas .dag-node').count(), 4);
  assert.equal(await page.inputValue('#estimand-measure'), 'rmst_difference');
  assert.equal(await page.inputValue('#estimand-horizon'), '28');
  await page.click('#power-run');
  await page.waitForSelector('#power-table');
  assert.match(await page.textContent('#power-method'), /RMST/);
  assert.equal(await page.locator('#power-table tbody tr').count(), 2, 'the cached worker calculates power offline');
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
