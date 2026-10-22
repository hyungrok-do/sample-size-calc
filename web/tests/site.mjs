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
for (const m of index.matchAll(/(?:src|href)="([^"#]+)"/g)) if (!/^https?:/.test(m[1])) reference(m[1], 'index.html');
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
const types = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.wasm': 'application/wasm' };
const prefix = '/sample-size-calc/';
const server = createServer((req, res) => {
  if (down) { req.socket.destroy(); return; }
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
async function ready(p) {
  try {
    await p.waitForFunction(() => /ready in/.test(document.getElementById('engine-status')?.textContent)
      && !document.getElementById('method').hidden, null, { timeout: 60000 });
  } catch (error) {
    const state = await p.evaluate(() => ({
      status: document.getElementById('engine-status')?.textContent, controlled: !!navigator.serviceWorker?.controller,
      title: document.title, error: document.getElementById('navigation-error')?.textContent,
    })).catch((e) => ({ page: e.message }));
    throw new Error(`the page did not become ready: ${JSON.stringify(state)}`);
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
  const context = await browser.newContext({ serviceWorkers: 'allow' });
  const page = await context.newPage();
  const errors = [];
  page.on('pageerror', (e) => errors.push(e.message));

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
  await page.goto(url);
  await ready(page);
  const isolated = await page.evaluate(() => self.crossOriginIsolated);
  const cores = await page.evaluate(() => navigator.hardwareConcurrency);
  const status = await page.textContent('#engine-status');
  if (kind === 'chromium' && cores > 1) assert.match(status, new RegExp(`${cores} threads`), 'second visit: every core');
  await calculates(page);
  console.log(`Site: second visit — ${isolated ? 'cross-origin isolated' : 'not isolated'}; ${status}`);

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
  assert.deepEqual(errors, []);
  await context.close();
  console.log(`Site: every check passed in ${kind}.`);
} finally {
  await browser.close();
  server.close();
}
