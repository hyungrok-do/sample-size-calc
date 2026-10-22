// The calculator page in a real browser with the real engine, served cross-origin isolated (the multi-core engine):
// each fix from the testers' reports, responsiveness with a 20,000-scenario answer, and the fallback to the
// one-thread engine.  One line per check; any failure fails the run.  SSC_BROWSER picks chromium (default),
// firefox or webkit.
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { dirname, resolve, extname, sep } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

const web = resolve(dirname(fileURLToPath(import.meta.url)), '..');
const types = { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.wasm': 'application/wasm' };
const isolation = { 'Cross-Origin-Opener-Policy': 'same-origin', 'Cross-Origin-Embedder-Policy': 'require-corp' };
const server = createServer(async (req, res) => {
  const pathname = new URL(req.url, 'http://localhost').pathname;
  const path = resolve(web, '.' + (pathname === '/' ? '/index.html' : pathname));
  if (!path.startsWith(web + sep)) { res.writeHead(403).end(); return; }
  try {
    const body = await readFile(path);
    res.writeHead(200, { 'Content-Type': types[extname(path)] ?? 'application/octet-stream', 'Cache-Control': 'no-store', ...isolation });
    res.end(body);
  } catch { res.writeHead(404).end(); }
});
await new Promise((done) => server.listen(0, '127.0.0.1', done));
const base = `http://127.0.0.1:${server.address().port}/`;

const kind = process.env.SSC_BROWSER ?? 'chromium';
const browser = await { chromium, firefox, webkit }[kind].launch(kind === 'chromium'
  ? { executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined, args: ['--disable-dev-shm-usage'] } : {});

let failed = 0;
async function check(name, fn) {
  try {
    const detail = await fn();
    console.log(`ok    ${name}${detail ? `  — ${detail}` : ''}`);
  } catch (error) {
    failed++;
    console.log(`FAIL  ${name}  — ${String(error.message).split('\n')[0]}`);
  }
}
const errors = [];
async function open(url, options = {}) {
  const context = await browser.newContext({ serviceWorkers: 'block', ...options });
  const page = await context.newPage();
  page.on('pageerror', (e) => errors.push(`${url}: ${e.message}`));
  await page.goto(url);
  return page;
}
const ready = (p, timeout = 60000) => p.waitForFunction(() => /ready in/.test(document.getElementById('engine-status').textContent)
  && !document.getElementById('method').hidden, null, { timeout });
const answer = (p) => p.waitForFunction(() => !!document.querySelector('#result .statement, #result .error'), null, { timeout: 120000 });
const idle = (p) => p.waitForFunction(() => document.getElementById('stop').hidden, null, { timeout: 120000 });
async function show(p, id) {
  await p.evaluate((h) => { location.hash = h; }, id);
  await p.waitForFunction((h) => document.querySelector('#method-list .method[aria-current]')?.dataset.id === h, id);
}

try {
  let p = await open(base);
  await ready(p);
  await check('the page opens with the multi-core engine', async () => {
    const status = await p.textContent('#engine-status');
    const cores = await p.evaluate(() => navigator.hardwareConcurrency);
    if (cores > 1) assert.match(status, /\d+ threads/);
    return status;
  });
  await check('dropout-adjusted sample sizes are shown', async () => {
    await p.fill('#in-dropout', '0.2'); await p.click('#calc'); await answer(p);
    assert.ok((await p.$$eval('#result h3', (h) => h.map((x) => x.textContent))).includes('Sample size allowing for dropout'));
    const cells = await p.$$eval('#result .table-wrap', (w) => w[1]?.innerText ?? '');
    assert.match(cells, /44/);
  });
  await check('an invalid field is marked, and the mark goes when it is edited', async () => {
    await p.fill('#in-sigma', 'abc'); await p.click('#calc'); await answer(p);
    assert.equal(await p.getAttribute('#in-sigma', 'aria-invalid'), 'true');
    await p.fill('#in-sigma', '1');
    assert.equal(await p.getAttribute('#in-sigma', 'aria-invalid'), null);
  });
  await check('alpha 5e-8 shows as 5e-08, as on the command line', async () => {
    await p.click('#solve label >> nth=1'); await p.fill('#in-alpha', '5e-8'); await p.fill('#in-dropout', '0');
    await p.click('#calc'); await answer(p);
    assert.match(await p.innerText('#result .table-wrap'), /5e-08/);
  });
  await check('focus stays on Calculate after a keyboard calculation', async () => {
    await p.focus('#calc'); await p.keyboard.press('Enter'); await answer(p); await p.waitForTimeout(200);
    assert.equal(await p.evaluate(() => document.activeElement?.id), 'calc');
  });
  await check('the answer is announced in a small live region, not the whole table', async () => {
    assert.ok((await p.textContent('#announce')).length > 20);
    assert.equal(await p.getAttribute('#result', 'aria-live'), null);
  });
  await check('help text expands on a click', async () => {
    await p.click('#in-alpha-help');
    assert.ok(await p.evaluate(() => document.getElementById('in-alpha-help').classList.contains('open')));
  });
  await check('clicking the method already open keeps its inputs', async () => {
    await p.fill('#in-sigma', '2'); await p.click('#method-list .method[data-id="means.two.t"]'); await p.waitForTimeout(300);
    assert.equal(await p.inputValue('#in-sigma'), '2');
  });
  await check('Back returns to the previous method', async () => {
    await p.click('#method-list .method[data-id="prop.two"]'); await p.waitForFunction(() => /Proportions/.test(document.getElementById('m-name').textContent));
    await p.click('#method-list .method[data-id="corr.one"]'); await p.waitForFunction(() => /Correlation/.test(document.getElementById('m-name').textContent));
    await p.goBack();
    await p.waitForFunction(() => /Two Proportions/.test(document.getElementById('m-name').textContent), null, { timeout: 5000 });
  });
  for (const [id, title] of [['gs.means', 'Boundaries'], ['prop.one.simon', ''], ['sim.multiarm', '']]) {
    await check(`the extra tables of ${id} are shown`, async () => {
      await show(p, id); await p.click('#calc'); await answer(p);
      const titles = await p.$$eval('#result h3', (h) => h.map((x) => x.textContent));
      assert.ok(titles.length > 0 && (!title || titles.includes(title)), titles.join(', '));
      return titles.join(', ');
    });
  }
  await check('no field asks for the decimal pad, which has no minus sign on iPhones', async () => {
    await show(p, 'gs.means');
    const modes = await p.$$eval('.field input[type=text]', (xs) => xs.map((x) => `${x.id}:${x.getAttribute('inputmode') || '-'}`));
    assert.ok(!modes.some((m) => m.endsWith(':decimal')) && modes.includes('in-gamma:-'), modes.join(' '));
  });
  await check('a Unicode minus sign and a no-break space are read', async () => {
    await show(p, 'means.two.t');
    await p.click('#solve label >> nth=1'); await p.fill('#in-mu1', '−0.5'); await p.fill('#in-alpha', '0.05 ');
    await p.click('#calc'); await answer(p);
    assert.ok(await p.$('#result .statement'), (await p.textContent('#result')).slice(0, 80));
  });
  await check('"10-20" asks for "10 to 20"', async () => {
    await p.fill('#in-mu1', '10-20'); await p.click('#calc'); await answer(p);
    assert.match(await p.textContent('#result'), /10 to 20/);
  });
  await check('an answer is not shown under a form opened since', async () => {
    await show(p, 'sim.means.robust');
    await p.click('#solve label >> nth=1'); await p.fill('#in-reps', '200000'); await p.click('#calc'); await p.waitForTimeout(400);
    await show(p, 'means.two.t'); await show(p, 'sim.means.robust');
    await idle(p); await p.waitForTimeout(300);
    assert.equal(await p.$('#result .statement'), null);
  });
  await check('search works during a calculation, and searches again after it', async () => {
    await p.click('#solve label >> nth=1'); await p.fill('#in-reps', '100000'); await p.click('#calc'); await p.waitForTimeout(300);
    await p.fill('#search', 'log-rank'); await p.waitForTimeout(400);
    const during = await p.$$eval('#method-list .method', (b) => b.length);
    await idle(p); await p.waitForTimeout(600);
    const note = await p.$eval('#method-list', (n) => n.querySelector('p')?.textContent || '');
    assert.ok(during > 0 && !/while the calculation/.test(note), `${during} during; note "${note}"`);
    await p.fill('#search', '');
  });
  await check('the speed test has its own Stop and keeps the answer\'s timing', async () => {
    await show(p, 'means.two.t'); await p.click('#calc'); await answer(p);
    const timing = await p.textContent('#timing');
    await p.click('#bench'); await p.waitForTimeout(300);
    assert.ok(await p.isVisible('#bench-stop'));
    await p.click('#bench-stop'); await p.waitForTimeout(500);
    assert.equal(await p.textContent('#timing'), timing);
    assert.match(await p.textContent('#bench-out'), /Stopped/);
  });
  await check('a 20,000-scenario answer never stalls the page for long', async () => {
    await show(p, 'means.two.t');
    await p.click('#solve label >> nth=1'); await p.fill('#in-n1', '2 to 10001'); await p.fill('#in-alpha', '0.05 0.01');
    await p.evaluate(() => {
      window.__gap = 0; let last = performance.now();
      const tick = (t) => { window.__gap = Math.max(window.__gap, t - last); last = t; if (!window.__stopTicks) requestAnimationFrame(tick); };
      requestAnimationFrame(tick);
    });
    await p.click('#calc'); await answer(p);
    const drawn = await p.waitForFunction(() => document.querySelectorAll('#result tbody tr').length >= 1000
      || !!document.querySelector('#result .error'), null, { timeout: 60000 }).then(() => true, () => false);
    const rows = await p.evaluate(() => document.querySelectorAll('#result tbody tr').length);
    const shown = (await p.textContent('#result')).replace(/\s+/g, ' ').slice(0, 160);
    assert.ok(drawn && rows >= 1000, `${rows} rows${drawn ? '' : ' after 60 s'}; status "${await p.textContent('#engine-status')}"; result "${shown}"`);
    await p.waitForTimeout(500);
    const gap = await p.evaluate(() => { window.__stopTicks = true; return window.__gap; });
    assert.ok(gap < 2000, `longest pause ${Math.round(gap)} ms`);
    assert.match(await p.textContent('#result'), /20,000/);
    return `longest pause ${Math.round(gap)} ms`;
  });
  await p.context().close();

  for (const hash of ['#%', '#100%', '#%E0%A4%A', '#nonexistent', '#PROP.TWO']) {
    await check(`the address ${hash} opens the first method and is replaced`, async () => {
      const q = await open(`${base}${hash}`); await ready(q);
      assert.equal(await q.evaluate(() => location.hash), '#means.two.t');
      await q.context().close();
    });
  }
  await check('a malformed address typed later keeps the method open', async () => {
    const q = await open(`${base}#prop.two`); await ready(q);
    await q.evaluate(() => { location.hash = '%'; }); await q.waitForTimeout(300);
    assert.equal(await q.evaluate(() => location.hash), '#prop.two');
    assert.equal(await q.isVisible('#method'), true);
    await q.context().close();
  });
  await check('a long error wraps at phone width, and fields use 16px text (no zoom on iPhones)', async () => {
    const q = await open(`${base}#means.two.t`, { viewport: { width: 320, height: 568 }, deviceScaleFactor: 2 }); await ready(q);
    await q.click('#solve label >> nth=1'); await q.fill('#in-n1', 'x'.repeat(300)); await q.click('#calc'); await answer(q);
    const [scroll, width] = await q.evaluate(() => [document.documentElement.scrollWidth, innerWidth]);
    assert.ok(scroll <= width, `${scroll} / ${width}`);
    const fonts = await q.$$eval('#search, .field input, .field select', (xs) => [...new Set(xs.map((x) => getComputedStyle(x).fontSize))]);
    assert.ok(fonts.every((f) => parseFloat(f) >= 16), fonts.join(' '));
    await q.context().close();
  });
  await check('a search typed while the engine loads is applied once it is ready', async () => {
    const context = await browser.newContext({ serviceWorkers: 'block' });
    await context.route('**/engine/*.wasm*', async (r) => { await new Promise((x) => setTimeout(x, 3000)); await r.continue(); });
    const q = await context.newPage();
    q.on('pageerror', (e) => errors.push(`slow load: ${e.message}`));
    await q.goto(`${base}#means.two.t`);
    await q.fill('#search', 'logrank'); await q.waitForTimeout(200); await q.fill('#search', ''); await q.fill('#search', 'logrank');
    await ready(q); await q.waitForTimeout(500);
    const n = await q.$$eval('#method-list .method', (b) => b.length);
    assert.ok(n > 0 && n < 10, `${n} methods listed`);
    await context.close();
  });

  // Workers' requests are routed reliably only in Chromium.
  if (kind === 'chromium') {
    async function fallBack(name, route, timeout) {
      await check(name, async () => {
        const context = await browser.newContext({ serviceWorkers: 'block' });
        await route(context);
        const q = await context.newPage();
        q.on('pageerror', (e) => errors.push(`${name}: ${e.message}`));
        const started = Date.now();
        await q.goto(`${base}#means.two.t`);
        await ready(q, timeout);
        assert.match(await q.textContent('#engine-status'), /1 thread/);
        await q.click('#calc'); await answer(q);
        assert.ok(await q.$('#result .statement'));
        await context.close();
        return `one thread after ${((Date.now() - started) / 1000).toFixed(1)} s`;
      });
    }
    await fallBack('without the multi-core engine file, the one-thread engine calculates', (c) =>
      c.route('**/engine/ssc_web_mt.wasm*', (r) => r.fulfill({ status: 404, body: '' })), 30000);
    const glue = await readFile(resolve(web, 'engine/ssc_web_mt.mjs'), 'utf8');
    const stuck = glue.replace('new Worker(new URL(import.meta.url)', 'new Worker(URL.createObjectURL(new Blob([""],{type:"text/javascript"}))');
    assert.notEqual(stuck, glue, 'the multi-core engine starts its threads from its own address (web/build.sh)');
    await fallBack('when the multi-core engine\'s threads never start, the one-thread engine takes over', (c) =>
      c.route('**/engine/ssc_web_mt.mjs*', (r) => r.fulfill({ body: stuck, contentType: 'text/javascript', headers: isolation })), 60000);
  }

  await check('no page errors', async () => { assert.deepEqual(errors, []); });
} finally {
  await browser.close();
  server.close();
}
console.log(failed ? `${failed} check(s) failed in ${kind}` : `Page: every check passed in ${kind}.`);
process.exitCode = failed ? 1 : 0;
