// Real browser -> Python transport -> Lean elaborator and kernel. No mocked solver.
import assert from 'node:assert/strict';
import { spawn } from 'node:child_process';
import { resolve } from 'node:path';
import { chromium } from 'playwright';
import { examples } from '../lean-native-examples.js';

const repo = resolve(new URL('../..', import.meta.url).pathname);
const server = spawn('python', ['scripts/serve_lean.py', '--port', '0'], { cwd: repo });
let startup = '';
const base = await new Promise((resolve, reject) => {
  const timeout = setTimeout(() => reject(Error('Native compiler service did not start')), 10000);
  server.stdout.on('data', chunk => {
    startup += chunk;
    const url = /http:\/\/127\.0\.0\.1:\d+/.exec(startup)?.[0];
    if (url) { clearTimeout(timeout); resolve(url); }
  });
  server.on('exit', code => { clearTimeout(timeout); reject(Error(`Native service exited ${code}`)); });
});
const browser = await chromium.launch({ executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined,
  args: ['--disable-dev-shm-usage', ...JSON.parse(process.env.SSC_CHROMIUM_ARGS || '[]')] });
const context = await browser.newContext({ viewport: { width: 1360, height: 1000 } });
const page = await context.newPage(), errors = [], violations = [];
page.on('pageerror', e => errors.push(e.message));
await page.exposeFunction('recordViolation', e => violations.push(e));
await page.addInitScript(() => document.addEventListener('securitypolicyviolation', e => window.recordViolation(e.violatedDirective)));
async function compile(source) {
  if (source) await page.locator('#lean-source').fill(source);
  await page.locator('#lean-compile').click();
  await page.waitForFunction(() => document.getElementById('lean-status').textContent === 'Lean compilation complete.' || !document.getElementById('lean-error').hidden, null, { timeout: 90000 });
  assert.equal(await page.locator('#lean-error').isVisible(), false, await page.locator('#lean-error').textContent());
}
try {
  await page.goto(`${base}/lean.html`);
  await page.waitForFunction(() => !document.getElementById('lean-compile').disabled);
  assert.equal(await page.locator('.workspace-tabs a').count(), 3);
  assert.match(await page.locator('#lean-source').inputValue(), /noncomputable def response/);
  await compile();
  assert.match(await page.locator('#lean-solution').textContent(), /x = \(Real.log 5 - 1\) \/ 2/);
  assert.match(await page.locator('#lean-proof-status').textContent(), /Lean kernel checked/);
  assert.match(await page.locator('#lean-derivation').textContent(), /exponential/);
  assert.match(await page.locator('#lean-theorem-type').textContent(), /response_solution/);
  assert.equal(await page.locator('#lean-cpp-details').getAttribute('open'), null);
  if (process.env.SCREENSHOTS) await page.screenshot({ path: resolve(process.env.SCREENSHOTS, 'native-lean-desktop.png'), fullPage: true });
  await compile(examples['custom-rule']);
  assert.equal(await page.locator('#lean-solution').textContent(), 'x = 3 / 2');
  assert.match(await page.locator('#lean-derivation').textContent(), /cubeLink_zero/);
  await compile(examples.numerical);
  assert.match(await page.locator('#lean-proof-status').textContent(), /unresolved/);
  assert.match(await page.locator('#lean-cpp').textContent(), /residual = \[=\]\(double x\)/);
  await compile(examples.parameters);
  assert.match(await page.locator('#lean-cpp-parameters').textContent(), /Runtime inputs: a, b, y/);
  assert.match(await page.locator('#lean-cpp-assumptions').textContent(), /a ≠ 0/);
  assert.match(await page.locator('#lean-cpp-assumptions').textContent(), /0 < y/);
  assert.match(await page.locator('#lean-cpp-command').textContent(), /'a=VALUE' 'b=VALUE' 'y=VALUE'/);
  await compile(examples.precision);
  assert.equal(await page.locator('#lean-solution').textContent(), 'n = variance / targetVariance');
  assert.match(await page.locator('#lean-cpp-command').textContent(), /'variance=VALUE' 'targetVariance=VALUE'/);
  await compile(examples.conditional);
  assert.match(await page.locator('#lean-cpp').textContent(), /all_reals/);
  assert.equal(await page.locator('#lean-cpp-guards').getAttribute('hidden'), '');
  await compile('solve one for x : ℝ where 2*x=6\nsolve two for y : ℝ where y^2=9');
  assert.equal(await page.locator('#lean-query option').count(), 2);
  await page.locator('#lean-query').selectOption('1');
  assert.match(await page.locator('#lean-solution').textContent(), /y = 3 ∨ y = -3/);
  await page.locator('#lean-source').fill('solve for x : ℝ where nonexistent x = 3');
  assert.equal(await page.locator('#lean-compiled').isVisible(), false);
  await page.locator('#lean-compile').click();
  await page.waitForFunction(() => !document.getElementById('lean-error').hidden, null, { timeout: 90000 });
  assert.match(await page.locator('#lean-error').textContent(), /unknown identifier/);
  await page.locator('#lean-source').fill('#eval IO.sleep 20000\nsolve for x : ℝ where x=1');
  await page.locator('#lean-compile').click();
  await page.locator('#lean-stop').click();
  await page.waitForFunction(() => document.getElementById('lean-status').textContent === 'Stopped.');
  await compile(examples.branches);
  await page.setViewportSize({ width: 390, height: 844 });
  assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth), false);
  if (process.env.SCREENSHOTS) await page.screenshot({ path: resolve(process.env.SCREENSHOTS, 'native-lean-mobile.png'), fullPage: true });
  await page.reload();
  assert.equal(await page.locator('#lean-source').inputValue(), examples.branches);
  assert.equal(await page.locator('#lean-compiled').isVisible(), false);
  // Static hosting must never pretend that a kernel has run.
  await context.route('**/api/lean/health', route => route.fulfill({ status: 404, body: '' }));
  await page.reload();
  await page.waitForFunction(() => !document.getElementById('lean-setup').hidden);
  assert.equal(await page.locator('#lean-compile').isDisabled(), true);
  assert.match(await page.locator('#lean-connection-status').textContent(), /unavailable/);
  assert.deepEqual(errors, []); assert.deepEqual(violations, []);
  console.log('Native Lean editor: actual kernel results, user rules, fallback, multiple queries, errors, cancellation, mobile, persistence, and static-host state passed.');
} finally {
  await context.close(); await browser.close(); server.kill('SIGINT');
}
