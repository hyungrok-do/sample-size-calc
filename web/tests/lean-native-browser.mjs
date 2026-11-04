// Real browser -> Python transport -> Lean elaborator and kernel. No mocked solver.
import assert from 'node:assert/strict';
import { spawn } from 'node:child_process';
import { resolve } from 'node:path';
import { chromium } from 'playwright';
import { examples } from '../lean-native-examples.js';
import { examples as studies } from '../lean-study-examples.js';

const repo = resolve(new URL('../..', import.meta.url).pathname);
const server = process.env.SSC_LEAN_TEST_URL ? null : spawn('python', ['scripts/serve_lean.py', '--port', '0'], { cwd: repo });
let startup = '';
const base = process.env.SSC_LEAN_TEST_URL || await new Promise((resolve, reject) => {
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
const context = await browser.newContext({ viewport: { width: 1360, height: 1000 },
  serviceWorkers: 'block',
  ignoreHTTPSErrors: process.env.SSC_LEAN_TEST_INSECURE_TLS === '1',
  httpCredentials: process.env.SSC_LEAN_TEST_PASSWORD ? {
    username: process.env.SSC_LEAN_TEST_USERNAME || 'researcher', password: process.env.SSC_LEAN_TEST_PASSWORD,
  } : undefined });
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
  if (process.env.SSC_LEAN_TEST_URL) assert.match(await page.locator('#lean-connection-status').textContent(), /Hosted Lean compiler connected/);
  assert.equal(await page.locator('.workspace-tabs a').count(), 3);
  assert.match(await page.locator('#lean-source').inputValue(), /gaussianContrast/);
  assert.equal(await page.locator('#lean-example option').count(), 5);
  await compile();
  assert.match(await page.locator('#lean-solution').textContent(), /n = 51 participants per arm/);
  assert.match(await page.locator('#lean-solution').textContent(), /Total participants: 102/);
  assert.match(await page.locator('#lean-study-numbers').textContent(), /79.8176%/);
  assert.match(await page.locator('#lean-derivation').textContent(), /Propagate the full covariance/);
  assert.match(await page.locator('#lean-theorem-type').textContent(), /variance study = 8 \/ 5/);
  if (process.env.SCREENSHOTS) await page.screenshot({ path: resolve(process.env.SCREENSHOTS, 'statistical-lean-desktop.png'), fullPage: true });
  await page.locator('#lean-example').selectOption('targeted-ate');
  await page.locator('#lean-load').click();
  assert.equal(await page.locator('#lean-compiled').isVisible(), false);
  await compile();
  assert.match(await page.locator('#lean-proof-status').textContent(), /Asymptotic Wald power/);
  assert.match(await page.locator('#lean-derivation').textContent(), /influence function/);
  assert.match(await page.locator('#lean-study-model').textContent(), /Effect-heterogeneity variance0.0006/);
  const downloadEvent = page.waitForEvent('download');
  await page.locator('#lean-cpp-details').evaluate(element => { element.open = true; });
  await page.locator('#lean-cpp-save').click();
  assert.equal((await downloadEvent).suggestedFilename(), 'study.cpp');
  await page.locator('#lean-cpp-details').evaluate(element => { element.open = false; });
  const threeQueries = studies['repeated-measures'] + '\nsolve power p100 using study at 100\nsolve effect e100 using study at 100';
  await compile(threeQueries);
  assert.equal(await page.locator('#lean-query option').count(), 3);
  await page.locator('#lean-query').selectOption('1');
  assert.match(await page.locator('#lean-solution').textContent(), /Power =/);
  await page.locator('#lean-query').selectOption('2');
  assert.match(await page.locator('#lean-solution').textContent(), /Detectable/);
  await compile(examples.composition);
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
  await compile(studies['repeated-measures']);
  await page.setViewportSize({ width: 390, height: 844 });
  assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth), false);
  if (process.env.SCREENSHOTS) await page.screenshot({ path: resolve(process.env.SCREENSHOTS, 'native-lean-mobile.png'), fullPage: true });
  await page.reload();
  assert.equal(await page.locator('#lean-source').inputValue(), studies['repeated-measures']);
  assert.equal(await page.locator('#lean-compiled').isVisible(), false);
  // A hosted outage must not tell users to install a local compiler.
  await context.route('**/api/lean/health', route => route.fulfill({ status: 503, contentType: 'application/json',
    body: JSON.stringify({ ready: false, execution: 'hosted', compiler: 'lean-equation-compiler-v1' }) }));
  await page.reload();
  await page.waitForFunction(() => document.getElementById('lean-connection-status').textContent.includes('temporarily unavailable'));
  assert.equal(await page.locator('#lean-setup').isVisible(), false);
  assert.equal(await page.locator('#lean-compile').isDisabled(), true);
  await context.unroute('**/api/lean/health');
  // Static hosting must never pretend that a kernel has run.
  await context.route('**/api/lean/health', route => route.fulfill({ status: 404, body: '' }));
  await page.reload();
  await page.waitForFunction(() => !document.getElementById('lean-setup').hidden);
  assert.equal(await page.locator('#lean-compile').isDisabled(), true);
  assert.match(await page.locator('#lean-connection-status').textContent(), /unavailable/);
  assert.deepEqual(errors, []); assert.deepEqual(violations, []);
  console.log('Native Lean editor: actual kernel results, user rules, fallback, multiple queries, errors, cancellation, mobile, persistence, and static-host state passed.');
} finally {
  await context.close(); await browser.close(); server?.kill('SIGINT');
}
