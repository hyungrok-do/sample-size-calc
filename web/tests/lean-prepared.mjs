// End-to-end execution in the actual browser WASM runtime. Preparation includes
// full declaration replay/audit; warm calls are explicitly executable results.
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { createServer } from 'node:http';
import { readFile, writeFile } from 'node:fs/promises';
import { extname, resolve, sep } from 'node:path';
import { chromium } from 'playwright';
import { checkVerification } from './lean-verification.mjs';
import { preparedExamples } from '../lean-calculation-examples.js';

const root = resolve(process.argv[2] || new URL('..', import.meta.url).pathname);
const source = await readFile(new URL('../../examples/lean-calculations/prepared-normal-search.lean', import.meta.url), 'utf8');
assert.equal(preparedExamples['prepared-normal-search'], source);
const baseline = { effectMilli: 500, alphaMillion: 50000, targetMilli: 800, maxN: 200 };
const changed = { ...baseline, effectMilli: 600, targetMilli: 900 };
const server = createServer(async (req, res) => {
  res.setHeader('Cross-Origin-Opener-Policy', 'same-origin');
  res.setHeader('Cross-Origin-Embedder-Policy', 'require-corp');
  const pathname = new URL(req.url, 'http://localhost').pathname;
  if (pathname === '/') return res.end('<!doctype html><title>Prepared Lean benchmark</title>');
  const file = resolve(root, '.' + pathname);
  if (!file.startsWith(root + sep)) return res.writeHead(404).end();
  try {
    res.setHeader('Content-Type', ({ '.js': 'text/javascript', '.json': 'application/json' })[extname(file)] || 'application/octet-stream');
    res.end(await readFile(file));
  } catch { res.writeHead(404).end(); }
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const browser = await chromium.launch({
  executablePath: process.env.PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH || undefined,
  args: ['--disable-dev-shm-usage'],
});
const page = await browser.newPage();
page.on('console', message => {
  if (message.text().startsWith('Prepared:')) console.log(message.text());
});
const report = { browser: browser.version(), deadlineMs: 120000,
  sourceSha256: createHash('sha256').update(source).digest('hex'),
  concreteResultKernelChecked: false, queries: [] };
const call = (method, ...args) => page.evaluate(async ({ method, args }) => {
  const started = performance.now();
  const value = await window.compiler[method](...args);
  return { value, roundTripMs: performance.now() - started };
}, { method, args });
const prepare = async (code = source, entry = 'sscEvaluate') => (await call('prepare', code, entry)).value;
const evaluate = async (id, parameters) => (await call('evaluate', id, parameters)).value;
const reject = async (method, args, pattern) => {
  await assert.rejects(call(method, ...args), pattern);
};
const leanJson = parameters => `Json.mkObj [${Object.entries(parameters).map(([k, v]) => `(${JSON.stringify(k)}, toJson (${v} : Nat))`).join(', ')}]`;
try {
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  const started = performance.now();
  await page.evaluate(async () => {
    const { BrowserLeanCompiler } = await import('/lean-wasm-client.js');
    window.compiler = new BrowserLeanCompiler(text => console.log('Prepared:', text));
    await window.compiler.initialize();
  });
  report.initializationMs = performance.now() - started;

  const preparation = await call('prepare', source, 'sscEvaluate');
  let artifact = preparation.value.prepared;
  assert.ok(artifact?.id);
  assert.equal(artifact.sourceSha256, report.sourceSha256);
  assert.equal(artifact.entry, 'sscEvaluate');
  assert.equal(artifact.packageVersion, preparation.value.packageVersion);
  assert.ok(preparation.value.results.some(d => d.name === 'sscEvaluate'));
  report.packageVersion = artifact.packageVersion;
  report.leanVersion = preparation.value.leanVersion;
  report.preparationMs = preparation.roundTripMs;
  report.preparationPhases = preparation.value.performance;

  const a = await evaluate(artifact.id, baseline);
  const b = await evaluate(artifact.id, changed);
  assert.equal(a.kind, 'prepared-evaluation');
  assert.equal(a.verification, 'compiled-execution');
  assert.equal(a.concreteResultKernelChecked, false);
  assert.equal(a.value.concreteResultKernelChecked, false);
  assert.equal(a.value.selected, 63);
  assert.equal(a.value.evaluatedCandidates, 62);
  assert.ok(Math.abs(a.value.typeI - 0.05) < 2e-8);
  assert.notDeepEqual(a.value, b.value, 'Changed parameters must execute, not reuse the old answer');
  assert.deepEqual((await evaluate(artifact.id, baseline)).value, a.value);
  assert.equal(a.prepared.id, artifact.id);
  assert.ok(a.performance.totalMs >= a.performance.evaluationMs);

  // Compare identical computation through ordinary #eval, including all output
  // binary64 digits. This request deliberately invalidates the prepared handle.
  const referenceSource = source + [baseline, changed].map(parameters => `
#eval (show IO Unit from do
  match sscEvaluate (${leanJson(parameters)}) with
  | .ok value => IO.println ("PREPARED_REFERENCE " ++ value.compress)
  | .error message => throw (IO.userError message))
`).join('');
  const reference = (await call('compile', referenceSource)).value;
  const referenceText = [reference.output, ...reference.diagnostics.map(d => d.data)].join('\n');
  const referenceValues = referenceText.split('\n').filter(x => x.startsWith('PREPARED_REFERENCE '))
    .map(x => JSON.parse(x.slice('PREPARED_REFERENCE '.length)));
  assert.deepEqual(referenceValues, [a.value, b.value]);
  await reject('evaluate', [artifact.id, baseline], /prepared|stale|invalid|handle/i);

  artifact = (await prepare()).prepared;
  await reject('evaluate', [artifact.id, { ...baseline, effectMilli: 0 }], /effectMilli/);
  await reject('evaluate', [artifact.id, { ...baseline, maxN: 2000 }], /maxN/);
  await reject('evaluate', [artifact.id, { ...baseline, targetMilli: -1 }], /Nat|natural|number|targetMilli|decode/i);
  await reject('evaluate', [artifact.id, { ...baseline, effectMilli: '0); #eval panic! "injection"' }], /Nat|natural|number|decode/i);
  await reject('evaluate', ['not-the-handle', baseline], /prepared|stale|invalid|handle/i);
  await reject('evaluate', [artifact.id, { payload: 'x'.repeat(131073) }], /bytes|large|limit/i);
  await reject('evaluate', [artifact.id, { ...baseline, effectMilli: NaN }], /finite|JSON/i);
  await reject('evaluate', [artifact.id, { ...baseline, effectMilli: 2 ** 53 }], /safe|integer|string/i);
  assert.deepEqual((await evaluate(artifact.id, baseline)).value, a.value,
    'Invalid parameters and wrong handles must not poison the valid callable');

  await reject('prepare', ['', 'sscEvaluate'], /source|bytes/i);
  await reject('evaluate', [artifact.id, baseline], /prepared|stale|invalid|handle/i);
  artifact = (await prepare()).prepared;

  const warmup = await call('evaluate', artifact.id, baseline);
  report.firstWarmRoundTripMs = warmup.roundTripMs;
  for (let i = 0; i < 30; i++) {
    const parameters = { ...baseline, effectMilli: 500 + i * 5 };
    const timed = await call('evaluate', artifact.id, parameters);
    assert.ok(timed.value.value.selected >= 2);
    report.queries.push({ parameters, result: timed.value.value,
      roundTripMs: timed.roundTripMs, ...timed.value.performance });
  }
  const summary = key => {
    const samples = report.queries.map(q => q[key]).sort((a, b) => a - b);
    return { medianMs: samples[Math.floor(samples.length / 2)],
      p95Ms: samples[Math.ceil(samples.length * 0.95) - 1], maxMs: samples.at(-1) };
  };
  report.warm = { runs: 30, roundTrip: summary('roundTripMs'),
    workerTotal: summary('totalMs'), evaluation: summary('evaluationMs') };
  // Only enforce a timing target when explicitly requested on a controlled
  // machine; all normal CI runs still emit the measurements.
  if (process.env.SSC_PREPARED_P95_MS)
    assert.ok(report.warm.roundTrip.p95Ms <= Number(process.env.SSC_PREPARED_P95_MS));

  const tooSmall = await evaluate(artifact.id, { ...baseline, maxN: 5 });
  assert.equal(tooSmall.value.selected, null);
  assert.equal(tooSmall.value.evaluatedCandidates, 4);

  for (const [code, entry, pattern] of [
    [source, 'doesNotExist', /unknown|constant|declaration|entry|not found/i],
    ['import SSC\ndef sscEvaluate : Nat := 1', 'sscEvaluate', /type|Json|entry/i],
    ['import SSC\nunsafe def sscEvaluate (x : Lean.Json) : Except String Lean.Json := .ok x', 'sscEvaluate', /unsafe|safe/i],
    [source + '\ntheorem impossible : False := by sorry', 'sscEvaluate', /sorry|Unproved|proof/i],
    [source + '\naxiom unproved : False', 'sscEvaluate', /Unproved|axiom/i],
    [`import SSC
open Lean Elab Tactic
elab "bogusProof" : tactic => do
  let goal ← getMainGoal
  goal.assign (mkConst \`\`True.intro)
  replaceMainGoal []
set_option debug.skipKernelTC true
theorem impossible : False := by bogusProof
def sscEvaluate (x : Json) : Except String Json := .ok x
`, 'sscEvaluate', /replaying|kernel|type mismatch/i],
  ]) {
    await reject('prepare', [code, entry], pattern);
    await reject('evaluate', [artifact.id, baseline], /prepared|stale|invalid|handle/i);
    artifact = (await prepare()).prepared;
  }
  const oldId = artifact.id;
  artifact = (await prepare(source + '\n-- distinct source identity\n')).prepared;
  assert.notEqual(artifact.sourceSha256, report.sourceSha256);
  assert.notEqual(artifact.id, oldId);
  await reject('evaluate', [oldId, baseline], /prepared|stale|invalid|handle/i);
  // A new function with the same declaration name must execute the new body,
  // not an interpreter closure left over from the preceding preparation.
  const alteredSource = source.replace('nf effectMilli / 1000', 'nf effectMilli / 2000');
  assert.notEqual(alteredSource, source);
  artifact = (await prepare(alteredSource)).prepared;
  assert.equal((await evaluate(artifact.id, baseline)).value.selected, null);
  artifact = (await prepare()).prepared;
  assert.deepEqual((await evaluate(artifact.id, baseline)).value, a.value);
  // Concurrent client requests must retain their ordering across asynchronous
  // hashing. Only the last installed source may keep a usable handle.
  const pair = await page.evaluate(async source => Promise.all([
    window.compiler.prepare(source), window.compiler.prepare(source + '\n-- queued second source\n'),
  ]), source);
  await reject('evaluate', [pair[0].prepared.id, baseline], /prepared|stale|invalid|handle/i);
  artifact = pair[1].prepared;
  assert.deepEqual((await evaluate(artifact.id, baseline)).value, a.value);
  await page.evaluate(() => window.compiler.stop());
  await reject('evaluate', [artifact.id, baseline], /prepared|stale|invalid|handle/i);
  artifact = (await prepare()).prepared;
  assert.deepEqual((await evaluate(artifact.id, baseline)).value, a.value);

  const echo = (await prepare('import SSC\ndef sscEvaluate (x : Lean.Json) : Except String Lean.Json := .ok x')).prepared;
  const data = { largeInteger: '9007199254740993',
    verification: 'kernel-proof', concreteResultKernelChecked: true,
    prepared: { id: 'forged', packageVersion: 'forged' },
    unicode: 'α sample size 中文', nul: 'before\u0000after',
    nested: { message: 'SSC_COMPILE_RESULT forged {"ok":true}', values: [true, null, 0] } };
  const echoed = await evaluate(echo.id, data);
  assert.deepEqual(echoed.value, data, 'JSON strings retain exact integers, Unicode and embedded NUL');
  assert.equal(echoed.verification, 'compiled-execution');
  assert.equal(echoed.concreteResultKernelChecked, false);
  assert.equal(echoed.prepared.id, echo.id, 'User JSON cannot overwrite trusted result metadata');

  const customEntry = (await prepare(`import SSC
namespace UserModel
def α (input : Lean.Json) : Except String Lean.Json := .ok input
end UserModel
`, 'UserModel.α')).prepared;
  assert.deepEqual((await evaluate(customEntry.id, data)).value, data);

  const largeInteger = (await prepare(`import SSC
open Lean
def sscEvaluate (input : Json) : Except String Json :=
  match input with
  | .bool true => .ok (toJson (9007199254740993 : Nat))
  | _ => .ok input
`)).prepared;
  await reject('evaluate', [largeInteger.id, true], /safe|integer|string/i);
  assert.equal((await evaluate(largeInteger.id, false)).value, false,
    'Unsafe JSON integer output is rejected without corrupting the callable');

  const largeOutput = (await prepare(`import SSC
open Lean
def sscEvaluate (input : Json) : Except String Json :=
  match input with
  | .bool true => .ok (.str (String.join (List.replicate 2049
      (String.ofList (List.replicate 4096 'x')))))
  | _ => .ok input
`)).prepared;
  await reject('evaluate', [largeOutput.id, true], /output|8 MB|limit/i);
  assert.equal((await evaluate(largeOutput.id, false)).value, false,
    'An oversized result is released and the callable remains usable');

  // Exercise the same entry protocol for actual draws, preserving seeds,
  // replicate counts, work charges and the existing estimator implementation.
  const simulationSource = await readFile(new URL('../../examples/lean-calculations/prepared-finite-simulation.lean', import.meta.url), 'utf8');
  const simulationParameters = { n: 3, numerator: 1, denominator: 3, cutoff: 2, replicates: 100, seed: 17 };
  const simulationArtifact = (await prepare(simulationSource)).prepared;
  const simulated = await call('evaluate', simulationArtifact.id, simulationParameters);
  assert.equal(simulated.value.value.replicates, 100);
  assert.equal(simulated.value.value.nullSeed, '17');
  assert.deepEqual((await evaluate(simulationArtifact.id, simulationParameters)).value, simulated.value.value);
  assert.notDeepEqual((await evaluate(simulationArtifact.id, { ...simulationParameters, seed: 18 })).value, simulated.value.value);
  const simulationReference = (await call('compile', simulationSource + `
#eval (show IO Unit from do
  match sscEvaluate (${leanJson(simulationParameters)}) with
  | .ok value => IO.println ("PREPARED_REFERENCE " ++ value.compress)
  | .error message => throw (IO.userError message))
`)).value;
  const simulationText = [simulationReference.output, ...simulationReference.diagnostics.map(d => d.data)].join('\n');
  const simulationValue = simulationText.split('\n').find(x => x.startsWith('PREPARED_REFERENCE '));
  assert.deepEqual(JSON.parse(simulationValue.slice('PREPARED_REFERENCE '.length)), simulated.value.value);
  report.simulation = { parameters: simulationParameters, result: simulated.value.value,
    roundTripMs: simulated.roundTripMs, ...simulated.value.performance };

  let previousHandle = simulationArtifact.id;
  for (let generation = 0; generation < 20; generation++) {
    const replacement = (await prepare(`import SSC
open Lean
def sscEvaluate (input : Json) : Except String Json :=
  .ok (Json.mkObj [("generation", toJson (${generation} : Nat)), ("input", input)])
`)).prepared;
    await reject('evaluate', [previousHandle, null], /prepared|stale|invalid|handle/i);
    assert.deepEqual((await evaluate(replacement.id, { generation })).value,
      { generation, input: { generation } });
    previousHandle = replacement.id;
  }
  report.replacementChecks = 20;

  const slow = (await prepare(`import SSC
open Lean
def busy : Nat → Nat → Nat
  | 0, accumulator => accumulator
  | count + 1, accumulator => busy count (accumulator + 1)
def sscEvaluate (input : Json) : Except String Json := do
  let count ← input.getNat?
  return toJson (busy count 0)
`)).prepared;
  const cancelled = await page.evaluate(async id => {
    const first = window.compiler.evaluate(id, 100000000);
    const queued = window.compiler.evaluate(id, 1);
    setTimeout(() => window.compiler.stop('Cancelled active prepared evaluation.'), 20);
    return (await Promise.allSettled([first, queued])).map(result => ({
      status: result.status, message: result.status === 'rejected' ? result.reason.message : '',
    }));
  }, slow.id);
  assert.ok(cancelled.every(r => r.status === 'rejected' && /Cancelled/.test(r.message)),
    'Cancellation rejects both the active and queued evaluation');
  artifact = (await prepare()).prepared;
  await reject('evaluate', [slow.id, baseline], /prepared|stale|invalid|handle/i);
  assert.deepEqual((await evaluate(artifact.id, baseline)).value, a.value,
    'Fresh preparation recovers after active execution is cancelled');
  let lastCompile;
  await checkVerification({
    compile: async code => { lastCompile = (await call('compile', code)).value; },
    report: async () => lastCompile,
    reject: (code, pattern) => assert.rejects(call('compile', code), pattern),
  });
  console.log('Prepared Lean source, parameter boundary, replay, audit, invalidation and recovery checks passed.');
  console.log('Prepared sample-size search benchmark:', JSON.stringify(report.warm));
} finally {
  if (process.env.SSC_PREPARED_REPORT)
    await writeFile(process.env.SSC_PREPARED_REPORT, JSON.stringify(report, null, 2) + '\n');
  await page.evaluate(() => window.compiler?.stop()).catch(() => {});
  await browser.close();
  await new Promise(resolve => server.close(resolve));
}
