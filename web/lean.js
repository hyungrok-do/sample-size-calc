import { examples } from './lean-native-examples.js?v=dev';

const $ = id => document.getElementById(`lean-${id}`);
const source = $('source'), storageKey = 'ssc-native-lean-v1';
let ready = false, busy = false, compiled = null, generation = 0, activeJob = null, controller = null;

function refresh() {
  $('compile').disabled = busy || !ready;
  $('stop').hidden = !busy;
  for (const id of ['source', 'example', 'load', 'import', 'query', 'reconnect', 'save', 'proof-save', 'cpp-save', 'result-save']) $(id).disabled = busy;
}
function invalidate() {
  ++generation; compiled = null; $('compiled').hidden = true; $('error').hidden = true;
  $('status').textContent = 'Source changed. Compile to update the derivation.';
  try { localStorage.setItem(storageKey, source.value); } catch {}
}
function download(name, body, type = 'text/plain') {
  const url = URL.createObjectURL(new Blob([body], { type })), a = document.createElement('a');
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
}
async function checkConnection() {
  ready = false; refresh();
  $('connection-status').textContent = 'Checking the native compiler…';
  try {
    const response = await fetch('api/lean/health', { cache: 'no-store', signal: AbortSignal.timeout(5000) });
    if (response.ok) {
      const health = await response.json();
      ready = health.ready === true && health.compiler === 'lean-equation-compiler-v1';
    }
  } catch {}
  $('connection-status').textContent = ready ? 'Native Lean compiler connected. Compilation runs on this computer.' : 'Native Lean compiler unavailable. The source editor is ready.';
  $('setup').hidden = ready; refresh();
}
function current() { return compiled?.results[Number($('query').value)]; }
function render() {
  const equation = current(); if (!equation) return;
  $('proof-status').textContent = equation.status === 'solved' ? 'Solution derived · Lean kernel checked' : 'Explicit solution unresolved · partial equivalence checked';
  $('solution').textContent = equation.solution;
  $('scope').textContent = equation.status === 'solved'
    ? 'This describes the solution set under the parameters and assumptions in the theorem below.'
    : 'Bounded proof search did not isolate the unknown. This does not show that an analytic solution is impossible.';
  $('derivation').replaceChildren();
  const visibleSteps = equation.steps.filter(step => step.before !== step.after);
  for (const step of visibleSteps) {
    const li = document.createElement('li'), label = document.createElement('strong'), expression = document.createElement('pre');
    label.textContent = step.rule.replace(/^DesignCompiler\.Equation\./, ''); label.title = step.rule;
    expression.textContent = `${step.before}\n↔\n${step.after}`;
    li.append(label, expression);
    if (step.premises.length) {
      const premises = document.createElement('p'); premises.textContent = `Proved premises: ${step.premises.join(', ')}`; li.append(premises);
    }
    $('derivation').append(li);
  }
  if (!visibleSteps.length) {
    const li = document.createElement('li');
    li.textContent = equation.status === 'solved' ? 'The equation is already explicit.' : 'No applicable transformation established an explicit solution.';
    $('derivation').append(li);
  }
  $('theorem-type').textContent = `${equation.name} :\n${equation.theoremType}`;
  $('proof-term').textContent = equation.proof;
  $('axioms').textContent = `Proof dependencies: ${equation.axioms.join(', ') || 'none'}. No sorry or user axioms. Explored ${equation.explored} search states.`;
  $('diagnostics').textContent = compiled.diagnostics;
  $('cpp-output').hidden = !equation.cpp;
  $('cpp').textContent = equation.cpp;
  $('cpp-status').textContent = equation.lowering.kind === 'analytic_roots'
    ? 'Evaluate the derived roots numerically in standalone C++17.'
    : equation.lowering.kind === 'bracketed_bisection'
      ? 'The elaborated equation can be evaluated in C++. Supply bounds to attempt numerical root finding.'
      : equation.lowering.reason;
  $('cpp-command').textContent = 'c++ -std=c++17 -O2 equation.cpp -o equation\n' +
    (equation.lowering.kind === 'bracketed_bisection' ? './equation LOWER UPPER' : './equation');
  $('compiled').hidden = false;
}

source.value = examples.composition;
try { const saved = localStorage.getItem(storageKey); if (saved && new TextEncoder().encode(saved).length <= 131072) source.value = saved; } catch {}
const requested = new URLSearchParams(location.search).get('example');
if (Object.hasOwn(examples, requested)) { source.value = examples[requested]; $('example').value = requested; }
source.addEventListener('input', invalidate);
$('load').addEventListener('click', () => { source.value = examples[$('example').value]; invalidate(); });
$('query').addEventListener('change', render);
$('reconnect').addEventListener('click', checkConnection);
$('save').addEventListener('click', () => download('equations.lean', source.value));
$('proof-save').addEventListener('click', () => { if (compiled) download('checked-equations.lean', compiled.module); });
$('cpp-save').addEventListener('click', () => { if (current()?.cpp) download('equation.cpp', current().cpp, 'text/x-c++src'); });
$('result-save').addEventListener('click', () => { if (compiled) download('lean-result.json', JSON.stringify(compiled, null, 2), 'application/json'); });
$('import').addEventListener('change', async () => {
  const file = $('import').files[0]; if (!file) return;
  const version = ++generation; busy = true; refresh();
  try {
    if (file.size > 131072) throw Error('Use a Lean file of at most 128 KB.');
    const text = await file.text();
    if (version === generation) { source.value = text; invalidate(); }
  } catch (error) { $('error').textContent = error.message; $('error').hidden = false; }
  finally { busy = false; $('import').value = ''; refresh(); }
});
$('compile').addEventListener('click', async () => {
  const version = ++generation;
  compiled = null; busy = true; $('compiled').hidden = true; $('error').hidden = true; refresh();
  $('status').textContent = 'Lean is elaborating the source, deriving the solution, and checking its proof…';
  const jobId = crypto.randomUUID(); activeJob = jobId; controller = new AbortController();
  try {
    const response = await fetch('api/lean/compile', { method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ source: source.value, jobId }), signal: controller.signal });
    const value = await response.json();
    if (!response.ok || !value.ok) throw Error(value.error || 'Native compilation failed.');
    if (version !== generation) return;
    compiled = value; $('query').replaceChildren();
    value.results.forEach((equation, index) => { const option = document.createElement('option'); option.value = index; option.textContent = equation.name; $('query').append(option); });
    render(); $('status').textContent = 'Lean compilation complete.';
  } catch (error) {
    if (version === generation) { $('error').textContent = error.message; $('error').hidden = false; $('status').textContent = 'Compilation failed.'; }
  } finally {
    if (version === generation) { busy = false; activeJob = null; controller = null; refresh(); }
  }
});
$('stop').addEventListener('click', async () => {
  ++generation; controller?.abort();
  const jobId = activeJob; activeJob = null;
  $('status').textContent = 'Stopping Lean…';
  try { if (jobId) await fetch(`api/lean/compile/${jobId}`, { method: 'DELETE', signal: AbortSignal.timeout(5000) }); }
  catch {}
  busy = false; refresh(); $('status').textContent = 'Stopped.';
});
window.addEventListener('pagehide', () => { if (activeJob) fetch(`api/lean/compile/${activeJob}`, { method: 'DELETE', keepalive: true }).catch(() => {}); });
checkConnection();
if ('serviceWorker' in navigator) navigator.serviceWorker.register('sw.js').catch(() => {});
