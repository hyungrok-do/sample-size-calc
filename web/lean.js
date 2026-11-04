import { BrowserLeanCompiler } from './lean-wasm-client.js?v=dev';
import { examples as proofs } from './lean-proof-examples.js?v=dev';
import { examples as calculations } from './lean-calculation-examples.js?v=dev';
import { examples as derivations } from './lean-derivation-examples.js?v=dev';
import { renderCalculations } from './lean-calculation-results.js?v=dev';
const examples = { ...proofs, ...calculations, ...derivations };

const $ = id => document.getElementById(`lean-${id}`);
const source = $('source'), storageKey = 'ssc-lean-proofs-v3';
const browserCompiler = new BrowserLeanCompiler(text => { $('connection-status').textContent = text; });
let ready = false, busy = false, compiled = null, generation = 0, connectionGeneration = 0;

function refresh() {
  $('compile').disabled = busy || !ready;
  $('stop').hidden = !busy;
  for (const id of ['source', 'example', 'load', 'import', 'query', 'reconnect', 'save', 'proof-save', 'result-save']) $(id).disabled = busy;
}
function invalidate() {
  ++generation; compiled = null; $('compiled').hidden = true; $('error').hidden = true;
  $('status').textContent = 'Source changed. Compile to check the derivation.';
  try { localStorage.setItem(storageKey, source.value); } catch {}
}
function download(name, body, type = 'text/plain') {
  const url = URL.createObjectURL(new Blob([body], { type })), a = document.createElement('a');
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
}
async function checkConnection() {
  const check = ++connectionGeneration;
  ready = false; refresh(); $('setup').hidden = true;
  try {
    await browserCompiler.initialize();
    if (check !== connectionGeneration) return;
    ready = true;
    $('connection-status').textContent = 'Lean compiler ready in this browser. Source, proof checking and evaluation stay here.';
  } catch (error) {
    if (check !== connectionGeneration) return;
    $('connection-status').textContent = error.message; $('setup').hidden = false;
  }
  refresh();
}
function current() { return compiled?.results[Number($('query').value)]; }
function render() {
  const declaration = current();
  renderCalculations($('calculations'), compiled.calculations || []);
  $('declarations').hidden = !declaration;
  $('evaluation').textContent = [compiled.output,
    ...compiled.diagnostics.filter(d => d.severity === 'information').map(d => d.data)].filter(Boolean).join('\n\n')
      || (compiled.calculations?.length ? 'Structured calculation output is shown above.' : 'No evaluation commands in this module.');
  $('diagnostics').textContent = compiled.diagnostics.filter(d => d.severity !== 'information')
    .map(d => `Line ${d.pos?.line || 1}: ${d.data}`).join('\n\n') || 'No warnings or errors.';
  $('compiled').hidden = false;
  if (!declaration) return;
  const theorem = declaration.kind === 'theorem';
  $('proof-status').textContent = theorem ? 'Theorem checked by the Lean kernel'
    : declaration.kind === 'unchecked definition' ? 'Executable definition · outside the proof guarantee'
    : 'Declaration checked by Lean';
  $('solution').textContent = `${declaration.name} :\n${declaration.statement}`;
  $('scope').textContent = theorem
    ? 'The theorem establishes exactly this proposition, conditional on its parameters and hypotheses.'
    : 'This declaration defines a type or value. A numerical claim needs its own theorem.';
  $('parameters').hidden = declaration.parameters.length === 0;
  $('hypotheses').textContent = declaration.parameters.map(p => `${p.hypothesis ? 'Hypothesis' : 'Parameter'} ${p.name} : ${p.type}`).join('\n');
  $('derivation').replaceChildren();
  for (const step of declaration.steps) {
    const li = document.createElement('li'), label = document.createElement('strong'), tactic = document.createElement('pre');
    label.textContent = `Source line ${step.line}`; tactic.textContent = step.tactic;
    const before = document.createElement('pre'), after = document.createElement('pre');
    before.textContent = `Before\n${step.before}`; after.textContent = `After\n${step.after}`;
    li.append(label, tactic, before, after); $('derivation').append(li);
  }
  if (!declaration.steps.length) {
    const li = document.createElement('li');
    li.textContent = theorem ? 'No source tactic trace is available. Inspect the checked proof term below.' : 'Inspect the definition and its dependencies below.';
    $('derivation').append(li);
  }
  if (compiled.traceTruncated) {
    const li = document.createElement('li');
    li.textContent = 'A source command reached the 256-step trace limit. Its displayed trace may be incomplete. The full theorems were checked; download the source to inspect the remaining steps.';
    $('derivation').append(li);
  }
  $('theorem-type').textContent = `${declaration.kind} ${declaration.name} : ${declaration.statement}`;
  $('proof-term').textContent = declaration.body || '(No defining expression for this declaration.)';
  $('axioms').textContent = `Foundational axioms: ${declaration.axioms.join(', ') || 'none'}. No unproved user axioms or sorry dependencies.`;
  $('dependencies').textContent = declaration.dependencies.join('\n') || 'No constant dependencies.';
}

source.value = examples['bernoulli-observations'];
try { const saved = localStorage.getItem(storageKey); if (saved && new TextEncoder().encode(saved).length <= 131072) source.value = saved; } catch {}
const requested = new URLSearchParams(location.search).get('example');
if (Object.hasOwn(examples, requested)) { source.value = examples[requested]; $('example').value = requested; }
source.addEventListener('input', invalidate);
$('load').addEventListener('click', () => { source.value = examples[$('example').value]; invalidate(); });
$('query').addEventListener('change', render);
$('reconnect').addEventListener('click', () => { browserCompiler.stop(); checkConnection(); });
$('save').addEventListener('click', () => download('source.lean', source.value));
$('proof-save').addEventListener('click', () => { if (compiled) download('checked-source.lean', compiled.module); });
$('result-save').addEventListener('click', () => { if (compiled) download('lean-report.json', JSON.stringify(compiled, null, 2), 'application/json'); });
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
  $('status').textContent = 'Elaborating Lean, checking proofs and evaluating source commands…';
  try {
    const value = await browserCompiler.compile(source.value);
    if (version !== generation) return;
    compiled = value; $('query').replaceChildren();
    value.results.forEach((declaration, index) => {
      const option = document.createElement('option'); option.value = index;
      option.textContent = `${declaration.kind} ${declaration.name}`; $('query').append(option);
    });
    const theorem = value.results.findLastIndex(d => d.kind === 'theorem');
    if (theorem >= 0) $('query').value = String(theorem);
    render(); $('status').textContent = 'Lean compilation complete.';
  } catch (error) {
    if (version === generation) {
      $('error').textContent = error.message; $('error').hidden = false;
      $('status').textContent = 'Not certified. Resolve the Lean errors or remaining proof goals.';
    }
  } finally {
    if (version === generation) { busy = false; refresh(); }
  }
});
$('stop').addEventListener('click', () => {
  ++generation; browserCompiler.stop(); ready = false; busy = false; refresh();
  $('status').textContent = 'Stopped.'; checkConnection();
});
window.addEventListener('pagehide', () => browserCompiler.stop());
checkConnection();
if ('serviceWorker' in navigator) navigator.serviceWorker.register('sw.js').catch(() => {});
