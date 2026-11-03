import { Engine } from './engine.js?v=dev';
import { examples } from './lean-examples.js?v=dev';

const $ = (id) => document.getElementById(`lean-${id}`);
const source = $('source'), target = $('target'), n = $('n');
let engine, compiled = null, result = null, busy = false, generation = 0;
const storageKey = 'ssc-lean-study-v1';
const titles = { sample_size: 'sample size', power: 'power', effect: 'detectable effect' };
const fmt = (value) => Number(value).toLocaleString(undefined, { maximumSignificantDigits: 7 });
const percent = (value) => `${(100 * value).toFixed(4).replace(/0+$/, '').replace(/\.$/, '')}%`;
function save() { try { localStorage.setItem(storageKey, JSON.stringify({ source: source.value, target: target.value, n: n.value })); } catch {} }
function refresh() {
  $('n-field').hidden = target.value === 'sample_size';
  $('n-label').textContent = /:\s*NormalInput\b/.test(source.value) ? 'Participants per arm' : 'Independent units';
  $('solve').textContent = `Solve ${titles[target.value]}`;
  $('solve').disabled = busy || !compiled;
}
function invalidate() {
  ++generation; compiled = null; result = null;
  $('compiled').hidden = true; $('result').hidden = true; $('prompt').hidden = false;
  $('prompt').textContent = 'Compile the study to enable solving.';
  $('error').hidden = true; $('status').textContent = 'Study changed. Compile to update the calculation.';
  refresh(); save();
}
function lock(value) {
  busy = value;
  for (const id of ['source', 'target', 'n', 'example', 'load', 'import', 'compile', 'save', 'cpp-save', 'result-save']) $(id).disabled = value;
  $('stop').hidden = !value; refresh();
}
function error(err) { $('error').textContent = err.message || String(err); $('error').hidden = false; }
function download(name, body, type) {
  const url = URL.createObjectURL(new Blob([body], { type })), a = document.createElement('a');
  a.href = url; a.download = name; a.click(); setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function tableRow(label, value) {
  const tr = document.createElement('tr'), th = document.createElement('th'), td = document.createElement('td');
  th.scope = 'row'; th.textContent = label; td.textContent = value; tr.append(th, td); $('result-table').append(tr);
}
function render(r) {
  $('prompt').hidden = true; $('result').hidden = false;
  const normal = r.method === 'normal_known';
  $('basis').textContent = normal ? 'Normal theory · known variance · two-sided test' : 'Influence-function variance · asymptotic planning approximation';
  $('answer-unit').textContent = '';
  if (r.solveFor === 'power') {
    $('answer').textContent = percent(r.power); $('answer-unit').textContent = 'power';
    $('result-detail').textContent = `At ${fmt(r.n)} ${normal ? 'participants per arm' : r.sizeUnit === 'cluster' ? 'independent clusters' : 'independent participants'}.`;
  } else if (r.solveFor === 'effect') {
    $('answer').textContent = fmt(r.detectableMagnitude); $('answer-unit').textContent = 'detectable effect magnitude';
    $('result-detail').textContent = `Alternative − null, in the outcome’s units. Holds the supplied variance fixed at n = ${fmt(r.n)}${normal ? ' per arm' : ''}.`;
  } else if (r.status === 'found') {
    $('answer').textContent = fmt(r.selectedN);
    $('answer-unit').textContent = normal ? 'participants per arm' : r.sizeUnit === 'cluster' ? 'independent clusters' : 'independent participants';
    $('result-detail').textContent = `${normal ? `${fmt(r.totalParticipants)} participants total. ` : ''}Calculated power: ${percent(r.power)}. Smallest qualifying integer in the specified range under floating-point evaluation.`;
  } else {
    $('answer').textContent = r.status === 'numerically_unresolved' ? 'Numerically unresolved' : 'No solution in range';
    $('result-detail').textContent = r.status === 'numerically_unresolved' ? `Power is too close to the target at n = ${r.search.unresolvedN} to select a size with the comparison guard.` : 'Increase maxN or revise the planning assumptions.';
  }
  $('result-table').replaceChildren();
  tableRow('Significance level', percent(r.alpha));
  if (r.solveFor !== 'power') tableRow('Target power', percent(r.targetPower));
  tableRow('Variance coefficient', fmt(r.varianceCoefficient));
  if (r.search?.predecessorN != null) tableRow(`Power at preceding n = ${r.search.predecessorN}`, percent(r.search.predecessorPower));
  if (r.solveFor === 'effect') { tableRow('Alternative above the null', fmt(r.alternativeAbove)); tableRow('Alternative below the null', fmt(r.alternativeBelow)); }
  $('assumptions').replaceChildren();
  for (const text of r.assumptions || []) { const li = document.createElement('li'); li.textContent = text; $('assumptions').append(li); }
}
source.value = examples.normal;
try {
  const previous = JSON.parse(localStorage.getItem(storageKey));
  if (previous && typeof previous.source === 'string' && previous.source.length <= 65536) {
    source.value = previous.source;
    if (Object.hasOwn(titles, previous.target)) target.value = previous.target;
    if (typeof previous.n === 'string') n.value = previous.n;
  }
} catch {}
refresh();
source.addEventListener('input', invalidate); target.addEventListener('change', invalidate); n.addEventListener('input', invalidate);
$('load').addEventListener('click', () => { source.value = examples[$('example').value]; n.value = $('example').value === 'influence' ? '126' : '63'; invalidate(); });
$('import').addEventListener('change', async () => {
  const file = $('import').files[0]; if (!file) return;
  const version = ++generation; lock(true);
  try {
    if (file.size > 65536) throw new Error('Use a Lean file of at most 64 KB.');
    const text = await file.text();
    if (version === generation) { source.value = text; lock(false); invalidate(); }
  } catch (err) { if (version === generation) { lock(false); error(err); } }
  finally { $('import').value = ''; }
});
$('save').addEventListener('click', () => download('study.lean', source.value, 'text/plain'));
$('cpp-save').addEventListener('click', () => { if (compiled) download('study.cpp', compiled.cpp, 'text/x-c++src'); });
$('result-save').addEventListener('click', () => { if (result) download('study-result.json', JSON.stringify({ source: source.value, plan: compiled.plan, result }, null, 2), 'application/json'); });
$('compile').addEventListener('click', async () => {
  const version = ++generation; compiled = null; result = null; $('result').hidden = true; $('compiled').hidden = true;
  $('error').hidden = true; $('prompt').hidden = false; lock(true); $('status').textContent = 'Starting the C++ compiler for Lean study declarations…';
  try {
    engine ??= new Engine();
    const request = { op: 'compile', source: source.value, solveFor: target.value };
    if (target.value !== 'sample_size') request.n = Number(n.value);
    const value = await engine.design(request);
    if (version !== generation) return;
    compiled = value; $('cpp').textContent = value.cpp; $('compiled').hidden = false;
    $('status').textContent = 'Study syntax and parameters checked. C++ generated.';
    $('prompt').textContent = `Ready to solve for ${titles[target.value]}.`;
  } catch (err) { if (version === generation) { error(err); $('status').textContent = 'Compilation failed.'; } }
  finally { if (version === generation) lock(false); }
});
$('solve').addEventListener('click', async () => {
  if (!compiled) return;
  const version = ++generation; lock(true); $('error').hidden = true; $('status').textContent = 'Solving in C++…';
  try {
    const value = await engine.design(compiled.plan);
    if (version !== generation) return;
    result = value; render(value); $('status').textContent = 'Calculation complete.';
  } catch (err) { if (version === generation) { error(err); $('status').textContent = 'Calculation failed.'; } }
  finally { if (version === generation) lock(false); }
});
$('stop').addEventListener('click', () => { ++generation; engine?.dispose(); engine = null; lock(false); invalidate(); $('status').textContent = 'Stopped.'; });
window.addEventListener('pagehide', () => engine?.dispose());
if ('serviceWorker' in navigator) navigator.serviceWorker.register('sw.js').catch(() => {});
