import { Engine } from './engine.js?v=dev';
import { examples, programs } from './lean-examples.js?v=dev';

const $ = (id) => document.getElementById(`lean-${id}`);
const source = $('source'), target = $('target'), n = $('n');
let engine, compiled = null, result = null, busy = false, generation = 0;
const storageKey = 'ssc-lean-study-v2';
const titles = { code: 'query', sample_size: 'sample size', power: 'power', effect: 'detectable effect' };
const fmt = (value) => Number(value).toLocaleString(undefined, { maximumSignificantDigits: 7 });
const percent = (value) => `${(100 * value).toFixed(4).replace(/0+$/, '').replace(/\.$/, '')}%`;
function save() { try { localStorage.setItem(storageKey, JSON.stringify({ source: source.value, target: target.value, n: n.value })); } catch {} }
function refresh() {
  const legacy = /:\s*(NormalInput|InfluenceInput)\b/.test(source.value);
  if (!legacy) target.value = 'code';
  else if (target.value === 'code') target.value = 'sample_size';
  target.disabled = busy || !legacy;
  $('n-field').hidden = target.value === 'sample_size' || target.value === 'code';
  $('n-label').textContent = /:\s*NormalInput\b/.test(source.value) ? 'Participants per arm' : 'Independent units';
  $('solve').textContent = target.value === 'code' ? 'Run query' : `Solve ${titles[target.value]}`;
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
  for (const id of ['source', 'target', 'n', 'example', 'load', 'import', 'compile', 'save', 'cpp-save', 'proof-save', 'result-save']) $(id).disabled = value;
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
  $('trace').hidden = true;
  if (r.kind === 'program') { renderProgram(r); return; }
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
function renderProgram(r) {
  const methods = { algebraic_linear: 'Algebraic linear solution', algebraic_quadratic: 'Algebraic quadratic solution', algebraic_constant: 'Constant equation', bracketed_bisection: 'Numerical equation solving', ascending_grid_search: 'Ascending grid search', expression_evaluation: 'Expression evaluation', monte_carlo: 'Monte Carlo simulation' };
  const probability = compiled.explanation?.probabilityQuantity || r.evaluation === 'monte_carlo';
  const display = value => probability ? percent(value) : fmt(value);
  $('basis').textContent = `${methods[r.method] || r.method}${r.evaluation === 'analytic_gaussian_probability' ? ' · derived Gaussian probability' : ''}`;
  $('answer-unit').textContent = '';
  $('result-detail').textContent = r.detail || r.scope || '';
  if (r.status === 'found') { $('answer').textContent = fmt(r.solution); $('answer-unit').textContent = r.variable; }
  else if (r.status === 'evaluated' || r.status === 'estimated') { $('answer').textContent = display(r.value); $('answer-unit').textContent = r.status === 'estimated' ? 'estimated probability' : 'value'; }
  else if (r.status === 'bounded') { $('answer').textContent = `${fmt(r.solutionLower)} – ${fmt(r.solutionUpper)}`; $('answer-unit').textContent = `${r.variable}: minimum lies in this grid range`; }
  else { $('answer').textContent = ({ not_found_in_domain: 'No solution in range', root_not_bracketed: 'Root not bracketed', all_values_satisfy: 'Every value satisfies', numerically_unresolved: 'Numerically unresolved', unresolved: 'Simulation unresolved' })[r.status] || r.status; }
  $('result-table').replaceChildren();
  if (r.variable) tableRow('Solve variable', r.variable);
  if (r.value != null) tableRow(probability ? 'Probability at reported candidate' : 'Expression value', display(r.value));
  if (r.variable) tableRow('Target', display(r.target));
  if (r.residual != null) tableRow('Equation residual', fmt(r.residual));
  if (r.bracketLower != null) tableRow('Numerical bracket', `${fmt(r.bracketLower)} to ${fmt(r.bracketUpper)}`);
  if (r.estimatedSolution != null) tableRow('Minimum from point estimates', fmt(r.estimatedSolution));
  if (r.lower != null) tableRow('Monte Carlo bounds', `${percent(r.lower)} to ${percent(r.upper)}`);
  if (r.confidence != null) tableRow('Simultaneous confidence level', percent(r.confidence));
  if (r.replicatesPerCandidate) tableRow('Replicates per candidate', fmt(r.replicatesPerCandidate));
  $('assumptions').replaceChildren();
  const notes = [compiled.explanation?.reason || '', r.intervalMethod || '', 'Results follow the formulas, random model, and decision rule written in the source. Numerical arithmetic is not a Lean proof certificate.'].filter(Boolean);
  for (const text of notes) { const li = document.createElement('li'); li.textContent = text; $('assumptions').append(li); }
  $('trace').hidden = !r.trace?.length; $('trace-rows').replaceChildren();
  for (const point of r.trace || []) { const tr = document.createElement('tr'); for (const value of [fmt(point.x), percent(point.estimate), percent(point.lower), percent(point.upper), fmt(point.replicates)]) { const td = document.createElement('td'); td.textContent = value; tr.append(td); } $('trace-rows').append(tr); }
}

source.value = programs.normalPower;
try {
  const previous = JSON.parse(localStorage.getItem(storageKey) || localStorage.getItem('ssc-lean-study-v1'));
  if (previous && typeof previous.source === 'string' && previous.source.length <= 131072) {
    source.value = previous.source;
    if (Object.hasOwn(titles, previous.target)) target.value = previous.target;
    if (typeof previous.n === 'string') n.value = previous.n;
  }
} catch {}
const requestedExample = new URLSearchParams(location.search).get('example');
if (Object.hasOwn(programs, requestedExample)) { source.value = programs[requestedExample]; $('example').value = requestedExample; }
refresh();
source.addEventListener('input', invalidate); target.addEventListener('change', invalidate); n.addEventListener('input', invalidate);
$('load').addEventListener('click', () => { source.value = programs[$('example').value] || examples[$('example').value]; n.value = $('example').value === 'influence' ? '126' : '63'; invalidate(); });
$('import').addEventListener('change', async () => {
  const file = $('import').files[0]; if (!file) return;
  const version = ++generation; lock(true);
  try {
    if (file.size > 131072) throw new Error('Use a Lean file of at most 128 KB.');
    const text = await file.text();
    if (version === generation) { source.value = text; lock(false); invalidate(); }
  } catch (err) { if (version === generation) { lock(false); error(err); } }
  finally { $('import').value = ''; }
});
$('save').addEventListener('click', () => download(target.value === 'code' ? 'study.study' : 'study.lean', source.value, 'text/plain'));
$('cpp-save').addEventListener('click', () => { if (compiled) download('study.cpp', compiled.cpp, 'text/x-c++src'); });
$('proof-save').addEventListener('click', () => { if (compiled?.derivation?.leanSource) download('derivation.lean', compiled.derivation.leanSource, 'text/plain'); });
$('result-save').addEventListener('click', () => { if (result) download('study-result.json', JSON.stringify({ source: source.value, plan: compiled.plan, derivation: compiled.derivation, result }, null, 2), 'application/json'); });
$('compile').addEventListener('click', async () => {
  const version = ++generation; compiled = null; result = null; $('result').hidden = true; $('compiled').hidden = true;
  $('error').hidden = true; $('prompt').hidden = false; lock(true); $('status').textContent = 'Parsing, checking types, and compiling the study…';
  try {
    engine ??= new Engine();
    const request = { op: 'compile', source: source.value, solveFor: target.value };
    if (target.value !== 'sample_size' && target.value !== 'code') request.n = Number(n.value);
    const value = await engine.design(request);
    if (version !== generation) return;
    compiled = value; $('cpp').textContent = value.cpp; $('compiled').hidden = false;
    $('status').textContent = 'Derivation ready. Calculation compiled.';
    const d = value.derivation;
    $('derivation-title').textContent = d?.kind === 'exact_algebraic' ? 'Exact algebraic derivation' : d?.kind === 'analytic_probability' ? 'Analytic probability → numerical solution' : d?.kind === 'simulation' ? 'Simulation method and uncertainty' : 'Solution method';
    $('derivation').replaceChildren();
    for (const step of d?.steps || [{ rule: 'Study formula', expression: value.language }]) {
      const li = document.createElement('li'), label = document.createElement('strong'), math = document.createElement('pre');
      label.textContent = step.rule; math.textContent = step.expression; li.append(label, math); $('derivation').append(li);
    }
    $('proof-status').textContent = d?.proofScope || 'The browser does not invoke the Lean proof kernel.';
    $('proof').hidden = !d?.leanSource; $('proof-source').textContent = d?.leanSource || '';
    $('method-steps').replaceChildren();
    const steps = value.explanation?.steps || ['Validate the study parameters.', 'Use the selected normal or influence-function formula.', 'Execute the search or inversion in the C++ runtime.'];
    for (const text of steps) { const li = document.createElement('li'); li.textContent = text; $('method-steps').append(li); }
    $('lowered').hidden = !value.explanation?.loweredExpression;
    $('expression').textContent = value.explanation?.loweredExpression || '';
    $('prompt').textContent = target.value === 'code' ? 'Ready to run the query written in your source.' : `Ready to solve for ${titles[target.value]}.`;
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
