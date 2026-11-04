import { examples as algebraExamples } from './lean-native-examples.js?v=dev';
import { examples as studyExamples } from './lean-study-examples.js?v=dev';
const examples = { ...algebraExamples, ...studyExamples };

const $ = id => document.getElementById(`lean-${id}`);
const source = $('source'), storageKey = 'ssc-native-lean-statistical-v2';
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
  $('study-detail').hidden = equation.kind !== 'statistical_design';
  if (equation.kind === 'statistical_design') { renderStudy(equation); return; }
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
  $('proof-heading').textContent = 'Proof term (pretty-printed)';
  $('proof-term').textContent = equation.proof;
  $('axioms').textContent = `Proof dependencies: ${equation.axioms.join(', ') || 'none'}. No sorry or user axioms. Explored ${equation.explored} search states.`;
  $('diagnostics').textContent = compiled.diagnostics;
  $('cpp-output').hidden = !equation.cpp;
  $('cpp').textContent = equation.cpp;
  $('cpp-status').textContent = equation.lowering.kind === 'analytic_roots'
    ? 'Evaluate the derived solution set numerically in standalone C++17. Conditional branches are preserved.'
    : equation.lowering.kind === 'bracketed_bisection'
      ? 'The elaborated equation can be evaluated in C++. Supply bounds to attempt numerical root finding.'
      : equation.lowering.reason;
  const parameters = equation.lowering.parameters || [];
  $('cpp-parameters').textContent = parameters.length
    ? `Runtime inputs: ${parameters.map(parameter => parameter.name).join(', ')}. Supply each as NAME=VALUE; order does not matter.`
    : '';
  $('cpp-assumptions').replaceChildren();
  for (const assumption of equation.lowering.assumptions || []) {
    const li = document.createElement('li'); li.textContent = assumption.statement;
    $('cpp-assumptions').append(li);
  }
  $('cpp-guards').hidden = !(equation.lowering.assumptions || []).length;
  const shellQuote = value => "'" + value.replaceAll("'", "'\\''") + "'";
  $('cpp-command').textContent = 'c++ -std=c++17 -O2 equation.cpp -o equation\n' +
    (equation.lowering.kind === 'bracketed_bisection' ? './equation LOWER UPPER' : './equation') +
    parameters.map(parameter => ' ' + shellQuote(`${parameter.name}=VALUE`)).join('');
  $('numerical-scope').textContent = 'Binary64 values and comparisons are approximate, including assumption checks near boundaries. The Lean theorem certifies the symbolic derivation. Numerical root finding needs a continuous, sign-changing bracket and does not prove uniqueness.';
  $('compiled').hidden = false;
}


function numeric(value) {
  if (typeof value === 'number') return value;
  const [a, b] = String(value).split('/').map(Number); return a / (b ?? 1);
}
function display(value, digits = 6) {
  return Number.isFinite(numeric(value)) ? Number(numeric(value).toPrecision(digits)).toLocaleString('en-US', { maximumSignificantDigits: digits }) : String(value);
}
function tableRow(table, label, value) {
  const row = document.createElement('tr'), th = document.createElement('th'), td = document.createElement('td');
  th.scope = 'row'; th.textContent = label; td.textContent = value; row.append(th, td); table.append(row);
}
function renderStudy(plan) {
  const result = plan.evaluation, mode = plan.solveFor;
  $('proof-status').textContent = plan.model.kind === 'gaussian_contrast'
    ? 'Known-covariance Gaussian power · exact model arithmetic checked in Lean'
    : 'Asymptotic Wald power · exact model arithmetic checked in Lean';
  const n = result.selectedN ?? result.n;
  $('solution').textContent = result.status === 'found'
    ? `n = ${n.toLocaleString()} ${plan.sizeUnit}\nTotal participants: ${result.totalParticipants.toLocaleString()}`
    : mode === 'power' ? `Power = ${display(100 * result.power)}% at n = ${n} ${plan.sizeUnit}`
    : mode === 'effect' ? `Detectable |effect| = ${display(result.detectableMagnitude)} at n = ${n} ${plan.sizeUnit}`
    : result.status === 'numerically_unresolved' ? `Power is too close to the target at n = ${result.unresolvedN} to select a size reliably.`
    : 'The target power was not reached in the specified search domain.';
  $('scope').textContent = `${plan.estimand}. ${plan.analysis}.`;
  const numbers = $('study-numbers'); numbers.replaceChildren();
  tableRow(numbers, 'Alternative effect relative to null', display(plan.statistics.effect));
  tableRow(numbers, 'Variance coefficient V', display(plan.statistics.variance));
  tableRow(numbers, 'Type I error / target power', `${display(100*numeric(plan.criterion.alpha))}% / ${display(100*numeric(plan.criterion.targetPower))}%`);
  tableRow(numbers, 'Critical value', display(result.criticalValue));
  for (const group of result.groups || []) tableRow(numbers, group.label, `${group.count.toLocaleString()} units · ${group.participants.toLocaleString()} participants`);
  if (result.power != null) tableRow(numbers, 'Power at this size', `${display(100*result.power)}%`);
  if (result.predecessorPower != null) tableRow(numbers, `Power at n = ${result.predecessorN}`, `${display(100*result.predecessorPower)}%`);
  if (result.dominantTailApproximation != null) tableRow(numbers, 'Dominant-tail approximation (continuous n)', display(result.dominantTailApproximation));
  if (mode === 'effect') tableRow(numbers, 'Power at detectable effect', `${display(100*result.powerAtUpper)}%`);
  const model = $('study-model'); model.replaceChildren();
  const data = document.createElement('pre'); data.className = 'command';
  if (plan.model.kind === 'gaussian_contrast') {
    data.textContent = `Components: ${plan.model.labels.join(' | ')}\nMeans: ${plan.model.means.map(v => display(v)).join(', ')}\nContrast: ${plan.model.contrast.map(v => display(v)).join(', ')}`;
    model.append(data);
    const covariance = document.createElement('table'), header = document.createElement('tr');
    for (const label of ['Covariance', ...plan.model.labels]) { const th = document.createElement('th'); th.textContent = label; header.append(th); }
    covariance.append(header);
    plan.model.covariance.forEach((values, i) => {
      const row = document.createElement('tr'), label = document.createElement('th'); label.textContent = plan.model.labels[i]; row.append(label);
      for (const value of values) { const td = document.createElement('td'); td.textContent = display(value); row.append(td); }
      covariance.append(row);
    }); model.append(covariance);
  } else if (plan.model.kind === 'efficient_ate') {
    const components = document.createElement('table');
    tableRow(components, 'Outcome-noise variance', display(plan.model.residualVariance));
    tableRow(components, 'Effect-heterogeneity variance', display(plan.model.heterogeneityVariance));
    const strata = document.createElement('table'), header = document.createElement('tr');
    for (const label of ['Stratum', 'Weight', 'Treatment probability', 'Mean Y₀', 'Mean Y₁', 'Variance Y₀', 'Variance Y₁']) {
      const th = document.createElement('th'); th.scope = 'col'; th.textContent = label; header.append(th);
    }
    strata.append(header);
    plan.model.strata.forEach((stratum, i) => {
      const row = document.createElement('tr'), label = document.createElement('th');
      label.scope = 'row'; label.textContent = String(i + 1); row.append(label);
      for (const key of ['weight', 'propensity', 'mean0', 'mean1', 'variance0', 'variance1']) {
        const td = document.createElement('td'); td.textContent = display(stratum[key]); row.append(td);
      }
      strata.append(row);
    });
    model.append(components, strata);
  } else {
    const description = document.createElement('p');
    description.textContent = plan.model.varianceSource; model.append(description);
  }
  $('study-assumptions').replaceChildren();
  for (const assumption of plan.assumptions) { const li = document.createElement('li'); li.textContent = assumption; $('study-assumptions').append(li); }
  $('derivation').replaceChildren();
  for (const step of plan.steps) {
    const li = document.createElement('li'), title = document.createElement('strong'), formula = document.createElement('pre'), detail = document.createElement('p');
    title.textContent = step.title; formula.textContent = step.formula; detail.textContent = step.detail;
    li.append(title, formula, detail); $('derivation').append(li);
  }
  if (result.status === 'found') {
    const li = document.createElement('li');
    li.textContent = `Selected n = ${n} after evaluating ${result.evaluated} admissible counts. ${result.predecessorPower != null ? `At n = ${result.predecessorN}, power was ${display(100*result.predecessorPower)}%.` : 'The selected size is the first count in the requested domain.'} Minimality refers to numerical evaluations in that domain.`;
    $('derivation').append(li);
  }
  $('theorem-type').textContent = plan.theoremType;
  $('proof-heading').textContent = 'Verification scope';
  $('proof-term').textContent = plan.verification;
  $('axioms').textContent = `Arithmetic proof dependencies: ${plan.axioms.join(', ') || 'none'}. The normal-tail identity is DesignCompiler.normal_rejection_formula.`;
  $('diagnostics').textContent = compiled.diagnostics;
  $('cpp-output').hidden = false; $('cpp').textContent = plan.cpp;
  $('cpp-status').textContent = 'The generated C++ has already calculated this result. Download it to reproduce the same power calculation and integer search.';
  $('cpp-parameters').textContent = 'The planning values are compiled from this Lean source. Change the study assumptions and compile again to evaluate a different design.';
  $('cpp-assumptions').replaceChildren(); $('cpp-guards').hidden = true;
  $('cpp-command').textContent = 'c++ -std=c++17 -O2 study.cpp -o study\n./study';
  $('numerical-scope').textContent = 'Power and sample-size comparisons use approximate binary64 arithmetic. A comparison within 1e-12 of target power is reported as unresolved. This guard is not a certified error bound; the integer minimum is numerical and restricted to the requested search domain.';
  $('compiled').hidden = false;
}

source.value = examples['repeated-measures'];
try { const saved = localStorage.getItem(storageKey); if (saved && new TextEncoder().encode(saved).length <= 131072) source.value = saved; } catch {}
const requested = new URLSearchParams(location.search).get('example');
if (Object.hasOwn(examples, requested)) { source.value = examples[requested]; $('example').value = requested; }
source.addEventListener('input', invalidate);
$('load').addEventListener('click', () => { source.value = examples[$('example').value]; invalidate(); });
$('query').addEventListener('change', render);
$('reconnect').addEventListener('click', checkConnection);
$('save').addEventListener('click', () => download('study.lean', source.value));
$('proof-save').addEventListener('click', () => { if (compiled) download('checked-study.lean', compiled.module); });
$('cpp-save').addEventListener('click', () => {
  const plan = current();
  if (plan?.cpp) download(plan.kind === 'statistical_design' ? 'study.cpp' : 'equation.cpp', plan.cpp, 'text/x-c++src');
});
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
  $('status').textContent = 'Compiling the Lean specification, checking its derivation, and evaluating the design…';
  const jobId = crypto.randomUUID(); activeJob = jobId; controller = new AbortController();
  try {
    const response = await fetch('api/lean/compile', { method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ source: source.value, jobId }), signal: controller.signal });
    const value = await response.json();
    if (!response.ok || !value.ok) throw Error(value.error || 'Native compilation failed.');
    if (version !== generation) return;
    compiled = value; $('query').replaceChildren();
    value.results.forEach((equation, index) => { const option = document.createElement('option'); option.value = index; option.textContent = equation.title ? `${equation.title} · ${equation.solveFor}` : equation.name; $('query').append(option); });
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
