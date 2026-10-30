// The power-analysis panel: what to detect (comparison, outcome, measure), the design, and the results, in one form.
// Its target is kept as an estimand object (estimand-model.js), which the C++ engine validates and analyzes.
import { DEFAULT_POWER, parsePowerSettings, powerCompatibility } from './power-model.js?v=dev';
import { MEASURES, availableMeasures, validateEstimand, describeEstimand, copyEstimand } from './estimand-model.js?v=dev';
import { Engine } from './engine.js?v=dev';

const POPULATION = 'Population represented by the data generating process';
const SHORT = { mean_difference: 'Mean difference', risk_difference: 'Risk difference', risk_ratio: 'Risk ratio', odds_ratio: 'Odds ratio',
  survival_difference: 'Survival difference', rmst_difference: 'RMST difference', coefficient: 'Coefficient', exponentiated_coefficient: 'exp(coefficient)' };

function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (value == null || value === false) continue;
    if (key.startsWith('on')) node.addEventListener(key.slice(2), value);
    else node.setAttribute(key, value === true ? '' : value);
  }
  node.append(...children.flat()); return node;
}
function field(title, input, help = '') {
  return el('div', { class: 'dag-field' }, el('label', { for: input.id }, title), input,
    ...(help ? [el('p', { class: 'small muted' }, help)] : []));
}
function options(input, list, value) {
  input.replaceChildren(...list.map(([key, label]) => el('option', { value: key }, label)));
  input.value = list.some(([key]) => key === value) ? value : list[0]?.[0] ?? '';
  return input;
}
const select = (id, list, value) => options(el('select', { id }), list, value);
const numeric = (id, value, attrs = {}) => el('input', { id, type: 'number', step: 'any', value: value ?? '', ...attrs });
const fmt = (value) => value == null || !Number.isFinite(value) ? '—' : Number(value.toPrecision(5)).toString();
const pct = (value) => `${(100 * value).toFixed(1)}%`;
function download(value, mime, filename) {
  const url = URL.createObjectURL(new Blob([value], { type: mime })), link = el('a', { href: url, download: filename });
  document.body.append(link); link.click(); link.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function csv(result) {
  const head = ['n', 'repetitions', 'power', 'mc_lower', 'mc_upper', 'mc_se', 'analysis_failures', 'mean_estimate', 'null_rejection', 'null_failures'];
  const rows = result.results.map((row) => [row.n, row.alternative.repetitions, row.alternative.probability,
    ...row.alternative.interval, row.alternative.mcSE, row.alternative.failures, row.alternative.meanEstimate,
    row.null?.probability ?? '', row.null?.failures ?? ''].join(','));
  return [head.join(','), ...rows].join('\n') + '\n';
}
const hasParents = (graph, id) => graph.edges.some((e) => e.target === id);
const hasChildren = (graph, id) => graph.edges.some((e) => e.source === id);
function below(graph, id) {
  const found = new Set([id]), queue = [id];
  for (let i = 0; i < queue.length; i++) for (const edge of graph.edges)
    if (edge.source === queue[i] && !found.has(edge.target)) { found.add(edge.target); queue.push(edge.target); }
  found.delete(id); return found;
}
const horizonFor = (node) => node.censoring?.mode === 'administrative' ? node.censoring.time : node.scale;
const levelsFor = (node) => node.kind === 'survival' ? [1, 2] : [0, 1];
// comparison levels a variable of this type can take (and two different ones)
const levelsFit = (node, a, b) => Number.isFinite(a) && Number.isFinite(b) && a !== b && [a, b].every((v) =>
  node.kind === 'binary' ? v === 0 || v === 1 : node.kind === 'count' ? Number.isInteger(v) && v >= 0 : node.kind === 'survival' ? v > 0 : true);

// A ready-to-run target for a graph: a binary variable that affects others, compared on a variable it affects.
export function defaultTarget(graph) {
  if (graph.nodes.length < 2) return null;
  const binary = graph.nodes.filter((n) => n.kind === 'binary');
  const exposure = binary.find((n) => !hasParents(graph, n.id) && hasChildren(graph, n.id))
    ?? binary.find((n) => hasChildren(graph, n.id)) ?? binary[0] ?? graph.nodes[0];
  const affected = below(graph, exposure.id), last = [...graph.nodes].reverse();
  const outcome = last.find((n) => affected.has(n.id) && !hasChildren(graph, n.id)) ?? last.find((n) => affected.has(n.id))
    ?? last.find((n) => n.id !== exposure.id);
  const measures = availableMeasures(graph, 'marginal', outcome.id);
  const measure = outcome.kind === 'binary' && measures.includes('risk_difference') ? 'risk_difference' : measures[0];
  const [reference, active] = levelsFor(exposure);
  const target = { version: 1, kind: 'marginal', population: POPULATION, definition: '', exposure: exposure.id, outcome: outcome.id,
    measure, reference, active, horizon: MEASURES[measure]?.horizon ? horizonFor(outcome) : null,
    nullValue: MEASURES[measure]?.nullValue ?? 0, alternative: null };
  try { validateEstimand(graph, target); return target; } catch { return null; }
}

export function mountPower(root, { getGraph, onTarget, onSettings }) {
  let worker = null, generation = 0, result = null, ready = false, assign = () => false;
  const $ = (id) => root.querySelector(`#${id}`);
  const idle = 'Each simulated study draws data from the graph, fits the analysis above and applies its test.';
  function stop() {
    generation++; worker?.dispose(); worker = null;
    if ($('power-run')) $('power-run').disabled = !ready;
    if ($('power-stop')) $('power-stop').hidden = true;
  }
  function invalidate(message = idle) {
    stop(); result = null;
    $('power-results')?.replaceChildren(el('p', { class: 'muted small' }, message));
    for (const id of ['power-json', 'power-csv']) if ($(id)) $(id).disabled = true;
    if ($('power-progress')) $('power-progress').textContent = '';
  }
  function display(answer) {
    const header = ['Total N', 'Power', '95% MC interval', 'Failed', ...(answer.nullScenario ? ['Null rejection'] : [])];
    const rows = answer.results.map((row) => el('tr', { class: row.n === answer.firstTestedN ? 'reached' : null },
      ...[row.n, row.alternative.valid ? pct(row.alternative.probability) : 'No valid analyses',
        `${pct(row.alternative.interval[0])}–${pct(row.alternative.interval[1])}`,
        `${row.alternative.failures} / ${row.alternative.repetitions}`,
        ...(row.null ? [row.null.valid ? pct(row.null.probability) : 'No valid analyses'] : [])]
        .map((value) => el('td', {}, String(value)))));
    const first = answer.results.find((row) => row.n === answer.firstTestedN);
    const nodes = [el('p', { class: 'dag-power-recommendation' }, first
      ? `N = ${first.n} is the first tested size with power ≥ ${pct(answer.power.targetPower)} (estimated ${pct(first.alternative.probability)}). Sizes between the tested values were not simulated.`
      : `No tested sample size reached ${pct(answer.power.targetPower)} power. Try larger sizes.`),
    el('div', { class: 'table-wrap', tabindex: 0, 'aria-label': 'Simulated power by sample size' },
      el('table', { id: 'power-table' }, el('thead', {}, el('tr', {}, header.map((title) => el('th', { scope: 'col' }, title)))), el('tbody', {}, rows)))];
    const truth = answer.alternative.truth, facts = [];
    facts.push(el('li', { id: 'power-truth' }, truth.exact ? `True effect in the graph: ${fmt(truth.estimate)} (model coefficient).`
      : `True effect in the graph: ${fmt(truth.estimate)} (Monte Carlo SE ${fmt(truth.mcSE)}).`));
    const parameter = answer.alternative.parameter;
    if (parameter) facts.push(el('li', { id: 'power-calibration' }, `${parameter.edge.replace('>', ' → ')} set to ${fmt(parameter.used)} (drawn ${fmt(parameter.drawn)}) to give the planning effect ${fmt(answer.alternative.target)}; the drawn graph is unchanged.`));
    facts.push(el('li', {}, `${answer.power.repetitions} studies per size · ${answer.power.side} test at α = ${answer.power.alpha} · intervals show simulation error only. ${answer.denominator}`));
    facts.push(el('li', { id: 'power-engine' }, `C++/WebAssembly · ${answer.engine.threads} thread${answer.engine.threads === 1 ? '' : 's'} · ${answer.engine.rng}`));
    nodes.push(el('ul', { class: 'dag-power-facts small muted' }, facts));
    const messages = [...answer.warnings];
    if (answer.results.some((row) => row.alternative.failures || row.null?.failures)) messages.push('Some studies could not be generated or analyzed. They count as non-rejections; check the failures before reading power.');
    if (!answer.results.some((row) => row.alternative.valid > 0)) messages.push('No analysis succeeded, so these results do not establish power for this design.');
    for (const message of messages) nodes.push(el('p', { class: 'error' }, message));
    const diagnostics = [];
    for (const row of answer.results) for (const name of ['alternative', 'null']) if (row[name]) {
      const item = row[name];
      diagnostics.push(el('li', {}, `N=${row.n}, ${name}: ${item.valid} valid; mean estimate ${fmt(item.meanEstimate)}; SD ${fmt(item.estimateSD)}; mean SE ${fmt(item.meanSE)}.`
        + (item.failures ? ` Failures: ${Object.entries(item.failureReasons).map(([reason, count]) => `${reason.replaceAll('_', ' ')} ${count}`).join(', ')}.` : '')));
    }
    if (answer.nullScenario) diagnostics.push(el('li', {}, `${answer.nullScenario.exactNull ? 'Null scenario' : 'Reference scenario'}: ${answer.nullScenario.construction} Its effect: ${fmt(answer.nullScenario.truth.estimate)}.`));
    nodes.push(el('details', {}, el('summary', {}, 'Diagnostics'), el('ul', {}, diagnostics)));
    $('power-results').replaceChildren(...nodes); $('power-json').disabled = false; $('power-csv').disabled = false;
  }

  function refresh(saved, savedSettings = null) {
    stop();
    const graph = getGraph();
    let settings;
    try { settings = parsePowerSettings(savedSettings ?? DEFAULT_POWER); } catch { settings = parsePowerSettings(); }
    let start = null;
    try { if (saved && saved.kind !== 'custom') start = validateEstimand(graph, saved); } catch { /* replaced below */ }
    start ??= defaultTarget(graph);
    if (start && start !== saved) onTarget(copyEstimand(start), 'auto');
    if (!start) {
      ready = false; assign = () => false;
      root.replaceChildren(el('div', { class: 'dag-card-heading' }, el('h3', {}, 'Power analysis')),
        el('p', { id: 'power-method', class: 'muted' }, 'Add at least two variables and an arrow to compare groups and calculate power.'));
      return;
    }
    const nodeList = graph.nodes.map((n) => [n.id, `${n.id} · ${n.label}`]);
    const kind = select('estimand-kind', [['marginal', 'Group difference'], ['coefficient', 'Regression coefficient']], start.kind);
    const exposure = el('select', { id: 'estimand-exposure' }), outcome = el('select', { id: 'estimand-outcome' });
    const measure = el('select', { id: 'estimand-measure' });
    const horizon = numeric('estimand-horizon', start.horizon, { min: 0 });
    const alternative = numeric('estimand-alternative', start.alternative, { placeholder: 'As drawn' });
    const active = numeric('estimand-active', start.active), reference = numeric('estimand-reference', start.reference);
    const nullValue = numeric('estimand-null', start.nullValue);
    const sizes = el('input', { id: 'power-sizes', type: 'text', value: settings.sampleSizes.join(', '), inputmode: 'numeric' });
    const repetitions = numeric('power-repetitions', settings.repetitions, { min: 10, max: 10000, step: 1 });
    const alpha = numeric('power-alpha', settings.alpha, { min: 0, max: 1 });
    const side = select('power-side', [['two-sided', 'Two-sided'], ['greater', 'One-sided, greater'], ['less', 'One-sided, less']], settings.side);
    const targetPower = numeric('power-target', settings.targetPower, { min: 0, max: 1 });
    const seed = numeric('power-seed', settings.seed, { min: 0, max: 4294967295, step: 1 });
    const population = numeric('power-population', settings.populationSize, { min: 1000, max: 50000, step: 1 });
    const edge = el('select', { id: 'power-edge' });
    const lower = numeric('power-lower', settings.calibrationLower), upper = numeric('power-upper', settings.calibrationUpper);
    const nullCheck = el('input', { id: 'power-null-check', type: 'checkbox', checked: settings.checkNull });
    const horizonField = field('Time horizon', horizon), levelsFields = [field('Active level', active), field('Reference level', reference)];
    const sizesField = field('Total sample sizes', sizes);
    sizesField.classList.add('dag-wide');
    const calibration = el('div', { class: 'dag-grid' }, field('Arrow set to the effect to detect', edge), field('Lowest value', lower), field('Highest value', upper));
    const formula = el('p', { id: 'power-estimand', class: 'dag-power-formula' }), method = el('p', { id: 'power-method' });
    const form = el('form', { id: 'power-form', novalidate: true },
      el('p', { class: 'dag-group' }, 'What to detect'),
      el('div', { class: 'dag-grid' }, field('Effect', kind), field('Compare', exposure), field('Outcome', outcome), field('Measure', measure),
        horizonField, field('Effect to detect', alternative)),
      el('p', { class: 'dag-group' }, 'Design'),
      el('div', { class: 'dag-grid' }, sizesField, field('Studies per size', repetitions), field('Significance level α', alpha),
        field('Test', side), field('Target power', targetPower)),
      el('details', { id: 'power-more' }, el('summary', {}, 'More settings'),
        el('div', { class: 'dag-grid' }, ...levelsFields, field('Null value', nullValue), field('Random seed', seed),
          field('Population draws for the true effect', population)),
        calibration,
        el('label', { class: 'dag-power-check', for: 'power-null-check' }, nullCheck, 'Also simulate the null scenario (type I error check)')),
      el('div', { class: 'dag-power-summary' }, method, formula),
      el('div', { class: 'dag-actions' }, el('button', { id: 'power-run', type: 'submit', class: 'primary' }, 'Calculate power'),
        el('button', { id: 'power-stop', type: 'button', hidden: true }, 'Stop'),
        el('button', { id: 'power-json', type: 'button', disabled: true }, 'Download report'),
        el('button', { id: 'power-csv', type: 'button', disabled: true }, 'Download CSV'),
        el('span', { id: 'power-progress', role: 'status', 'aria-live': 'polite', class: 'small muted' })));
    root.replaceChildren(el('div', { class: 'dag-card-heading' }, el('h3', {}, 'Power analysis'),
      el('span', { class: 'small muted' }, 'Repeated simulated studies from the graph above')),
    el('div', { class: 'dag-power-layout' }, form, el('div', { id: 'power-results', class: 'dag-power-results', 'aria-live': 'polite' })));

    const num = (input, optional = false) => input.value.trim() === '' ? (optional ? null : NaN) : Number(input.value);
    const node = (id) => graph.nodes.find((n) => n.id === id);
    function choices(keep = {}) {
      if (kind.value === 'coefficient') {
        options(outcome, nodeList.filter(([id]) => hasParents(graph, id)), keep.outcome ?? outcome.value);
        options(exposure, graph.edges.filter((e) => e.target === outcome.value).map((e) => [e.source, `${e.source} · ${node(e.source).label}`]), keep.exposure ?? exposure.value);
      } else {
        options(exposure, nodeList, keep.exposure ?? exposure.value);
        options(outcome, nodeList, keep.outcome ?? outcome.value);
      }
      // for a 0/1 outcome the mean difference is the risk difference: offer it once
      const measures = availableMeasures(graph, kind.value, outcome.value)
        .filter((key) => !(key === 'mean_difference' && node(outcome.value)?.kind === 'binary'));
      options(measure, measures.map((key) => [key, SHORT[key] ?? MEASURES[key].label]), keep.measure ?? measure.value);
      exposure.previousElementSibling.textContent = kind.value === 'coefficient' ? 'Coefficient of' : 'Compare';
      for (const item of levelsFields) item.hidden = kind.value !== 'marginal';
      horizonField.hidden = !MEASURES[measure.value]?.horizon;
      if (!horizonField.hidden && horizon.value.trim() === '' && node(outcome.value)) horizon.value = horizonFor(node(outcome.value));
      if (kind.value === 'marginal' && node(exposure.value) && !levelsFit(node(exposure.value), num(reference), num(active)))
        [reference.value, active.value] = levelsFor(node(exposure.value));
    }
    function measureChanged() { nullValue.value = MEASURES[measure.value]?.nullValue ?? 0; alternative.value = ''; }
    const target = () => ({ version: 1, kind: kind.value, population: POPULATION, definition: '', exposure: exposure.value,
      outcome: outcome.value, measure: measure.value, reference: num(reference), active: num(active),
      horizon: MEASURES[measure.value]?.horizon ? num(horizon) : null, nullValue: num(nullValue), alternative: num(alternative, true) });
    const readSettings = () => parsePowerSettings({ sampleSizes: sizes.value, repetitions: num(repetitions), alpha: num(alpha), side: side.value,
      seed: num(seed), targetPower: num(targetPower), populationSize: num(population), calibrationEdge: edge.value,
      calibrationLower: num(lower), calibrationUpper: num(upper), checkNull: nullCheck.checked });
    // Reads the form, keeps what is valid, and says what the run would do or why it cannot run.  mode: 'auto' (the
    // form as built), 'edit' (typing) or 'commit' (a finished change, which undo can take back).
    function sync(mode) {
      const value = target();
      let problem = null, plan = null;
      try { validateEstimand(graph, value); onTarget(copyEstimand(value), mode); plan = powerCompatibility(graph, value); }
      catch (error) { problem = error.message; }
      const edges = plan?.supported ? plan.calibrationEdges : [];
      const wanted = edge.value || settings.calibrationEdge;
      options(edge, edges.map((e) => [`${e.source}>${e.target}`, `${e.source} → ${e.target} (drawn ${fmt(e.coefficient)})`]), wanted);
      edge.disabled = !edges.length || value.kind === 'coefficient';
      lower.closest('.dag-field').hidden = upper.closest('.dag-field').hidden = value.kind === 'coefficient';
      calibration.hidden = value.alternative === null && value.nullValue === (MEASURES[value.measure]?.nullValue ?? 0);
      try { settings = readSettings(); onSettings(settings); } catch (error) { problem ??= error.message; }
      if (!problem && plan && !plan.supported) problem = plan.reason;
      ready = !problem;
      method.className = problem ? 'error' : '';
      method.textContent = problem ?? `Analysis: ${plan.label}.`;
      formula.textContent = problem ? '' : `${describeEstimand(graph, value).formula}; H₀: θ = ${value.nullValue}.`;
      formula.hidden = Boolean(problem);
      if (!worker) $('power-run').disabled = !ready;
    }
    // a role set from the graph: the other role swaps in if it was this variable; a coefficient target becomes a
    // group difference when the variable cannot take part in it
    assign = (role, id) => {
      const other = role === 'exposure' ? 'outcome' : 'exposure', current = { exposure: exposure.value, outcome: outcome.value };
      if (kind.value === 'coefficient' && (role === 'outcome' ? !hasParents(graph, id)
        : !graph.edges.some((e) => e.source === id && e.target === outcome.value))) kind.value = 'marginal';
      const was = measure.value;
      choices({ [role]: id, [other]: current[other] === id ? current[role] : current[other] });
      if (measure.value !== was) measureChanged();
      invalidate('Settings changed. Calculate again to update the results.'); sync('commit');
      return true;
    };
    kind.addEventListener('change', () => { choices(); measureChanged(); });
    outcome.addEventListener('change', () => { const was = measure.value; choices(); if (measure.value !== was) measureChanged(); });
    exposure.addEventListener('change', () => choices());
    measure.addEventListener('change', () => { choices(); measureChanged(); });
    form.addEventListener('input', () => { invalidate('Settings changed. Calculate again to update the results.'); sync('edit'); });
    form.addEventListener('change', () => { invalidate('Settings changed. Calculate again to update the results.'); sync('commit'); });
    form.addEventListener('submit', async (event) => {
      event.preventDefault(); invalidate(); sync('commit');
      if (!ready) return;
      try {
        const graph = getGraph(), estimand = copyEstimand(target()), configuration = readSettings();
        const token = generation, instance = new Engine();
        worker = instance; $('power-run').disabled = true; $('power-stop').hidden = false;
        $('power-progress').textContent = 'Starting the engine…';
        const answer = await instance.dag({ op: 'power', graph, estimand, settings: configuration }, (info) => {
          if (token !== generation) return;
          $('power-progress').textContent = info.message ?? `${info.completed} / ${info.total} studies · N = ${info.n} · ${info.scenario}`;
        });
        if (token !== generation) return;
        instance.dispose(); worker = null; $('power-run').disabled = false; $('power-stop').hidden = true;
        result = answer; display(result); $('power-progress').textContent = '';
      } catch (error) { if (!error.stopped) { invalidate(); $('power-results').replaceChildren(el('p', { class: 'error', role: 'alert' }, error.message)); } }
    });
    $('power-stop').addEventListener('click', () => invalidate('Power calculation stopped.'));
    $('power-json').addEventListener('click', () => { if (result) download(JSON.stringify(result, null, 2) + '\n', 'application/json', 'power-report.json'); });
    $('power-csv').addEventListener('click', () => { if (result) download(csv(result), 'text/csv', 'power.csv'); });
    choices({ exposure: start.exposure, outcome: start.outcome, measure: start.measure });
    invalidate(); sync('auto');
  }
  return { refresh, stop, assign: (role, id) => assign(role, id) };
}
