import { sizeUnit, scopeOf } from './dag-study.js?v=dev';
import { BASES, parseAnalysis } from './dag-functions.js?v=dev';
// The power-analysis panel: what to detect (comparison, outcome, measure), the design, and the results, in one form.
// Its target is kept as an estimand object (estimand-model.js), which the C++ engine validates and analyzes.
import { DEFAULT_POWER, parsePowerSettings, powerCompatibility } from './power-model.js?v=dev';
import { MEASURES, availableMeasures, validateEstimand, describeEstimand, copyEstimand } from './estimand-model.js?v=dev';
import { Engine } from './engine.js?v=dev';
import { questionSummary } from './dag-workflow.js?v=dev';
import { estimateSize, nextSizes, rowChecks, effectText, summaryText } from './dag-results.js?v=dev';

const POPULATION = 'Population represented by the data generating process';
const SHORT = { mean_difference: 'Mean difference', risk_difference: 'Risk difference', risk_ratio: 'Risk ratio', odds_ratio: 'Odds ratio',
  survival_difference: 'Survival difference', rmst_difference: 'RMST difference', coefficient: 'Coefficient', exponentiated_coefficient: 'exp(coefficient)', overall_effect: 'Overall effect', nonlinearity: 'Nonlinearity', interaction: 'Interaction' };

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

// Map engine results to the display; the estimated size, when there is one, is drawn as an estimate.
function powerPlot(answer, estimate) {
  const rows = answer.results.filter((row) => row.alternative.valid > 0);
  if (!rows.length) return [];
  const svg = (tag, attrs, ...children) => {
    const node = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
    node.append(...children); return node;
  };
  const min = Math.min(...rows.map((r) => r.n)), max = Math.max(...rows.map((r) => r.n));
  // doubling sizes sit evenly on a log axis; a narrow range reads better on a linear one
  const log = max / min >= 4, scale = (n) => log ? Math.log(n) : n;
  const x = (n) => min === max ? 350 : 60 + (scale(n) - scale(min)) / (scale(max) - scale(min)) * 580, y = (p) => 225 - p * 190;
  const plot = svg('svg', { id: 'power-plot', viewBox: '0 0 700 285', role: 'img', 'aria-labelledby': 'power-plot-title power-plot-desc' },
    svg('title', { id: 'power-plot-title' }, 'Estimated power at each tested sample size'),
    svg('desc', { id: 'power-plot-desc' }, 'Dots show simulated power, bars show 95% Monte Carlo intervals, and the dashed line shows target power. Exact values are in the table below.'));
  for (const p of [0, .25, .5, .75, 1]) plot.append(svg('line', { x1: 60, y1: y(p), x2: 640, y2: y(p), class: 'power-plot-grid' }),
    svg('text', { x: 48, y: y(p) + 4, 'text-anchor': 'end' }, `${p * 100}%`));
  plot.append(svg('line', { x1: 60, y1: y(answer.power.targetPower), x2: 640, y2: y(answer.power.targetPower), class: 'power-plot-target' }),
    svg('text', { x: 640, y: 18, 'text-anchor': 'end', class: 'power-plot-target-label' }, `Dashed line: target ${pct(answer.power.targetPower)}`));
  if (estimate.status === 'between') plot.append(svg('line', { x1: x(estimate.n), y1: 30, x2: x(estimate.n), y2: 225, class: 'power-plot-estimate' }),
    svg('text', { x: x(estimate.n) + 6, y: 40, class: 'power-plot-estimate-label' }, `≈ ${estimate.n} (estimate)`));
  rows.forEach((row, i) => {
    const px = x(row.n), a = row.alternative, lo = y(a.interval[0]), hi = y(a.interval[1]);
    plot.append(svg('path', { d: `M ${px} ${lo} V ${hi} M ${px - 5} ${lo} H ${px + 5} M ${px - 5} ${hi} H ${px + 5}`, class: 'power-plot-interval' }),
      svg('circle', { cx: px, cy: y(a.probability), r: 5, class: 'power-plot-dot' }, svg('title', {}, `${row.n}: ${pct(a.probability)} power, ${pct(a.interval[0])} to ${pct(a.interval[1])}`)));
    // Dense grids retain their exact sizes in the table and hover text.
    if (rows.length <= 8 || i === 0 || i === rows.length - 1) plot.append(svg('text', { x: px, y: 247, 'text-anchor': 'middle' }, String(row.n)));
  });
  plot.append(svg('text', { x: 350, y: 275, 'text-anchor': 'middle' }, (answer.graph.study ? `Number of ${answer.sizeUnit}` : 'Total independent observations') + (log ? ' (log scale)' : '')));
  return [el('div', { class: 'power-plot-wrap' }, plot), el('p', { class: 'small muted' }, 'Dots: estimated power · bars: 95% Monte Carlo intervals. Only the displayed sample sizes were simulated.')];
}
function download(value, mime, filename) {
  const url = URL.createObjectURL(new Blob([value], { type: mime })), link = el('a', { href: url, download: filename });
  document.body.append(link); link.click(); link.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function csv(result) {
  const head = ['n', 'repetitions', 'power', 'mc_lower', 'mc_upper', 'mc_se', 'analysis_failures', 'mean_estimate', 'bias', 'ci_coverage_valid', 'null_rejection', 'null_failures'];
  const rows = result.results.map((row) => [row.n, row.alternative.repetitions, row.alternative.probability,
    ...row.alternative.interval, row.alternative.mcSE, row.alternative.failures, row.alternative.meanEstimate, row.alternative.bias ?? '', row.alternative.coverage ?? '',
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
// The learners carrying weight, largest first: "gbm 0.62, kernel 0.21, gam 0.17".
const learnerWeights = (learners, model) => learners.names.map((name, i) => [name, learners[model][i]]).filter(([, w]) => w >= .005)
  .sort((a, b) => b[1] - a[1]).map(([name, w]) => `${name} ${w.toFixed(2)}`).join(', ');
export function defaultTarget(graph) {
  if (graph.nodes.length < 2) return null;
  const binary = graph.nodes.filter((n) => n.kind === 'binary' && n.measured !== false && (!graph.study || ['subject', 'cluster'].includes(scopeOf(n))));
  const exposure = binary.find((n) => !hasParents(graph, n.id) && hasChildren(graph, n.id))
    ?? binary.find((n) => hasChildren(graph, n.id)) ?? binary[0] ?? graph.nodes.find((n) => n.measured !== false);
  if (!exposure) return null;
  const affected = below(graph, exposure.id), last = [...graph.nodes].reverse().filter((n) => n.measured !== false);
  const outcome = last.find((n) => affected.has(n.id) && !hasChildren(graph, n.id)) ?? last.find((n) => affected.has(n.id))
    ?? last.find((n) => n.id !== exposure.id);
  if (!outcome) return null;
  const measures = availableMeasures(graph, 'marginal', outcome.id);
  const measure = outcome.kind === 'binary' && measures.includes('risk_difference') ? 'risk_difference' : measures[0];
  const [reference, active] = levelsFor(exposure);
  const target = { version: 1, kind: 'marginal', population: POPULATION, definition: '', exposure: exposure.id, outcome: outcome.id,
    measure, reference, active, horizon: MEASURES[measure]?.horizon ? horizonFor(outcome) : null,
    nullValue: MEASURES[measure]?.nullValue ?? 0, alternative: null };
  try { validateEstimand(graph, target); return target; } catch { return null; }
}

export function mountPower(root, { getGraph, onTarget, onSettings, onRefresh = () => {}, onRun = () => {}, engine = null }) {
  // background requests share the page's warm engine when it offers one
  const borrow = engine?.borrow ?? (() => new Engine()), give = engine?.give ?? ((instance) => instance.dispose());
  let worker = null, generation = 0, result = null, ready = false, assign = () => false;
  // the effect the model implies for the current question, computed in the background without a run
  const hint = { timer: 0, token: 0, key: '', text: '', number: '' };
  function implied(graph, value, show) {
    clearTimeout(hint.timer);
    const name = (SHORT[value.measure] ?? value.measure).toLowerCase(), digits = (x) => Number(x.toPrecision(3)).toString();
    if (value.kind === 'joint') return show('');
    if (value.kind === 'coefficient') {
      const edge = graph.edges.find((e) => e.source === value.exposure && e.target === value.outcome);
      const c = edge && fmt(value.measure === 'exponentiated_coefficient' ? Math.exp(edge.coefficient) : edge.coefficient);
      return show(edge ? `${name} ${c}` : '', edge ? c : '');
    }
    const key = JSON.stringify([graph, { ...value, alternative: null, nullValue: null }]);
    if (key === hint.key) return show(hint.text, hint.number);
    show('', '');
    hint.timer = setTimeout(async () => {
      const token = ++hint.token, instance = borrow();
      try {
        const truth = await instance.dag({ op: 'population', graph, estimand: { ...value, alternative: null }, options: { n: graph.study ? 400 : 5000, seed: 1 } });
        if (token !== hint.token) return;
        const approximate = truth.mcSE > 1e-9 * Math.max(1, Math.abs(truth.estimate));
        hint.key = key; hint.number = Number.isFinite(truth.estimate) ? `${approximate ? '≈ ' : ''}${digits(truth.estimate)}` : '';
        hint.text = hint.number ? `${name} ${hint.number}` : '';
        show(hint.text, hint.number);
      } catch { if (token === hint.token) show('', ''); }
      finally { give(instance); }
    }, 400);
  }
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
  // simulates a new set of sizes from the results: the sizes field takes them and the form runs again
  function simulateSizes(list) {
    const input = $('power-sizes');
    if (!input || !list.length) return;
    input.value = list.join(', ');
    input.dispatchEvent(new Event('change', { bubbles: true }));
    $('power-form').requestSubmit();
  }
  function headline(answer, estimate, doubtful) {
    const unit = answer.graph.study ? answer.sizeUnit : 'in total', target = pct(answer.power.targetPower);
    const truth = answer.alternative.truth, null0 = answer.estimand.nullValue;
    // refining the size of a test that does not hold its level would only refine a misleading number
    const action = (label, list) => list.length && !doubtful ? [el('button', { type: 'button', class: 'dag-power-next', 'data-sizes': list.join(','),
      onclick: () => simulateSizes(list) }, `${label}: ${list.join(', ')}`)] : [];
    const caution = doubtful ? ' — not reliable yet' : '';
    if (answer.estimand.kind !== 'joint' && truth && Number.isFinite(truth.estimate)
      && Math.abs(truth.estimate - null0) <= Math.max(2 * (truth.mcSE ?? 0), 1e-12 * Math.max(1, Math.abs(null0))))
      return [el('p', { class: 'dag-power-answer' }, 'Your model has no effect to detect.'),
        el('p', { class: 'small' }, 'The true effect in the graph equals the null value, so no sample size gives more than the false-positive rate. Change the effect in the Model or the Question step.')];
    if (estimate.status === 'between') {
      const { lower, upper } = estimate;
      return [el('p', { class: 'dag-power-answer' }, `About ${estimate.n} ${unit} for ${target} power${caution}`),
        el('p', { class: 'small' }, `Estimated between the tested sizes ${lower.n} (${pct(lower.alternative.probability)}) and ${upper.n} (${pct(upper.alternative.probability)}); ${estimate.n} itself has not been simulated. ${upper.n} is the smallest tested size that reaches the target.`),
        ...action('Simulate sizes around the estimate', nextSizes(estimate))];
    }
    if (estimate.status === 'below') return [el('p', { class: 'dag-power-answer' }, `${estimate.upper.n} ${unit} already gives ${pct(estimate.upper.alternative.probability)} power${caution}`),
      el('p', { class: 'small' }, doubtful ? 'Read the check below before relying on this size.' : `Every tested size reaches ${target}; a smaller study may be enough.`), ...action('Try smaller sizes', nextSizes(estimate))];
    if (estimate.status === 'above') return [el('p', { class: 'dag-power-answer' }, `Not reached: ${estimate.lower.n} ${unit} gives ${pct(estimate.lower.alternative.probability)} power${caution}`),
      el('p', { class: 'small' }, doubtful ? `No tested size reaches ${target}, and the check below questions even that.` : `No tested size reaches ${target}.`), ...action('Try larger sizes', nextSizes(estimate))];
    return [el('p', { class: 'dag-power-answer' }, 'No analysis succeeded'), el('p', { class: 'small' }, 'See the failures below before changing the sample sizes.')];
  }
  function display(answer) {
    const estimate = estimateSize(answer), checks = new Map(answer.results.map((row) => [row.n, rowChecks(answer, row)]));
    const header = [answer.graph.study ? `Number of ${answer.sizeUnit}` : 'Total N', ...(answer.graph.study ? ['Subjects', 'Measurement rows'] : []), 'Power', '95% MC interval', 'Failed', ...(answer.nullScenario ? ['Null rejection'] : []), 'Checks'];
    const rows = answer.results.map((row) => {
      const found = checks.get(row.n);
      return el('tr', { class: [row.n === answer.firstTestedN ? 'reached' : '', found.length ? 'flagged' : ''].join(' ').trim() || null },
        ...[row.n, ...(answer.graph.study ? [row.subjects, row.measurementRows] : []), row.alternative.valid ? pct(row.alternative.probability) : 'No valid analyses',
          `${pct(row.alternative.interval[0])}–${pct(row.alternative.interval[1])}`,
          `${row.alternative.failures} / ${row.alternative.repetitions}`,
          ...(row.null ? [row.null.valid ? `${pct(row.null.probability)} (${pct(row.null.interval[0])}–${pct(row.null.interval[1])})` : 'No valid analyses'] : [])]
          .map((value) => el('td', {}, String(value))),
        el('td', { class: found.length ? 'dag-check-flag' : 'dag-check-ok' }, found.length ? `⚠ ${found.map((c) => c.short).join(' · ')}` : 'OK'));
    });
    // a size whose test does not hold its level or whose estimates are off does not support its power
    const relevant = estimate.status === 'between' ? [estimate.lower, estimate.upper] : estimate.status === 'below' ? [estimate.upper] : estimate.status === 'above' ? [estimate.lower] : [];
    const doubtful = relevant.flatMap((row) => checks.get(row.n).map((c) => `at ${row.n}, ${c.text}`));
    const nodes = [el('div', { class: `dag-power-recommendation${doubtful.length ? ' doubtful' : ''}`, id: 'power-answer' }, ...headline(answer, estimate, doubtful.length > 0))];
    if (doubtful.length) nodes.push(el('p', { class: 'dag-power-warning', role: 'note' }, `Check before using this size: ${doubtful.join('; ')}. Power from a test that over-rejects or from biased estimates overstates what a study would achieve; choose a different analysis in the Question step or inspect Diagnostics.`));
    nodes.push(...powerPlot(answer, estimate), el('div', { class: 'table-wrap', tabindex: 0, 'aria-label': 'Simulated power by sample size' },
      el('table', { id: 'power-table' }, el('thead', {}, el('tr', {}, header.map((title) => el('th', { scope: 'col' }, title)))), el('tbody', {}, rows))));
    const truth = answer.alternative.truth, facts = [];
    facts.push(el('li', { id: 'power-truth' }, answer.estimand.kind === 'joint' ? `Generating predictor-block coefficients: ${truth.coefficients?.map(fmt).join(', ') || 'none'}. The hypothesis concerns the specified regression basis.`
      : `True effect in the graph: ${effectText(answer)}${truth.exact ? ' (model coefficient)' : ''}.`));
    const parameter = answer.alternative.parameter;
    if (parameter) facts.push(el('li', { id: 'power-calibration' }, `${parameter.edge.replace('>', ' → ')} set to ${fmt(parameter.used)} (drawn ${fmt(parameter.drawn)}) to give the planning effect ${fmt(answer.alternative.target)}; the drawn graph is unchanged.`));
    facts.push(el('li', {}, `${answer.power.repetitions} studies per size · ${answer.power.side} test at α = ${answer.power.alpha} · intervals show simulation error only. ${answer.denominator}`));
    nodes.push(el('ul', { class: 'dag-power-facts small muted' }, facts));
    if (answer.assumptions && !answer.graph.study) nodes.push(el('p', { class: 'dag-power-note' }, answer.assumptions));
    // the analysis's own notes read as notes; anything else the run found, as a warning
    const notes = answer.analysis.notes ?? [];
    for (const note of notes) nodes.push(el('p', { class: 'dag-power-note' }, note));
    const messages = answer.warnings.filter((message) => !notes.includes(message));
    if (answer.results.some((row) => row.alternative.failures || row.null?.failures)) messages.push('Some studies could not be generated or analyzed. They count as non-rejections; check the failures before reading power.');
    if (!answer.results.some((row) => row.alternative.valid > 0)) messages.push('No analysis succeeded, so these results do not establish power for this design.');
    for (const message of messages) nodes.push(el('p', { class: 'error' }, message));
    // a paragraph to paste into a protocol, from the numbers shown
    const text = summaryText(answer, questionSummary(answer.graph, answer.estimand));
    const copied = el('span', { class: 'small muted', role: 'status', 'aria-live': 'polite' });
    nodes.push(el('details', { class: 'dag-power-summary-text' }, el('summary', {}, 'Text for a protocol or grant'),
      el('p', { id: 'power-summary' }, text),
      el('button', { type: 'button', id: 'power-copy', onclick: async () => {
        try { await navigator.clipboard.writeText(text); copied.textContent = 'Copied.'; }
        catch { const range = document.createRange(); range.selectNodeContents($('power-summary')); getSelection().removeAllRanges(); getSelection().addRange(range); copied.textContent = 'Selected · press Ctrl+C or ⌘C to copy.'; }
      } }, 'Copy text'), copied));
    const diagnostics = [];
    for (const row of answer.results) for (const name of ['alternative', 'null']) if (row[name]) {
      const item = row[name];
      diagnostics.push(el('li', {}, `N=${row.n}, ${name}: ${item.valid} valid.` + (answer.estimand.kind === 'joint' ? '' : ` Mean estimate ${fmt(item.meanEstimate)}; bias ${fmt(item.bias)}; SD ${fmt(item.estimateSD)}; mean SE ${fmt(item.meanSE)}; ${pct(1 - answer.power.alpha)} CI coverage ${item.coverage == null ? '—' : pct(item.coverage)} among valid fits.`)
        + (item.failures ? ` Failures: ${Object.entries(item.failureReasons).map(([reason, count]) => `${reason.replaceAll('_', ' ')} ${count}`).join(', ')}.` : '')));
    }
    for (const row of answer.results) for (const name of ['alternative', 'null']) {
      const d = row[name]?.targetedDiagnostics;
      if (d) diagnostics.push(el('li', {}, `N=${row.n}, ${name}, successful targeted fits: largest inverse weight ${fmt(d.maxInverseWeight)}; ${d.propensityTruncated} subject propensities truncated; ${d.initialPredictionsClipped} initial outcome predictions clipped across replicates.`
        + (d.learners ? ` Mean Super Learner weights, outcome: ${learnerWeights(d.learners, 'outcome')}; propensity: ${learnerWeights(d.learners, 'propensity')}.` : '')));
    }
    if (answer.nullScenario) diagnostics.push(el('li', {}, `${answer.nullScenario.exactNull ? 'Null scenario' : 'Reference scenario'}: ${answer.nullScenario.construction} Its effect: ${fmt(answer.nullScenario.truth.estimate)}.`));
    diagnostics.push(el('li', { id: 'power-engine' }, `Computed by the C++/WebAssembly engine · ${answer.engine.threads} thread${answer.engine.threads === 1 ? '' : 's'} · ${answer.engine.rng}.`));
    nodes.push(el('details', {}, el('summary', {}, 'Diagnostics'), el('ul', {}, diagnostics)));
    $('power-results').replaceChildren(...nodes); $('power-json').disabled = false; $('power-csv').disabled = false;
  }

  function refresh(saved, savedSettings = null) {
    stop(); clearTimeout(hint.timer); hint.token++;
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
        el('p', { id: 'power-method', class: 'muted' }, 'Add at least two variables and an arrow to compare groups and calculate power.'),
        el('button', { type: 'button', 'data-goto': 'model' }, 'Go to Model →'));
      onRefresh();
      return;
    }
    const nodeList = graph.nodes.map((n) => [n.id, `${n.id} · ${n.label}`]);
    const kind = select('estimand-kind', graph.study ? [['marginal', 'Trajectory contrast']] : [['marginal', 'Group / dose contrast'], ['coefficient', 'Regression coefficient'], ['joint', 'Joint test (overall, nonlinearity, interaction)']], start.kind);
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
    const nullMode = select('power-null-mode', [['no_effect', 'Remove all exposure effects'], ['average_effect', 'Calibrate the contrast; retain other effects']], settings.nullMode);
    const currentAnalysis = settings.analysis;
    const analysisMethod = select('analysis-method', graph.study ? [['auto', 'Automatic'], ['mmrm', 'MMRM'], ['lmm', 'Linear mixed model'], ['trajectory', 'Trajectory means · robust inference'], ['tmle', 'TMLE · complete outcomes']] : [['auto', 'Automatic'], ['gformula', 'Parametric g-computation'], ['aipw', 'AIPW'], ['tmle', 'TMLE · binary / continuous / count']], currentAnalysis.method);
    if (!analysisMethod.value) analysisMethod.value = 'auto';
    const trajectory = select('estimand-trajectory', [['final', 'Selected visit'], ['average', 'Average across scheduled visits'], ['change', 'Change from first to last visit'], ['slope', 'Difference in linear trend'], ['custom', 'Custom visit weights']], start.repeated?.contrast ?? 'final');
    const visit = select('estimand-visit', (graph.study?.times ?? []).map((t, i) => [String(i), `Time ${t}`]), String(start.repeated?.visit ?? (graph.study?.times.length ?? 1) - 1));
    const weights = el('input', { id: 'estimand-weights', value: (start.repeated?.weights ?? graph.study?.times.map(() => 0) ?? []).join(', ') });
    const visitField = field('Visit to compare', visit), weightsField = field('Weights in visit order', weights);
    const trajectoryPanel = el('div', { class: 'dag-grid dag-trajectory-grid', hidden: !graph.study }, field('Trajectory target', trajectory), visitField, weightsField);
    const contrastFields = () => { visitField.hidden = trajectory.value !== 'final'; weightsField.hidden = trajectory.value !== 'custom'; };
    trajectory.addEventListener('change', contrastFields); contrastFields();
    const timeBasis = select('analysis-time-basis', [['categorical', 'Separate mean at every visit'], ['linear', 'Linear time trend']], currentAnalysis.timeBasis ?? 'categorical');
    const randomSlope = el('input', { id: 'analysis-random-slope', type: 'checkbox', checked: currentAnalysis.randomSlope });
    const analysisBasis = select('analysis-basis', [...BASES, ...(currentAnalysis.basis === 'custom' ? [['custom', 'Custom imported terms']] : [])], currentAnalysis.basis);
    const propensityBasis = select('analysis-propensity', [['same', 'Use the analysis basis'], ...BASES, ...(currentAnalysis.propensity?.basis === 'custom' ? [['custom', 'Custom imported terms']] : [])], currentAnalysis.propensity?.basis ?? 'same');
    const folds = numeric('analysis-folds', currentAnalysis.folds, { min: 1, max: 10, step: 1 });
    const learner = select('analysis-learner', [['glm', 'Parametric: the terms chosen below'], ['superlearner', 'Super Learner: a cross-validated ensemble, no functional form']], currentAnalysis.learner ?? 'glm');
    const propensityBound = numeric('analysis-propensity-bound', currentAnalysis.propensityBound ?? 0, { min: 0, max: .499, step: .01 });
    const outcomeBounds = el('input', { id: 'analysis-outcome-bounds', type: 'text', value: (currentAnalysis.outcomeBounds ?? []).join(', '), placeholder: 'Unbounded, or e.g. 0, 100' });
    const targetedOptions = el('div', { class: 'dag-grid' }, field('Propensity bound (0 = none)', propensityBound, 'Explicitly limit treatment probabilities to [b, 1−b]. This may add bias; inspect diagnostics.'),
      field('Known outcome-summary bounds (TMLE)', outcomeBounds, 'Optional lower, upper bounds for a continuous summary. Binary bounds are automatic. Values outside the bounds fail the fit.'));
    const knots = el('input', { id: 'analysis-knots', type: 'text', value: (currentAnalysis.knots ?? [-1, 0, 1]).join(', ') });
    const propensityKnots = el('input', { id: 'analysis-propensity-knots', type: 'text', value: (currentAnalysis.propensity?.knots ?? [-1, 0, 1]).join(', ') });
    const interactions = el('input', { id: 'analysis-interactions', type: 'checkbox', checked: currentAnalysis.interactions });
    const propensityInteractions = el('input', { id: 'analysis-propensity-interactions', type: 'checkbox', checked: currentAnalysis.propensity?.interactions });
    const overrides = new Map();
    const overrideRows = graph.nodes.filter((n) => hasParents(graph, n.id)).map((n) => {
      const cfg = currentAnalysis.models?.[n.id];
      const choice = select(`analysis-model-${n.id}`, [['default', 'Use common settings'], ...BASES, ...(cfg?.basis === 'custom' ? [['custom', 'Custom imported terms']] : [])], cfg?.basis ?? 'default');
      const modelKnots = el('input', { id: `analysis-model-knots-${n.id}`, type: 'text', value: (cfg?.knots ?? [-1, 0, 1]).join(', ') });
      const modelInteractions = el('input', { id: `analysis-model-interactions-${n.id}`, type: 'checkbox', checked: cfg?.interactions });
      overrides.set(n.id, { choice, knots: modelKnots, interactions: modelInteractions, cfg });
      const k = field('Spline knots', modelKnots), inter = el('label', {}, modelInteractions, 'Include pairwise interactions');
      const visibility = () => { k.hidden = choice.value !== 'spline'; inter.hidden = ['default', 'match', 'custom'].includes(choice.value); };
      choice.addEventListener('change', visibility); visibility();
      return el('div', { class: 'dag-grid' }, field(`Model for ${n.id}`, choice), k, inter);
    });
    const learnerField = field('AIPW / TMLE nuisance models', learner, 'The Super Learner weighs a mean, GLMs, additive splines, boosted trees and kernel ridge regression by 5-fold cross-validation; it uses the variables of the terms below, not their form. Automatic analysis chooses it, with TMLE and 3-fold cross-fitting, when a treatment or outcome model has a relationship of unknown form.');
    const knotsField = field('Analysis spline knots', knots, 'Fixed knots in the original predictor units; includes linear terms.'),
      propensityKnotsField = field('Propensity spline knots', propensityKnots);
    function analysisFields() {
      if (graph.study) { timeBasis.disabled = randomSlope.disabled = analysisMethod.value !== 'lmm'; nuisancePanel.hidden = analysisMethod.value !== 'tmle'; trajectoryAnalysis.hidden = analysisMethod.value === 'tmle'; }
      targetedOptions.hidden = !['aipw', 'tmle'].includes(analysisMethod.value);
      learnerField.hidden = !['aipw', 'tmle'].includes(analysisMethod.value);
      if (learnerField.hidden) learner.value = 'glm';
      knotsField.hidden = analysisBasis.value !== 'spline'; propensityKnotsField.hidden = propensityBasis.value !== 'spline';
      folds.disabled = !['auto', 'aipw', 'tmle'].includes(analysisMethod.value);
      if (folds.disabled) folds.value = 1;
      if (kind.value !== 'marginal') { analysisMethod.value = 'auto'; folds.value = 1; }
      analysisMethod.disabled = kind.value !== 'marginal';
      nullMode.disabled = kind.value !== 'marginal';
    }
    const parseKnots = (input) => input.value.trim().split(/[\s,]+/).filter(Boolean).map(Number);
    const readAnalysis = () => {
      if (graph.study && analysisMethod.value !== 'tmle') return parseAnalysis({ method: analysisMethod.value, basis: 'match', folds: 1, timeBasis: timeBasis.value, randomSlope: randomSlope.checked });
      const config = { ...(analysisBasis.value === 'custom' ? currentAnalysis : {}), method: analysisMethod.value,
        basis: analysisBasis.value, folds: Number(folds.value), interactions: interactions.checked,
        ...(analysisBasis.value === 'spline' ? { knots: parseKnots(knots) } : {}) };
      delete config.models; delete config.propensity; delete config.outcomeBounds;
      if (['aipw', 'tmle'].includes(analysisMethod.value)) config.propensityBound = Number(propensityBound.value);
      if (['aipw', 'tmle'].includes(analysisMethod.value) && learner.value === 'superlearner') config.learner = 'superlearner';
      if (analysisMethod.value === 'tmle' && outcomeBounds.value.trim()) config.outcomeBounds = parseKnots(outcomeBounds);
      if (propensityBasis.value !== 'same') config.propensity = { ...(propensityBasis.value === 'custom' ? currentAnalysis.propensity : {}), basis: propensityBasis.value,
        interactions: propensityInteractions.checked, ...(propensityBasis.value === 'spline' ? { knots: parseKnots(propensityKnots) } : {}) };
      const models = Object.fromEntries([...overrides].filter(([, v]) => v.choice.value !== 'default').map(([id, v]) => [id,
        { ...(v.choice.value === 'custom' ? v.cfg : {}), basis: v.choice.value, interactions: v.interactions.checked,
          ...(v.choice.value === 'spline' ? { knots: parseKnots(v.knots) } : {}) }]));
      if (Object.keys(models).length) config.models = models;
      return parseAnalysis(config);
    };
    const nuisancePanel = el('div', {},
      el('p', { class: 'small muted' }, 'These settings change the fitted analysis, independently of the generating graph. Cross-fitting applies to AIPW and TMLE. For repeated data, adjustments use baseline covariates and folds keep entire subjects or clusters together.'),
      learnerField,
      el('div', { class: 'dag-grid' }, field('Outcome / mediator terms', analysisBasis), field('Cross-fitting folds (1 = none)', folds), knotsField),
      el('label', {}, interactions, 'Include pairwise interactions with linear, quadratic or spline terms'),
      el('div', { class: 'dag-grid' }, field('Propensity terms (AIPW / TMLE)', propensityBasis), propensityKnotsField),
      el('label', {}, propensityInteractions, 'Include pairwise propensity interactions'),
      el('details', {}, el('summary', {}, 'Override individual models'), overrideRows), targetedOptions);
    const trajectoryAnalysis = el('div', { hidden: !graph.study },
      el('p', { class: 'small muted' }, 'MMRM uses an unstructured covariance. LMM uses the time and random-effect choices below. Trajectory means use robust inference at the subject or cluster level.'),
      el('div', { class: 'dag-grid dag-trajectory-grid' }, field('LMM time model', timeBasis)), el('label', {}, randomSlope, 'LMM subject random slope (in addition to intercept)'));
    const analysisPanel = el('details', { id: 'power-analysis-models' }, el('summary', {}, 'Analysis details and assumptions'), nuisancePanel, trajectoryAnalysis);
    for (const input of [analysisMethod, analysisBasis, propensityBasis]) input.addEventListener('change', analysisFields);
    const horizonField = field('Time horizon', horizon), levelsFields = [field('Active level', active), field('Reference level', reference)];
    const sizesField = field(graph.study ? `Number of ${sizeUnit(graph)}` : 'Total sample sizes', sizes, 'Separate the sizes to test with commas. These are totals across both groups.');
    sizesField.classList.add('dag-wide');
    const calibration = el('div', { class: 'dag-grid' }, field('Arrow set to the effect to detect', edge), field('Lowest value', lower), field('Highest value', upper));
    const formula = el('p', { id: 'power-estimand', class: 'dag-power-formula' }), method = el('p', { id: 'power-method' });
    const impliedText = el('span', { id: 'estimand-implied' });
    const adjustNote = el('p', { id: 'analysis-adjust-note', class: 'dag-adjust-note small', hidden: true });
    let impliedValue = '';
    const notes = el('ul', { id: 'power-notes', class: 'dag-power-notes', hidden: true });
    const question = el('p', { id: 'power-question-summary', class: 'dag-question-summary' });
    const review = el('p', { id: 'power-review-question' }), reviewMethod = el('p', { id: 'power-review-method', class: 'small' }),
      reviewEffect = el('p', { id: 'power-review-effect', class: 'small' });
    const heading = (step, id, title, help) => el('div', { class: 'dag-section-heading' }, el('p', { class: 'dag-eyebrow' }, step),
      el('h3', { id, tabindex: -1 }, title), el('p', {}, help));
    const questionPane = el('section', { 'data-dag-step': 'question' },
      heading('03 / Question', 'dag-question-title', 'What do you want to detect?', 'Choose the outcome, comparison and analysis for every simulated study.'), question,
      el('div', { class: 'dag-grid dag-question-grid' }, field('Type of question', kind), field('Compare', exposure), field('Outcome', outcome), field('Measure', measure)),
      el('div', { class: 'dag-grid dag-question-grid' }, ...levelsFields, horizonField,
        el('div', { class: 'dag-field' }, el('label', { for: 'estimand-alternative' }, 'Effect to detect'), alternative,
          el('p', { class: 'small muted' }, 'Leave blank to use the effect implied by your model', impliedText, '.'))),
      trajectoryPanel,
      el('div', { class: 'dag-analysis-choice' }, field('Analysis method', analysisMethod), method), adjustNote, analysisPanel,
      el('details', { class: 'dag-question-details' }, el('summary', {}, 'Formal target and assumptions'),
        field('Null value', nullValue), formula, notes));
    const powerPane = el('section', { 'data-dag-step': 'power' },
      heading('04 / Power', 'dag-power-title', 'How many do you need?', 'Test a range of sample sizes, then compare their power with your target.'),
      el('div', { class: 'dag-run-review' }, el('span', { class: 'dag-eyebrow' }, 'Your comparison'), review, reviewEffect, reviewMethod,
        el('button', { type: 'button', 'data-goto': 'question' }, 'Edit question')),
      el('div', { class: 'dag-grid dag-design-grid' }, sizesField,
        field('Target power', targetPower, '0.80 means an 80% chance of detecting the specified effect.'),
        field('Simulated studies per size', repetitions, 'More studies reduce Monte Carlo error and take longer.'),
        field('Significance level α', alpha), field('Test', side)),
      el('details', { id: 'power-more' }, el('summary', {}, 'Simulation settings and null checks'),
        el('div', { class: 'dag-grid' }, field('Random seed', seed), field('Population draws for the true effect', population)),
        calibration, field('Reference / null construction', nullMode),
        el('label', { class: 'dag-power-check', for: 'power-null-check' }, nullCheck, 'Also simulate the null scenario (type I error check)')),
      ...(graph.study ? [el('p', { class: 'small muted' }, 'Power is conditional on your generating model, dependence and missingness assumptions. Monte Carlo intervals measure simulation error, not uncertainty about those assumptions.')] : []),
      el('div', { class: 'dag-actions' }, el('button', { id: 'power-run', type: 'submit', class: 'primary' }, 'Calculate power'),
        el('button', { id: 'power-stop', type: 'button', hidden: true }, 'Stop'),
        el('button', { id: 'power-json', type: 'button', disabled: true }, 'Download report'),
        el('button', { id: 'power-csv', type: 'button', disabled: true }, 'Download CSV'),
        el('span', { id: 'power-progress', role: 'status', 'aria-live': 'polite', class: 'small muted' })));
    const form = el('form', { id: 'power-form', novalidate: true }, questionPane, powerPane);
    root.replaceChildren(form, el('div', { id: 'power-results', class: 'dag-power-results', 'data-dag-step': 'power', 'aria-live': 'polite' }));

    const num = (input, optional = false) => input.value.trim() === '' ? (optional ? null : NaN) : Number(input.value);
    const node = (id) => graph.nodes.find((n) => n.id === id);
    function choices(keep = {}) {
      if (['coefficient', 'joint'].includes(kind.value)) {
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
      alternative.closest('.dag-field').hidden = kind.value === 'joint';
      nullValue.disabled = side.disabled = kind.value === 'joint';
      if (kind.value === 'joint') { alternative.value = ''; nullValue.value = 0; side.value = 'two-sided'; }
      analysisFields();
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
      horizon: MEASURES[measure.value]?.horizon ? num(horizon) : null, nullValue: num(nullValue), alternative: num(alternative, true),
      ...(graph.study ? { repeated: { contrast: trajectory.value, ...(trajectory.value === 'final' ? { visit: Number(visit.value) } : {}), ...(trajectory.value === 'custom' ? { weights: weights.value.trim().split(/[\s,]+/).filter(Boolean).map(Number) } : {}) } } : {}) });
    const readSettings = () => parsePowerSettings({ sampleSizes: sizes.value, repetitions: num(repetitions), alpha: num(alpha), side: side.value,
      seed: num(seed), targetPower: num(targetPower), populationSize: num(population), calibrationEdge: edge.value,
      calibrationLower: num(lower), calibrationUpper: num(upper), checkNull: nullCheck.checked, nullMode: nullMode.value, analysis: readAnalysis() });
    // Reads the form, keeps what is valid, and says what the run would do or why it cannot run.  mode: 'auto' (the
    // form as built), 'edit' (typing) or 'commit' (a finished change, which undo can take back).
    function sync(mode) {
      const value = target();
      let problem = null, plan = null;
      try { validateEstimand(graph, value); onTarget(copyEstimand(value), mode === 'auto' ? 'auto' : 'edit'); plan = powerCompatibility(graph, value, readAnalysis()); }
      catch (error) { problem = error.message; }
      const edges = plan?.supported ? plan.calibrationEdges : [];
      const wanted = edge.value || settings.calibrationEdge;
      options(edge, edges.map((e) => [`${e.source}>${e.target}`, `${e.source} → ${e.target} (drawn ${fmt(e.coefficient)})`]), wanted);
      edge.disabled = !edges.length || value.kind === 'coefficient';
      lower.closest('.dag-field').hidden = upper.closest('.dag-field').hidden = value.kind === 'coefficient';
      calibration.hidden = value.kind === 'joint' || (value.alternative === null && value.nullValue === (MEASURES[value.measure]?.nullValue ?? 0)
        && !(value.kind === 'marginal' && nullCheck.checked && nullMode.value === 'average_effect'));
      try { settings = readSettings(); onSettings(settings, mode); } catch (error) { problem ??= error.message; }
      if (!problem && plan && !plan.supported) problem = plan.reason;
      ready = !problem;
      method.className = problem ? 'error' : '';
      method.textContent = problem ?? `Analysis: ${plan.label}.`;
      question.textContent = review.textContent = questionSummary(graph, value);
      const renderEffect = () => { reviewEffect.textContent = value.kind === 'joint' ? 'Effect to detect: the predictor block specified in your model.'
        : `Effect to detect: ${value.alternative == null ? `as implied by your model${impliedValue ? ` (${impliedValue})` : ''}` : value.alternative}. Null value: ${value.nullValue}.`; };
      renderEffect();
      if (!problem) implied(graph, value, (text, number) => {
        impliedValue = text; impliedText.textContent = text ? `: ${text}` : '';
        alternative.placeholder = number ? `As drawn (${number})` : 'As drawn'; renderEffect();
      });
      else { impliedText.textContent = ''; alternative.placeholder = 'As drawn'; }
      // an unadjusted test leaves out the outcome's other causes; say so where the analysis is chosen
      const unadjusted = !problem && value.kind === 'marginal' && ['welch', 'risk_score', 'risk_wald', 'risk_ratio', 'odds_ratio', 'km_survival', 'km_rmst'].includes(plan.method);
      const affected = below(graph, value.exposure);
      // only causes a study records can be adjusted for
      const recorded = (id) => graph.nodes.find((n) => n.id === id)?.measured !== false;
      const prognostic = unadjusted ? graph.edges.filter((e) => e.target === value.outcome && e.source !== value.exposure && !affected.has(e.source) && recorded(e.source)).map((e) => e.source) : [];
      adjustNote.hidden = !prognostic.length;
      adjustNote.textContent = prognostic.length ? `This test is unadjusted: ${prognostic.join(', ')} ${prognostic.length === 1 ? 'also causes' : 'also cause'} ${value.outcome} but ${prognostic.length === 1 ? 'is' : 'are'} not used. Adjusting for ${prognostic.length === 1 ? 'it' : 'them'} (choose g-computation or AIPW) usually needs fewer participants.` : '';
      reviewMethod.textContent = method.textContent;
      reviewMethod.className = problem ? 'error' : 'small';
      formula.textContent = problem ? '' : describeEstimand(graph, value).formula + (value.kind === 'joint' ? '' : `. H₀: θ = ${value.nullValue}.`);
      formula.hidden = Boolean(problem);
      notes.replaceChildren(...(problem ? [] : plan.notes).map((text) => el('li', {}, text)));
      notes.hidden = !notes.children.length;
      if (!worker) $('power-run').disabled = !ready;
    }
    // a role set from the graph: the other role swaps in if it was this variable; a coefficient target becomes a
    // group difference when the variable cannot take part in it
    assign = (role, id) => {
      const other = role === 'exposure' ? 'outcome' : 'exposure', current = { exposure: exposure.value, outcome: outcome.value };
      if (['coefficient', 'joint'].includes(kind.value) && (role === 'outcome' ? !hasParents(graph, id)
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
      onRun();
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
    onRefresh();
  }
  // The design changed while the panel is hidden: a run of the old design stops and its results go, and the
  // graph keeps a runnable comparison to show; the form itself is rebuilt when it is next shown.
  function outdated(saved) {
    clearTimeout(hint.timer); hint.token++;
    if (worker || result) invalidate('The design changed. Calculate again to update the results.');
    const graph = getGraph();
    let start = null;
    try { if (saved && saved.kind !== 'custom') start = validateEstimand(graph, saved); } catch { /* replaced below */ }
    start ??= defaultTarget(graph);
    if (start && start !== saved) onTarget(copyEstimand(start), 'auto');
  }
  return { refresh, outdated, stop, assign: (role, id) => assign(role, id) };
}

