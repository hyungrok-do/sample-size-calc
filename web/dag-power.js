import { DEFAULT_POWER, parsePowerSettings, powerCompatibility } from './power-model.js?v=dev';
import { describeEstimand } from './estimand-model.js?v=dev';
import { Engine } from './engine.js?v=dev';

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
function select(id, options, value) {
  const input = el('select', { id }, options.map(([key, label]) => el('option', { value: key }, label)));
  input.value = value; return input;
}
const numeric = (id, value, attrs = {}) => el('input', { id, type: 'number', step: 'any', value, ...attrs });
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

export function mountPower(root, { getDesign, onSettings }) {
  let worker = null, generation = 0, draft = false, result = null, activePlan = null;
  const $ = (id) => root.querySelector(`#${id}`);
  function stop() {
    generation++; worker?.dispose(); worker = null;
    if ($('power-run')) $('power-run').disabled = draft || !activePlan?.supported;
    if ($('power-stop')) $('power-stop').hidden = true;
  }
  function invalidate(message = 'Run the selected analysis across simulated studies to estimate power.') {
    stop(); result = null;
    if ($('power-results')) $('power-results').replaceChildren(el('p', { class: 'muted' }, message));
    if ($('power-json')) $('power-json').disabled = true;
    if ($('power-csv')) $('power-csv').disabled = true;
    if ($('power-progress')) $('power-progress').textContent = '';
  }
  function setDraft(value) {
    draft = value;
    invalidate(value ? 'The estimand has unsaved changes. Save it above before calculating power.' : undefined);
  }
  function display(answer) {
    const header = ['Total N', 'Power', '95% MC interval', 'Failed studies', ...(answer.nullScenario ? ['Null rejection', 'Null failures'] : [])];
    const rows = answer.results.map((row) => el('tr', {},
      ...[row.n, row.alternative.valid ? pct(row.alternative.probability) : 'No valid analyses',
        `${pct(row.alternative.interval[0])}–${pct(row.alternative.interval[1])}`,
        `${row.alternative.failures} / ${row.alternative.repetitions}`,
        ...(row.null ? [row.null.valid ? pct(row.null.probability) : 'No valid analyses', `${row.null.failures} / ${row.null.repetitions}`] : [])]
        .map((value) => el('td', {}, String(value)))));
    const truth = answer.alternative.truth;
    const truthText = truth.exact ? `DGP estimand: ${fmt(truth.estimate)} (exact model coefficient).`
      : `Verified DGP estimand: ${fmt(truth.estimate)}; population Monte Carlo SE ${fmt(truth.mcSE)}${truth.interval ? `, 95% interval [${truth.interval.map(fmt).join(', ')}]` : ''}.`;
    const nodes = [el('h4', {}, 'Power for the saved estimand'),
      el('p', { class: 'dag-power-target' }, describeEstimand(answer.graph, answer.estimand).formula),
      el('p', {}, `${answer.analysis.label}. H₀: θ = ${answer.estimand.nullValue}; ${answer.power.side}; α = ${answer.power.alpha}.`),
      el('p', { id: 'power-truth' }, truthText)];
    nodes.push(el('p', { id: 'power-engine', class: 'small muted' }, `C++/WebAssembly · ${answer.engine.threads} simulation thread${answer.engine.threads === 1 ? '' : 's'} · ${answer.engine.rng}`));
    if (answer.alternative.parameter) {
      const parameter = answer.alternative.parameter;
      nodes.push(el('p', { id: 'power-calibration' }, `Run DGP: ${parameter.edge.replace('>', ' → ')} = ${fmt(parameter.used)} (drawn value ${fmt(parameter.drawn)}), calibrated to planning effect ${fmt(answer.alternative.target)}. The graph is retained as drawn; the calibrated run graph is included in the report.`));
    } else nodes.push(el('p', {}, answer.alternative.target === null ? 'The alternative uses the graph parameters as drawn.' : 'The drawn DGP already matches the planning effect within calibration tolerance.'));
    nodes.push(el('div', { class: 'table-wrap', tabindex: 0, 'aria-label': 'Simulated power by sample size' },
      el('table', { id: 'power-table' }, el('thead', {}, el('tr', {}, header.map((title) => el('th', { scope: 'col' }, title)))), el('tbody', {}, rows))),
      el('p', { class: 'small muted' }, `${answer.power.repetitions} studies per sample size and scenario. Intervals describe Monte Carlo uncertainty, not uncertainty in the assumed DGP. ${answer.denominator}`));
    const validRows = answer.results.filter((row) => row.alternative.valid > 0);
    nodes.push(el('p', { class: 'dag-power-recommendation' }, answer.firstTestedN === null
      ? `No tested sample size reached the ${pct(answer.power.targetPower)} target by its point estimate.`
      : `First tested sample size with estimated power at least ${pct(answer.power.targetPower)}: N = ${answer.firstTestedN}. Review its Monte Carlo interval; sizes between the tested values have not been evaluated.`));
    if (answer.nullScenario) nodes.push(el('details', {}, el('summary', {}, answer.nullScenario.exactNull ? 'Null scenario and type I error check' : 'Calibrated reference scenario'),
      el('p', {}, answer.nullScenario.construction), el('p', {}, `Reference DGP estimand: ${fmt(answer.nullScenario.truth.estimate)}; population MC SE ${fmt(answer.nullScenario.truth.mcSE)}.`)));
    const diagnostics = [];
    for (const row of answer.results) for (const name of ['alternative', 'null']) if (row[name]) {
      const item = row[name];
      diagnostics.push(el('li', {}, `N=${row.n}, ${name}: ${item.valid} valid analyses; mean estimate ${fmt(item.meanEstimate)}; empirical SD ${fmt(item.estimateSD)}; mean SE on the test scale ${fmt(item.meanSE)}.`
        + (item.failures ? ` Failures: ${Object.entries(item.failureReasons).map(([reason, count]) => `${reason.replaceAll('_', ' ')} ${count}`).join(', ')}.` : '')));
    }
    nodes.push(el('details', {}, el('summary', {}, 'Analysis diagnostics'), el('ul', {}, diagnostics)));
    const messages = [...answer.warnings];
    if (answer.results.some((row) => row.alternative.failures || row.null?.failures)) messages.push('Some studies could not be generated or analyzed. They remain in the denominator as non-rejections; inspect the failure counts before interpreting power.');
    if (!validRows.length) messages.push('No alternative-scenario analysis succeeded. These results do not establish usable power for this design.');
    for (const message of messages) nodes.push(el('p', { class: 'error' }, message));
    $('power-results').replaceChildren(...nodes); $('power-json').disabled = false; $('power-csv').disabled = false;
  }
  function refresh(savedSettings = null) {
    draft = false; invalidate();
    const design = getDesign(), plan = powerCompatibility(design.graph, design.estimand); activePlan = plan;
    let settings;
    try { settings = parsePowerSettings(savedSettings ?? DEFAULT_POWER); } catch { settings = parsePowerSettings(); }
    const edges = plan.supported ? plan.calibrationEdges : [];
    const edgeValue = edges.some((edge) => `${edge.source}>${edge.target}` === settings.calibrationEdge) ? settings.calibrationEdge : edges[0] ? `${edges[0].source}>${edges[0].target}` : '';
    const sizes = el('input', { id: 'power-sizes', type: 'text', value: settings.sampleSizes.join(', '), inputmode: 'numeric' });
    const repetitions = numeric('power-repetitions', settings.repetitions, { min: 10, max: 10000, step: 1 });
    const alpha = numeric('power-alpha', settings.alpha, { min: 0, max: 1 });
    const side = select('power-side', [['two-sided', 'Two-sided'], ['greater', 'Greater than null'], ['less', 'Less than null']], settings.side);
    const seed = numeric('power-seed', settings.seed, { min: 0, max: 4294967295, step: 1 });
    const target = numeric('power-target', settings.targetPower, { min: 0, max: 1 });
    const population = numeric('power-population', settings.populationSize, { min: 1000, max: 50000, step: 1 });
    const edge = select('power-edge', edges.map((item) => [`${item.source}>${item.target}`, `${item.source} → ${item.target} (drawn ${fmt(item.coefficient)})`]), edgeValue);
    const lower = numeric('power-lower', settings.calibrationLower), upper = numeric('power-upper', settings.calibrationUpper);
    const nullCheck = el('input', { id: 'power-null-check', type: 'checkbox', checked: settings.checkNull });
    const needsCalibration = design.estimand?.alternative !== null && design.estimand?.alternative !== undefined;
    edge.disabled = !edges.length || design.estimand?.kind === 'coefficient';
    const form = el('form', { id: 'power-form', novalidate: true },
      el('div', { class: 'dag-estimand-grid' }, field('Total sample sizes', sizes, 'Comma-separated sample sizes for repeated simulation.'), field('Simulated studies per sample size', repetitions)),
      el('div', { class: 'dag-estimand-grid' }, field('Significance level α', alpha), field('Rejection rule', side)),
      el('div', { class: 'dag-estimand-grid' }, field('Random seed', seed), field('Target power', target)),
      el('div', { class: 'dag-power-calibration-fields', hidden: !needsCalibration && design.estimand?.nullValue === (design.estimand?.measure?.includes('ratio') ? 1 : 0) },
        field('Arrow to calibrate to the planning effect', edge, 'Only arrows on causal paths to the selected outcome are eligible. A conditional-coefficient target uses its own arrow.'),
        el('div', { class: 'dag-estimand-grid', hidden: design.estimand?.kind === 'coefficient' }, field('Lower coefficient bound', lower), field('Upper coefficient bound', upper))),
      el('details', {}, el('summary', {}, 'Population calibration precision'), field('Population observations for effect calibration', population,
        'Paired intervention simulations evaluate the estimand. An independent stream verifies the calibrated effect.')),
      el('label', { class: 'dag-power-check', for: 'power-null-check' }, nullCheck, 'Also simulate a reference/null scenario'),
      el('div', { class: 'dag-actions' }, el('button', { id: 'power-run', type: 'submit', class: 'primary', disabled: !plan.supported }, 'Calculate power'),
        el('button', { id: 'power-stop', type: 'button', hidden: true }, 'Stop power calculation'),
        el('button', { id: 'power-json', type: 'button', disabled: true }, 'Download report'),
        el('button', { id: 'power-csv', type: 'button', disabled: true }, 'Download power CSV')));
    const method = plan.supported ? el('div', { class: 'dag-power-method' },
      el('p', { class: 'eyebrow' }, 'Analysis selected from the saved estimand'), el('p', { id: 'power-method' }, plan.label),
      el('p', { class: 'small' }, plan.formula),
      el('p', { class: 'small muted' }, needsCalibration ? `The alternative will be calibrated to θ = ${design.estimand.alternative} before simulating studies.` : 'Power uses the effect implied by the graph as drawn.'))
      : el('p', { id: 'power-method', class: 'muted' }, plan.reason);
    root.replaceChildren(el('h3', {}, 'Power for this estimand'), method, form,
      el('p', { id: 'power-progress', role: 'status', 'aria-live': 'polite', class: 'small muted' }),
      el('div', { id: 'power-results' }, el('p', { class: 'muted' }, 'Each replicate will generate data from the DAG, estimate the selected target, and apply its rejection rule.')));
    const readNumber = (input) => input.value.trim() === '' ? NaN : Number(input.value);
    const read = () => parsePowerSettings({ sampleSizes: sizes.value, repetitions: readNumber(repetitions), alpha: readNumber(alpha), side: side.value,
      seed: readNumber(seed), targetPower: readNumber(target), populationSize: readNumber(population), calibrationEdge: edge.value,
      calibrationLower: readNumber(lower), calibrationUpper: readNumber(upper), checkNull: nullCheck.checked });
    const edited = () => { invalidate('Power settings changed. Calculate again to update the result.'); try { onSettings(read()); } catch { /* the run button reports invalid inputs */ } };
    form.addEventListener('input', edited); form.addEventListener('change', edited);
    form.addEventListener('submit', async (event) => {
      event.preventDefault(); invalidate();
      if (draft) return;
      try {
        const current = getDesign(), configuration = read();
        const compatible = powerCompatibility(current.graph, current.estimand);
        if (!compatible.supported) throw new Error(compatible.reason);
        onSettings(configuration);
        const token = generation, instance = new Engine();
        worker = instance; $('power-run').disabled = true; $('power-stop').hidden = false;
        $('power-progress').textContent = 'Starting the C++/WebAssembly engine…';
        const answer = await instance.dag({ op: 'power', graph: current.graph, estimand: current.estimand, settings: configuration }, (info) => {
          if (token !== generation) return;
          $('power-progress').textContent = info.message ?? `${info.completed} / ${info.total} studies · N=${info.n} · ${info.scenario} · ${info.threads} C++ thread${info.threads === 1 ? '' : 's'}`;
        });
        if (token !== generation) return;
        instance.dispose(); worker = null; $('power-run').disabled = false; $('power-stop').hidden = true;
        result = answer; display(result); $('power-progress').textContent = 'Power calculation complete for the saved estimand.';
      } catch (error) { if (!error.stopped) { invalidate(); $('power-results').replaceChildren(el('p', { class: 'error', role: 'alert' }, error.message)); } }
    });
    $('power-stop').addEventListener('click', () => { invalidate('Power calculation stopped.'); $('power-progress').textContent = ''; });
    $('power-json').addEventListener('click', () => { if (result) download(JSON.stringify(result, null, 2) + '\n', 'application/json', 'estimand-power-report.json'); });
    $('power-csv').addEventListener('click', () => { if (result) download(csv(result), 'text/csv', 'estimand-power.csv'); });
    onSettings({ ...settings, calibrationEdge: edgeValue });
  }
  return { refresh, setDraft, stop };
}
