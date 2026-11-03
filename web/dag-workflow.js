// Presentation only. Changing sections never changes the design or stops a simulation.
const STEPS = ['study', 'model', 'question', 'power'];
const NAMES = { study: 'Study', model: 'Model', question: 'Question', power: 'Power' };
const TITLES = { study: 'dag-study-title', model: 'dag-model-step-title', question: 'dag-question-title', power: 'dag-power-title' };
const VIEW_KEY = 'ssc-dag-view-v1';

export function studySummary(graph) {
  if (!graph.nodes.length) return 'No design yet · choose a starting point below or build your own model.';
  const s = graph.study;
  return `${graph.nodes.length} variables · ${s ? `${s.times.length} scheduled visit${s.times.length === 1 ? '' : 's'} · sample size counts ${s.clusterSizes.length ? 'clusters' : 'subjects'}` : 'independent observations · sample size counts observations'}`;
}

export function questionSummary(graph, target) {
  if (!target) return 'Choose a compared variable and an outcome to define your question.';
  const name = (id) => graph.nodes.find((n) => n.id === id)?.label || id;
  const outcome = name(target.outcome), exposure = name(target.exposure);
  if (target.kind === 'coefficient') return `Estimate the ${target.measure === 'exponentiated_coefficient' ? 'exponentiated ' : ''}coefficient of ${exposure} in the model for ${outcome}.`;
  if (target.kind === 'joint') return `Test ${target.measure === 'nonlinearity' ? 'nonlinearity' : 'the overall effect'} of ${exposure} in the model for ${outcome}.`;
  const measures = { mean_difference: 'mean difference', risk_difference: 'risk difference', risk_ratio: 'risk ratio', odds_ratio: 'odds ratio', survival_difference: 'survival probability difference', rmst_difference: 'restricted mean survival time difference' };
  let response = outcome, when = '';
  if (graph.study) {
    const r = target.repeated ?? { contrast: 'final' }, times = graph.study.times;
    if (r.contrast === 'change') response = `the change in ${outcome} from time ${times[0]} to ${times.at(-1)}`;
    else if (r.contrast === 'slope') response = `the linear trend in ${outcome} over time`;
    else when = r.contrast === 'average' ? ' averaged across scheduled visits'
      : r.contrast === 'custom' ? ' using your custom visit weights' : ` at time ${times[r.visit ?? times.length - 1]}`;
  } else if (target.horizon != null) when = ` ${target.measure === 'rmst_difference' ? 'through' : 'at'} time ${target.horizon}`;
  return `Compare ${response}${when} if ${exposure} were set to ${target.active} versus ${target.reference}, using the ${measures[target.measure] ?? target.measure.replaceAll('_', ' ')}.`;
}

export function mountWorkflow(root, { onChange = () => {} } = {}) {
  let step = 'study', all = false;
  try {
    const saved = JSON.parse(sessionStorage.getItem(VIEW_KEY));
    if (STEPS.includes(saved?.step)) step = saved.step;
    all = saved?.all === true;
  } catch { /* A blocked or stale view preference does not affect the study. */ }
  const toggle = root.querySelector('#dag-show-all');
  const back = root.querySelector('#dag-step-back'), next = root.querySelector('#dag-step-next');
  function refresh() {
    root.dataset.view = all ? 'all' : step;
    toggle.checked = all;
    for (const panel of root.querySelectorAll('[data-dag-step]')) panel.hidden = !all && !panel.dataset.dagStep.split(' ').includes(step);
    for (const button of root.querySelectorAll('[data-step]')) {
      if (!all && button.dataset.step === step) button.setAttribute('aria-current', 'step');
      else button.removeAttribute('aria-current');
    }
    const index = STEPS.indexOf(step);
    back.hidden = all || index === 0;
    next.hidden = all || index === STEPS.length - 1;
    back.textContent = `← Back to ${NAMES[STEPS[index - 1]] ?? 'Study'}`;
    next.textContent = `Continue to ${NAMES[STEPS[index + 1]] ?? 'Power'} →`;
    root.querySelector('#dag-step-position').textContent = all ? 'All sections · your design is saved in this browser' : `Step ${index + 1} of 4 · you can jump to any step`;
    try { sessionStorage.setItem(VIEW_KEY, JSON.stringify({ step, all })); } catch { /* optional preference */ }
    onChange(all ? 'all' : step);
  }
  function go(value, focus = true) {
    if (!STEPS.includes(value)) return;
    step = value; all = false; refresh();
    if (focus) {
      const heading = root.querySelector(`#${TITLES[step]}`) ?? root.querySelector(`[data-step="${step}"]`);
      heading?.focus({ preventScroll: true });
      root.querySelector('.dag-steps').scrollIntoView({ block: 'start' });
    }
  }
  root.addEventListener('click', (event) => {
    const button = event.target.closest('[data-step], [data-goto]');
    if (button) go(button.dataset.step ?? button.dataset.goto);
  });
  toggle.addEventListener('change', () => { all = toggle.checked; refresh(); });
  back.addEventListener('click', () => go(STEPS[STEPS.indexOf(step) - 1]));
  next.addEventListener('click', () => go(STEPS[STEPS.indexOf(step) + 1]));
  refresh();
  return { refresh, go, update(graph) { root.querySelector('#dag-design-summary').textContent = studySummary(graph); refresh(); } };
}
