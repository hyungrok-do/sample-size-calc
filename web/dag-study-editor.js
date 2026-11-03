import { parseStudy, scopeOf, repeatedKeys, sizeUnit } from './dag-study.js?v=dev';
const el = (tag, attrs = {}, ...children) => {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== false && v != null) {
    if (k.startsWith('on')) n.addEventListener(k.slice(2), v); else n.setAttribute(k, v === true ? '' : v);
  }
  n.append(...children.flat()); return n;
};
const number = (id, value) => el('input', { id, type: 'number', step: 'any', value });
const field = (name, input, help = '') => el('div', { class: 'dag-field' }, el('label', { for: input.id }, name), input, ...(help ? [el('p', { class: 'small muted' }, help)] : []));
const choice = (id, options, value) => { const n = el('select', { id }, options.map(([v, label]) => el('option', { value: v }, label))); n.value = value; return n; };
const readNumber = (input) => { if (!input.value.trim() || !Number.isFinite(Number(input.value))) throw new Error('Enter finite numeric settings.'); return Number(input.value); };
const numbers = (input) => input.value.trim().split(/[\s,]+/).filter(Boolean).map(Number);

// The study's schedule, cluster sizes and missed visits. The Clusters and Repeated visits buttons add or remove the
// levels themselves; this shows only the settings of the levels the study has.
export function mountStudy(root, graph, onChange, onError) {
  root.hidden = !graph.study;
  if (!graph.study) { root.replaceChildren(); return; }
  const s = parseStudy(graph.study), visits = s.times.length > 1, clustered = s.clusterSizes.length > 0;
  const times = el('input', { id: 'study-times', value: s.times.join(', ') });
  const clusters = el('input', { id: 'study-clusters', value: s.clusterSizes.join(', ') });
  const missing = choice('study-missing', [['none', 'All visits observed'], ['dropout', 'Monotone dropout'], ['intermittent', 'Intermittent missing visits']], s.observation.mode);
  const history = choice('study-history', [['', 'None'], ...graph.nodes.filter((n) => scopeOf(n) === 'visit').map((n) => [n.id, n.label])], s.observation.historyVariable);
  const treatment = choice('study-treatment', [['', 'None'], ...graph.nodes.filter((n) => ['subject', 'cluster'].includes(scopeOf(n))).map((n) => [n.id, n.label])], s.observation.treatmentVariable);
  const intercept = number('study-drop-intercept', s.observation.intercept), past = number('study-past-effect', s.observation.historyCoefficient),
    current = number('study-current-effect', s.observation.currentCoefficient), arm = number('study-arm-effect', s.observation.treatmentCoefficient);
  const observation = el('div', {}, el('p', { class: 'small muted' }, 'The first visit is observed. Later missing-visit probabilities use a logistic model. Outcomes continue to exist when a visit is missed.'),
    el('div', { class: 'dag-grid' }, field('Missing-visit log-odds intercept', intercept), field('History variable', history), field('Last observed value coefficient', past),
      field('Unobserved current value coefficient (MNAR)', current), field('Treatment variable', treatment), field('Treatment coefficient', arm)));
  const missingField = field('Observation process', missing);
  const visibility = () => { observation.hidden = !visits || missing.value === 'none'; };
  missing.addEventListener('change', visibility); visibility();
  const form = el('form', { novalidate: true }, el('div', { class: 'dag-grid' },
    ...(visits ? [field('Visit times', times, 'Times may be irregular.')] : []),
    ...(clustered ? [field('Subjects per cluster', clusters, 'One size for equal clusters, or a repeating pattern such as 8, 12, 20.')] : []),
    ...(visits ? [missingField] : [])), observation, el('button', { type: 'submit' }, 'Apply'));
  form.addEventListener('submit', (event) => {
    event.preventDefault();
    try {
      const next = structuredClone(graph);
      next.study = parseStudy({ times: visits ? numbers(times) : s.times, clusterSizes: clustered ? numbers(clusters) : [], observation: { mode: visits ? missing.value : 'none', intercept: readNumber(intercept),
        historyVariable: history.value, historyCoefficient: readNumber(past), currentCoefficient: readNumber(current), treatmentVariable: treatment.value, treatmentCoefficient: readNumber(arm) } });
      onChange(next);
    } catch (e) { onError(e.message); }
  });
  const parts = [...(visits ? [`visits at ${s.times.join(', ')}`] : []), ...(clustered ? [`${s.clusterSizes.join(', ')} subjects per cluster`] : []),
    ...(visits ? [{ none: 'all visits observed', dropout: 'monotone dropout', intermittent: 'intermittent missed visits' }[s.observation.mode]] : [])];
  root.replaceChildren(el('summary', {}, 'Study settings', el('span', { class: 'small muted' }, ` · ${parts.join(' · ')} · size counts ${sizeUnit(graph)}`)), form);
}

// onEdit: called when a lag row is added or removed, so the variable's settings apply at once
export function nodeStudyFields(graph, node, onEdit = () => {}) {
  if (!graph.study) return { element: el('div'), read: () => ({}) };
  const scope = choice('dag-node-scope', [['subject', 'Once per subject'], ['cluster', 'Once per cluster'], ['visit', 'At every visit'], ['time', 'Visit time (from schedule)']], scopeOf(node));
  const labels = ['Subject intercept SD', 'Subject slope SD', 'Intercept–slope correlation', 'Cluster intercept SD', 'Residual correlation per time unit'];
  const inputs = repeatedKeys.map((k) => number(`repeat-${k}`, node.repeated?.[k] ?? 0));
  const repeated = el('details', { id: 'dag-node-effects', class: 'dag-dependence' }, el('summary', {}, 'Shared effects and correlation'), el('p', { class: 'small muted' }, 'Random effects act on the variable’s predictor scale. Slopes use time since the first visit. Residual correlation is available for Normal and log-normal errors.'),
    el('div', { class: 'dag-grid' }, inputs.map((input, i) => field(labels[i], input))));
  const lags = el('div'), rows = [];
  function addLag(raw = { source: node.id, lag: 1, coefficient: .3, initial: 0 }) {
    const source = choice(`lag-source-${rows.length}`, graph.nodes.map((n) => [n.id, n.id]), raw.source), lag = number(`lag-length-${rows.length}`, raw.lag),
      coefficient = number(`lag-coef-${rows.length}`, raw.coefficient), initial = number(`lag-initial-${rows.length}`, raw.initial ?? 0);
    const row = { source, lag, coefficient, initial, deleted: false };
    const element = el('div', { class: 'dag-grid' }, field('Previous variable', source), field('Visits back', lag), field('Lag coefficient', coefficient), field('Value before first visit', initial),
      el('button', { type: 'button', onclick: () => { row.deleted = true; element.remove(); onEdit(); } }, 'Remove lag'));
    rows.push(row); lags.append(element);
  }
  for (const l of node.lags ?? []) addLag(l);
  const temporal = el('details', { id: 'dag-node-lags' }, el('summary', {}, 'Lagged relationships'), el('p', { class: 'small muted' }, 'Earlier visits may affect later visits. These are separate from shared subject effects.'), lags,
    el('button', { type: 'button', onclick: () => { addLag(); onEdit(); } }, 'Add lag'));
  const visibility = () => { repeated.hidden = temporal.hidden = scope.value !== 'visit'; };
  scope.addEventListener('change', visibility); visibility();
  return { element: el('div', { class: 'dag-node-study' }, field('Measured at', scope), repeated, temporal), read: () => ({ scope: scope.value,
    repeated: Object.fromEntries(repeatedKeys.map((k, i) => [k, scope.value === 'visit' ? readNumber(inputs[i]) : 0])),
    lags: scope.value === 'visit' ? rows.filter((r) => !r.deleted).map((r) => ({ source: r.source.value, lag: readNumber(r.lag), coefficient: readNumber(r.coefficient), initial: readNumber(r.initial) })) : [] }) };
}

export function mountTimeline(root, graph) {
  root.hidden = !graph.study;
  if (!graph.study) { root.replaceChildren(); return; }
  const svg = (tag, attrs = {}, ...children) => {
    const n = document.createElementNS('http://www.w3.org/2000/svg', tag);
    for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, v);
    n.append(...children); return n;
  };
  const times = graph.study.times, width = 220 * (times.length + 1), height = 72 * graph.nodes.length + 60;
  const picture = svg('svg', { width, height, viewBox: `0 0 ${width} ${height}`, role: 'img', 'aria-label': 'DAG expanded across scheduled visits' });
  picture.append(svg('defs', {}, svg('marker', { id: 'timeline-arrow', viewBox: '0 0 10 10', refX: 9, refY: 5, markerWidth: 5, markerHeight: 5, orient: 'auto' }, svg('path', { d: 'M 0 0 L 10 5 L 0 10 z' }))));
  const byId = new Map(graph.nodes.map((n, i) => [n.id, { node: n, row: i }]));
  const point = (id, visit) => {
    const { node, row } = byId.get(id), repeated = ['visit', 'time'].includes(scopeOf(node));
    return [20 + (repeated ? visit + 1 : 0) * 220, 45 + row * 72];
  };
  picture.append(svg('text', { x: 20, y: 22 }, 'Shared across visits'));
  times.forEach((time, i) => picture.append(svg('text', { x: 240 + i * 220, y: 22 }, `Time ${time}`)));
  const links = svg('g', { class: 'timeline-links' }), boxes = svg('g', { class: 'timeline-boxes' });
  function link(source, fromVisit, target, toVisit, label, lag = false) {
    const [x, y] = point(source, fromVisit), [x2, y2] = point(target, toVisit);
    const path = svg('path', { d: `M ${x + 85} ${y + 38} C ${x + 85} ${y + 60}, ${x2 + 85} ${y2 - 20}, ${x2 + 85} ${y2}`, 'marker-end': 'url(#timeline-arrow)', ...(lag ? { 'stroke-dasharray': '5 4' } : {}) });
    path.append(svg('title', {}, label)); links.append(path);
  }
  for (const edge of graph.edges) {
    const repeated = ['visit', 'time'].includes(scopeOf(byId.get(edge.target).node));
    for (let v = 0; v < (repeated ? times.length : 1); v++) link(edge.source, v, edge.target, v, `${edge.source} → ${edge.target}: ${edge.coefficient}`);
  }
  for (const n of graph.nodes) for (const l of n.lags ?? []) for (let v = l.lag; v < times.length; v++) link(l.source, v - l.lag, n.id, v, `${l.source}, ${l.lag} visits earlier → ${n.id}: ${l.coefficient}`, true);
  for (const n of graph.nodes) {
    const repeated = ['visit', 'time'].includes(scopeOf(n));
    for (let v = 0; v < (repeated ? times.length : 1); v++) {
      const [x, y] = point(n.id, v);
      boxes.append(svg('rect', { x, y, width: 170, height: 38, rx: 6 }), svg('text', { x: x + 10, y: y + 24 }, repeated ? `${n.id} [${times[v]}]` : `${n.id} · ${scopeOf(n)}`));
    }
  }
  picture.append(links, boxes);
  root.replaceChildren(el('summary', {}, 'Expand the DAG across visits'), el('p', { class: 'small muted' }, 'Shared variables appear once. Visit variables repeat on the schedule. Dashed arrows are lagged relationships; hover an arrow for its coefficient. Random-effect settings remain on the variable.'), el('div', { class: 'dag-scroll' }, picture));
}
