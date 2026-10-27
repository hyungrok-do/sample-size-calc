import { addEdge, removeEdge, removeNode, makeNode, exampleGraph, emptyGraph, renameNode, validateGraph,
  equation, dataCsv, KINDS, DISTRIBUTIONS, FAMILIES, MAX_NODES } from './dag-model.js?v=dev';
import { validateEstimand, renameEstimand, exportDesign, parseDesign } from './estimand-model.js?v=dev';
import { mountEstimand } from './dag-estimand.js?v=dev';

const SVG = 'http://www.w3.org/2000/svg';
const STORE = 'ssc-dag-v1';
const number = (x) => x === null ? '—' : Number.isInteger(x) ? String(x) : x.toFixed(3);
const kindNames = { continuous: 'Continuous', binary: 'Binary', count: 'Count', survival: 'Survival time' };
const parameterLabels = { sd: 'Standard deviation', width: 'Uniform width', shape: 'Shape', scale: 'Baseline time scale' };

function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (value == null || value === false) continue;
    if (key.startsWith('on')) node.addEventListener(key.slice(2), value);
    else node.setAttribute(key, value === true ? '' : value);
  }
  node.append(...children.flat().filter((c) => c != null));
  return node;
}
function svg(tag, attrs = {}, ...children) {
  const node = document.createElementNS(SVG, tag);
  for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
  node.append(...children);
  return node;
}
function inputNumber(input, label) {
  if (input.value.trim() === '' || !Number.isFinite(Number(input.value))) throw new Error(`${label}: enter a finite number.`);
  return Number(input.value);
}
function field(label, input, hint = '') {
  const wrap = el('div', { class: 'dag-field' }, el('label', { for: input.id }, label), input);
  if (hint) wrap.append(el('p', { class: 'small muted' }, hint));
  return wrap;
}
function select(id, options, value) {
  const node = el('select', { id }, options.map(([v, label]) => el('option', { value: v }, label)));
  node.value = value;
  return node;
}
const numeric = (id, value, attrs = {}) => el('input', { id, type: 'number', step: 'any', value, ...attrs });
function download(content, type, name) {
  const url = URL.createObjectURL(new Blob([content], { type }));
  const a = el('a', { href: url, download: name });
  document.body.append(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

export function mountDag(root) {
  const $ = (id) => root.querySelector(`#${id}`);
  let graph = emptyGraph(), selected = null, connectFrom = null, result = null, worker = null;
  let estimand = null, estimandEditor;
  let history = [], future = [], storageNote = '';
  try {
    const saved = localStorage.getItem(STORE);
    if (saved) { const design = parseDesign(saved); graph = design.graph; estimand = design.estimand; }
  } catch { storageNote = 'The saved graph could not be restored. An empty canvas is open.'; }
  selected = graph.nodes.at(-1)?.id;
  const canvas = $('dag-canvas');
  function status(message, error = false) {
    $('dag-status').textContent = message;
    $('dag-status').className = error ? 'error' : 'dag-message';
  }
  function persist() {
    try { localStorage.setItem(STORE, JSON.stringify(exportDesign(graph, estimand))); }
    catch { status('Browser storage is unavailable. Export the model to keep your work.'); }
  }
  function clearResult() {
    worker?.terminate(); worker = null; result = null;
    $('dag-simulate').disabled = !graph.nodes.length;
    $('dag-cancel').hidden = true;
    $('dag-download-csv').disabled = true;
    $('dag-results').replaceChildren(el('p', { class: 'muted' }, 'Simulate data to inspect the generated variables.'));
  }
  function remember() {
    history.push(JSON.stringify(exportDesign(graph, estimand)));
    if (history.length > 50) history.shift();
    future = [];
  }
  function apply(next, message, planning = {}) {
    try {
      validateGraph(next);
      let nextEstimand = 'estimand' in planning ? planning.estimand : estimand;
      if (planning.renameFrom) nextEstimand = renameEstimand(nextEstimand, planning.renameFrom, planning.renameTo);
      let note = '';
      if (nextEstimand) {
        try { validateEstimand(next, nextEstimand); }
        catch { nextEstimand = null; note = ' The previous estimand is no longer compatible; define it again below. Undo restores the previous design.'; }
      }
      remember(); graph = next; estimand = nextEstimand;
      if (!graph.nodes.some((n) => n.id === selected)) selected = graph.nodes.at(-1)?.id;
      connectFrom = null;
      clearResult(); persist(); render();
      status((message ?? 'Model updated. Simulate again to see its data.') + note);
      return true;
    } catch (error) { status(error.message, true); return false; }
  }
  function updateSelection(id) {
    selected = id;
    draw(); inspector();
  }
  function connectTo(id) {
    if (!connectFrom) {
      connectFrom = id; draw();
      status(`Arrow from ${id}: choose the destination node, or press Escape to cancel.`);
    } else {
      try {
        const source = connectFrom;
        if (apply(addEdge(graph, source, id), `Added ${source} → ${id}. Edit its coefficient in the destination node.`)) updateSelection(id);
      } catch (error) { status(error.message, true); }
    }
  }
  function draw() {
    const defs = svg('defs', {}, svg('marker', { id: 'dag-arrow', viewBox: '0 0 10 10', refX: 9, refY: 5,
      markerWidth: 7, markerHeight: 7, orient: 'auto-start-reverse' }, svg('path', { d: 'M 0 0 L 10 5 L 0 10 z' })));
    const lines = svg('g', { class: 'dag-arrows' });
    const map = new Map(graph.nodes.map((n) => [n.id, n]));
    for (const edge of graph.edges) {
      const source = map.get(edge.source), target = map.get(edge.target);
      const dx = target.x - source.x, dy = target.y - source.y;
      const vertical = Math.abs(dy) > Math.abs(dx) * 0.6;
      const x1 = source.x + (vertical ? 85 : dx >= 0 ? 170 : 0);
      const y1 = source.y + (vertical ? dy >= 0 ? 84 : 0 : 42);
      const x2 = target.x + (vertical ? 85 : dx >= 0 ? 0 : 170);
      const y2 = target.y + (vertical ? dy >= 0 ? 0 : 84 : 42);
      const bend = Math.max(30, Math.abs(vertical ? y2 - y1 : x2 - x1) * 0.4);
      const d = vertical
        ? `M ${x1} ${y1} C ${x1} ${y1 + Math.sign(dy) * bend}, ${x2} ${y2 - Math.sign(dy) * bend}, ${x2} ${y2}`
        : `M ${x1} ${y1} C ${x1 + Math.sign(dx) * bend} ${y1}, ${x2 - Math.sign(dx) * bend} ${y2}, ${x2} ${y2}`;
      const group = svg('g', { class: edge.target === selected ? 'dag-edge selected' : 'dag-edge',
        role: 'button', tabindex: 0, 'aria-label': `Edit ${edge.source} to ${edge.target}, coefficient ${edge.coefficient}` },
      svg('path', { d, 'marker-end': 'url(#dag-arrow)' }),
      svg('text', { x: (x1 + x2) / 2, y: (y1 + y2) / 2 - 8, 'text-anchor': 'middle' }, String(Number(edge.coefficient.toPrecision(4)))));
      const pick = () => { updateSelection(edge.target); document.getElementById(`dag-edge-${edge.source}`)?.focus(); };
      group.addEventListener('click', pick);
      group.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); pick(); } });
      lines.append(group);
    }
    const nodes = svg('g', {});
    for (const node of graph.nodes) {
      const group = svg('g', { class: `dag-node ${node.kind}${node.id === selected ? ' selected' : ''}${node.id === connectFrom ? ' connecting' : ''}`,
        transform: `translate(${node.x} ${node.y})`, tabindex: 0, role: 'button',
        'aria-label': `${node.id}, ${node.label}, ${kindNames[node.kind]}. Enter to edit. Arrow keys to move.` },
      svg('rect', { width: 170, height: 84, rx: 12 }),
      svg('text', { x: 14, y: 22, class: 'dag-node-role' }, node.kind.toUpperCase()),
      svg('text', { x: 14, y: 45, class: 'dag-node-title' }, `${node.id} · ${node.label.length > 15 ? node.label.slice(0, 14) + '…' : node.label}`),
      svg('text', { x: 14, y: 66, class: 'dag-node-family' }, DISTRIBUTIONS[node.family].label
        + (node.kind === 'survival' && node.censoring.mode !== 'none'
          ? ` · ${node.censoring.mode === 'variable' ? 'censor: ' + node.censoring.variable : 'cutoff: ' + node.censoring.time}` : '')));
      const port = svg('g', { class: 'dag-port', role: 'button', tabindex: 0, 'aria-label': `Start an arrow from ${node.id}` },
        svg('circle', { cx: 170, cy: 42, r: 10 }), svg('text', { x: 170, y: 46, 'text-anchor': 'middle' }, '+'));
      port.addEventListener('pointerdown', (e) => e.stopPropagation());
      port.addEventListener('pointerup', (e) => e.stopPropagation());
      port.addEventListener('click', (e) => { e.stopPropagation(); connectTo(node.id); });
      port.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' || e.key === ' ') { e.stopPropagation(); e.preventDefault(); connectTo(node.id); }
      });
      group.append(port);
      let drag = null;
      group.addEventListener('pointerdown', (event) => {
        if (event.button !== 0 || connectFrom) return;
        const p = point(event);
        drag = { start: p, x: node.x, y: node.y, moved: false };
        group.setPointerCapture(event.pointerId);
      });
      group.addEventListener('pointermove', (event) => {
        if (!drag) return;
        const p = point(event), dx = p.x - drag.start.x, dy = p.y - drag.start.y;
        if (!drag.moved && Math.hypot(dx, dy) < 4) return;
        if (!drag.moved) { remember(); drag.moved = true; }
        node.x = Math.max(0, Math.min(920, Math.round(drag.x + dx)));
        node.y = Math.max(0, Math.min(440, Math.round(drag.y + dy)));
        group.setAttribute('transform', `translate(${node.x} ${node.y})`);
      });
      const finish = () => {
        const moved = drag?.moved;
        drag = null;
        if (moved) { persist(); draw(); buttons(); }
        return moved;
      };
      group.addEventListener('pointercancel', finish);
      group.addEventListener('pointerup', (event) => {
        if (finish()) { event.stopPropagation(); return; }
        if (connectFrom) connectTo(node.id);
        else updateSelection(node.id);
      });
      group.addEventListener('keydown', (event) => {
        if (event.target !== group) return;
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault();
          if (connectFrom) connectTo(node.id); else { updateSelection(node.id); $('dag-node-label').focus(); }
        }
        const delta = { ArrowLeft: [-10, 0], ArrowRight: [10, 0], ArrowUp: [0, -10], ArrowDown: [0, 10] }[event.key];
        if (delta) {
          event.preventDefault(); remember();
          node.x = Math.max(0, Math.min(920, node.x + delta[0]));
          node.y = Math.max(0, Math.min(440, node.y + delta[1]));
          persist(); draw(); buttons();
          [...canvas.querySelectorAll('.dag-node')].find((g) => g.getAttribute('aria-label').startsWith(`${node.id},`))?.focus();
        }
      });
      nodes.append(group);
    }
    if (!graph.nodes.length) nodes.append(svg('text', { x: 550, y: 260, 'text-anchor': 'middle', class: 'dag-empty' }, 'Add a variable to begin'));
    canvas.replaceChildren(svg('title', {}, 'Data generating process'),
      svg('desc', {}, 'Select a node to edit. Drag nodes or use arrow keys to move. Use plus ports or the arrow form to connect nodes.'), defs, lines, nodes);
    $('dag-node-count').textContent = `${graph.nodes.length} variable${graph.nodes.length === 1 ? '' : 's'} · ${graph.edges.length} arrow${graph.edges.length === 1 ? '' : 's'}`;
    $('dag-connect').querySelector('button').disabled = graph.nodes.length < 2;
    const options = graph.nodes.map((n) => [n.id, `${n.id} · ${n.label}`]);
    for (const key of ['dag-source', 'dag-target', 'dag-node-select']) {
      const previous = $(key).value;
      $(key).replaceChildren(...options.map(([value, label]) => el('option', { value }, label)));
      $(key).value = key === 'dag-node-select' ? selected : graph.nodes.some((n) => n.id === previous) ? previous : graph.nodes[0]?.id ?? '';
      $(key).disabled = !graph.nodes.length;
    }
    if ($('dag-source').value === $('dag-target').value && graph.nodes.length > 1)
      $('dag-target').value = graph.nodes.find((n) => n.id !== $('dag-source').value).id;
  }
  function point(event) {
    const p = canvas.createSVGPoint(); p.x = event.clientX; p.y = event.clientY;
    return p.matrixTransform(canvas.getScreenCTM().inverse());
  }
  function inspector() {
    const node = graph.nodes.find((n) => n.id === selected);
    const box = $('dag-inspector'); box.replaceChildren();
    if (!node) { box.append(el('p', { class: 'muted small' }, 'Add a variable, then choose its type and distribution here.')); return; }
    box.append(el('div', { class: 'dag-inspector-heading' }, el('h3', {}, `Edit ${node.id}`),
      el('span', { class: 'small muted' }, kindNames[node.kind])));
    const id = el('input', { id: 'dag-node-id', type: 'text', maxlength: 32, value: node.id, required: true });
    const label = el('input', { id: 'dag-node-label', type: 'text', maxlength: 80, value: node.label, required: true });
    const kind = select('dag-node-kind', KINDS.map((k) => [k, kindNames[k]]), node.kind);
    const family = select('dag-node-family', [], node.family);
    const intercept = numeric('dag-node-intercept', node.intercept);
    const probability = numeric('dag-node-probability', 1 / (1 + Math.exp(-node.intercept)), { min: 0, max: 1 });
    const inputs = Object.fromEntries(['sd', 'width', 'shape', 'scale'].map((key) => [key, numeric(`dag-node-${key}`, node[key], { min: 0 })]));
    const censorMode = select('dag-censor-mode', [['none', 'None · fully observed'], ['administrative', 'Administrative cutoff'], ['variable', 'Use another variable']], node.censoring.mode);
    const censorTime = numeric('dag-censor-time', node.censoring.time, { min: 0 });
    const censorVariable = select('dag-censor-variable', graph.nodes.filter((n) => n.id !== node.id).map((n) => [n.id, `${n.id} · ${n.label}`]), node.censoring.variable);
    const form = el('form', { novalidate: true }, field('Variable ID', id), field('Label', label), field('Variable type', kind), field('Distribution', family));
    const parameters = el('div'), observation = el('div');
    function observationFields() {
      observation.replaceChildren();
      if (kind.value !== 'survival') return;
      observation.append(field('Right-censoring observation rule', censorMode,
        'The variable itself is the latent event time. Observed time and event indicator are additional output columns.'));
      if (censorMode.value === 'administrative') observation.append(field('Administrative cutoff', censorTime));
      if (censorMode.value === 'variable') observation.append(field('Censoring time variable', censorVariable,
        'Its generated values must be nonnegative. It need not be a parent of the event-time variable.'));
    }
    function params() {
      const spec = DISTRIBUTIONS[family.value];
      parameters.replaceChildren(family.value === 'bernoulli'
        ? field('Success probability when all parents equal zero', probability, 'The conditional predictor uses its log odds plus incoming parent terms. Use a probability strictly between 0 and 1.')
        : field(`Intercept · ${spec.link}`, intercept,
          'The predictor is evaluated when all parents equal zero. Incoming arrows add terms in each parent’s original units.'));
      for (const key of spec.parameters) parameters.append(field(
        key === 'sd' && family.value === 'lognormal' ? 'Standard deviation of log values' : parameterLabels[key], inputs[key]));
      if (kind.value === 'survival') parameters.append(el('p', { class: 'small muted' },
        'Weibull/exponential proportional hazards: a parent coefficient is a log hazard ratio. The baseline scale is the time scale at predictor zero.'));
      observationFields();
    }
    function families(preferred) {
      const options = FAMILIES.filter((f) => DISTRIBUTIONS[f].kind === kind.value);
      family.replaceChildren(...options.map((f) => el('option', { value: f }, DISTRIBUTIONS[f].label)));
      family.value = options.includes(preferred) ? preferred : options[0]; params();
    }
    kind.addEventListener('change', () => families(family.value));
    family.addEventListener('change', params); censorMode.addEventListener('change', observationFields);
    families(node.family);
    form.append(parameters, observation, el('button', { type: 'submit', class: 'primary' }, 'Apply variable settings'));
    form.addEventListener('submit', (event) => {
      event.preventDefault();
      try {
        let baseline = inputNumber(intercept, 'Intercept');
        if (family.value === 'bernoulli') {
          const p = inputNumber(probability, 'Success probability');
          if (!(p > 0 && p < 1)) throw new Error('Success probability must be strictly between 0 and 1.');
          baseline = Math.log(p) - Math.log1p(-p);
        }
        const changed = { ...node, label: label.value.trim(), kind: kind.value, family: family.value,
          intercept: baseline,
          censoring: kind.value === 'survival' ? { mode: censorMode.value,
            time: censorMode.value === 'administrative' ? inputNumber(censorTime, 'Administrative cutoff') : 20,
            variable: censorMode.value === 'variable' ? censorVariable.value : '' } : { mode: 'none', time: 20, variable: '' } };
        for (const key of DISTRIBUTIONS[family.value].parameters) changed[key] = inputNumber(inputs[key], parameterLabels[key]);
        const next = renameNode({ ...graph, nodes: graph.nodes.map((n) => n.id === selected ? changed : n) }, node.id, id.value.trim());
        if (apply(next, undefined, { renameFrom: node.id, renameTo: id.value.trim() })) updateSelection(id.value.trim());
      } catch (error) { status(error.message, true); }
    });
    box.append(form);
    const edges = graph.edges.filter((e) => e.target === node.id);
    box.append(el('h4', {}, 'Incoming effects'), el('p', { class: 'small muted' },
      `Coefficients on the ${DISTRIBUTIONS[node.family].link} scale. With interactions, a main coefficient applies when its interacting parents equal zero.`));
    for (const edge of edges) {
      const input = numeric(`dag-edge-${edge.source}`, edge.coefficient, { 'aria-label': `${edge.source} to ${node.id} coefficient` });
      input.addEventListener('change', () => {
        try { apply({ ...graph, edges: graph.edges.map((e) => e === edge ? { ...e, coefficient: inputNumber(input, 'Arrow coefficient') } : e) }); }
        catch (error) { status(error.message, true); }
      });
      box.append(el('div', { class: 'dag-effect-row' }, el('label', { for: input.id }, `${edge.source} → ${node.id}`), input,
        el('button', { type: 'button', class: 'secondary', 'aria-label': `Remove ${edge.source} to ${node.id}`,
          onclick: () => apply(removeEdge(graph, edge.source, edge.target), 'Arrow and its associated interactions removed.') }, '×')));
    }
    if (!edges.length) box.append(el('p', { class: 'small muted' }, 'No parents. This variable is generated independently.'));
    box.append(el('h4', {}, 'Parent interactions'));
    for (const [index, term] of node.interactions.entries()) {
      const input = numeric(`dag-term-${index}`, term.coefficient, { 'aria-label': `${term.left} times ${term.right} coefficient` });
      input.addEventListener('change', () => {
        try { apply({ ...graph, nodes: graph.nodes.map((n) => n.id !== selected ? n : { ...n,
          interactions: n.interactions.map((t, i) => i === index ? { ...t, coefficient: inputNumber(input, 'Interaction coefficient') } : t) }) }); }
        catch (error) { status(error.message, true); }
      });
      box.append(el('div', { class: 'dag-effect-row' }, el('label', { for: input.id }, `${term.left} × ${term.right}`), input,
        el('button', { type: 'button', class: 'secondary', 'aria-label': `Remove ${term.left} times ${term.right}`, onclick: () => apply({ ...graph,
          nodes: graph.nodes.map((n) => n.id !== selected ? n : { ...n, interactions: n.interactions.filter((_, i) => i !== index) }) }) }, '×')));
    }
    if (edges.length >= 2) {
      const choices = edges.map((e) => [e.source, e.source]);
      const left = select('dag-interaction-left', choices, edges[0].source);
      const right = select('dag-interaction-right', choices, edges[1].source);
      const coefficient = numeric('dag-interaction-value', 0.2);
      const interaction = el('form', { class: 'dag-interaction-form', novalidate: true },
        field('Parent 1', left), field('Parent 2', right), field('Coefficient', coefficient),
        el('button', { type: 'submit', class: 'secondary' }, 'Add interaction'));
      interaction.addEventListener('submit', (event) => {
        event.preventDefault();
        try { apply({ ...graph, nodes: graph.nodes.map((n) => n.id !== selected ? n : { ...n, interactions: [...n.interactions,
          { left: left.value, right: right.value, coefficient: inputNumber(coefficient, 'Interaction coefficient') }] }) }); }
        catch (error) { status(error.message, true); }
      });
      box.append(interaction);
    } else box.append(el('p', { class: 'small muted' }, 'Connect at least two parents to add an interaction.'));
    box.append(el('button', { type: 'button', class: 'secondary dag-delete', onclick: () => apply(removeNode(graph, selected), 'Variable and its arrows removed. Any observation rule using it was reset. Undo restores them.') }, 'Delete variable'));
  }
  function buttons() { $('dag-undo').disabled = !history.length; $('dag-redo').disabled = !future.length; }
  function render() {
    draw(); inspector(); buttons();
    estimandEditor?.refresh(estimand);
    const { order, byId } = validateGraph(graph);
    $('dag-equations').replaceChildren(...order.map((id) => el('li', {}, el('code', {}, equation(graph, byId.get(id))))));
    $('dag-order').textContent = order.length
      ? `Generation order: ${order.join(' → ')}. Exogenous draws are independent across variables and rows.`
      : 'Add variables to define the conditional distributions and generation order.';
  }
  function showResult() {
    const table = (head, rows, label) => el('div', { class: 'table-wrap', tabindex: 0, 'aria-label': label },
      el('table', {}, el('thead', {}, el('tr', {}, head.map((v) => el('th', { scope: 'col' }, v)))),
        el('tbody', {}, rows.map((row) => el('tr', {}, row.map((v) => el('td', {}, v)))))));
    const ids = Object.keys(result.columns);
    const survivalRows = result.survival.map((s) => [s.id, String(s.events), String(s.censored), `${(100 * s.censored / result.n).toFixed(1)}%`]);
    $('dag-results').replaceChildren(
      el('p', { class: 'dag-result-heading' }, `One simulated dataset · N = ${result.n} · seed ${result.seed}`),
      el('p', { class: 'small muted' }, 'Descriptive summaries of generated variables. Power analysis will use a separately specified analysis and repeated simulations.'),
      ...(survivalRows.length ? [el('h4', {}, 'Survival observations'),
        table(['Variable', 'Events', 'Censored', 'Censored %'], survivalRows, 'Survival observation summaries'),
        el('p', { class: 'small muted' }, 'For each survival variable T: T is the latent event time, T_time is observed follow-up, and T_event is 1 for an event and 0 for censoring.')] : []),
      el('h4', {}, 'All variables and observation columns'),
      table(['Variable', 'Mean', 'SD', 'Minimum', 'Maximum'], result.summaries.map((s) =>
        [s.id, number(s.mean), number(s.sd), number(s.min), number(s.max)]), 'Variable summaries'),
      el('h4', {}, 'First 10 observations'),
      table(ids, Array.from({ length: Math.min(10, result.n) }, (_, i) =>
        ids.map((id) => number(result.columns[id][i]))), 'Generated data preview'));
    $('dag-download-csv').disabled = false;
  }
  $('dag-add').addEventListener('click', () => {
    if (graph.nodes.length >= MAX_NODES) { status(`Use at most ${MAX_NODES} variables.`, true); return; }
    let id = 'V1', k = 1;
    while (graph.nodes.some((n) => n.id === id)) id = `V${++k}`;
    const count = graph.nodes.length;
    const node = makeNode(id, 'continuous', 30 + (count % 4) * 260, 30 + (Math.floor(count / 4) % 3) * 150);
    node.label = 'Variable';
    if (apply({ ...graph, nodes: [...graph.nodes, node] }, `Added ${id}. Choose its type and distribution in the inspector.`)) updateSelection(id);
  });
  $('dag-connect').addEventListener('submit', (event) => {
    event.preventDefault();
    try {
      const source = $('dag-source').value, target = $('dag-target').value;
      if (apply(addEdge(graph, source, target, inputNumber($('dag-coefficient'), 'Arrow coefficient')), `Added ${source} → ${target}.`)) updateSelection(target);
    } catch (error) { status(error.message, true); }
  });
  $('dag-node-select').addEventListener('change', (e) => updateSelection(e.target.value));
  root.addEventListener('keydown', (event) => {
    if (event.key === 'Escape') { connectFrom = null; draw(); status('Arrow connection cancelled.'); }
  });
  function timeTravel(from, to) {
    if (!from.length) return;
    to.push(JSON.stringify(exportDesign(graph, estimand)));
    const design = parseDesign(from.pop()); graph = design.graph; estimand = design.estimand;
    selected = graph.nodes.some((n) => n.id === selected) ? selected : graph.nodes.at(-1)?.id;
    connectFrom = null; clearResult(); persist(); render(); status('Model restored.');
  }
  $('dag-undo').addEventListener('click', () => timeTravel(history, future));
  $('dag-redo').addEventListener('click', () => timeTravel(future, history));
  $('dag-example-load').addEventListener('click', () => apply(exampleGraph($('dag-example').value), 'Example loaded. All numbers are editable example assumptions.'));
  $('dag-export').addEventListener('click', () => download(JSON.stringify(exportDesign(graph, estimand), null, 2) + '\n', 'application/json', 'simulation-design.json'));
  $('dag-import').addEventListener('change', async (event) => {
    const file = event.target.files?.[0];
    if (!file) return;
    try {
      if (file.size > 1e6) throw new Error('Use a model JSON file smaller than 1 MB.');
      const design = parseDesign(await file.text());
      apply(design.graph, 'Simulation design imported.', { estimand: design.estimand });
    } catch (error) { status(`Could not import model: ${error.message}`, true); }
    event.target.value = '';
  });
  $('dag-simulation').addEventListener('submit', (event) => {
    event.preventDefault(); clearResult();
    try {
      validateGraph(graph);
      const n = inputNumber($('dag-n'), 'Sample size'), seed = inputNumber($('dag-seed'), 'Seed');
      worker = new Worker(new URL('./dag-worker.js?v=dev', import.meta.url), { type: 'module' });
      $('dag-simulate').disabled = true; $('dag-cancel').hidden = false;
      status('Generating data…');
      worker.onmessage = ({ data }) => {
        worker.terminate(); worker = null; $('dag-simulate').disabled = false; $('dag-cancel').hidden = true;
        if (data.error) status(data.error, true);
        else { result = data.result; showResult(); status('Data generated. Inspect summaries and download the full dataset.'); }
      };
      worker.onerror = (error) => { clearResult(); status(`The simulator could not start: ${error.message}`, true); };
      worker.postMessage({ graph, options: { n, seed } });
    } catch (error) { clearResult(); status(error.message, true); }
  });
  $('dag-cancel').addEventListener('click', () => { clearResult(); status('Simulation stopped.'); });
  for (const id of ['dag-n', 'dag-seed']) $(id).addEventListener('input', () => {
    clearResult(); status('Simulation settings changed. Run again to update the data.');
  });
  $('dag-download-csv').addEventListener('click', () => {
    if (result) download(dataCsv(graph, result), 'text/csv', `data-n${result.n}-seed${result.seed}.csv`);
  });
  estimandEditor = mountEstimand($('dag-estimand'), {
    getGraph: () => graph,
    onSave(value) {
      remember(); estimand = value; persist(); buttons(); estimandEditor.refresh(estimand);
      status('Estimand saved with this design. Next, specify the analysis and test for power evaluation.');
    },
  });
  render(); clearResult();
  status(storageNote || (graph.nodes.length
    ? 'Design restored. Edit a variable or the estimand to continue planning.'
    : 'Add a variable to start, or load an example. Each variable has its own type and distribution.'));
}
