import { SCENARIOS, scenarioDesign } from './dag-scenarios.js?v=dev';
import { DEFAULT_POWER } from './power-settings.js?v=dev';
import { mountStudy, nodeStudyFields, mountTimeline } from './dag-study-editor.js?v=dev';
import { scopeOf, sizeUnit, parseStudy } from './dag-study.js?v=dev';
import { reconcileAnalysis } from './dag-functions.js?v=dev';
import { mountFunctions } from './dag-function-editor.js?v=dev';
import { addEdge, removeEdge, removeNode, makeNode, emptyGraph, renameNode, validateGraph,
  equation, dataCsv, KINDS, DISTRIBUTIONS, FAMILIES, MAX_NODES, parameterValue } from './dag-model.js?v=dev';
import { NODE_W, NODE_H, CANVAS_LIMIT, ASPECT, plateLayout, plateAt, changeScope, setStructure, edgePath, loopPath, layout,
  freeSpot, connectTargets, inferRoles, contentBounds, fitView, zoomView, viewBox } from './dag-canvas.js?v=dev';
import { validateEstimand, renameEstimand, exportDesign, parseDesign } from './estimand-model.js?v=dev';
import { mountPower } from './dag-power.js?v=dev';
import { studySummary } from './dag-workflow.js?v=dev';
import { calibrationEdges, powerCompatibility } from './power-model.js?v=dev';
import { sketch, findPaths, classifyPath, pathText } from './dag-view.js?v=dev';
import { Engine } from './engine.js?v=dev';

const SVG = 'http://www.w3.org/2000/svg';
const STORE = 'ssc-dag-v1';
const number = (x) => x === null ? '—' : Number.isInteger(x) ? String(x) : x.toFixed(3);
const kindNames = { continuous: 'Continuous', binary: 'Binary', count: 'Count', survival: 'Survival time' };
const parameterLabels = { sd: 'Standard deviation', width: 'Uniform width', shape: 'Shape', scale: 'Baseline time scale',
  df: 'Degrees of freedom', precision: 'Precision φ (variance = μ(1 − μ) / (1 + φ))', overdispersion: 'Overdispersion α (variance = μ + αμ²)' };
const PREVIEW_N = 2000;  // observations drawn to sketch each variable's distribution
const SKETCH_W = 220;  // width of the selected variable's sketch in its settings
// how an arrow marks a response-curve term of its source in the child
const TERM_MARK = { square: 'x²', log: 'log', threshold: 'step', hinge: 'hinge', saturation: 'tanh', spline: 'spline' };
const short = (x) => String(Number(x.toPrecision(3))).replace('-', '−');
// What a node's box says about its distribution: a root's values, or the scale its incoming arrows act on.
function summary(graph, node) {
  if (graph.study && scopeOf(node) === 'time') return `Visit times · ${graph.study.times.join(', ')}`;
  const i = node.intercept, e = Math.exp(i), family = DISTRIBUTIONS[node.family].label;
  const param = (key) => short(parameterValue(node, key));  // every family's text is built, so read parameters safely
  if (graph.edges.some((edge) => edge.target === node.id)) {
    const rest = { gaussian: `SD ${short(node.sd)}`, uniform: `width ${short(node.width)}`, lognormal: 'log scale', gamma: 'log mean',
      student: `scale ${short(node.sd)}, ${param('df')} df`, beta: 'log odds of the mean', bernoulli: 'log-odds', poisson: 'log mean',
      negbin: `log mean, α ${param('overdispersion')}`, weibull: 'log HR', exponential: 'log HR' }[node.family];
    return `${family} · ${rest}`;
  }
  return `${family} · ${{ gaussian: `mean ${short(i)}, SD ${short(node.sd)}`, uniform: `center ${short(i)}, width ${short(node.width)}`,
    lognormal: `median ${short(e)}`, gamma: `mean ${short(e)}, shape ${short(node.shape)}`, bernoulli: `P(1) = ${short(1 / (1 + Math.exp(-i)))}`,
    student: `center ${short(i)}, scale ${short(node.sd)}, ${param('df')} df`, beta: `mean ${short(1 / (1 + Math.exp(-i)))}, precision ${param('precision')}`,
    poisson: `mean ${short(e)}`, negbin: `mean ${short(e)}, α ${param('overdispersion')}`,
    weibull: `k ${short(node.shape)}, scale ${short(node.scale)}`, exponential: `scale ${short(node.scale)}` }[node.family]}`;
}
// The second line of a box: the distribution in a few characters, with a root's parameters.
function compact(graph, node) {
  if (graph.study && scopeOf(node) === 'time') return 'visit time';
  if (node.measured === false) return 'unmeasured';
  const name = DISTRIBUTIONS[node.family].label;
  if (graph.edges.some((edge) => edge.target === node.id)) return name;
  const i = node.intercept;
  return { gaussian: `N(${short(i)}, ${short(node.sd)}²)`, uniform: `U(${short(i - node.width / 2)}, ${short(i + node.width / 2)})`,
    bernoulli: `P(1) = ${short(1 / (1 + Math.exp(-i)))}`, poisson: `Poisson(${short(Math.exp(i))})` }[node.family] ?? name;
}
// What the canvas says about a variable's inferred role for the comparison
const roleHelp = { confounder: 'Affects both the compared variable and the outcome: adjust for it', mediator: 'On a path from the compared variable to the outcome',
  collider: 'Affected by both the compared variable and the outcome: do not adjust for it', instrument: 'Affects the outcome only through the compared variable',
  prognostic: 'Affects the outcome only: adjusting for it can add precision' };
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
  let estimand = null, powerEditor, powerSettings = null;
  // pending: the design before target edits that have not yet become an undo step
  let history = [], future = [], pending = null, storageNote = '';
  // a sample from the current graph, drawn in the boxes; `key` is the model it came from (positions aside)
  const preview = { columns: null, key: '', token: 0, timer: 0, sketches: new Map() };
  // Starting an engine loads and compiles the WebAssembly module, which costs more than a preview itself, so
  // previews share one engine and the page closes it a moment after the last answer.
  const background = { engine: null, active: 0, idle: 0 };
  function borrowEngine() {
    clearTimeout(background.idle); background.active++;
    return background.engine ??= new Engine();
  }
  function returnEngine(instance) {
    background.active = Math.max(0, background.active - 1);
    if (background.active || background.engine !== instance) return;
    background.idle = setTimeout(() => { background.engine?.dispose(); background.engine = null; }, 1500);
  }
  const modelKey = () => JSON.stringify([graph.nodes.map(({ x, y, ...rest }) => rest), graph.edges, graph.study]);
  let menu = null;
  // the canvas: what part of the drawing is in view, the selected arrow, the plates and random effects last drawn,
  // the open non-causal paths' arrows, a rubber-band arrow being drawn, and what is being dragged
  let view = null, selectedEdge = null, drawn = { plates: [], effects: [] }, biased = new Set(), pointer = null, gesture = null;
  try {
    const saved = localStorage.getItem(STORE);
    if (saved) { const design = parseDesign(saved); graph = design.graph; estimand = design.estimand; powerSettings = design.power ?? null; }
    if (estimand?.kind === 'custom') estimand = null;
  } catch { storageNote = 'The saved graph could not be restored. An empty canvas is open.'; }
  selected = graph.nodes.at(-1)?.id;
  const canvas = $('dag-canvas');
  function status(message, error = false) {
    $('dag-status').textContent = message;
    $('dag-status').className = error ? 'error' : 'dag-message';
  }
  function persist() {
    try { localStorage.setItem(STORE, JSON.stringify(exportDesign(graph, estimand, powerSettings))); }
    catch { status('Browser storage is unavailable. Export the model to keep your work.'); }
  }
  function clearResult() {
    worker?.dispose(); worker = null; result = null;
    $('dag-simulate').disabled = !graph.nodes.length;
    $('dag-cancel').hidden = true;
    $('dag-download-csv').disabled = true;
    $('dag-results').replaceChildren();
  }
  function remember() {
    const current = JSON.stringify(exportDesign(graph, estimand, powerSettings));
    if (pending !== null && pending !== current) history.push(pending);
    history.push(current);
    while (history.length > 50) history.shift();
    future = []; pending = null;
  }
  // planning.remembered: the undo step was already taken (a drag that also moved a variable to another level)
  function apply(next, message, planning = {}) {
    try {
      validateGraph(next);
      let nextEstimand = 'estimand' in planning ? planning.estimand : estimand;
      if (planning.renameFrom) nextEstimand = renameEstimand(nextEstimand, planning.renameFrom, planning.renameTo);
      let note = '';
      if (nextEstimand) {
        try { validateEstimand(next, nextEstimand); }
        catch { nextEstimand = null; note = ' The power analysis no longer fit the graph and was reset.'; }
      }
      if (!planning.remembered) remember();
      graph = next; estimand = nextEstimand;
      if ('power' in planning) powerSettings = planning.power;
      if (planning.renameFrom && powerSettings?.calibrationEdge) {
        powerSettings = { ...powerSettings, calibrationEdge: powerSettings.calibrationEdge.split('>').map((id) => id === planning.renameFrom ? planning.renameTo : id).join('>') };
      }
      if (powerSettings?.analysis) powerSettings = { ...powerSettings, analysis: reconcileAnalysis(powerSettings.analysis, graph, planning.renameFrom, planning.renameTo) };
      if (!graph.nodes.some((n) => n.id === selected)) selected = graph.nodes.at(-1)?.id;
      if (selectedEdge && !graph.edges.some((e) => `${e.source}>${e.target}` === selectedEdge)) selectedEdge = null;
      connectFrom = null;
      clearResult(); persist();
      if (planning.fit) view = null;
      render();
      status((message ?? 'Model updated. Simulate again to see its data.') + note);
      return true;
    } catch (error) { status(error.message, true); return false; }
  }
  const nodeElement = (id) => [...canvas.querySelectorAll('.dag-node')].find((g) => g.dataset.node === id);
  const edgeElement = (key) => canvas.querySelector(`[data-edge="${key}"]`);
  function updateSelection(id) {
    selected = id;
    draw(); inspector();
  }
  function selectEdge(key) {
    selectedEdge = key; selected = key.split('>')[1];
    draw(); inspector();
    edgeElement(key)?.focus();
  }
  function connectTo(id) {
    if (!connectFrom) {
      connectFrom = id; selectedEdge = null; draw();
      status(`Arrow from ${id}: click or drag to the variable it affects. Escape cancels.`);
    } else {
      const source = connectFrom;
      const check = connectTargets(graph, source).get(id);
      if (check && !check.ok) { status(check.reason, true); return; }
      try {
        if (apply(addEdge(graph, source, id), `Added ${source} → ${id}. Click its number to change the effect.`)) selectEdge(`${source}>${id}`);
      } catch (error) { status(error.message, true); }
    }
  }
  function cancelConnect(message = 'Arrow connection cancelled.') {
    if (!connectFrom && gesture?.type !== 'connect') return;
    connectFrom = null; gesture = null; draw(); status(message);
  }
  function deleteSelection(target) {
    const edgeKey = target?.closest?.('.dag-edge')?.dataset.edge ?? (target?.closest?.('.dag-node') ? null : selectedEdge);
    if (edgeKey) {
      const [source, edgeTarget] = edgeKey.split('>');
      apply(removeEdge(graph, source, edgeTarget), `Removed ${source} → ${edgeTarget} and any interactions using it. Undo restores it.`);
      return;
    }
    const id = target?.closest?.('.dag-node')?.dataset.node ?? selected;
    if (id) apply(removeNode(graph, id), `${id} and its arrows removed. Undo restores them.`);
  }
  // The visible part of the drawing.  Fitting shows every variable and plate; the view otherwise stays where the
  // person left it.
  function fit() { view = fitView(contentBounds(graph, drawn.plates, drawn.effects)); canvas.setAttribute('viewBox', viewBox(view)); }
  function zoom(factor, at) { view = zoomView(view, factor, at); canvas.setAttribute('viewBox', viewBox(view)); }
  function reveal(id) {
    const n = graph.nodes.find((v) => v.id === id);
    if (!n || !view) return;
    const h = view.w * ASPECT;
    if (n.x < view.x || n.y < view.y || n.x + NODE_W > view.x + view.w || n.y + NODE_H > view.y + h) {
      view = { ...view, x: n.x + NODE_W / 2 - view.w / 2, y: n.y + NODE_H / 2 - h / 2 };
      canvas.setAttribute('viewBox', viewBox(view));
    }
  }
  function addVariable(kind, at) {
    if (graph.nodes.length >= MAX_NODES) { status(`Use at most ${MAX_NODES} variables.`, true); return; }
    const prefix = { continuous: 'X', binary: 'B', count: 'N', survival: 'T' }[kind];
    let k = 1, id = `${prefix}${k}`;
    while (graph.nodes.some((n) => n.id === id)) id = `${prefix}${++k}`;
    const centre = at ?? (view ? { x: view.x + view.w / 2 - NODE_W / 2, y: view.y + view.w * ASPECT / 2 - NODE_H / 2 } : { x: 40, y: 40 });
    const spot = at ? { x: Math.max(0, Math.round(at.x)), y: Math.max(0, Math.round(at.y)) } : freeSpot(graph, centre);
    let node = { ...makeNode(id, kind, spot.x, spot.y), label: { continuous: 'Variable', binary: 'Binary variable', count: 'Count', survival: 'Event time' }[kind] };
    let level = '';
    if (graph.study) {
      level = plateAt(drawn.plates, { x: spot.x + NODE_W / 2, y: spot.y + NODE_H / 2 }) ?? 'subject';
      if (kind === 'survival' && level === 'visit') level = 'subject';
      node = { ...node, scope: level, repeated: { subjectSD: 0, slopeSD: 0, slopeCorrelation: 0, clusterSD: 0, residualCorrelation: 0 }, lags: [] };
    }
    if (apply({ ...graph, nodes: [...graph.nodes, node] }, `Added ${id}${level ? `, measured once per ${level}` : ''}. Choose its distribution on the right; drag its + to connect it.`)) {
      updateSelection(id); reveal(id);
    }
  }
  function draw() {
    const W = NODE_W, H = NODE_H;
    drawn = plateLayout(graph);
    const found = analysisPaths();
    biased = found?.open ?? new Set();
    if (!view) fit();
    canvas.setAttribute('viewBox', viewBox(view));
    const marker = (id) => svg('marker', { id, viewBox: '0 0 10 10', refX: 9, refY: 5, markerWidth: 13, markerHeight: 13,
      orient: 'auto-start-reverse', markerUnits: 'userSpaceOnUse' }, svg('path', { d: 'M 0 0 L 10 5 L 0 10 z' }));
    const defs = svg('defs', {}, marker('dag-arrow'), marker('dag-arrow-effect'), marker('dag-arrow-lag'), marker('dag-arrow-bias'));
    const map = new Map(graph.nodes.map((n) => [n.id, n]));
    // the arrows that carry the effect the power analysis tests
    const effect = new Set((estimand ? calibrationEdges(graph, estimand) : []).map((e) => `${e.source}>${e.target}`));
    const connecting = gesture?.type === 'connect' ? gesture.source : connectFrom;
    const targets = connecting ? connectTargets(graph, connecting) : null;
    const roles = estimand ? inferRoles(graph, estimand.exposure, estimand.outcome) : new Map();
    // Plates: the study's levels as nested boxes.  A box dragged into a plate is measured at that level; a plate's
    // label edits the visit schedule or the cluster sizes.
    const plates = svg('g', { class: 'dag-plates' });
    for (const p of drawn.plates) {
      const editable = p.level !== 'subject';
      const label = svg('g', { class: `dag-plate-label${editable ? ' editable' : ''}`, ...(editable ? { role: 'button', tabindex: 0,
        'aria-label': p.level === 'visit' ? `Visit times ${graph.study.times.join(', ')}. Enter to edit.` : `Subjects per cluster ${graph.study.clusterSizes.join(', ')}. Enter to edit.` } : {}) },
      svg('text', { x: p.x + 16, y: p.y + 21 }, p.label));
      if (editable) {
        const edit = () => editPlate(p.level, label);
        label.addEventListener('click', (e) => { e.stopPropagation(); edit(); });
        label.addEventListener('pointerdown', (e) => e.stopPropagation());
        label.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); edit(); } });
      }
      plates.append(svg('g', { class: `dag-plate ${p.level}${p.empty ? ' empty' : ''}`, 'data-plate': p.level },
        svg('rect', { class: 'dag-plate-box', x: p.x, y: p.y, width: p.w, height: p.h, rx: 16 }),
        label, svg('text', { class: 'dag-plate-count', x: p.x + p.w - 14, y: p.y + p.h - 10, 'text-anchor': 'end' }, p.multiplicity),
        ...(p.empty ? [svg('text', { class: 'dag-plate-hint', x: p.x + p.w / 2, y: p.y + p.h / 2 + 14, 'text-anchor': 'middle' },
          p.level === 'visit' ? 'Drag a variable here to measure it at every visit' : 'Drag a variable here')] : [])));
    }
    const lines = svg('g', { class: 'dag-arrows' });
    for (const edge of graph.edges) {
      const key = `${edge.source}>${edge.target}`, { d, mid } = edgePath(map.get(edge.source), map.get(edge.target));
      // the arrow's linear coefficient, and a mark for each response-curve term of the source in the child
      const shapes = [...new Set((map.get(edge.target).functions ?? []).filter((t) => t.left === edge.source || t.right === edge.source)
        .map((t) => t.right ? '×' : TERM_MARK[t.transform?.kind] ?? 'f')
        .concat(map.get(edge.target).unknown?.inputs.includes(edge.source) ? ['?'] : []))];
      const number = short(edge.coefficient) + (shapes.length ? ` + ${shapes.join(', ')}` : ''), width = 7.6 * number.length + 12;
      const isSelected = key === selectedEdge, kind = effect.has(key) ? 'effect' : biased.has(key) ? 'bias' : '';
      const group = svg('g', { class: `dag-edge${edge.target === selected && !selectedEdge ? ' selected' : ''}${isSelected ? ' chosen' : ''}${kind ? ` ${kind}` : ''}`,
        'data-edge': key, role: 'button', tabindex: 0,
        'aria-label': `${edge.source} to ${edge.target}, coefficient ${edge.coefficient}. Enter or click the number to change it; Delete removes the arrow.` },
      svg('path', { class: 'dag-edge-hit', d }),
      svg('path', { class: 'dag-edge-line', d, 'marker-end': `url(#dag-arrow${kind ? `-${kind}` : ''})` }),
      // the number sits on a backing so it reads over the curve and a click anywhere on it edits it
      svg('g', { class: 'dag-edge-label' },
        svg('rect', { x: mid[0] - width / 2, y: mid[1] - 22, width, height: 19, rx: 4 }),
        svg('text', { x: mid[0], y: mid[1] - 8, 'text-anchor': 'middle' }, number)));
      if (isSelected) {
        const remove = svg('g', { class: 'dag-edge-remove', role: 'button', tabindex: 0, 'aria-label': `Remove the arrow ${edge.source} to ${edge.target}` },
          svg('circle', { cx: mid[0] + width / 2 + 13, cy: mid[1] - 12.5, r: 10 }), svg('text', { x: mid[0] + width / 2 + 13, y: mid[1] - 8, 'text-anchor': 'middle' }, '×'));
        const go = (e) => { e.stopPropagation(); e.preventDefault(); deleteSelection(group); };
        remove.addEventListener('click', go);
        remove.addEventListener('pointerdown', (e) => e.stopPropagation());
        remove.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') go(e); });
        group.append(remove);
      }
      group.addEventListener('pointerdown', (e) => e.stopPropagation());
      group.addEventListener('click', (e) => {
        if (e.target.closest('.dag-edge-remove')) return;
        const label = e.target.closest('.dag-edge-label');
        if (label) editCoefficient(edge, label);
        else selectEdge(key);
      });
      group.addEventListener('keydown', (e) => {
        if (e.target !== group) return;
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); editCoefficient(edge, group.querySelector('.dag-edge-label')); }
        if (e.key === 'Delete' || e.key === 'Backspace') { e.preventDefault(); deleteSelection(group); }
      });
      group.addEventListener('focus', () => { if (selectedEdge !== key) { selectedEdge = key; group.classList.add('chosen'); } });
      lines.append(group);
    }
    // Lagged relationships (a variable at an earlier visit): dashed, beside any same-visit arrow, a loop for the
    // variable's own past.  Their settings are on the child.
    const lags = svg('g', { class: 'dag-lags' });
    if (graph.study) for (const node of graph.nodes) for (const l of node.lags ?? []) {
      const source = map.get(l.source);
      if (!source) continue;
      const same = graph.edges.some((e) => e.source === l.source && e.target === node.id);
      const { d, mid } = l.source === node.id ? loopPath(node) : edgePath(source, node, same ? 46 : 0);
      const text = `lag ${l.lag} · ${short(l.coefficient)}`, width = 6.4 * text.length + 12;
      const group = svg('g', { class: 'dag-lag', role: 'button', tabindex: 0,
        'aria-label': `${l.source} ${l.lag} visit${l.lag === 1 ? '' : 's'} earlier affects ${node.id}, coefficient ${l.coefficient}. Enter to edit.` },
      svg('path', { d, 'marker-end': 'url(#dag-arrow-lag)' }),
      svg('rect', { x: mid[0] - width / 2, y: mid[1] - 10, width, height: 18, rx: 9 }),
      svg('text', { x: mid[0], y: mid[1] + 3, 'text-anchor': 'middle' }, text));
      const open = () => openSection(node.id, 'dag-node-lags');
      group.addEventListener('pointerdown', (e) => e.stopPropagation());
      group.addEventListener('click', open);
      group.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(); } });
      lags.append(group);
    }
    // Shared random effects as small latent circles on the plate they vary over.
    const effects = svg('g', { class: 'dag-effects' });
    for (const re of drawn.effects) {
      const target = map.get(re.node);
      const group = svg('g', { class: `dag-random ${re.level}`, role: 'button', tabindex: 0, 'aria-label': `${re.text}. Enter to edit.` },
        svg('title', {}, re.text),
        svg('path', { d: `M ${re.cx} ${re.cy + 14} L ${target.x + W / 2 + (re.cx - target.x - W / 2) * 0.35} ${target.y - 2}`, 'marker-end': 'url(#dag-arrow-lag)' }),
        svg('circle', { cx: re.cx, cy: re.cy, r: 14 }),
        svg('text', { x: re.cx, y: re.cy + 4.5, 'text-anchor': 'middle' }, re.symbol));
      const open = () => openSection(re.node, 'dag-node-effects', `repeat-${re.key}`);
      group.addEventListener('pointerdown', (e) => e.stopPropagation());
      group.addEventListener('click', open);
      group.addEventListener('keydown', (e) => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(); } });
      effects.append(group);
    }
    const nodes = svg('g', {});
    for (const node of graph.nodes) {
      const target = targets?.get(node.id);
      const r = node.repeated ?? {};
      const sketchKey = `${node.id}|${node.kind}`, values = preview.columns?.[node.id];
      if (values && !preview.sketches.has(sketchKey)) preview.sketches.set(sketchKey, sketch(node.kind, values, SKETCH_W, 30));
      const sketched = values && preview.sketches.get(sketchKey);
      const inferred = roles.get(node.id), shownRole = ['confounder', 'mediator', 'collider'].includes(inferred) ? inferred : null;
      // A box says what the variable is and how it is drawn, in two lines; everything else is in its tooltip,
      // the variable settings and the equations.
      const tip = [`${node.id} · ${node.label}${node.measured === false ? ' (unmeasured: no analysis can use it)' : ''}`, `${kindNames[node.kind]}: ${summary(graph, node)}`,
        ...(node.kind === 'survival' && node.censoring.mode !== 'none' ? [node.censoring.mode === 'variable' ? `Censored by ${node.censoring.variable}` : `Follow-up ends at ${short(node.censoring.time)}`] : []),
        ...(r.residualCorrelation ? [`Residual correlation ${short(r.residualCorrelation)} per unit time`] : []),
        ...(node.unknown ? [`Unknown form: f(${node.unknown.inputs.join(', ')}), SD ${short(node.unknown.size)}, smoothness ${short(node.unknown.smoothness)}`] : []),
        ...(inferred ? [`${inferred[0].toUpperCase()}${inferred.slice(1)}: ${roleHelp[inferred].toLowerCase()}`] : []),
        ...(sketched ? [sketched.title] : [])].join('\n');
      const group = svg('g', { class: `dag-node ${node.kind}${node.measured === false ? ' unmeasured' : ''}${node.id === selected ? ' selected' : ''}${node.id === connecting ? ' connecting' : ''}`
        + `${target && node.id !== connecting ? target.ok ? ' target-ok' : ' target-bad' : ''}`
        + `${!graph.edges.some((e) => e.source === node.id || e.target === node.id) && graph.nodes.length > 1 ? ' isolated' : ''}`,
        'data-node': node.id, transform: `translate(${node.x} ${node.y})`, tabindex: 0, role: 'button',
        'aria-label': `${node.id}, ${node.label}, ${kindNames[node.kind]}${graph.study ? `, measured ${{ cluster: 'once per cluster', subject: 'once per subject', visit: 'at every visit', time: 'as the visit time' }[scopeOf(node)]}` : ''}. Enter to edit, arrow keys to move, C to compare it, O to make it the outcome, Delete to remove it.` },
      svg('title', {}, target && !target.ok && node.id !== connecting ? target.reason : tip),
      svg('rect', { class: 'dag-box', width: W, height: H, rx: 10 }),
      svg('rect', { class: 'dag-stripe', width: 4, height: H - 22, x: 0, y: 11, rx: 2 }),
      svg('text', { x: 14, y: 25, class: 'dag-node-title' }, node.label.length > 17 ? node.label.slice(0, 16) + '…' : node.label),
      svg('text', { x: 14, y: 43, class: 'dag-node-family' }, `${node.id} · ${compact(graph, node)}`));
      const role = estimand?.exposure === node.id ? (estimand.kind === 'coefficient' ? 'Predictor' : 'Compared')
        : estimand?.outcome === node.id ? 'Outcome' : null;
      if (role) group.append(svg('g', { class: 'dag-node-tag' }, svg('rect', { x: 10, y: -9, width: 66, height: 17, rx: 8.5 }),
        svg('text', { x: 43, y: 3.5, 'text-anchor': 'middle' }, role)));
      if (shownRole) {
        const width = 6.2 * shownRole.length + 14;
        group.append(svg('g', { class: `dag-node-inferred ${shownRole}` },
          svg('rect', { x: W - width - 8, y: -8, width, height: 16, rx: 8 }), svg('text', { x: W - 8 - width / 2, y: 3.5, 'text-anchor': 'middle' }, shownRole)));
      }
      const port = svg('g', { class: 'dag-port', role: 'button', tabindex: 0, 'aria-label': `Start an arrow from ${node.id}` },
        svg('title', {}, `Drag to another variable to draw an arrow from ${node.id}`),
        svg('circle', { cx: W, cy: H / 2, r: 9 }), svg('text', { x: W, y: H / 2 + 4.5, 'text-anchor': 'middle' }, '+'));
      port.addEventListener('pointerdown', (e) => {
        e.stopPropagation();
        if (e.button !== 0) return;
        gesture = { type: 'connect', source: node.id, start: point(e), moved: false };
      });
      // released on the port itself: a click (which starts click-to-connect) or a drag that came back (cancelled)
      port.addEventListener('pointerup', (e) => {
        e.stopPropagation();
        if (gesture?.type !== 'connect') return;
        if (gesture.moved) { gesture = { suppressClick: true }; pointer = null; draw(); status('Arrow cancelled: release over another variable to connect them.'); }
        else gesture = null;
      });
      port.addEventListener('click', (e) => {
        e.stopPropagation();
        if (gesture?.suppressClick) { gesture = null; return; }
        connectTo(node.id);
      });
      port.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' || e.key === ' ') { e.stopPropagation(); e.preventDefault(); connectTo(node.id); }
      });
      group.append(port);
      let drag = null;
      group.addEventListener('pointerdown', (event) => {
        event.stopPropagation();
        if (event.button !== 0 || connectFrom) return;
        const p = point(event);
        drag = { start: p, x: node.x, y: node.y, moved: false, plates: drawn.plates };
        group.setPointerCapture(event.pointerId);
      });
      group.addEventListener('pointermove', (event) => {
        if (!drag) return;
        const p = point(event), dx = p.x - drag.start.x, dy = p.y - drag.start.y;
        if (!drag.moved && Math.hypot(dx, dy) < 4) return;
        if (!drag.moved) { remember(); drag.moved = true; group.classList.add('dragging'); }
        node.x = Math.max(0, Math.min(CANVAS_LIMIT.x, Math.round(drag.x + dx)));
        node.y = Math.max(0, Math.min(CANVAS_LIMIT.y, Math.round(drag.y + dy)));
        group.setAttribute('transform', `translate(${node.x} ${node.y})`);
        if (graph.study) {
          const level = plateAt(drag.plates, { x: node.x + W / 2, y: node.y + H / 2 });
          for (const g of canvas.querySelectorAll('.dag-plate')) g.classList.toggle('drop', g.dataset.plate === level && level !== levelName(node));
        }
      });
      const finish = () => {
        const moved = drag?.moved, plates = drag?.plates;
        drag = null;
        if (!moved) return false;
        selected = node.id; selectedEdge = null;  // the box just moved is the one being worked on
        const level = graph.study ? plateAt(plates, { x: node.x + W / 2, y: node.y + H / 2 }) : null;
        if (level && level !== levelName(node) && scopeOf(node) !== 'time') {
          const word = { cluster: 'once per cluster', subject: 'once per subject', visit: 'at every visit' }[level];
          if (!apply(changeScope(graph, node.id, level), `${node.id} is now measured ${word}. Undo restores it.`, { remembered: true })) { persist(); updateSelection(node.id); }
        } else { persist(); updateSelection(node.id); }
        buttons();
        return true;
      };
      group.addEventListener('pointercancel', finish);
      group.addEventListener('pointerup', (event) => {
        if (gesture?.type === 'connect') return;  // an arrow dragged onto this box: the canvas connects it
        if (finish()) { event.stopPropagation(); return; }
        if (connectFrom) connectTo(node.id);
        else { selectedEdge = null; updateSelection(node.id); }
      });
      group.addEventListener('dblclick', (event) => { event.stopPropagation(); updateSelection(node.id); $('dag-node-label')?.focus(); });
      group.addEventListener('contextmenu', (event) => { event.preventDefault(); openMenu(node, event.clientX, event.clientY); });
      group.addEventListener('keydown', (event) => {
        if (event.target !== group) return;
        const role = { c: 'exposure', o: 'outcome' }[event.key.toLowerCase()];
        if (role && !event.ctrlKey && !event.metaKey && !event.altKey) {
          event.preventDefault(); assignRole(role, node.id);
          nodeElement(node.id)?.focus();
          return;
        }
        if (event.key === 'ContextMenu' || (event.shiftKey && event.key === 'F10')) {
          event.preventDefault(); const box = group.getBoundingClientRect(); openMenu(node, box.left + 20, box.bottom - 10); return;
        }
        if (event.key === 'Delete' || event.key === 'Backspace') { event.preventDefault(); deleteSelection(group); return; }
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault();
          if (connectFrom) connectTo(node.id); else { updateSelection(node.id); $('dag-node-label').focus(); }
        }
        const delta = { ArrowLeft: [-10, 0], ArrowRight: [10, 0], ArrowUp: [0, -10], ArrowDown: [0, 10] }[event.key];
        if (delta) {
          event.preventDefault(); remember();
          node.x = Math.max(0, Math.min(CANVAS_LIMIT.x, node.x + delta[0]));
          node.y = Math.max(0, Math.min(CANVAS_LIMIT.y, node.y + delta[1]));
          persist(); draw(); buttons();
          nodeElement(node.id)?.focus();
        }
      });
      nodes.append(group);
    }
    // the arrow being drawn, from the source's port to the pointer
    const overlay = svg('g', { class: 'dag-overlay' });
    if (connecting && pointer && map.has(connecting)) {
      const s = map.get(connecting);
      overlay.append(svg('path', { class: 'dag-rubber', d: `M ${s.x + W} ${s.y + H / 2} L ${pointer.x} ${pointer.y}`, 'marker-end': 'url(#dag-arrow)' }));
    }
    if (!graph.nodes.length) nodes.append(svg('text', { x: view.x + view.w / 2, y: view.y + view.w * ASPECT / 2, 'text-anchor': 'middle', class: 'dag-empty' },
      'Choose a starting point or add a variable above, or double-click here'));
    canvas.replaceChildren(svg('title', {}, 'Data generating process'),
      svg('desc', {}, 'Select a box to edit it. Drag boxes to move them, or into a plate to change how often they are measured. Drag from a box’s plus port to another box to draw an arrow.'),
      defs, plates, effects, lines, lags, nodes, overlay);
    $('dag-node-count').replaceChildren(studySummary(graph),
      ...(effect.size ? [el('span', { class: 'dag-legend-effect' }, 'tested effect')] : []),
      ...(biased.size ? [el('span', { class: 'dag-legend-bias' }, 'open non-causal path')] : []));
    $('dag-node-select').replaceChildren(...graph.nodes.map((n) => el('option', { value: n.id }, `${n.id} · ${n.label}`)));
    $('dag-node-select').value = selected ?? '';
    $('dag-node-select').disabled = !graph.nodes.length;
    $('dag-tidy').disabled = graph.nodes.length < 2;
    $('dag-clear').disabled = !graph.nodes.length;
    $('dag-level-clusters').setAttribute('aria-pressed', String(Boolean(graph.study?.clusterSizes?.length)));
    $('dag-level-visits').setAttribute('aria-pressed', String((graph.study?.times?.length ?? 0) > 1));
    for (const [role, id] of [['exposure', 'dag-role-exposure'], ['outcome', 'dag-role-outcome']])
      $(id)?.setAttribute('aria-pressed', String(estimand?.[role] === selected));
    paths(found);
  }
  const levelName = (node) => ({ cluster: 'cluster', subject: 'subject', visit: 'visit', time: 'visit' })[scopeOf(node)];
  // The paths between the compared variable and the outcome, and whether the planned analysis leaves each open.
  // Also returns the arrows on open non-causal paths, which the canvas marks.
  // the paths depend on the model's structure, the question and the analysis, not on where boxes are drawn: the
  // last answer is kept until one of those changes
  let pathsMemo = { key: null, value: null };
  function analysisPaths() {
    const key = JSON.stringify([modelKey(), estimand, powerSettings?.analysis ?? null]);
    if (key !== pathsMemo.key) pathsMemo = { key, value: computePaths() };
    return pathsMemo.value;
  }
  function computePaths() {
    const a = estimand?.exposure, y = estimand?.outcome;
    if (graph.study || !a || !y || !graph.nodes.some((n) => n.id === a) || !graph.nodes.some((n) => n.id === y)) return null;
    const plan = powerCompatibility(graph, estimand, powerSettings?.analysis);
    // the g-formula redraws the variables it refits (mediators and the outcome), so only the others are held fixed
    const terms = !plan.supported ? [] : plan.method === 'gformula' ? plan.models.flatMap((m) => m.terms) : plan.terms;
    const redrawn = new Set([a, ...(plan.models ?? []).map((m) => m.node)]);
    const adjusted = new Set(terms.flatMap((t) => [t.left, t.right]).filter((v) => v && !redrawn.has(v)));
    const { paths: list, truncated } = findPaths(graph, a, y);
    const rows = list.map((path) => ({ path, c: classifyPath(graph, path, adjusted),
      keys: path.map((s) => s.forward ? `${s.from}>${s.to}` : `${s.to}>${s.from}`) }));
    const open = new Set(rows.filter(({ c }) => c.open && c.kind !== 'causal').flatMap(({ keys }) => keys));
    return { a, y, adjusted, rows, truncated, open, plan };
  }
  function paths(found = analysisPaths()) {
    const box = $('dag-paths');
    issues(found);
    if (!found) { box.replaceChildren(); return; }
    const { a, y, adjusted, rows, truncated } = found;
    const mark = (keys) => { for (const g of canvas.querySelectorAll('.dag-edge')) g.classList.toggle('path', keys.includes(g.dataset.edge)); };
    const items = rows.map(({ path, c, keys }) => {
      const what = c.kind === 'causal' ? 'causal' : c.kind === 'backdoor' ? 'back-door' : 'non-causal';
      const why = c.open ? (c.kind === 'causal' ? 'open: carries effect' : 'open: biases the estimate')
        : c.collider ? `blocked at collider ${c.by}` : `blocked: ${c.by} is adjusted for`;
      const row = el('li', { class: `${c.kind}${c.open ? ' open' : ''}`, tabindex: 0 },
        el('span', { class: 'dag-path-text' }, pathText(path)), el('span', { class: 'dag-path-kind' }, `${what} · ${why}`));
      const show = (on) => mark(on ? keys : []);
      row.addEventListener('mouseenter', () => show(true)); row.addEventListener('mouseleave', () => show(false));
      row.addEventListener('focus', () => show(true)); row.addEventListener('blur', () => show(false));
      return { row, c };
    });
    const count = items.filter(({ c }) => c.open && c.kind !== 'causal').length;
    const heading = `Paths from ${a} to ${y}` + (adjusted.size ? ` · analysis adjusts for ${[...adjusted].join(', ')}` : '');
    box.replaceChildren(el('h4', {}, heading),
      items.length ? el('ul', {}, items.map(({ row }) => row)) : el('p', { class: 'small muted' }, `No path connects ${a} and ${y}: the effect is zero.`),
      ...(truncated ? [el('p', { class: 'small muted' }, `The first ${rows.length} paths are shown.`)] : []),
      ...(count ? [el('p', { class: 'small error' }, `${count} open non-causal path${count === 1 ? '' : 's'}: the estimate would mix confounding into the effect.`)] : []));
  }
  // What needs attention, shown over the canvas: no comparison chosen yet, a question the planned analysis cannot
  // answer, open non-causal paths.
  function issues(found) {
    const list = [];
    if (graph.nodes.length >= 2 && (!estimand?.exposure || !estimand?.outcome))
      list.push(['info', 'Choose the comparison: select a box and press C for the compared variable, O for the outcome (or right-click it).']);
    else if (estimand) {
      const plan = found?.plan ?? powerCompatibility(graph, estimand, powerSettings?.analysis);
      if (!plan.supported) list.push(['warn', plan.reason]);
      const open = found ? found.rows.filter(({ c }) => c.open && c.kind !== 'causal').length : 0;
      if (open) list.push(['warn', `${open} open non-causal path${open === 1 ? '' : 's'} from ${found.a} to ${found.y} (orange): the estimate would mix in confounding. See the paths below the canvas.`]);
    }
    const isolated = graph.nodes.length > 1 ? graph.nodes.filter((n) => !graph.edges.some((e) => e.source === n.id || e.target === n.id)
      && !(n.lags ?? []).length && scopeOf(n) !== 'time') : [];
    if (isolated.length && list.length < 2) list.push(['info', `${isolated.map((n) => n.id).join(', ')} ${isolated.length === 1 ? 'has' : 'have'} no arrows: generated independently of everything else.`]);
    $('dag-issues').replaceChildren(...list.map(([kind, text]) => el('p', { class: `dag-issue ${kind}` }, text)));
    $('dag-issues').hidden = !list.length;
  }
  function assignRole(role, id) {
    // a role set from the graph goes through the power form, which keeps the question consistent
    if (!powerEditor?.assign(role, id)) status('The power analysis needs at least two variables.', true);
    else status(`${id} is now the ${role === 'exposure' ? 'compared variable' : 'outcome'} in the power analysis.`);
  }
  function closeMenu() { menu?.remove(); menu = null; }
  function openMenu(node, x, y) {
    closeMenu();
    const item = (text, action) => el('button', { type: 'button', role: 'menuitem', onclick: () => { closeMenu(); action(); } }, text);
    const level = (scope, text) => graph.study && scopeOf(node) !== 'time' && scopeOf(node) !== scope && (scope !== 'cluster' || graph.study.clusterSizes?.length)
      ? [item(text, () => apply(changeScope(graph, node.id, scope), `${node.id} is now measured ${text.toLowerCase().replace('measure ', '')}. Undo restores it.`))] : [];
    menu = el('div', { class: 'dag-menu', role: 'menu', 'aria-label': `${node.id} actions` },
      item('Compare this variable', () => assignRole('exposure', node.id)),
      item('Use as the outcome', () => assignRole('outcome', node.id)),
      item('Edit', () => { updateSelection(node.id); $('dag-node-label').focus(); }),
      item('Draw an arrow from it', () => connectTo(node.id)),
      ...level('cluster', 'Measure once per cluster'), ...level('subject', 'Measure once per subject'), ...level('visit', 'Measure at every visit'),
      ...(graph.study && scopeOf(node) === 'time' ? [] : [item(node.measured === false ? 'Mark as measured' : 'Mark as unmeasured', () => {
        const { measured: _, ...rest } = node;
        apply({ ...graph, nodes: graph.nodes.map((n) => n.id === node.id ? (node.measured === false ? rest : { ...rest, measured: false }) : n) },
          node.measured === false ? `${node.id} is measured: analyses can use it.` : `${node.id} is unmeasured: it is simulated, but no analysis can use it.`);
      })]),
      item('Delete', () => apply(removeNode(graph, node.id), `${node.id} and its arrows removed. Undo restores them.`)));
    document.body.append(menu);
    menu.style.left = `${Math.min(x, innerWidth - menu.offsetWidth - 8) + scrollX}px`;
    menu.style.top = `${Math.min(y, innerHeight - menu.offsetHeight - 8) + scrollY}px`;
    menu.querySelector('button').focus();
    menu.addEventListener('keydown', (e) => {
      const items = [...menu.querySelectorAll('button')], at = items.indexOf(document.activeElement);
      if (e.key === 'Escape') { closeMenu(); e.stopPropagation(); }
      if (e.key === 'ArrowDown' || e.key === 'ArrowUp') { e.preventDefault(); items[(at + (e.key === 'ArrowDown' ? 1 : items.length - 1)) % items.length].focus(); }
    });
  }
  // A small field over something drawn on the canvas (an arrow's number, a plate's label), saved on Enter or when
  // it loses focus, cancelled by Escape.
  function inlineEdit(anchor, { value, label, width = 72, type = 'number', onSave, refocus }) {
    const host = canvas.parentElement, at = anchor.getBoundingClientRect(), base = host.getBoundingClientRect();
    const input = el('input', { type, ...(type === 'number' ? { step: 'any' } : {}), value, class: 'dag-inline', 'aria-label': label });
    input.style.width = `${width}px`;
    host.append(input);
    input.style.left = `${at.left - base.left + host.scrollLeft + Math.max(0, at.width / 2 - width / 2)}px`;
    input.style.top = `${at.top - base.top + host.scrollTop - 6}px`;
    let done = false;
    const finish = (save) => {
      if (done) return; done = true; input.remove();
      if (save && input.value.trim() !== '' && input.value !== String(value)) onSave(input.value);
      refocus?.()?.focus();
    };
    input.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') { e.preventDefault(); finish(true); }
      if (e.key === 'Escape') { e.preventDefault(); e.stopPropagation(); finish(false); }
    });
    input.addEventListener('blur', () => finish(true));
    input.focus(); input.select();
  }
  // A coefficient changed where it is drawn: a small field over the arrow's number.
  function editCoefficient(edge, label) {
    inlineEdit(label, { value: edge.coefficient, label: `Coefficient of ${edge.source} → ${edge.target}`,
      refocus: () => edgeElement(`${edge.source}>${edge.target}`),
      onSave(text) {
        try { apply({ ...graph, edges: graph.edges.map((e) => e === edge ? { ...e, coefficient: inputNumber({ value: text }, 'Arrow coefficient') } : e) },
          `${edge.source} → ${edge.target} set to ${text}.`); }
        catch (error) { status(error.message, true); }
      } });
  }
  // The visit schedule or the cluster sizes, edited on the plate's label.
  function editPlate(level, label) {
    const key = level === 'visit' ? 'times' : 'clusterSizes', current = graph.study[key].join(', ');
    inlineEdit(label, { value: current, type: 'text', width: 180,
      label: level === 'visit' ? 'Visit times, separated by commas' : 'Subjects per cluster: one size, or a repeating pattern',
      refocus: () => canvas.querySelector(`.dag-plate.${level} .dag-plate-label`),
      onSave(text) {
        const values = text.trim().split(/[\s,]+/).filter(Boolean).map(Number);
        try {
          const study = parseStudy({ ...graph.study, [key]: values });
          apply({ ...graph, study }, level === 'visit' ? `Visits at ${study.times.join(', ')}.` : `Subjects per cluster: ${study.clusterSizes.join(', ')}.`);
        } catch (error) { status(error.message, true); }
      } });
  }
  // Selects a variable and opens one part of its settings (its random effects or lags), focusing a field there.
  function openSection(id, section, focus) {
    selectedEdge = null; updateSelection(id);
    const details = $(section);
    if (details) { details.hidden = false; details.open = true; (focus ? $(focus) : details.querySelector('summary'))?.focus(); details.scrollIntoView?.({ block: 'nearest' }); }
  }
  // Draws a fresh sample from the model for the boxes, a moment after the last change.
  function schedulePreview() {
    clearTimeout(preview.timer);
    preview.timer = setTimeout(async () => {
      const token = ++preview.token, key = modelKey();
      if (!graph.nodes.length) { preview.columns = null; preview.sketches.clear(); draw(); return; }
      const instance = borrowEngine();
      try {
        const data = await instance.dag({ op: 'simulate', graph, options: { n: graph.study ? Math.max(1, Math.floor(PREVIEW_N / (graph.study.times.length * Math.max(1, ...(graph.study.clusterSizes ?? []))))) : PREVIEW_N, seed: 1 } });
        if (token !== preview.token) return;
        preview.columns = data.columns; preview.key = key; preview.sketches.clear();
      } catch (error) {
        if (token !== preview.token || error.stopped) return;
        preview.columns = null; preview.key = key; preview.sketches.clear();
      } finally {
        returnEngine(instance);
      }
      draw(); showSketch();
    }, 250);
  }
  function point(event) {
    const p = canvas.createSVGPoint(); p.x = event.clientX; p.y = event.clientY;
    return p.matrixTransform(canvas.getScreenCTM().inverse());
  }
  // The selected variable's simulated distribution, under its settings; redrawn when a new preview arrives.
  function showSketch() {
    const figure = $('dag-node-sketch'), node = graph.nodes.find((n) => n.id === selected);
    if (!figure || !node) return;
    const key = `${node.id}|${node.kind}`, values = preview.columns?.[node.id];
    if (values && !preview.sketches.has(key)) preview.sketches.set(key, sketch(node.kind, values, SKETCH_W, 30));
    const s = values && preview.sketches.get(key);
    if (!s) {
      figure.replaceChildren(el('figcaption', { class: 'small muted' }, graph.study && scopeOf(node) === 'time' ? 'Takes the visit times.'
        : preview.key === modelKey() && !preview.columns ? 'Its values cannot be simulated until the model is valid.' : 'Simulating its distribution…'));
      return;
    }
    figure.classList.toggle('stale', preview.key !== modelKey());
    figure.replaceChildren(svg('svg', { viewBox: `0 -2 ${SKETCH_W} 36`, width: SKETCH_W, height: 36, role: 'img', 'aria-label': s.title },
      svg('line', { x1: 0, x2: SKETCH_W, y1: 30.5, y2: 30.5 }),
      ...(s.bars ? [svg('path', { class: 'bars', d: s.bars })] : []), ...(s.line ? [svg('path', { class: 'curve', d: s.line })] : [])),
    el('figcaption', { class: 'small muted' }, s.title));
  }
  function inspector() {
    const node = graph.nodes.find((n) => n.id === selected);
    const box = $('dag-inspector'); box.replaceChildren();
    if (!node) { box.append(el('p', { class: 'muted small' }, 'Add a variable, then choose its type and distribution here.')); return; }
    box.append(el('div', { class: 'dag-inspector-heading' }, el('h3', {}, node.id),
      el('span', { class: 'small muted' }, kindNames[node.kind] + (graph.study ? ` · ${scopeOf(node)}` : ''))),
    el('div', { class: 'dag-roles', role: 'group', 'aria-label': 'Role in the power analysis' },
      el('button', { id: 'dag-role-exposure', type: 'button', 'aria-pressed': String(estimand?.exposure === node.id), onclick: () => assignRole('exposure', node.id) }, 'Compared'),
      el('button', { id: 'dag-role-outcome', type: 'button', 'aria-pressed': String(estimand?.outcome === node.id), onclick: () => assignRole('outcome', node.id) }, 'Outcome')),
    el('figure', { id: 'dag-node-sketch', class: 'dag-sketch' }));
    showSketch();
    const id = el('input', { id: 'dag-node-id', type: 'text', maxlength: 32, value: node.id, required: true });
    const label = el('input', { id: 'dag-node-label', type: 'text', maxlength: 80, value: node.label, required: true });
    const kind = select('dag-node-kind', KINDS.map((k) => [k, kindNames[k]]), node.kind);
    const family = select('dag-node-family', [], node.family);
    const intercept = numeric('dag-node-intercept', node.intercept);
    // shown to six digits; the intercept is kept exactly unless the probability is edited
    const shownProbability = String(Number((1 / (1 + Math.exp(-node.intercept))).toPrecision(6)));
    const probability = numeric('dag-node-probability', shownProbability, { min: 0, max: 1 });
    const inputs = Object.fromEntries(['sd', 'width', 'shape', 'scale', 'df', 'precision', 'overdispersion'].map((key) => [key, numeric(`dag-node-${key}`, parameterValue(node, key), { min: 0 })]));
    const measured = el('input', { id: 'dag-node-measured', type: 'checkbox', checked: node.measured !== false });
    const censorMode = select('dag-censor-mode', [['none', 'None · fully observed'], ['administrative', 'Administrative cutoff'], ['variable', 'Use another variable']], node.censoring.mode);
    const censorTime = numeric('dag-censor-time', node.censoring.time, { min: 0 });
    const censorVariable = select('dag-censor-variable', graph.nodes.filter((n) => n.id !== node.id).map((n) => [n.id, `${n.id} · ${n.label}`]), node.censoring.variable);
    const form = el('form', { novalidate: true }, el('div', { class: 'dag-pair' }, field('ID', id), field('Label', label)),
      el('div', { class: 'dag-pair' }, field('Type', kind), field('Distribution', family)),
      ...(graph.study && scopeOf(node) === 'time' ? [] : [el('div', { class: 'dag-field dag-check' }, el('label', {}, measured, ' Measured in the study'),
        el('p', { class: 'small muted' }, 'Clear it for a variable that is simulated but that no analysis can use: an unmeasured confounder, the true value behind an error-prone measurement, or a hidden class.'))]));
    const parameters = el('div'), observation = el('div');
    function observationFields() {
      observation.replaceChildren();
      if (kind.value !== 'survival') return;
      observation.append(field('Censoring', censorMode, 'Adds observed-time and event columns.'));
      if (censorMode.value === 'administrative') observation.append(field('Administrative cutoff', censorTime));
      if (censorMode.value === 'variable') observation.append(field('Censoring time variable', censorVariable, 'Its values must be nonnegative.'));
    }
    function params() {
      const spec = DISTRIBUTIONS[family.value];
      parameters.replaceChildren(family.value === 'bernoulli'
        ? field('Baseline probability (before parent terms)', probability, 'Arrows add to its log odds.')
        : field(`Intercept (${spec.link})`, intercept, 'Baseline predictor; arrows and response terms add to it.'));
      for (const key of spec.parameters) parameters.append(field(
        key === 'sd' && family.value === 'lognormal' ? 'Standard deviation of log values' : parameterLabels[key], inputs[key]));
      if (kind.value === 'survival') parameters.append(el('p', { class: 'small muted' }, 'Arrow coefficients are log hazard ratios.'));
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
    const studyFields = nodeStudyFields(graph, node, () => setTimeout(() => save(true)));
    form.append(studyFields.element, parameters, observation, el('p', { class: 'small muted dag-autosave' }, 'Changes apply as you make them; press Enter to apply a typed value at once.'));
    form.addEventListener('keydown', (event) => { if (event.key === 'Enter' && event.target.matches('input')) { event.preventDefault(); save(false); } });
    // Every change applies at once (after focus has moved, so the field the person went to keeps the focus).
    form.addEventListener('change', (event) => { if (event.target.closest('.dag-field, .dag-grid, details')) setTimeout(() => save(true)); });
    form.addEventListener('submit', (event) => { event.preventDefault(); save(false); });
    function save(auto) {
      if (!form.isConnected) return;
      const focus = auto ? document.activeElement?.id : null;
      try {
        let baseline = inputNumber(intercept, 'Intercept');
        if (family.value === 'bernoulli' && !(node.family === 'bernoulli' && probability.value === shownProbability)) {
          const p = inputNumber(probability, 'Success probability');
          if (!(p > 0 && p < 1)) throw new Error('Success probability must be strictly between 0 and 1.');
          baseline = Math.log(p) - Math.log1p(-p);
        }
        const { measured: _, ...rest } = node;
        const changed = { ...rest, ...(measured.checked ? {} : { measured: false }), ...studyFields.read(), label: label.value.trim(), kind: kind.value, family: family.value,
          intercept: baseline,
          censoring: kind.value === 'survival' ? { mode: censorMode.value,
            time: censorMode.value === 'administrative' ? inputNumber(censorTime, 'Administrative cutoff') : 20,
            variable: censorMode.value === 'variable' ? censorVariable.value : '' } : { mode: 'none', time: 20, variable: '' } };
        for (const key of DISTRIBUTIONS[family.value].parameters) changed[key] = inputNumber(inputs[key], parameterLabels[key]);
        const next = renameNode({ ...graph, nodes: graph.nodes.map((n) => n.id === selected ? changed : n) }, node.id, id.value.trim());
        if (JSON.stringify(next) === JSON.stringify(graph)) return;
        // select the (possibly renamed) variable before the page redraws, so it draws once
        const before = selected; selected = id.value.trim();
        if (apply(next, `${id.value.trim()} updated.`, { renameFrom: node.id, renameTo: id.value.trim() })) { if (focus) $(focus)?.focus(); }
        else selected = before;
      } catch (error) { status(error.message, true); }
    }
    box.append(form);
    const edges = graph.edges.filter((e) => e.target === node.id);
    box.append(el('h4', {}, 'Incoming arrows'), el('p', { class: 'small muted' }, `Coefficients on the ${DISTRIBUTIONS[node.family].link} scale.`));
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
    if (!edges.length) box.append(el('p', { class: 'small muted' }, 'None: generated independently.'));
    // possible parents: not the variable itself, not one already pointing at it, and nothing it affects (no cycles)
    const affected = new Set([node.id]);
    for (const id of affected) for (const e of graph.edges) if (e.source === id) affected.add(e.target);
    const candidates = graph.nodes.filter((n) => !affected.has(n.id) && !edges.some((e) => e.source === n.id));
    if (candidates.length) {
      const source = select('dag-source', candidates.map((n) => [n.id, `${n.id} · ${n.label}`]), candidates[0].id);
      const coefficient = numeric('dag-coefficient', 0.5, { 'aria-label': `Coefficient of the new arrow into ${node.id}` });
      const add = el('form', { class: 'dag-add-arrow', novalidate: true }, el('label', { for: 'dag-source' }, `New arrow into ${node.id} from`),
        source, coefficient, el('button', { type: 'submit', class: 'secondary', 'aria-label': 'Add arrow' }, 'Add'));
      add.addEventListener('submit', (event) => {
        event.preventDefault();
        try { apply(addEdge(graph, source.value, node.id, inputNumber(coefficient, 'Arrow coefficient')), `Added ${source.value} → ${node.id}.`); }
        catch (error) { status(error.message, true); }
      });
      box.append(add);
    }
    if (edges.length >= 2 || node.interactions.length) box.append(el('h4', {}, 'Interactions'));
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
    }
    // A relationship of unknown form: parents that enter a smooth function no analysis is told about.
    const parentIds = graph.edges.filter((e) => e.target === node.id).map((e) => e.source);
    if (parentIds.length && !(graph.study && scopeOf(node) === 'time')) {
      const u = node.unknown;
      // a new relationship gets a shape number no other variable uses, so two variables never share a shape
      const freshSeed = 1 + Math.max(0, ...graph.nodes.map((n) => n.unknown?.seed ?? 0));
      const checks = parentIds.map((pid) => el('input', { type: 'checkbox', value: pid, id: `dag-unknown-input-${pid}`, checked: Boolean(u?.inputs.includes(pid)) }));
      const size = numeric('dag-unknown-size', u?.size ?? 0.5, { min: 0 }), smoothness = numeric('dag-unknown-smoothness', u?.smoothness ?? 2, { min: 0 });
      const saveUnknown = (seed = u?.seed ?? freshSeed) => {
        if (!unknownBox.isConnected) return;
        const focus = document.activeElement?.id;
        try {
          const inputs = checks.filter((c) => c.checked).map((c) => c.value);
          if (!inputs.length && !u) return;  // no parent chosen yet: keep what was typed
          const value = inputs.length ? { inputs, size: inputNumber(size, 'Size'), smoothness: inputNumber(smoothness, 'Smoothness'), seed } : null;
          if (JSON.stringify(value) === JSON.stringify(u ?? null)) return;
          const { unknown: _, ...rest } = node;
          if (apply({ ...graph, nodes: graph.nodes.map((n) => n.id === node.id ? (value ? { ...rest, unknown: value } : rest) : n) },
            value ? `${node.id} depends on ${inputs.join(', ')} through a relationship of unknown form; no analysis is told its shape.` : `${node.id} has no relationship of unknown form.`)
            && focus) $(focus)?.focus();
        } catch (error) { status(error.message, true); }
      };
      const unknownBox = el('details', { id: 'dag-unknown', class: 'dag-unknown', open: Boolean(u) },
        el('summary', {}, u ? `Unknown form: f(${u.inputs.join(', ')})` : 'Relationship of unknown form'),
        el('p', { class: 'small muted' }, 'Adds a smooth random function of the chosen parents to this variable, which no analysis is told about. Size is its SD over the population on this variable’s predictor scale; smoothness is roughly how many SDs of a parent one bend spans (1 is wiggly, 3 gentle). Its arrows show “?”.'),
        el('div', { class: 'dag-check-row', role: 'group', 'aria-label': 'Parents in the relationship of unknown form' },
          ...checks.map((c) => el('label', { for: c.id }, c, ` ${c.value}`))),
        el('div', { class: 'dag-pair' }, field('Size (SD)', size), field('Smoothness', smoothness)),
        ...(u ? [el('button', { type: 'button', class: 'secondary', id: 'dag-unknown-redraw', onclick: () => saveUnknown(freshSeed) }, `Draw another shape (now #${u.seed ?? 0})`)] : []));
      // A ticked parent saves at once; a typed size or smoothness on Enter or when the focus leaves the section, so a
      // click on a checkbox or the button inside it is never lost to the settings being rebuilt under it.
      unknownBox.addEventListener('change', (event) => { if (event.target.type === 'checkbox') setTimeout(() => saveUnknown()); });
      unknownBox.addEventListener('keydown', (event) => {
        if (event.key === 'Enter' && event.target.matches('input:not([type="checkbox"])')) { event.preventDefault(); saveUnknown(); }
      });
      unknownBox.addEventListener('focusout', (event) => { if (!unknownBox.contains(event.relatedTarget)) setTimeout(() => saveUnknown()); });
      box.append(unknownBox);
    }
    mountFunctions(box, { graph, node, onChange: (next) => apply(next), onError: (message) => status(message, true) });
    box.append(el('button', { type: 'button', class: 'secondary dag-delete', onclick: () => apply(removeNode(graph, selected), 'Variable and its arrows removed. Any observation rule using it was reset. Undo restores them.') }, 'Delete variable'));
  }
  function buttons() { $('dag-undo').disabled = !history.length; $('dag-redo').disabled = !future.length; }
  function refreshPower() { powerEditor?.refresh(estimand, powerSettings); }
  function render() {
    draw(); inspector(); buttons();
    mountTimeline($('dag-timeline'), graph);
    mountStudy($('dag-study'), graph, (next) => apply(next, 'Study structure updated.', { fit: true }), (message) => status(message, true));
    $('dag-n').previousElementSibling.textContent = graph.study ? `Number of ${sizeUnit(graph)}` : 'Observations';
    refreshPower();
    if (preview.key !== modelKey()) schedulePreview();
    const { order, byId } = validateGraph(graph);
    $('dag-equations').replaceChildren(...order.map((id) => el('li', {}, el('code', {}, equation(graph, byId.get(id))))));
    $('dag-order').textContent = order.length
      ? graph.study ? `Within-visit generation order: ${order.join(' → ')}. Variables retain their cluster or subject values across visits; lagged links use earlier visits.` : `Generation order: ${order.join(' → ')}. Exogenous draws are independent across variables and rows.`
      : 'Add variables to define the conditional distributions and generation order.';
  }
  function showResult() {
    const table = (head, rows, label) => el('div', { class: 'table-wrap', tabindex: 0, 'aria-label': label },
      el('table', {}, el('thead', {}, el('tr', {}, head.map((v) => el('th', { scope: 'col' }, v)))),
        el('tbody', {}, rows.map((row) => el('tr', {}, row.map((v) => el('td', {}, v)))))));
    const survivalRows = result.survival.map((s) => [s.id, String(s.events), String(s.censored), `${(100 * s.censored / result.n).toFixed(1)}%`]);
    $('dag-results').replaceChildren(
      el('p', { class: 'dag-result-heading' }, graph.study ? `${result.units} ${result.sizeUnit} · ${result.subjects} subjects · ${result.n} measurement rows · seed ${result.seed}` : `N = ${result.n} · seed ${result.seed}`),
      ...(graph.study ? [el('p', { class: 'small muted' }, 'Summaries describe the generated values, including unobserved outcomes. In CSV, *_recorded is blank when missing and *_observed identifies observed measurements.')] : []),
      table(['Variable', 'Mean', 'SD', 'Minimum', 'Maximum'], result.summaries.map((s) =>
        [s.id, number(s.mean), number(s.sd), number(s.min), number(s.max)]), 'Variable summaries'),
      ...(survivalRows.length ? [table(['Survival variable', 'Events', 'Censored', 'Censored %'], survivalRows, 'Survival observation summaries'),
        el('p', { class: 'small muted' }, 'T is the latent event time, T_time the observed follow-up, T_event 1 for an event and 0 if censored.')] : []));
    $('dag-download-csv').disabled = false;
  }
  // Adding a variable: click a type to place it in view, or drag the type onto the canvas (into a plate, to
  // measure it at that level).  Double-clicking empty canvas adds a continuous variable there.
  for (const button of root.querySelectorAll('[data-add-kind]')) {
    let press = null;
    button.addEventListener('pointerdown', (e) => { if (e.button === 0) { press = { x: e.clientX, y: e.clientY, moved: false }; button.setPointerCapture(e.pointerId); } });
    button.addEventListener('pointermove', (e) => {
      if (press && Math.hypot(e.clientX - press.x, e.clientY - press.y) > 6) { press.moved = true; canvas.classList.add('dropping'); }
    });
    button.addEventListener('pointerup', (e) => {
      const moved = press?.moved; press = null; canvas.classList.remove('dropping');
      if (!moved) return;
      const box = canvas.getBoundingClientRect();
      if (e.clientX >= box.left && e.clientX <= box.right && e.clientY >= box.top && e.clientY <= box.bottom) {
        const p = point(e);
        addVariable(button.dataset.addKind, { x: p.x - NODE_W / 2, y: p.y - NODE_H / 2 });
      }
      button.dataset.dragged = '1';
    });
    button.addEventListener('pointercancel', () => { press = null; canvas.classList.remove('dropping'); });
    button.addEventListener('click', () => {
      if (button.dataset.dragged) { delete button.dataset.dragged; return; }
      addVariable(button.dataset.addKind);
    });
  }
  // Start over on an empty canvas: no variables, arrows, study levels or question; the analysis settings return to
  // their defaults. One undo brings the whole model back.
  $('dag-clear').addEventListener('click', () => {
    if (!graph.nodes.length) return;
    const blank = scenarioDesign('empty'), count = graph.nodes.length;
    if (apply(blank.graph, `Cleared ${count} variable${count === 1 ? '' : 's'} and their arrows. Undo restores the model.`, { estimand: blank.estimand, power: blank.power, fit: true }))
      $('dag-add').focus();
  });
  $('dag-tidy').addEventListener('click', () => {
    remember(); graph = layout(graph); persist(); view = null; draw(); buttons();
    status(graph.study ? 'Arranged by level (clusters, subjects, visits), then generation order. Undo restores the old layout.'
      : 'Arranged left to right in generation order. Undo restores the old layout.');
  });
  // The study's levels, switched on the canvas: each adds or removes a plate.
  function toggleLevel(which) {
    const clusters = Boolean(graph.study?.clusterSizes?.length), visits = (graph.study?.times?.length ?? 0) > 1;
    const next = setStructure(graph, { clusters: which === 'clusters' ? !clusters : clusters, visits: which === 'visits' ? !visits : visits });
    const on = which === 'clusters' ? !clusters : !visits;
    const message = which === 'clusters'
      ? on ? 'Clusters added: the sample size now counts clusters. Drag a variable into the cluster plate to measure it once per cluster; click the plate’s label to set subjects per cluster.'
        : 'Clusters removed. Cluster-level variables are now measured once per subject.'
      : on ? 'Repeated visits added. Drag the outcome into the visit plate to measure it at every visit; click the plate’s label to set the visit times.'
        : 'Repeated visits removed. Visit-level variables are now measured once per subject.';
    if (apply(next, message, { fit: true })) {
      const layered = JSON.stringify(graph.nodes.map((n) => [n.x, n.y]));
      // the new plates need room: lay the graph out again when the plates would overlap other boxes
      const { plates } = plateLayout(graph);
      const crowded = plates.some((p) => graph.nodes.some((n) => !p.members.includes(n.id) && p.level !== 'cluster'
        && (p.level === 'visit' ? true : scopeOf(n) === 'cluster') && n.x + NODE_W > p.x && n.x < p.x + p.w && n.y + NODE_H > p.y && n.y < p.y + p.h));
      if (crowded) { graph = layout(graph); persist(); view = null; draw(); }
      if (JSON.stringify(graph.nodes.map((n) => [n.x, n.y])) !== layered) status(`${message} The boxes were rearranged to make room.`);
    }
  }
  $('dag-level-clusters').addEventListener('click', () => toggleLevel('clusters'));
  $('dag-level-visits').addEventListener('click', () => toggleLevel('visits'));
  $('dag-zoom-in').addEventListener('click', () => zoom(1.25));
  $('dag-zoom-out').addEventListener('click', () => zoom(1 / 1.25));
  $('dag-zoom-fit').addEventListener('click', fit);
  $('dag-help-toggle').addEventListener('click', () => {
    const open = $('dag-shortcuts').hidden;
    $('dag-shortcuts').hidden = !open; $('dag-help-toggle').setAttribute('aria-expanded', String(open));
  });
  // The canvas background: drag to pan, Ctrl/⌘ + wheel (or a pinch) to zoom, double-click to add a variable, click
  // to clear the arrow selection or cancel an arrow being drawn.
  // Only the empty canvas (and the plates' backgrounds) pans. The pointer is captured once a pan has really begun, so
  // a plain click on an arrow, a number or a plate label still reaches it.
  const handles = '.dag-node, .dag-edge, .dag-lag, .dag-random, .dag-plate-label, .dag-port';
  canvas.addEventListener('pointerdown', (e) => {
    if (e.button !== 0 || gesture || e.target.closest?.(handles)) return;
    gesture = { type: 'pan', x: e.clientX, y: e.clientY, view: { ...view }, moved: false, pointerId: e.pointerId };
  });
  canvas.addEventListener('pointermove', (e) => {
    if (connectFrom || gesture?.type === 'connect') {
      pointer = point(e);
      if (gesture?.type === 'connect' && !gesture.moved && Math.hypot(pointer.x - gesture.start.x, pointer.y - gesture.start.y) > 6) {
        gesture.moved = true; status(`Arrow from ${gesture.source}: release over the variable it affects.`); draw();
      }
      // only the rubber band moves: redrawing the boxes under the pointer would replace the one about to be clicked
      const s = graph.nodes.find((n) => n.id === (gesture?.source ?? connectFrom));
      if (!s || !(connectFrom || gesture.moved)) return;
      let rubber = canvas.querySelector('.dag-rubber');
      if (!rubber) {
        rubber = svg('path', { class: 'dag-rubber', 'marker-end': 'url(#dag-arrow)' });
        canvas.querySelector('.dag-overlay')?.append(rubber);
      }
      rubber.setAttribute('d', `M ${s.x + NODE_W} ${s.y + NODE_H / 2} L ${pointer.x} ${pointer.y}`);
      return;
    }
    if (gesture?.type !== 'pan') return;
    const scale = 1 / (canvas.getScreenCTM()?.a || 1), dx = e.clientX - gesture.x, dy = e.clientY - gesture.y;
    if (!gesture.moved && Math.hypot(dx, dy) < 4) return;
    if (!gesture.moved) { try { canvas.setPointerCapture(gesture.pointerId); } catch { /* the pointer is gone */ } }
    gesture.moved = true; canvas.classList.add('panning');
    view = { ...view, x: gesture.view.x - dx * scale, y: gesture.view.y - dy * scale };
    canvas.setAttribute('viewBox', viewBox(view));
  });
  const endGesture = (e) => {
    const g = gesture;
    if (!g) return;
    canvas.classList.remove('panning');
    if (g.type === 'connect') {
      pointer = null;
      if (!g.moved) { gesture = null; return; }  // a click on the port: its click handler starts click-to-connect
      gesture = { suppressClick: true };
      const under = document.elementFromPoint(e.clientX, e.clientY)?.closest?.('.dag-node');
      const target = under?.dataset.node;
      if (target && target !== g.source) { connectFrom = g.source; connectTo(target); if (connectFrom) { connectFrom = null; draw(); } }
      else { draw(); status('Arrow cancelled: release over another variable to connect them.'); }
      setTimeout(() => { if (gesture?.suppressClick) gesture = null; });
      return;
    }
    gesture = null;
    if (g.type === 'pan' && !g.moved && e.type === 'pointerup') {
      if (connectFrom) cancelConnect();
      else if (selectedEdge) { selectedEdge = null; draw(); }
    }
  };
  canvas.addEventListener('pointerup', endGesture);
  // an arrow dragged off the canvas and released there is cancelled
  document.addEventListener('pointerup', (e) => { if (gesture?.type === 'connect' && !canvas.contains(e.target)) { gesture = null; pointer = null; draw(); } });
  canvas.addEventListener('pointercancel', endGesture);
  canvas.addEventListener('pointerleave', () => { if (connectFrom && !gesture) { pointer = null; canvas.querySelector('.dag-rubber')?.remove(); } });
  canvas.addEventListener('dblclick', (e) => {
    if (e.target.closest('.dag-node, .dag-edge, .dag-lag, .dag-random, .dag-plate-label')) return;
    const p = point(e);
    addVariable('continuous', { x: p.x - NODE_W / 2, y: p.y - NODE_H / 2 });
  });
  canvas.addEventListener('wheel', (e) => {
    if (!e.ctrlKey && !e.metaKey) return;
    e.preventDefault();
    zoom(Math.exp(-e.deltaY * (e.deltaMode === 1 ? 0.05 : 0.0025)), point(e));
  }, { passive: false });
  $('dag-node-select').addEventListener('change', (e) => { selectedEdge = null; updateSelection(e.target.value); });
  // Shortcuts: undo and redo anywhere outside a text field; on the canvas, F fits, + and − zoom, ? lists them.
  root.addEventListener('keydown', (event) => {
    if (event.key === 'Escape') {
      closeMenu();
      if (connectFrom || gesture?.type === 'connect') cancelConnect();
      if (!$('dag-shortcuts').hidden) { $('dag-shortcuts').hidden = true; $('dag-help-toggle').setAttribute('aria-expanded', 'false'); }
      return;
    }
    const typing = event.target.closest?.('input, select, textarea, [contenteditable]');
    if (typing || event.ctrlKey || event.metaKey || event.altKey || !event.target.closest?.('.dag-scroll')?.contains(canvas)) return;
    if (event.key === 'f' || event.key === 'F') { event.preventDefault(); fit(); }
    else if (event.key === '+' || event.key === '=') { event.preventDefault(); zoom(1.25); }
    else if (event.key === '-' || event.key === '_') { event.preventDefault(); zoom(1 / 1.25); }
    else if (event.key === '?') { event.preventDefault(); $('dag-help-toggle').click(); }
  });
  // Undo and redo also work when focus has left the page's controls (after a deleted box, focus falls to the body).
  document.addEventListener('keydown', (event) => {
    if (event.key === 'Escape' && event.target === document.body && (connectFrom || gesture?.type === 'connect')) { cancelConnect(); return; }
    const mod = event.ctrlKey || event.metaKey, key = event.key.toLowerCase();
    if (!mod || event.altKey || (key !== 'z' && key !== 'y') || !root.isConnected) return;
    if (event.target !== document.body && !root.contains(event.target)) return;
    if (event.target.closest?.('input, select, textarea, [contenteditable]')) return;
    event.preventDefault();
    $(key === 'y' || event.shiftKey ? 'dag-redo' : 'dag-undo').click();
  });
  function timeTravel(from, to) {
    if (!from.length) return;
    to.push(JSON.stringify(exportDesign(graph, estimand, powerSettings)));
    const design = parseDesign(from.pop()); graph = design.graph; estimand = design.estimand; powerSettings = design.power ?? null; pending = null;
    selected = graph.nodes.some((n) => n.id === selected) ? selected : graph.nodes.at(-1)?.id;
    connectFrom = null; selectedEdge = null; clearResult(); persist(); render(); status('Model restored.');
  }
  $('dag-undo').addEventListener('click', () => timeTravel(history, future));
  $('dag-redo').addEventListener('click', () => timeTravel(future, history));
  function loadExample(name) {
    const next = scenarioDesign(name);
    const label = SCENARIOS.find((s) => s.id === name)?.label ?? 'Scenario';
    if (!apply(next.graph, `${label} loaded with its question and analysis. All values are editable example assumptions.`, { estimand: next.estimand, power: next.power, fit: true })) return;
    $('dag-example').value = name; describeScenario();
  }
  // every starting point, grouped by kind of study; loading one replaces the design (undo restores it)
  function describeScenario() {
    const s = SCENARIOS.find((s) => s.id === $('dag-example').value);
    $('dag-scenario-description').textContent = s ? `${s.description} ${s.assumption}` : '';
  }
  $('dag-example').replaceChildren(...[...new Set(SCENARIOS.map((s) => s.category))].map((category) => el('optgroup', { label: category },
    SCENARIOS.filter((s) => s.category === category).map((s) => new Option(s.label, s.id)))));
  $('dag-example').value = 'trial'; describeScenario();
  $('dag-example').addEventListener('change', describeScenario);
  $('dag-example-load').addEventListener('click', () => loadExample($('dag-example').value));
  $('dag-export').addEventListener('click', () => download(JSON.stringify(exportDesign(graph, estimand, powerSettings), null, 2) + '\n', 'application/json', 'simulation-design.json'));
  $('dag-import').addEventListener('change', async (event) => {
    const file = event.target.files?.[0];
    if (!file) return;
    try {
      if (file.size > 5e6) throw new Error('Use a model or power report JSON file smaller than 5 MB.');
      const design = parseDesign(await file.text());
      apply(design.graph, 'Simulation design imported.', { estimand: design.estimand, power: design.power ?? null, fit: true });
    } catch (error) { status(`Could not import model: ${error.message}`, true); }
    event.target.value = '';
  });
  $('dag-simulation').addEventListener('submit', async (event) => {
    event.preventDefault(); clearResult();
    try {
      validateGraph(graph);
      const n = inputNumber($('dag-n'), 'Sample size'), seed = inputNumber($('dag-seed'), 'Seed');
      const instance = new Engine(); worker = instance;
      $('dag-simulate').disabled = true; $('dag-cancel').hidden = false;
      status('Generating data in C++/WebAssembly…');
      const data = await instance.dag({ op: 'simulate', graph, options: { n, seed } });
      if (worker !== instance) return;
      instance.dispose(); worker = null; $('dag-simulate').disabled = false; $('dag-cancel').hidden = true;
      result = data; showResult(); status('Dataset generated. Download it as CSV for the full data.');
    } catch (error) { if (!error.stopped) { clearResult(); status(error.message, true); } }
  });
  $('dag-cancel').addEventListener('click', () => { clearResult(); status('Simulation stopped.'); });
  for (const id of ['dag-n', 'dag-seed']) $(id).addEventListener('input', () => {
    clearResult(); status('Simulation settings changed. Run again to update the data.');
  });
  $('dag-download-csv').addEventListener('click', () => {
    if (result) download(dataCsv(graph, result), 'text/csv', `data-n${result.n}-seed${result.seed}.csv`);
  });
  powerEditor = mountPower($('dag-power'), {
    getGraph: () => graph,
    engine: { borrow: borrowEngine, give: returnEngine },
    onTarget(value, mode) {
      if (JSON.stringify(value) !== JSON.stringify(estimand)) {
        if (mode !== 'auto') pending ??= JSON.stringify(exportDesign(graph, estimand, powerSettings));
        // the canvas shows the roles, the tested arrows and the analysis's paths; an effect size or null value
        // typed in the form changes none of them
        const drawn = (e) => e && JSON.stringify([e.kind, e.exposure, e.outcome, e.measure, e.reference, e.active, e.horizon, e.repeated]);
        const redraw = drawn(value) !== drawn(estimand);
        estimand = value; persist();
        if (redraw) draw();
      }
      if (mode === 'commit' && pending !== null) {
        if (pending !== JSON.stringify(exportDesign(graph, estimand, powerSettings))) {
          history.push(pending); if (history.length > 50) history.shift(); future = [];
        }
        pending = null; buttons();
      }
    },
    onSettings(settings, mode) {
      if (JSON.stringify(settings) !== JSON.stringify(powerSettings)) {
        if (mode !== 'auto') pending ??= JSON.stringify(exportDesign(graph, estimand, powerSettings));
        // only the analysis changes what the canvas shows (the paths it leaves open); sizes and α do not
        const redraw = JSON.stringify(settings.analysis) !== JSON.stringify(powerSettings?.analysis);
        powerSettings = settings; persist();
        if (redraw) draw();
      }
      if (mode === 'commit' && pending !== null) {
        if (pending !== JSON.stringify(exportDesign(graph, estimand, powerSettings))) {
          history.push(pending); if (history.length > 50) history.shift(); future = [];
        }
        pending = null; buttons();
      }
    },
  });
  document.addEventListener('pointerdown', (event) => { if (menu && !menu.contains(event.target)) closeMenu(); });
  render(); clearResult();
  status(storageNote || (graph.nodes.length
    ? 'Design restored from this browser.'
    : 'Choose a starting point and press Load, or add a variable.'));
}
