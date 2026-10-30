import { addEdge, removeEdge, removeNode, makeNode, exampleGraph, emptyGraph, renameNode, validateGraph,
  equation, dataCsv, topologicalOrder, KINDS, DISTRIBUTIONS, FAMILIES, MAX_NODES } from './dag-model.js?v=dev';
import { validateEstimand, renameEstimand, exportDesign, parseDesign } from './estimand-model.js?v=dev';
import { mountPower } from './dag-power.js?v=dev';
import { calibrationEdges, powerCompatibility } from './power-model.js?v=dev';
import { sketch, findPaths, classifyPath, pathText } from './dag-view.js?v=dev';
import { Engine } from './engine.js?v=dev';

const SVG = 'http://www.w3.org/2000/svg';
const STORE = 'ssc-dag-v1';
const number = (x) => x === null ? '—' : Number.isInteger(x) ? String(x) : x.toFixed(3);
const kindNames = { continuous: 'Continuous', binary: 'Binary', count: 'Count', survival: 'Survival time' };
const parameterLabels = { sd: 'Standard deviation', width: 'Uniform width', shape: 'Shape', scale: 'Baseline time scale' };
const W = 180, H = 104;  // node box; positions run 0–920 × 0–440 in a 1100 × 564 drawing
const PREVIEW_N = 2000;  // observations drawn to sketch each variable's distribution
const short = (x) => String(Number(x.toPrecision(3))).replace('-', '−');
// What a node's box says about its distribution: a root's values, or the scale its incoming arrows act on.
function summary(graph, node) {
  const i = node.intercept, e = Math.exp(i), family = DISTRIBUTIONS[node.family].label;
  if (graph.edges.some((edge) => edge.target === node.id)) {
    const rest = { gaussian: `SD ${short(node.sd)}`, uniform: `width ${short(node.width)}`, lognormal: 'log scale', gamma: 'log mean',
      bernoulli: 'log-odds', poisson: 'log mean', weibull: 'log HR', exponential: 'log HR' }[node.family];
    return `${family} · ${rest}`;
  }
  return `${family} · ${{ gaussian: `mean ${short(i)}, SD ${short(node.sd)}`, uniform: `center ${short(i)}, width ${short(node.width)}`,
    lognormal: `median ${short(e)}`, gamma: `mean ${short(e)}, shape ${short(node.shape)}`, bernoulli: `P(1) = ${short(1 / (1 + Math.exp(-i)))}`,
    poisson: `mean ${short(e)}`, weibull: `k ${short(node.shape)}, scale ${short(node.scale)}`, exponential: `scale ${short(node.scale)}` }[node.family]}`;
}
// Columns by generation depth, left to right; within a column, near the parents.
function tidy(graph) {
  const depth = new Map(), y = new Map();
  for (const id of topologicalOrder(graph))
    depth.set(id, Math.max(0, ...graph.edges.filter((e) => e.target === id).map((e) => depth.get(e.source) + 1)));
  const columns = Math.max(...depth.values()) + 1, dx = columns > 1 ? Math.min(280, 920 / (columns - 1)) : 0;
  const x0 = (920 - dx * (columns - 1)) / 2;
  for (let c = 0; c < columns; c++) {
    const parentsY = (id) => { const ys = graph.edges.filter((e) => e.target === id).map((e) => y.get(e.source)); return ys.length ? ys.reduce((a, b) => a + b) / ys.length : 0; };
    const column = graph.nodes.filter((n) => depth.get(n.id) === c).sort((a, b) => c ? parentsY(a.id) - parentsY(b.id) : a.y - b.y);
    const dy = column.length > 1 ? Math.min(150, 440 / (column.length - 1)) : 0, y0 = (440 - dy * (column.length - 1)) / 2;
    column.forEach((n, k) => y.set(n.id, Math.round(y0 + k * dy)));
  }
  return { ...graph, nodes: graph.nodes.map((n) => ({ ...n, x: Math.round(x0 + depth.get(n.id) * dx), y: y.get(n.id) })) };
}

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
  const preview = { columns: null, key: '', token: 0, engine: null, timer: 0 };
  const modelKey = () => JSON.stringify([graph.nodes.map(({ x, y, ...rest }) => rest), graph.edges]);
  let menu = null;
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
      remember(); graph = next; estimand = nextEstimand;
      if ('power' in planning) powerSettings = planning.power;
      if (planning.renameFrom && powerSettings?.calibrationEdge) {
        powerSettings = { ...powerSettings, calibrationEdge: powerSettings.calibrationEdge.split('>').map((id) => id === planning.renameFrom ? planning.renameTo : id).join('>') };
      }
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
    const marker = (id) => svg('marker', { id, viewBox: '0 0 10 10', refX: 9, refY: 5, markerWidth: 13, markerHeight: 13,
      orient: 'auto-start-reverse', markerUnits: 'userSpaceOnUse' }, svg('path', { d: 'M 0 0 L 10 5 L 0 10 z' }));
    const defs = svg('defs', {}, marker('dag-arrow'), marker('dag-arrow-effect'));
    const lines = svg('g', { class: 'dag-arrows' });
    const map = new Map(graph.nodes.map((n) => [n.id, n]));
    // the arrows that carry the effect the power analysis tests
    const effect = new Set((estimand ? calibrationEdges(graph, estimand) : []).map((e) => `${e.source}>${e.target}`));
    for (const edge of graph.edges) {
      const source = map.get(edge.source), target = map.get(edge.target);
      const dx = target.x - source.x, dy = target.y - source.y;
      const vertical = Math.abs(dy) > Math.abs(dx) * 0.6;
      const x1 = source.x + (vertical ? W / 2 : dx >= 0 ? W : 0);
      const y1 = source.y + (vertical ? dy >= 0 ? H : 0 : H / 2);
      const x2 = target.x + (vertical ? W / 2 : dx >= 0 ? 0 : W);
      const y2 = target.y + (vertical ? dy >= 0 ? 0 : H : H / 2);
      const bend = Math.max(30, Math.abs(vertical ? y2 - y1 : x2 - x1) * 0.4);
      const [c1, c2] = vertical ? [[x1, y1 + Math.sign(dy) * bend], [x2, y2 - Math.sign(dy) * bend]]
        : [[x1 + Math.sign(dx) * bend, y1], [x2 - Math.sign(dx) * bend, y2]];
      const d = `M ${x1} ${y1} C ${c1[0]} ${c1[1]}, ${c2[0]} ${c2[1]}, ${x2} ${y2}`;
      const mid = [0, 1].map((k) => ([x1, y1][k] + 3 * c1[k] + 3 * c2[k] + [x2, y2][k]) / 8);  // the curve at t = 1/2
      const key = `${edge.source}>${edge.target}`;
      const group = svg('g', { class: `dag-edge${edge.target === selected ? ' selected' : ''}${effect.has(key) ? ' effect' : ''}`,
        'data-edge': key, role: 'button', tabindex: 0,
        'aria-label': `${edge.source} to ${edge.target}, coefficient ${edge.coefficient}. Enter or click the number to change it.` },
      svg('path', { d, 'marker-end': `url(#dag-arrow${effect.has(key) ? '-effect' : ''})` }),
      svg('text', { x: mid[0], y: mid[1] - 8, 'text-anchor': 'middle', class: 'dag-edge-label' }, short(edge.coefficient)));
      group.addEventListener('click', (e) => {
        if (e.target.classList.contains('dag-edge-label')) editCoefficient(edge, e.target);
        else updateSelection(edge.target);
      });
      group.addEventListener('keydown', (e) => {
        if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); editCoefficient(edge, group.querySelector('text')); }
      });
      lines.append(group);
    }
    const nodes = svg('g', {});
    for (const node of graph.nodes) {
      const group = svg('g', { class: `dag-node ${node.kind}${node.id === selected ? ' selected' : ''}${node.id === connectFrom ? ' connecting' : ''}`,
        transform: `translate(${node.x} ${node.y})`, tabindex: 0, role: 'button',
        'aria-label': `${node.id}, ${node.label}, ${kindNames[node.kind]}. Enter to edit, arrow keys to move, C to compare it, O to make it the outcome.` },
      svg('rect', { class: 'dag-box', width: W, height: H, rx: 12 }),
      svg('text', { x: 14, y: 22, class: 'dag-node-role' }, node.kind.toUpperCase()
        + (node.kind === 'survival' && node.censoring.mode !== 'none'
          ? ` · ${node.censoring.mode === 'variable' ? 'censored by ' + node.censoring.variable : 'cutoff ' + short(node.censoring.time)}` : '')),
      svg('text', { x: 14, y: 45, class: 'dag-node-title' }, `${node.id} · ${node.label.length > 16 ? node.label.slice(0, 15) + '…' : node.label}`),
      svg('text', { x: 14, y: 66, class: 'dag-node-family' }, summary(graph, node)));
      const role = estimand?.exposure === node.id ? (estimand.kind === 'coefficient' ? 'Predictor' : 'Compared')
        : estimand?.outcome === node.id ? 'Outcome' : null;
      const values = preview.columns?.[node.id];
      const drawn = values && sketch(node.kind, values, values && (node.kind === 'binary' || node.kind === 'survival') ? 92 : W - 28, 22);
      if (drawn) group.append(svg('g', { class: `dag-spark${preview.key === modelKey() ? '' : ' stale'}`, transform: `translate(14 ${H - 30})` },
        svg('title', {}, drawn.title), svg('line', { x1: 0, x2: W - 28, y1: 22.5, y2: 22.5 }),
        ...(drawn.bars ? [svg('path', { class: 'bars', d: drawn.bars })] : []),
        ...(drawn.line ? [svg('path', { class: 'curve', d: drawn.line })] : []),
        ...(drawn.label ? [svg('text', { x: 100, y: 17 }, drawn.label)] : [])));
      if (role) group.append(svg('g', { class: 'dag-node-tag' }, svg('rect', { x: 10, y: -11, width: 70, height: 20, rx: 10 }),
        svg('text', { x: 45, y: 3, 'text-anchor': 'middle' }, role)));
      const port = svg('g', { class: 'dag-port', role: 'button', tabindex: 0, 'aria-label': `Start an arrow from ${node.id}` },
        svg('circle', { cx: W, cy: H / 2, r: 10 }), svg('text', { x: W, y: H / 2 + 4, 'text-anchor': 'middle' }, '+'));
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
        node.x = Math.max(0, Math.min(1100 - W, Math.round(drag.x + dx)));
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
      group.addEventListener('contextmenu', (event) => { event.preventDefault(); openMenu(node, event.clientX, event.clientY); });
      group.addEventListener('keydown', (event) => {
        if (event.target !== group) return;
        const role = { c: 'exposure', o: 'outcome' }[event.key.toLowerCase()];
        if (role && !event.ctrlKey && !event.metaKey && !event.altKey) {
          event.preventDefault(); assignRole(role, node.id);
          [...canvas.querySelectorAll('.dag-node')].find((g) => g.getAttribute('aria-label').startsWith(`${node.id},`))?.focus();
          return;
        }
        if (event.key === 'ContextMenu' || (event.shiftKey && event.key === 'F10')) {
          event.preventDefault(); const box = group.getBoundingClientRect(); openMenu(node, box.left + 20, box.bottom - 10); return;
        }
        if (event.key === 'Enter' || event.key === ' ') {
          event.preventDefault();
          if (connectFrom) connectTo(node.id); else { updateSelection(node.id); $('dag-node-label').focus(); }
        }
        const delta = { ArrowLeft: [-10, 0], ArrowRight: [10, 0], ArrowUp: [0, -10], ArrowDown: [0, 10] }[event.key];
        if (delta) {
          event.preventDefault(); remember();
          node.x = Math.max(0, Math.min(1100 - W, node.x + delta[0]));
          node.y = Math.max(0, Math.min(440, node.y + delta[1]));
          persist(); draw(); buttons();
          [...canvas.querySelectorAll('.dag-node')].find((g) => g.getAttribute('aria-label').startsWith(`${node.id},`))?.focus();
        }
      });
      nodes.append(group);
    }
    if (!graph.nodes.length) nodes.append(svg('text', { x: 550, y: 260, 'text-anchor': 'middle', class: 'dag-empty' }, 'Add a variable to begin'));
    canvas.replaceChildren(svg('title', {}, 'Data generating process'),
      svg('desc', {}, 'Select a node to edit. Drag nodes or use arrow keys to move. Use plus ports, or Add arrow in the variable settings, to connect nodes.'), defs, lines, nodes);
    $('dag-node-count').replaceChildren(`${graph.nodes.length} variable${graph.nodes.length === 1 ? '' : 's'} · ${graph.edges.length} arrow${graph.edges.length === 1 ? '' : 's'}`,
      ...(effect.size ? [el('span', { class: 'dag-legend-effect' }, 'tested effect')] : []));
    $('dag-node-select').replaceChildren(...graph.nodes.map((n) => el('option', { value: n.id }, `${n.id} · ${n.label}`)));
    $('dag-node-select').value = selected ?? '';
    $('dag-node-select').disabled = !graph.nodes.length;
    $('dag-tidy').disabled = graph.nodes.length < 2;
    for (const [role, id] of [['exposure', 'dag-role-exposure'], ['outcome', 'dag-role-outcome']])
      $(id)?.setAttribute('aria-pressed', String(estimand?.[role] === selected));
    paths();
  }
  // The paths between the compared variable and the outcome, and whether the planned analysis leaves each open.
  function paths() {
    const box = $('dag-paths');
    const a = estimand?.exposure, y = estimand?.outcome;
    const mark = (keys) => { for (const g of canvas.querySelectorAll('.dag-edge')) g.classList.toggle('path', keys.includes(g.dataset.edge)); };
    if (!a || !y || !graph.nodes.some((n) => n.id === a) || !graph.nodes.some((n) => n.id === y)) { box.replaceChildren(); return; }
    const plan = powerCompatibility(graph, estimand);
    const adjusted = new Set((plan.supported ? plan.terms : []).flatMap((t) => [t.left, t.right]).filter((v) => v && v !== a));
    const { paths: list, truncated } = findPaths(graph, a, y);
    const rows = list.map((path) => {
      const c = classifyPath(graph, path, adjusted);
      const what = c.kind === 'causal' ? 'causal' : c.kind === 'backdoor' ? 'back-door' : 'non-causal';
      const why = c.open ? (c.kind === 'causal' ? 'open: carries effect' : 'open: biases the estimate')
        : c.collider ? `blocked at collider ${c.by}` : `blocked: ${c.by} is adjusted for`;
      const keys = path.map((s) => s.forward ? `${s.from}>${s.to}` : `${s.to}>${s.from}`);
      const row = el('li', { class: `${c.kind}${c.open ? ' open' : ''}`, tabindex: 0 },
        el('span', { class: 'dag-path-text' }, pathText(path)), el('span', { class: 'dag-path-kind' }, `${what} · ${why}`));
      const show = (on) => mark(on ? keys : []);
      row.addEventListener('mouseenter', () => show(true)); row.addEventListener('mouseleave', () => show(false));
      row.addEventListener('focus', () => show(true)); row.addEventListener('blur', () => show(false));
      return { row, c };
    });
    const biased = rows.filter(({ c }) => c.open && c.kind !== 'causal').length;
    const heading = `Paths from ${a} to ${y}` + (adjusted.size ? ` · analysis adjusts for ${[...adjusted].join(', ')}` : '');
    box.replaceChildren(el('h4', {}, heading),
      rows.length ? el('ul', {}, rows.map(({ row }) => row)) : el('p', { class: 'small muted' }, `No path connects ${a} and ${y}: the effect is zero.`),
      ...(truncated ? [el('p', { class: 'small muted' }, `The first ${list.length} paths are shown.`)] : []),
      ...(biased ? [el('p', { class: 'small error' }, `${biased} open non-causal path${biased === 1 ? '' : 's'}: the estimate would mix confounding into the effect.`)] : []));
  }
  function assignRole(role, id) {
    if (!powerEditor?.assign(role, id)) status('The power analysis needs at least two variables.', true);
    else status(`${id} is now the ${role === 'exposure' ? 'compared variable' : 'outcome'} in the power analysis.`);
  }
  function closeMenu() { menu?.remove(); menu = null; }
  function openMenu(node, x, y) {
    closeMenu();
    const item = (text, action) => el('button', { type: 'button', role: 'menuitem', onclick: () => { closeMenu(); action(); } }, text);
    menu = el('div', { class: 'dag-menu', role: 'menu', 'aria-label': `${node.id} actions` },
      item('Compare this variable', () => assignRole('exposure', node.id)),
      item('Use as the outcome', () => assignRole('outcome', node.id)),
      item('Edit', () => { updateSelection(node.id); $('dag-node-label').focus(); }),
      item('Draw an arrow from it', () => connectTo(node.id)),
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
  // A coefficient changed where it is drawn: a small field over the arrow's number.
  function editCoefficient(edge, label) {
    const host = canvas.parentElement, at = label.getBoundingClientRect(), base = host.getBoundingClientRect();
    const input = el('input', { type: 'number', step: 'any', value: edge.coefficient, class: 'dag-inline',
      'aria-label': `Coefficient of ${edge.source} → ${edge.target}` });
    host.append(input);
    input.style.left = `${at.left - base.left + host.scrollLeft + at.width / 2 - 36}px`;
    input.style.top = `${at.top - base.top + host.scrollTop - 6}px`;
    let done = false;
    const finish = (save) => {
      if (done) return; done = true; input.remove();
      if (save && input.value.trim() !== '' && Number(input.value) !== edge.coefficient) {
        try { apply({ ...graph, edges: graph.edges.map((e) => e === edge ? { ...e, coefficient: inputNumber(input, 'Arrow coefficient') } : e) },
          `${edge.source} → ${edge.target} set to ${input.value}.`); }
        catch (error) { status(error.message, true); }
      }
      canvas.querySelector(`[data-edge="${edge.source}>${edge.target}"]`)?.focus();
    };
    input.addEventListener('keydown', (e) => {
      if (e.key === 'Enter') { e.preventDefault(); finish(true); }
      if (e.key === 'Escape') { e.preventDefault(); e.stopPropagation(); finish(false); }
    });
    input.addEventListener('blur', () => finish(true));
    input.focus(); input.select();
  }
  // Draws a fresh sample from the model for the boxes, a moment after the last change.
  function schedulePreview() {
    clearTimeout(preview.timer);
    preview.timer = setTimeout(async () => {
      const token = ++preview.token, key = modelKey();
      preview.engine?.dispose();
      if (!graph.nodes.length) { preview.columns = null; preview.engine = null; draw(); return; }
      const instance = new Engine(); preview.engine = instance;
      try {
        const data = await instance.dag({ op: 'simulate', graph, options: { n: PREVIEW_N, seed: 1 } });
        if (token !== preview.token) return;
        preview.columns = data.columns; preview.key = key;
      } catch (error) {
        if (token !== preview.token || error.stopped) return;
        preview.columns = null;
      } finally {
        if (preview.engine === instance) { instance.dispose(); preview.engine = null; }
      }
      draw();
    }, 250);
  }
  function point(event) {
    const p = canvas.createSVGPoint(); p.x = event.clientX; p.y = event.clientY;
    return p.matrixTransform(canvas.getScreenCTM().inverse());
  }
  function inspector() {
    const node = graph.nodes.find((n) => n.id === selected);
    const box = $('dag-inspector'); box.replaceChildren();
    if (!node) { box.append(el('p', { class: 'muted small' }, 'Add a variable, then choose its type and distribution here.')); return; }
    box.append(el('div', { class: 'dag-inspector-heading' }, el('h3', {}, node.id),
      el('span', { class: 'small muted' }, kindNames[node.kind])),
    el('div', { class: 'dag-roles', role: 'group', 'aria-label': 'Role in the power analysis' },
      el('button', { id: 'dag-role-exposure', type: 'button', 'aria-pressed': String(estimand?.exposure === node.id), onclick: () => assignRole('exposure', node.id) }, 'Compared'),
      el('button', { id: 'dag-role-outcome', type: 'button', 'aria-pressed': String(estimand?.outcome === node.id), onclick: () => assignRole('outcome', node.id) }, 'Outcome')));
    const id = el('input', { id: 'dag-node-id', type: 'text', maxlength: 32, value: node.id, required: true });
    const label = el('input', { id: 'dag-node-label', type: 'text', maxlength: 80, value: node.label, required: true });
    const kind = select('dag-node-kind', KINDS.map((k) => [k, kindNames[k]]), node.kind);
    const family = select('dag-node-family', [], node.family);
    const intercept = numeric('dag-node-intercept', node.intercept);
    // shown to six digits; the intercept is kept exactly unless the probability is edited
    const shownProbability = String(Number((1 / (1 + Math.exp(-node.intercept))).toPrecision(6)));
    const probability = numeric('dag-node-probability', shownProbability, { min: 0, max: 1 });
    const inputs = Object.fromEntries(['sd', 'width', 'shape', 'scale'].map((key) => [key, numeric(`dag-node-${key}`, node[key], { min: 0 })]));
    const censorMode = select('dag-censor-mode', [['none', 'None · fully observed'], ['administrative', 'Administrative cutoff'], ['variable', 'Use another variable']], node.censoring.mode);
    const censorTime = numeric('dag-censor-time', node.censoring.time, { min: 0 });
    const censorVariable = select('dag-censor-variable', graph.nodes.filter((n) => n.id !== node.id).map((n) => [n.id, `${n.id} · ${n.label}`]), node.censoring.variable);
    const form = el('form', { novalidate: true }, el('div', { class: 'dag-pair' }, field('ID', id), field('Label', label)),
      el('div', { class: 'dag-pair' }, field('Type', kind), field('Distribution', family)));
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
        ? field('Probability when parents are 0', probability, 'Arrows add to its log odds.')
        : field(`Intercept (${spec.link})`, intercept, 'Value when parents are 0; arrows add to it.'));
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
    form.append(parameters, observation, el('button', { type: 'submit', class: 'primary' }, 'Apply'));
    form.addEventListener('submit', (event) => {
      event.preventDefault();
      try {
        let baseline = inputNumber(intercept, 'Intercept');
        if (family.value === 'bernoulli' && !(node.family === 'bernoulli' && probability.value === shownProbability)) {
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
    box.append(el('button', { type: 'button', class: 'secondary dag-delete', onclick: () => apply(removeNode(graph, selected), 'Variable and its arrows removed. Any observation rule using it was reset. Undo restores them.') }, 'Delete variable'));
  }
  function buttons() { $('dag-undo').disabled = !history.length; $('dag-redo').disabled = !future.length; }
  function render() {
    draw(); inspector(); buttons();
    powerEditor?.refresh(estimand, powerSettings);
    if (preview.key !== modelKey()) schedulePreview();
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
    const survivalRows = result.survival.map((s) => [s.id, String(s.events), String(s.censored), `${(100 * s.censored / result.n).toFixed(1)}%`]);
    $('dag-results').replaceChildren(
      el('p', { class: 'dag-result-heading' }, `N = ${result.n} · seed ${result.seed}`),
      table(['Variable', 'Mean', 'SD', 'Minimum', 'Maximum'], result.summaries.map((s) =>
        [s.id, number(s.mean), number(s.sd), number(s.min), number(s.max)]), 'Variable summaries'),
      ...(survivalRows.length ? [table(['Survival variable', 'Events', 'Censored', 'Censored %'], survivalRows, 'Survival observation summaries'),
        el('p', { class: 'small muted' }, 'T is the latent event time, T_time the observed follow-up, T_event 1 for an event and 0 if censored.')] : []));
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
  $('dag-tidy').addEventListener('click', () => {
    remember(); graph = tidy(graph); persist(); draw(); buttons();
    status('Arranged left to right in generation order.');
  });
  $('dag-node-select').addEventListener('change', (e) => updateSelection(e.target.value));
  root.addEventListener('keydown', (event) => {
    if (event.key === 'Escape') { closeMenu(); connectFrom = null; draw(); status('Arrow connection cancelled.'); }
  });
  function timeTravel(from, to) {
    if (!from.length) return;
    to.push(JSON.stringify(exportDesign(graph, estimand, powerSettings)));
    const design = parseDesign(from.pop()); graph = design.graph; estimand = design.estimand; powerSettings = design.power ?? null; pending = null;
    selected = graph.nodes.some((n) => n.id === selected) ? selected : graph.nodes.at(-1)?.id;
    connectFrom = null; clearResult(); persist(); render(); status('Model restored.');
  }
  $('dag-undo').addEventListener('click', () => timeTravel(history, future));
  $('dag-redo').addEventListener('click', () => timeTravel(future, history));
  $('dag-example-load').addEventListener('click', () => apply(exampleGraph($('dag-example').value), 'Example loaded. All numbers are editable example assumptions.'));
  $('dag-export').addEventListener('click', () => download(JSON.stringify(exportDesign(graph, estimand, powerSettings), null, 2) + '\n', 'application/json', 'simulation-design.json'));
  $('dag-import').addEventListener('change', async (event) => {
    const file = event.target.files?.[0];
    if (!file) return;
    try {
      if (file.size > 5e6) throw new Error('Use a model or power report JSON file smaller than 5 MB.');
      const design = parseDesign(await file.text());
      apply(design.graph, 'Simulation design imported.', { estimand: design.estimand, power: design.power ?? null });
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
    onTarget(value, mode) {
      const before = JSON.stringify(exportDesign(graph, estimand, powerSettings));
      if (JSON.stringify(value) !== JSON.stringify(estimand)) {
        if (mode !== 'auto') pending ??= before;
        estimand = value; persist(); draw();
      }
      if (mode === 'commit' && pending !== null) {
        if (pending !== JSON.stringify(exportDesign(graph, estimand, powerSettings))) {
          history.push(pending); if (history.length > 50) history.shift(); future = [];
        }
        pending = null; buttons();
      }
    },
    onSettings(settings) { powerSettings = settings; persist(); },
  });
  document.addEventListener('pointerdown', (event) => { if (menu && !menu.contains(event.target)) closeMenu(); });
  render(); clearResult();
  status(storageNote || (graph.nodes.length
    ? 'Design restored from this browser.'
    : 'Add a variable or load an example to start.'));
}
