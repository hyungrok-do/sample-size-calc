import { TRANSFORMS, parseTransform, termText, transformValue } from './dag-functions.js?v=dev';

const el = (tag, attrs = {}, ...children) => {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) if (k.startsWith('on')) n.addEventListener(k.slice(2), v); else n.setAttribute(k, v);
  n.append(...children.flat()); return n;
};
const field = (label, input) => el('label', { class: 'dag-field' }, el('span', {}, label), input);
const numeric = (id, value) => el('input', { id, type: 'number', step: 'any', value });
const readNumber = (input) => {
  if (!input.value.trim() || !Number.isFinite(Number(input.value))) throw new Error('Enter a finite number for each response parameter.');
  return Number(input.value);
};
const select = (id, list, value) => {
  const n = el('select', { id }, list.map(([v, name]) => el('option', { value: v }, name))); n.value = value; return n;
};
function transformation(prefix, initial = 'square') {
  const kind = select(`${prefix}-kind`, TRANSFORMS, initial), params = el('div', { class: 'dag-pair' });
  const center = numeric(`${prefix}-center`, 0), scale = numeric(`${prefix}-scale`, 1), knot = numeric(`${prefix}-knot`, 0);
  const knots = el('input', { id: `${prefix}-knots`, value: '-1, 0, 1', type: 'text' }), index = numeric(`${prefix}-index`, 1);
  function render() {
    params.replaceChildren(...(['square', 'saturation'].includes(kind.value) ? [field('Center', center), field('Scale', scale)]
      : ['threshold', 'hinge'].includes(kind.value) ? [field('Threshold / knot', knot)]
      : kind.value === 'spline' ? [field('Knots (original units)', knots), field('Spline component (1 to knots − 2)', index)] : []));
  }
  kind.addEventListener('change', render); render();
  return { element: el('div', {}, field('Function', kind), params), read: () => parseTransform({ kind: kind.value,
    center: readNumber(center), scale: readNumber(scale), knot: readNumber(knot),
    knots: knots.value.trim().split(/[\s,]+/).filter(Boolean).map(Number), index: readNumber(index) - 1 }) };
}
export function mountFunctions(box, { graph, node, onChange, onError }) {
  const parents = graph.edges.filter((e) => e.target === node.id).map((e) => e.source);
  if (!parents.length) return;
  const change = (functions) => onChange({ ...graph, nodes: graph.nodes.map((n) => n.id === node.id ? { ...n, functions } : n) });
  const functions = node.functions ?? [];
  box.append(el('h4', {}, 'Response curves'), el('p', { class: 'small muted' },
    'Add terms to the linear arrows. Spline components and transformed interactions can be combined.'));
  for (const [i, t] of functions.entries()) {
    const coefficient = numeric(`dag-function-value-${i}`, t.coefficient);
    coefficient.setAttribute('aria-label', `${termText(t)} coefficient`);
    coefficient.addEventListener('change', () => { try { change(functions.map((f, j) => j === i ? { ...f, coefficient: readNumber(coefficient) } : f)); } catch (e) { onError(e.message); } });
    box.append(el('div', { class: 'dag-effect-row' }, el('label', { for: coefficient.id }, termText(t)), coefficient,
      el('button', { type: 'button', class: 'secondary', 'aria-label': `Remove response term ${i + 1}`, onclick: () => change(functions.filter((_, j) => j !== i)) }, '×')));
  }
  const left = select('dag-function-left', parents.map((p) => [p, p]), parents[0]);
  const right = select('dag-function-right', [['', 'None'], ...parents.map((p) => [p, p])], '');
  const transform = transformation('dag-function'), rightTransform = transformation('dag-function-right', 'identity');
  rightTransform.element.hidden = true;
  right.addEventListener('change', () => { rightTransform.element.hidden = !right.value; });
  const coefficient = numeric('dag-function-coefficient', .3);
  const form = el('form', { id: 'dag-function-form', novalidate: '' }, field('Parent', left), transform.element,
    field('Multiply by another parent', right), rightTransform.element, field('Coefficient', coefficient),
    el('button', { type: 'submit', class: 'secondary' }, 'Add response term'));
  form.addEventListener('submit', (event) => {
    event.preventDefault();
    try { change([...functions, { left: left.value, transform: transform.read(), coefficient: readNumber(coefficient),
      ...(right.value ? { right: right.value, rightTransform: rightTransform.read() } : {}) }]); } catch (e) { onError(e.message); }
  });
  box.append(el('details', { id: 'dag-function-add' }, el('summary', {}, 'Add a curve or transformed interaction'), form));
  // Deterministic conditional response preview; simulated data still come only from C++.
  const variable = select('dag-curve-variable', parents.map((p) => [p, p]), parents[0]);
  const lower = numeric('dag-curve-lower', -3), upper = numeric('dag-curve-upper', 3);
  const fixed = Object.fromEntries(parents.map((p) => [p, numeric(`dag-curve-fixed-${p}`, 0)]));
  const picture = el('div', { id: 'dag-curve-preview', 'aria-live': 'polite' });
  const preview = el('details', { id: 'dag-curve-panel' }, el('summary', {}, 'Preview the conditional response'),
    field('Vary parent', variable), el('div', { class: 'dag-pair' }, field('From', lower), field('To', upper)),
    el('details', {}, el('summary', {}, 'Values of other parents'), parents.map((p) => field(p, fixed[p]))), picture);
  function draw() {
    try {
      const lo = readNumber(lower), hi = readNumber(upper);
      if (!(hi > lo)) throw new Error('Preview range must be increasing.');
      const values = Object.fromEntries(parents.map((p) => [p, readNumber(fixed[p])]));
      const response = [];
      for (let i = 0; i <= 160; i++) {
        values[variable.value] = lo + (hi - lo) * i / 160;
        let eta = node.intercept;
        for (const e of graph.edges) if (e.target === node.id) eta += e.coefficient * values[e.source];
        for (const t of node.interactions) eta += t.coefficient * values[t.left] * values[t.right];
        for (const t of functions) eta += t.coefficient * transformValue(values[t.left], t.transform) * (t.right ? transformValue(values[t.right], t.rightTransform) : 1);
        response.push(node.family === 'bernoulli' ? 1 / (1 + Math.exp(-eta)) : ['gaussian', 'uniform'].includes(node.family) ? eta : Math.exp(eta + (node.family === 'lognormal' ? node.sd ** 2 / 2 : 0)));
      }
      if (response.some((x) => !Number.isFinite(x))) throw new Error('The response overflows within this preview range.');
      const min = Math.min(...response), max = Math.max(...response), span = max - min || 1;
      const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
      svg.setAttribute('viewBox', '0 0 320 150'); svg.setAttribute('role', 'img');
      svg.setAttribute('aria-label', `${node.id} conditional ${node.kind === 'survival' ? 'hazard multiplier' : 'mean'} as ${variable.value} varies from ${lo} to ${hi}`);
      const path = document.createElementNS(svg.namespaceURI, 'path');
      path.setAttribute('d', response.map((v, i) => `${i ? 'L' : 'M'}${40 + i * 260 / 160},${120 - (v - min) * 100 / span}`).join(' '));
      path.setAttribute('fill', 'none'); path.setAttribute('stroke', 'currentColor'); path.setAttribute('stroke-width', '2'); svg.append(path);
      for (const [x, y, value] of [[4, 25, max], [4, 120, min], [40, 145, lo], [275, 145, hi]]) {
        const text = document.createElementNS(svg.namespaceURI, 'text'); text.setAttribute('x', x); text.setAttribute('y', y); text.setAttribute('font-size', 11); text.setAttribute('fill', 'currentColor'); text.textContent = Number(value.toPrecision(3)); svg.append(text);
      }
      picture.replaceChildren(el('p', { class: 'small muted' }, `${node.kind === 'survival' ? 'Conditional hazard multiplier' : node.kind === 'binary' ? 'Conditional probability' : 'Conditional mean'}; other parents fixed at the values above.`), svg);
    } catch (e) { picture.replaceChildren(el('p', { class: 'error' }, e.message)); }
  }
  // drawn when opened (and as its inputs change), not on every rebuild of the variable's settings
  preview.addEventListener('input', draw); preview.addEventListener('change', draw);
  // the summary's click (or Enter) opens it: draw then, so the curve is there as the panel opens
  preview.querySelector('summary').addEventListener('click', () => { if (!preview.open) draw(); });
  preview.addEventListener('toggle', () => { if (preview.open && !picture.childElementCount) draw(); }); box.append(preview);
}
