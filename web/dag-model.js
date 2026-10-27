// Executable distribution DAG. Variables have types and distributions, without prescribed scientific roles.
// Parent effects and interactions are explicit terms on the child's conditional distribution scale. No eval.
export const DAG_VERSION = 1;
export const MAX_NODES = 40;
export const MAX_ROWS = 50000;
export const KINDS = ['continuous', 'binary', 'count', 'survival'];
export const DISTRIBUTIONS = {
  gaussian: { kind: 'continuous', label: 'Normal', link: 'mean', parameters: ['sd'] },
  uniform: { kind: 'continuous', label: 'Uniform', link: 'center', parameters: ['width'] },
  lognormal: { kind: 'continuous', label: 'Log-normal', link: 'log median', parameters: ['sd'] },
  gamma: { kind: 'continuous', label: 'Gamma', link: 'log mean', parameters: ['shape'] },
  bernoulli: { kind: 'binary', label: 'Bernoulli', link: 'log odds', parameters: [] },
  poisson: { kind: 'count', label: 'Poisson', link: 'log mean', parameters: [] },
  weibull: { kind: 'survival', label: 'Weibull', link: 'log hazard multiplier', parameters: ['shape', 'scale'] },
  exponential: { kind: 'survival', label: 'Exponential', link: 'log hazard multiplier', parameters: ['scale'] },
};
export const FAMILIES = Object.keys(DISTRIBUTIONS);
const DEFAULT_FAMILY = { continuous: 'gaussian', binary: 'bernoulli', count: 'poisson', survival: 'weibull' };
const DEFAULT_PARAMETERS = { sd: 1, width: 2, shape: 2, scale: 10 };
const idPattern = /^[A-Za-z][A-Za-z0-9_]{0,31}$/;
const finite = (x) => typeof x === 'number' && Number.isFinite(x);
const fail = (message) => { throw new Error(message); };

export function topologicalOrder(graph) {
  const nodes = graph.nodes ?? [], edges = graph.edges ?? [];
  const indegree = new Map(nodes.map((node) => [node.id, 0]));
  const children = new Map(nodes.map((node) => [node.id, []]));
  for (const edge of edges) {
    if (!indegree.has(edge.source) || !indegree.has(edge.target)) fail('An arrow refers to a missing variable.');
    indegree.set(edge.target, indegree.get(edge.target) + 1);
    children.get(edge.source).push(edge.target);
  }
  const queue = nodes.filter((node) => indegree.get(node.id) === 0).map((node) => node.id);
  for (let i = 0; i < queue.length; i++) for (const child of children.get(queue[i])) {
    indegree.set(child, indegree.get(child) - 1);
    if (indegree.get(child) === 0) queue.push(child);
  }
  if (queue.length !== nodes.length) fail('That arrow creates a cycle. A DAG must have no directed cycles.');
  return queue;
}
export function validateGraph(graph) {
  if (!graph || graph.version !== DAG_VERSION) fail('Unsupported DAG file version.');
  if (!Array.isArray(graph.nodes) || graph.nodes.length > MAX_NODES) fail(`Use at most ${MAX_NODES} variables.`);
  if (!Array.isArray(graph.edges) || graph.edges.length > MAX_NODES * (MAX_NODES - 1) / 2) fail('The arrow list is invalid or too large.');
  const ids = new Set();
  for (const node of graph.nodes) {
    if (!node || typeof node.id !== 'string' || !idPattern.test(node.id) || ids.has(node.id))
      fail('Variable IDs must be unique letters, digits and underscores, starting with a letter (up to 32 characters).');
    ids.add(node.id);
    if (typeof node.label !== 'string' || !node.label.trim() || node.label.length > 80) fail('Each variable needs a label of 1–80 characters.');
    if (!KINDS.includes(node.kind) || !FAMILIES.includes(node.family) || DISTRIBUTIONS[node.family].kind !== node.kind)
      fail(`${node.id}: choose a distribution compatible with the variable type.`);
    if (!finite(node.intercept)) fail(`${node.label}: enter a finite intercept.`);
    for (const param of DISTRIBUTIONS[node.family].parameters)
      if (!finite(node[param]) || node[param] <= 0) fail(`${node.label}: ${param} must be positive and finite.`);
    if (!finite(node.x) || !finite(node.y) || node.x < 0 || node.x > 920 || node.y < 0 || node.y > 440)
      fail(`${node.label}: invalid canvas position.`);
    if (!Array.isArray(node.interactions) || node.interactions.length > MAX_NODES * MAX_NODES) fail(`${node.label}: invalid interaction list.`);
    if (node.kind === 'survival') {
      const censor = node.censoring;
      if (!censor || !['none', 'administrative', 'variable'].includes(censor.mode)) fail(`${node.label}: choose a censoring observation rule.`);
      if (censor.mode === 'administrative' && (!finite(censor.time) || censor.time <= 0)) fail(`${node.label}: the administrative cutoff must be positive.`);
    }
  }
  const pairs = new Set();
  for (const edge of graph.edges) {
    if (!edge || !ids.has(edge.source) || !ids.has(edge.target)) fail('An arrow refers to a missing variable.');
    if (edge.source === edge.target) fail('A variable cannot point to itself.');
    const pair = `${edge.source}>${edge.target}`;
    if (pairs.has(pair)) fail('These variables are already connected in that direction.');
    pairs.add(pair);
    if (!finite(edge.coefficient)) fail('Every arrow needs a finite coefficient.');
  }
  const order = topologicalOrder(graph);
  const byId = new Map(graph.nodes.map((node) => [node.id, node]));
  const parents = new Map(graph.nodes.map((node) => [node.id, graph.edges.filter((edge) => edge.target === node.id).map((edge) => edge.source)]));
  const outputs = new Set(ids);
  for (const node of graph.nodes) {
    const p = parents.get(node.id), terms = new Set();
    for (const term of node.interactions) {
      if (!term || !p.includes(term.left) || !p.includes(term.right) || term.left === term.right || !finite(term.coefficient))
        fail(`${node.label}: interactions require two different parents and a finite coefficient.`);
      const pair = [term.left, term.right].sort().join('*');
      if (terms.has(pair)) fail(`${node.label}: that interaction already exists.`);
      terms.add(pair);
    }
    if (node.kind === 'survival') {
      if (node.censoring.mode === 'variable') {
        const other = byId.get(node.censoring.variable);
        if (!other || other.id === node.id)
          fail(`${node.label}: select a different variable for censoring.`);
      }
      for (const name of [`${node.id}_time`, `${node.id}_event`]) {
        if (outputs.has(name)) fail(`The derived survival column ${name} conflicts with a variable ID.`);
        outputs.add(name);
      }
    }
  }
  return { order, parents, byId };
}
export function parseGraph(input) {
  const raw = typeof input === 'string' ? JSON.parse(input) : input;
  validateGraph(raw);
  return { version: DAG_VERSION,
    nodes: raw.nodes.map((n) => ({ id: n.id, label: n.label, kind: n.kind, family: n.family, intercept: n.intercept,
      ...Object.fromEntries(Object.entries(DEFAULT_PARAMETERS).map(([key, value]) => [key, finite(n[key]) && n[key] > 0 ? n[key] : value])),
      x: n.x, y: n.y,
      censoring: n.kind === 'survival' ? { mode: n.censoring.mode,
        time: n.censoring.mode === 'administrative' ? n.censoring.time : 20,
        variable: n.censoring.mode === 'variable' ? n.censoring.variable : '' } : { mode: 'none', time: 20, variable: '' },
      interactions: n.interactions.map((t) => ({ left: t.left, right: t.right, coefficient: t.coefficient })) })),
    edges: raw.edges.map((e) => ({ source: e.source, target: e.target, coefficient: e.coefficient })) };
}
export function addEdge(graph, source, target, coefficient = 0.5) {
  const next = { ...graph, edges: [...graph.edges, { source, target, coefficient }] };
  validateGraph(next); return next;
}
export function removeNode(graph, id) {
  return { ...graph, nodes: graph.nodes.filter((n) => n.id !== id).map((n) => ({ ...n,
    censoring: n.censoring?.mode === 'variable' && n.censoring.variable === id ? { mode: 'none', time: 20, variable: '' } : n.censoring,
    interactions: n.interactions.filter((t) => t.left !== id && t.right !== id) })),
  edges: graph.edges.filter((e) => e.source !== id && e.target !== id) };
}
export function removeEdge(graph, source, target) {
  return { ...graph, edges: graph.edges.filter((e) => e.source !== source || e.target !== target),
    nodes: graph.nodes.map((n) => n.id !== target ? n : { ...n, interactions: n.interactions.filter((t) => t.left !== source && t.right !== source) }) };
}
export function renameNode(graph, oldId, newId) {
  const next = { ...graph, nodes: graph.nodes.map((n) => ({ ...n, id: n.id === oldId ? newId : n.id,
    censoring: n.censoring?.mode === 'variable' && n.censoring.variable === oldId ? { ...n.censoring, variable: newId } : n.censoring,
    interactions: n.interactions.map((t) => ({ ...t, left: t.left === oldId ? newId : t.left, right: t.right === oldId ? newId : t.right })) })),
    edges: graph.edges.map((e) => ({ ...e, source: e.source === oldId ? newId : e.source, target: e.target === oldId ? newId : e.target })) };
  validateGraph(next); return next;
}
export function makeNode(id, kind = 'continuous', x = 30, y = 30) {
  return { id, label: id, kind, family: DEFAULT_FAMILY[kind], intercept: 0, ...DEFAULT_PARAMETERS,
    censoring: { mode: 'none', time: 20, variable: '' }, interactions: [], x, y };
}
export function emptyGraph() { return { version: DAG_VERSION, nodes: [], edges: [] }; }
export function exampleGraph(name = 'mediation') {
  if (name === 'empty') return emptyGraph();
  const x = { ...makeNode('X', 'continuous', 30, 60), label: 'Covariate' };
  const a = { ...makeNode('A', 'binary', 30, 300), label: 'Treatment' };
  const m = { ...makeNode('M', 'continuous', 430, 300), label: 'Mediator' };
  const y = { ...makeNode('Y', name === 'binary' ? 'binary' : 'continuous', 830, 60), label: 'Outcome' };
  if (name === 'binary') {
    y.intercept = -1;
    return { version: DAG_VERSION, nodes: [x, a, y], edges: [
      { source: 'X', target: 'Y', coefficient: 0.6 }, { source: 'A', target: 'Y', coefficient: -0.5 }] };
  }
  if (name === 'survival') {
    const t = { ...makeNode('T', 'survival', 830, 60), label: 'Event time', shape: 1.6,
      censoring: { mode: 'variable', time: 20, variable: 'C' } };
    const c = { ...makeNode('C', 'survival', 830, 300), label: 'Censoring time', family: 'exponential', scale: 18 };
    return { version: DAG_VERSION, nodes: [x, a, t, c], edges: [
      { source: 'X', target: 'T', coefficient: 0.35 }, { source: 'A', target: 'T', coefficient: Math.log(0.7) }] };
  }
  return { version: DAG_VERSION, nodes: [x, a, m, y], edges: [
    { source: 'X', target: 'M', coefficient: 0.4 }, { source: 'X', target: 'Y', coefficient: 0.4 },
    { source: 'A', target: 'M', coefficient: 0.8 }, { source: 'M', target: 'Y', coefficient: 0.6 },
    { source: 'A', target: 'Y', coefficient: 0.3 }] };
}
export function equation(graph, node) {
  let predictor = String(node.intercept);
  for (const term of [...graph.edges.filter((e) => e.target === node.id).map((e) => ({ value: e.coefficient, name: e.source })),
    ...node.interactions.map((t) => ({ value: t.coefficient, name: `${t.left} × ${t.right}` }))])
    predictor += ` ${term.value < 0 ? '−' : '+'} ${Math.abs(term.value)} × ${term.name}`;
  if (node.family === 'gaussian') return `${node.id} = ${predictor} + ε, ε ~ Normal(0, SD = ${node.sd})`;
  if (node.family === 'uniform') return `${node.id} ~ Uniform(η − ${node.width / 2}, η + ${node.width / 2}), η = ${predictor}`;
  if (node.family === 'lognormal') return `log(${node.id}) = ${predictor} + ε, ε ~ Normal(0, SD = ${node.sd})`;
  if (node.family === 'gamma') return `${node.id} ~ Gamma(shape = ${node.shape}, scale = μ / ${node.shape}), log(μ) = ${predictor}`;
  if (node.family === 'bernoulli') return `${node.id} ~ Bernoulli(p), logit(p) = ${predictor}`;
  if (node.family === 'poisson') return `${node.id} ~ Poisson(μ), log(μ) = ${predictor}`;
  const shape = node.family === 'exponential' ? 1 : node.shape;
  const obs = node.censoring.mode === 'none' ? '∞' : node.censoring.mode === 'administrative' ? node.censoring.time : node.censoring.variable;
  return `S_${node.id}(t | parents) = exp(−(t / ${node.scale})^${shape} × exp(η)), η = ${predictor}. ${node.id}_time = min(${node.id}, ${obs}), ${node.id}_event = 1{${node.id} ≤ ${obs}}`;
}
export function dataCsv(_graph, result) {
  const ids = Object.keys(result.columns);
  return [ids.join(','), ...Array.from({ length: result.n }, (_, i) => ids.map((id) => result.columns[id][i]).join(','))].join('\n') + '\n';
}
