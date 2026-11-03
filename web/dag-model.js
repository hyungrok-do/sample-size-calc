import { validateStudy, parseStudy, updateStudyReferences } from './dag-study.js?v=dev';
import { parseFunction, termText, termKey } from './dag-functions.js?v=dev';
import { CANVAS_LIMIT } from './dag-canvas.js?v=dev';
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
  student: { kind: 'continuous', label: 'Student t', link: 'center', parameters: ['sd', 'df'] },
  beta: { kind: 'continuous', label: 'Beta (proportion)', link: 'log odds of the mean', parameters: ['precision'] },
  bernoulli: { kind: 'binary', label: 'Bernoulli', link: 'log odds', parameters: [] },
  poisson: { kind: 'count', label: 'Poisson', link: 'log mean', parameters: [] },
  negbin: { kind: 'count', label: 'Negative binomial', link: 'log mean', parameters: ['overdispersion'] },
  weibull: { kind: 'survival', label: 'Weibull', link: 'log hazard multiplier', parameters: ['shape', 'scale'] },
  exponential: { kind: 'survival', label: 'Exponential', link: 'log hazard multiplier', parameters: ['scale'] },
};
export const FAMILIES = Object.keys(DISTRIBUTIONS);
const DEFAULT_FAMILY = { continuous: 'gaussian', binary: 'bernoulli', count: 'poisson', survival: 'weibull' };
const DEFAULT_PARAMETERS = { sd: 1, width: 2, shape: 2, scale: 10 };
// Kept on a variable only while its distribution uses them (older files never have them).
const EXTRA_PARAMETERS = { df: 5, precision: 10, overdispersion: 0.5 };
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
    if (node.family === 'student' && !(node.df > 1)) fail(`${node.label}: a Student t variable needs more than 1 degree of freedom, so that its mean exists.`);
    if (node.measured != null && typeof node.measured !== 'boolean') fail(`${node.label}: measured must be true or false.`);
    if (node.unknown != null) {
      const u = node.unknown, parents = new Set((graph.edges ?? []).filter((e) => e.target === node.id).map((e) => e.source));
      if (typeof u !== 'object' || !Array.isArray(u.inputs) || u.inputs.length < 1 || u.inputs.length > 8 || new Set(u.inputs).size !== u.inputs.length)
        fail(`${node.label}: a relationship of unknown form needs 1–8 different parents.`);
      if (u.inputs.some((id) => !parents.has(id))) fail(`${node.label}: a relationship of unknown form can use only the variable's parents.`);
      if (!(finite(u.size) && u.size > 0 && u.size <= 50) || !(finite(u.smoothness) && u.smoothness > 0 && u.smoothness <= 100))
        fail(`${node.label}: the unknown relationship needs a size in (0, 50] and a smoothness in (0, 100].`);
      if (!Number.isInteger(u.seed ?? 0) || (u.seed ?? 0) < 0 || (u.seed ?? 0) > 1e9) fail(`${node.label}: the unknown relationship's shape number must be a whole number.`);
    }
    if (!finite(node.x) || !finite(node.y) || node.x < 0 || node.x > CANVAS_LIMIT.x || node.y < 0 || node.y > CANVAS_LIMIT.y)
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
    if (node.functions != null) {
      if (!Array.isArray(node.functions) || node.functions.length > 120) fail('Use at most 120 additional response terms per variable.');
      const seen = new Set([...p.map((left) => termKey({ left })), ...node.interactions.map(termKey)]);
      for (const raw of node.functions) {
        const term = parseFunction(raw);
        if (!p.includes(term.left) || term.right && !p.includes(term.right)) fail('Response terms must use existing parents of the child.');
        const key = termKey(term);
        if (seen.has(key)) fail('Duplicate response term.');
        seen.add(key);
      }
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
  validateStudy(graph);
  return { order, parents, byId };
}
export function parseGraph(input) {
  const raw = typeof input === 'string' ? JSON.parse(input) : input;
  validateGraph(raw);
  return { version: DAG_VERSION, ...(raw.study ? { study: parseStudy(raw.study) } : {}),
    nodes: raw.nodes.map((n) => ({ id: n.id, label: n.label, kind: n.kind, family: n.family, intercept: n.intercept,
      ...Object.fromEntries(Object.entries(DEFAULT_PARAMETERS).map(([key, value]) => [key, finite(n[key]) && n[key] > 0 ? n[key] : value])),
      ...Object.fromEntries(Object.entries(EXTRA_PARAMETERS).filter(([key]) => DISTRIBUTIONS[n.family].parameters.includes(key))
        .map(([key, value]) => [key, finite(n[key]) && n[key] > 0 ? n[key] : value])),
      ...(n.measured === false ? { measured: false } : {}),
      ...(n.unknown ? { unknown: { inputs: [...n.unknown.inputs], size: n.unknown.size, smoothness: n.unknown.smoothness, seed: n.unknown.seed ?? 0 } } : {}),
      x: n.x, y: n.y,
      ...(raw.study ? { scope: n.scope ?? 'subject', repeated: Object.fromEntries(['subjectSD', 'slopeSD', 'slopeCorrelation', 'clusterSD', 'residualCorrelation'].map((k) => [k, n.repeated?.[k] ?? 0])),
        lags: (n.lags ?? []).map((l) => ({ source: l.source, lag: l.lag, coefficient: l.coefficient, initial: l.initial ?? 0 })) } : {}),
      censoring: n.kind === 'survival' ? { mode: n.censoring.mode,
        time: n.censoring.mode === 'administrative' ? n.censoring.time : 20,
        variable: n.censoring.mode === 'variable' ? n.censoring.variable : '' } : { mode: 'none', time: 20, variable: '' },
      ...(n.functions?.length ? { functions: n.functions.map((t) => parseFunction(t)) } : {}),
      interactions: n.interactions.map((t) => ({ left: t.left, right: t.right, coefficient: t.coefficient })) })),
    edges: raw.edges.map((e) => ({ source: e.source, target: e.target, coefficient: e.coefficient })) };
}
export function addEdge(graph, source, target, coefficient = 0.5) {
  const next = { ...graph, edges: [...graph.edges, { source, target, coefficient }] };
  validateGraph(next); return next;
}
// A variable without `id` among its unknown relationship's inputs (and without the relationship when none remain).
function withoutInput(n, id) {
  if (!n.unknown?.inputs.includes(id)) return {};
  const inputs = n.unknown.inputs.filter((x) => x !== id);
  return { unknown: inputs.length ? { ...n.unknown, inputs } : undefined };
}
export function removeNode(graph, id) {
  return { ...graph, ...(graph.study ? { study: updateStudyReferences(graph.study, id, null) } : {}), nodes: graph.nodes.filter((n) => n.id !== id).map((n) => ({ ...n,
    ...(n.lags ? { lags: n.lags.filter((l) => l.source !== id) } : {}),
    censoring: n.censoring?.mode === 'variable' && n.censoring.variable === id ? { mode: 'none', time: 20, variable: '' } : n.censoring,
    ...(n.functions ? { functions: n.functions.filter((t) => t.left !== id && t.right !== id) } : {}),
    ...withoutInput(n, id),
    interactions: n.interactions.filter((t) => t.left !== id && t.right !== id) })),
  edges: graph.edges.filter((e) => e.source !== id && e.target !== id) };
}
export function removeEdge(graph, source, target) {
  return { ...graph, edges: graph.edges.filter((e) => e.source !== source || e.target !== target),
    nodes: graph.nodes.map((n) => n.id !== target ? n : { ...n, ...(n.functions ? { functions: n.functions.filter((t) => t.left !== source && t.right !== source) } : {}), ...withoutInput(n, source), interactions: n.interactions.filter((t) => t.left !== source && t.right !== source) }) };
}
export function renameNode(graph, oldId, newId) {
  const next = { ...graph, ...(graph.study ? { study: updateStudyReferences(graph.study, oldId, newId) } : {}), nodes: graph.nodes.map((n) => ({ ...n, id: n.id === oldId ? newId : n.id,
    ...(n.lags ? { lags: n.lags.map((l) => ({ ...l, source: l.source === oldId ? newId : l.source })) } : {}),
    censoring: n.censoring?.mode === 'variable' && n.censoring.variable === oldId ? { ...n.censoring, variable: newId } : n.censoring,
    ...(n.functions ? { functions: n.functions.map((t) => ({ ...t, left: t.left === oldId ? newId : t.left, ...(t.right ? { right: t.right === oldId ? newId : t.right } : {}) })) } : {}),
    ...(n.unknown ? { unknown: { ...n.unknown, inputs: n.unknown.inputs.map((x) => x === oldId ? newId : x) } } : {}),
    interactions: n.interactions.map((t) => ({ ...t, left: t.left === oldId ? newId : t.left, right: t.right === oldId ? newId : t.right })) })),
    edges: graph.edges.map((e) => ({ ...e, source: e.source === oldId ? newId : e.source, target: e.target === oldId ? newId : e.target })) };
  validateGraph(next); return next;
}
// A parameter's value for a distribution that uses it: the variable's own, else the default.
export const parameterValue = (node, key) => finite(node[key]) && node[key] > 0 ? node[key] : (DEFAULT_PARAMETERS[key] ?? EXTRA_PARAMETERS[key]);
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
  if (name === 'trial') return { version: DAG_VERSION, nodes: [{ ...a, x: 240, y: 180 }, { ...y, x: 680, y: 180 }],
    edges: [{ source: 'A', target: 'Y', coefficient: .5 }] };
  if (['repeated', 'clustered-repeated', 'repeated-binary'].includes(name)) {
    const time = { ...makeNode('Time', 'continuous', 400, 300), label: 'Visit time', scope: 'time' };
    a.scope = name === 'clustered-repeated' ? 'cluster' : 'subject'; x.scope = 'subject';
    y.scope = 'visit'; y.repeated = { subjectSD: .7, slopeSD: .04, slopeCorrelation: .2, clusterSD: name === 'clustered-repeated' ? .4 : 0, residualCorrelation: 0 };
    if (name === 'repeated-binary') { y.kind = 'binary'; y.family = 'bernoulli'; y.intercept = -1; }
    y.interactions = [{ left: 'A', right: 'Time', coefficient: .06 }];
    return { version: DAG_VERSION, study: parseStudy({ times: [0, 3, 6, 12], clusterSizes: name === 'clustered-repeated' ? [8, 12, 16] : [] }), nodes: [x, a, time, y],
      edges: [{ source: 'X', target: 'Y', coefficient: .4 }, { source: 'A', target: 'Y', coefficient: 0 }, { source: 'Time', target: 'Y', coefficient: .1 }] };
  }
  if (name === 'nonlinear-dose') {
    const dose = { ...makeNode('A', 'continuous', 30, 180), label: 'Dose' };
    y.functions = [{ left: 'A', transform: { kind: 'square' }, coefficient: .8 }];
    return { version: DAG_VERSION, nodes: [dose, y], edges: [{ source: 'A', target: 'Y', coefficient: 0 }] };
  }
  if (name === 'nonlinear-confounding') {
    x.family = 'uniform'; x.width = 4;
    a.intercept = -.5;
    a.functions = [{ left: 'X', transform: { kind: 'square' }, coefficient: .3 }];
    y.functions = [{ left: 'X', transform: { kind: 'square' }, coefficient: .7 }];
    return { version: DAG_VERSION, nodes: [x, a, y], edges: [
      { source: 'X', target: 'A', coefficient: .2 }, { source: 'X', target: 'Y', coefficient: .2 }, { source: 'A', target: 'Y', coefficient: .4 }] };
  }
  if (name === 'heterogeneous-null') {
    y.functions = [{ left: 'X', transform: { kind: 'square' }, right: 'A', coefficient: .5 }];
    return { version: DAG_VERSION, nodes: [x, a, y], edges: [
      { source: 'X', target: 'Y', coefficient: .1 }, { source: 'A', target: 'Y', coefficient: -.5 }] };
  }
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
const shown = (x) => String(Number(x.toPrecision(6)));  // equations show six significant digits
export function equation(graph, node) {
  if (graph.study && node.scope === 'time') return `${node.id} = visit time from [${graph.study.times.join(', ')}]`;

  let predictor = shown(node.intercept);
  for (const term of [...graph.edges.filter((e) => e.target === node.id).map((e) => ({ value: e.coefficient, name: e.source })),
    ...node.interactions.map((t) => ({ value: t.coefficient, name: `${t.left} × ${t.right}` })),
    ...(node.functions ?? []).map((t) => ({ value: t.coefficient, name: termText(t) }))])
    predictor += ` ${term.value < 0 ? '−' : '+'} ${shown(Math.abs(term.value))} × ${term.name}`;
  if (node.unknown) predictor += ` + f(${node.unknown.inputs.join(', ')}) [unknown form, SD ${node.unknown.size}]`;
  if (graph.study) {
    for (const l of node.lags ?? []) predictor += ` + ${l.coefficient} × ${l.source}[visit − ${l.lag}] (initial ${l.initial ?? 0})`;
    const r = node.repeated ?? {};
    if (r.subjectSD) predictor += ` + subject intercept (SD ${r.subjectSD})`;
    if (r.slopeSD) predictor += ` + subject slope × elapsed time (SD ${r.slopeSD})`;
    if (r.clusterSD) predictor += ` + cluster intercept (SD ${r.clusterSD})`;
  }
  const hidden = node.measured === false ? ' (unmeasured)' : '';
  if (node.family === 'gaussian') return `${node.id} = ${predictor} + ε, ε ~ Normal(0, SD = ${node.sd})${hidden}`;
  if (node.family === 'uniform') return `${node.id} ~ Uniform(η − ${node.width / 2}, η + ${node.width / 2}), η = ${predictor}${hidden}`;
  if (node.family === 'lognormal') return `log(${node.id}) = ${predictor} + ε, ε ~ Normal(0, SD = ${node.sd})${hidden}`;
  if (node.family === 'gamma') return `${node.id} ~ Gamma(shape = ${node.shape}, scale = μ / ${node.shape}), log(μ) = ${predictor}${hidden}`;
  if (node.family === 'student') return `${node.id} = ${predictor} + ${node.sd} × t, t ~ Student t(${node.df} df)${hidden}`;
  if (node.family === 'beta') return `${node.id} ~ Beta(μφ, (1 − μ)φ), φ = ${node.precision}, logit(μ) = ${predictor}${hidden}`;
  if (node.family === 'bernoulli') return `${node.id} ~ Bernoulli(p), logit(p) = ${predictor}${hidden}`;
  if (node.family === 'poisson') return `${node.id} ~ Poisson(μ), log(μ) = ${predictor}${hidden}`;
  if (node.family === 'negbin') return `${node.id} ~ NegBin(μ, Var = μ + ${node.overdispersion} μ²), log(μ) = ${predictor}${hidden}`;
  const shape = node.family === 'exponential' ? 1 : node.shape;
  const obs = node.censoring.mode === 'none' ? '∞' : node.censoring.mode === 'administrative' ? node.censoring.time : node.censoring.variable;
  return `S_${node.id}(t | parents) = exp(−(t / ${node.scale})^${shape} × exp(η)), η = ${predictor}. ${node.id}_time = min(${node.id}, ${obs}), ${node.id}_event = 1{${node.id} ≤ ${obs}}${hidden}`;
}
export function dataCsv(_graph, result) {
  const ids = Object.keys(result.columns);
  return [ids.join(','), ...Array.from({ length: result.n }, (_, i) => ids.map((id) => result.columns[id][i]).join(','))].join('\n') + '\n';
}

