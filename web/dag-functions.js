// Explicit basis functions shared by the editor, equations and analysis preview.
// The C++ engine independently validates and evaluates every term and fit.
export const TRANSFORMS = [['identity', 'Linear'], ['square', 'Quadratic'], ['threshold', 'Threshold'],
  ['hinge', 'Hinge'], ['saturation', 'Saturation (tanh)'], ['spline', 'Restricted cubic spline']];
export const BASES = [['match', 'Match the generating terms'], ['linear', 'Linear main effects'],
  ['quadratic', 'Linear + quadratic terms'], ['spline', 'Restricted cubic splines']];
const fail = (s) => { throw new Error(s); };
const finite = (x) => typeof x === 'number' && Number.isFinite(x);
export function parseTransform(raw = null) {
  if (raw == null) return { kind: 'identity' };
  if (typeof raw !== 'object' || Array.isArray(raw) || !TRANSFORMS.some(([k]) => k === raw.kind)) fail('Unknown response transformation.');
  const { kind } = raw, t = { kind };
  for (const [key, fallback] of [['center', 0], ['scale', 1], ['knot', 0]]) {
    const v = Object.hasOwn(raw, key) ? raw[key] : fallback;
    if (!finite(v) || key === 'scale' && v <= 0) fail('Transformation parameters must be finite and scale must be positive.');
    if (['square', 'saturation'].includes(kind) && ['center', 'scale'].includes(key) || ['threshold', 'hinge'].includes(kind) && key === 'knot') t[key] = v;
  }
  if (kind === 'spline') {
    const k = raw.knots;
    if (!Array.isArray(k) || k.length < 3 || k.length > 7 || k.some((x, i) => !finite(x) || i && x <= k[i - 1]))
      fail('A restricted cubic spline needs 3–7 finite, strictly increasing knots.');
    t.knots = [...k]; t.index = Object.hasOwn(raw, 'index') ? raw.index : 0;
    if (!Number.isInteger(t.index) || t.index < 0 || t.index > k.length - 3) fail('Choose a valid spline basis index.');
  }
  return t;
}
export function transformValue(x, raw) {
  const t = parseTransform(raw);
  if (t.kind === 'identity') return x;
  if (t.kind === 'square') return ((x - t.center) / t.scale) ** 2;
  if (t.kind === 'saturation') return Math.tanh((x - t.center) / t.scale);
  if (t.kind === 'threshold') return Number(x > t.knot);
  if (t.kind === 'hinge') return Math.max(0, x - t.knot);
  const k = t.knots, a = k[t.index], b = k.at(-2), c = k.at(-1), width = c - k[0], z = Math.min(x, c);
  const u = Math.max(0, z - a) / width, v = Math.max(0, z - b) / width, w = Math.max(0, z - c) / width;
  const r = (c - a) / (c - b), s = (b - a) / (c - b);
  return width * (u ** 3 - r * v ** 3 + s * w ** 3) + (x > c ? 3 * (u * u - r * v * v) * (x - c) : 0);
}
export function transformText(id, raw) {
  const t = parseTransform(raw);
  if (t.kind === 'identity') return id;
  if (['square', 'saturation'].includes(t.kind)) return `${t.kind === 'square' ? 'square' : 'tanh'}((${id} - ${t.center}) / ${t.scale})`;
  if (t.kind === 'threshold') return `I(${id} > ${t.knot})`;
  if (t.kind === 'hinge') return `max(0, ${id} - ${t.knot})`;
  return `rcs${t.index + 1}(${id}; [${t.knots.join(',')}])`;
}
export const termText = (t) => transformText(t.left, t.transform) + (t.right ? ` × ${transformText(t.right, t.rightTransform)}` : '');
export const nonlinearIn = (t, id) => t.left === id && (t.transform?.kind ?? 'identity') !== 'identity'
  || t.right === id && (t.rightTransform?.kind ?? 'identity') !== 'identity' || t.left === id && t.right === id;
export function parseFunction(raw, coefficient = true) {
  if (!raw || typeof raw.left !== 'string' || raw.right != null && typeof raw.right !== 'string') fail('A response term needs parent identifiers.');
  const t = { left: raw.left }, left = parseTransform(raw.transform), right = parseTransform(raw.rightTransform);
  if (raw.right) t.right = raw.right;
  if (left.kind !== 'identity') t.transform = left;
  if (t.right && right.kind !== 'identity') t.rightTransform = right;
  if (coefficient) { if (!finite(raw.coefficient)) fail('Response coefficients must be finite.'); t.coefficient = raw.coefficient; }
  return t;
}
export function termKey(t) {
  const a = `${t.left}:${JSON.stringify(parseTransform(t.transform))}`;
  return t.right ? [a, `${t.right}:${JSON.stringify(parseTransform(t.rightTransform))}`].sort().join('*') : a;
}
export const nodeTerms = (graph, node) => [...graph.edges.filter((e) => e.target === node.id).map((e) => ({ left: e.source })),
  ...node.interactions.map(({ left, right }) => ({ left, right })), ...(node.functions ?? []).map((t) => parseFunction(t, false))];
function parseBasis(raw = {}) {
  if (!raw || typeof raw !== 'object' || Array.isArray(raw)) fail('An analysis model must be an object.');
  const basis = raw.basis ?? 'match';
  if (![...BASES.map(([x]) => x), 'custom'].includes(basis)) fail('Choose matching, linear, quadratic, spline or custom analysis terms.');
  const out = { basis };
  if (raw.interactions != null) { if (typeof raw.interactions !== 'boolean') fail('Analysis interactions must be true or false.'); out.interactions = raw.interactions; }
  if (basis === 'spline') out.knots = parseTransform({ kind: 'spline', knots: raw.knots ?? [-1, 0, 1] }).knots;
  if (basis === 'custom') {
    if (!Array.isArray(raw.terms) || raw.terms.length > 120) fail('A custom analysis needs an explicit list of at most 120 terms.');
    out.terms = raw.terms.map((t) => parseFunction(t, false));
  }
  return out;
}
export function parseAnalysis(raw = {}) {
  const out = { ...parseBasis(raw), method: raw.method ?? 'auto', folds: raw.folds ?? 1 };
  if (!['auto', 'gformula', 'aipw', 'tmle', 'mmrm', 'lmm', 'trajectory'].includes(out.method)) fail('Choose automatic analysis, g-computation, AIPW or TMLE.');
  if (!Number.isInteger(out.folds) || out.folds < 1 || out.folds > 10) fail('Use 1–10 cross-fitting folds.');
  out.propensityBound = raw.propensityBound ?? 0;
  if (!Number.isFinite(out.propensityBound) || out.propensityBound < 0 || out.propensityBound >= .5) fail('The propensity bound must be in [0, 0.5).');
  if (raw.outcomeBounds != null) {
    if (!Array.isArray(raw.outcomeBounds) || raw.outcomeBounds.length !== 2 || !raw.outcomeBounds.every(Number.isFinite) || raw.outcomeBounds[0] >= raw.outcomeBounds[1]) fail('Supply finite, increasing prespecified outcome bounds.');
    out.outcomeBounds = [...raw.outcomeBounds];
  }
  if (raw.timeBasis != null) { if (!['categorical', 'linear'].includes(raw.timeBasis)) fail('Choose categorical or linear analysis time.'); out.timeBasis = raw.timeBasis; }
  if (raw.randomSlope != null) { if (typeof raw.randomSlope !== 'boolean') fail('Random slope must be true or false.'); out.randomSlope = raw.randomSlope; }
  if (raw.propensity) out.propensity = parseBasis(raw.propensity);
  if (raw.models) {
    if (typeof raw.models !== 'object' || Array.isArray(raw.models) || Object.keys(raw.models).length > 40) fail('Invalid analysis model overrides.');
    out.models = Object.fromEntries(Object.entries(raw.models).map(([id, v]) => [id, parseBasis(v)]));
  }
  return out;
}
export function analysisTerms(graph, terms, spec, id, propensity = false) {
  const cfg = spec.models?.[id] ?? (propensity && spec.propensity ? spec.propensity : spec);
  if (cfg.basis === 'match') return terms;
  let out = [];
  if (cfg.basis === 'custom') out = cfg.terms;
  else {
    // C++ sorts by the graph's node index.
    const chosen = new Set(terms.flatMap((t) => [t.left, ...(t.right ? [t.right] : [])]));
    const variables = graph.nodes.filter((n) => chosen.has(n.id));
    for (const n of variables) {
      out.push({ left: n.id });
      if (n.kind === 'binary') continue;
      if (cfg.basis === 'quadratic') out.push({ left: n.id, transform: { kind: 'square', center: 0, scale: 1 } });
      if (cfg.basis === 'spline') for (let index = 0; index < cfg.knots.length - 2; index++) out.push({ left: n.id, transform: { kind: 'spline', knots: cfg.knots, index } });
    }
    if (cfg.interactions) variables.forEach((a, i) => variables.slice(i + 1).forEach((b) => out.push({ left: a.id, right: b.id })));
  }
  if (out.length > 120) fail('An analysis model supports at most 120 terms.');
  if (out.some((t) => !graph.nodes.some((n) => n.id === t.left) || t.right && !graph.nodes.some((n) => n.id === t.right))) fail('An analysis term refers to a missing variable.');
  return out;
}

// Preserve independently specified analysis terms when graph identifiers change.
export function reconcileAnalysis(raw, graph, from = null, to = null) {
  const spec = parseAnalysis(raw), ids = new Set(graph.nodes.map((n) => n.id));
  const rename = (id) => id === from ? to : id;
  function model(cfg) {
    if (cfg.basis !== 'custom') return cfg;
    return { ...cfg, terms: cfg.terms.map((t) => ({ ...t, left: rename(t.left), ...(t.right ? { right: rename(t.right) } : {}) }))
      .filter((t) => ids.has(t.left) && (!t.right || ids.has(t.right))) };
  }
  const out = { ...model(spec) };
  if (out.propensity) out.propensity = model(out.propensity);
  if (out.models) out.models = Object.fromEntries(Object.entries(out.models).map(([id, cfg]) => [rename(id), model(cfg)]).filter(([id]) => ids.has(id)));
  return out;
}
