// The DAG canvas's geometry and reading of the graph: where boxes, arrows, plates and shared random effects are drawn,
// how the view zooms and fits, how the graph is laid out, which arrows may be drawn, what role each variable plays for
// the compared variable and the outcome, and how the study's levels (clusters, subjects, visits) are switched on and
// off.  Plain functions of data, no DOM, so they are tested under Node (tests/canvas.test.mjs).

export const NODE_W = 150, NODE_H = 56;
// Positions are kept nonnegative and within this area (the saved design's x and y).  Earlier versions allowed
// 0–920 × 0–440, which stays valid.
export const CANVAS_LIMIT = { x: 3600, y: 2400 };
export const ASPECT = 580 / 1100;  // the drawing's height over its width
const PAD = 22, LABEL = 30, EFFECT_GAP = 56, EMPTY_W = NODE_W + 40, EMPTY_H = NODE_H;

const scopeOf = (node) => node.scope ?? 'subject';
// cluster 0, subject 1, visit 2; visit time is measured at the visits
export const levelOf = (graph, node) => !graph.study ? 1 : { cluster: 0, subject: 1, visit: 2, time: 2 }[scopeOf(node)] ?? 1;
const LEVELS = ['cluster', 'subject', 'visit'];

function bounds(boxes) {
  if (!boxes.length) return null;
  const x = Math.min(...boxes.map((b) => b.x)), y = Math.min(...boxes.map((b) => b.y));
  return { x, y, w: Math.max(...boxes.map((b) => b.x + b.w)) - x, h: Math.max(...boxes.map((b) => b.y + b.h)) - y };
}
const nodeBox = (n) => ({ x: n.x, y: n.y, w: NODE_W, h: NODE_H });
const around = (box, extraTop = 0) => ({ x: box.x - PAD, y: box.y - PAD - LABEL - extraTop, w: box.w + 2 * PAD, h: box.h + 2 * PAD + LABEL + extraTop });
const short = (x) => String(Number(Number(x).toPrecision(3))).replace('-', '−');
const plural = (n, word) => `${n} ${word}${n === 1 ? '' : 's'}`;

// The shared random effects drawn as small circles (plate notation): a subject's intercept and slope just outside the
// visit plate, a cluster's intercept just outside the subject plate, each above the variable it shifts.
function effectsOf(node) {
  const r = node.repeated ?? {}, list = [];
  if (r.subjectSD) list.push({ key: 'subjectSD', level: 'subject', symbol: 'u₀', text: `Subject random intercept on ${node.id}, SD ${short(r.subjectSD)}` });
  if (r.slopeSD) list.push({ key: 'slopeSD', level: 'subject', symbol: 'u₁', text: `Subject random slope on ${node.id}, SD ${short(r.slopeSD)}${r.slopeCorrelation ? `, correlation with the intercept ${short(r.slopeCorrelation)}` : ''}` });
  if (r.clusterSD) list.push({ key: 'clusterSD', level: 'cluster', symbol: 'v', text: `Cluster random intercept on ${node.id}, SD ${short(r.clusterSD)}` });
  return list;
}

// Nested plates for a structured study, outermost first: clusters (when sizes are set) ⊃ subjects ⊃ visits.  A
// plate holds the boxes of the variables measured at its level and the plate inside it; an empty plate is placed to
// the right as somewhere to drop a variable.  Also returns the random-effect circles, positioned on the plates.
export function plateLayout(graph) {
  if (!graph.study) return { plates: [], effects: [] };
  const members = Object.fromEntries(LEVELS.map((l) => [l, graph.nodes.filter((n) => LEVELS[levelOf(graph, n)] === l)]));
  const all = bounds(graph.nodes.map(nodeBox)) ?? { x: 40, y: 140, w: NODE_W, h: NODE_H };
  let right = all.x + all.w;
  const empty = () => { const box = { x: right + 60, y: all.y, w: EMPTY_W, h: EMPTY_H }; right = box.x + box.w; return box; };
  const times = graph.study.times, clusters = graph.study.clusterSizes ?? [];
  const plates = [], effects = [];
  const visitNodes = members.visit;
  let inner = null;
  if (visitNodes.length || times.length > 1) {
    const content = bounds(visitNodes.map(nodeBox)) ?? empty();
    inner = { level: 'visit', ...around(content), members: visitNodes.map((n) => n.id), empty: !visitNodes.length,
      label: `Every visit · t = ${times.join(', ')}`, multiplicity: plural(times.length, 'visit') };
    plates.push(inner);
  }
  const subjectEffects = visitNodes.flatMap((n) => effectsOf(n).filter((e) => e.level === 'subject').map((e, k, list) => ({ ...e, node: n.id,
    cx: n.x + NODE_W / 2 + (k - (list.length - 1) / 2) * 34, cy: inner.y - EFFECT_GAP / 2 })));
  const subjectContent = bounds([...members.subject.map(nodeBox), ...(inner ? [inner] : []),
    ...subjectEffects.map((e) => ({ x: e.cx - 14, y: e.cy - 14, w: 28, h: 28 }))]) ?? empty();
  const subject = { level: 'subject', ...around(subjectContent), members: members.subject.map((n) => n.id),
    empty: !members.subject.length && !inner, label: clusters.length ? `Every subject · ${clusters.join(', ')} per cluster` : 'Every subject',
    multiplicity: clusters.length ? 'subjects in each cluster' : 'sample size counts subjects' };
  plates.unshift(subject);
  effects.push(...subjectEffects);
  if (clusters.length) {
    const clusterEffects = visitNodes.flatMap((n) => effectsOf(n).filter((e) => e.level === 'cluster').map((e) => ({ ...e, node: n.id,
      cx: n.x + NODE_W / 2, cy: subject.y - EFFECT_GAP / 2 })));
    const content = bounds([...members.cluster.map(nodeBox), subject, ...clusterEffects.map((e) => ({ x: e.cx - 14, y: e.cy - 14, w: 28, h: 28 }))]);
    plates.unshift({ level: 'cluster', ...around(content), members: members.cluster.map((n) => n.id), empty: !members.cluster.length,
      label: `Every cluster · ${clusters.join(', ')} subjects`, multiplicity: 'sample size counts clusters' });
    effects.push(...clusterEffects);
  }
  return { plates, effects };
}

// The level a point falls in: the innermost plate containing it, or null outside every plate.
export function plateAt(plates, point) {
  let found = null;
  for (const p of plates) if (point.x >= p.x && point.x <= p.x + p.w && point.y >= p.y && point.y <= p.y + p.h) found = p;
  return found?.level ?? null;
}

// The graph with one variable measured at another level, its repeated settings and lags cleared when it leaves the
// visits.  Visit time stays where it is.  The caller validates the result.
export function changeScope(graph, id, level) {
  return { ...graph, nodes: graph.nodes.map((n) => {
    if (n.id !== id || scopeOf(n) === 'time' || scopeOf(n) === level) return n;
    const zero = { subjectSD: 0, slopeSD: 0, slopeCorrelation: 0, clusterSD: 0, residualCorrelation: 0 };
    return { ...n, scope: level, ...(level !== 'visit' ? { repeated: zero, lags: [] } : {}) };
  }) };
}

// The study's levels switched on or off.  Visits on keeps the schedule (or starts 0, 3, 6, 12); off keeps the first
// time and moves visit variables to the subject.  Clusters on keeps the sizes (or starts 10 per cluster); off moves
// cluster variables to the subject.  Both off is independent observations.
export function setStructure(graph, { clusters, visits }) {
  const zero = { subjectSD: 0, slopeSD: 0, slopeCorrelation: 0, clusterSD: 0, residualCorrelation: 0 };
  if (!clusters && !visits) {
    const { study, ...rest } = graph;
    return { ...rest, nodes: graph.nodes.map(({ scope, repeated, lags, ...n }) => n) };
  }
  const old = graph.study, oldTimes = old?.times ?? [];
  const times = visits ? (oldTimes.length > 1 ? [...oldTimes] : [0, 3, 6, 12]) : [oldTimes[0] ?? 0];
  const clusterSizes = clusters ? (old?.clusterSizes?.length ? [...old.clusterSizes] : [10]) : [];
  const observation = visits && old?.observation ? { ...old.observation } : { mode: 'none', intercept: -3, historyVariable: '', historyCoefficient: 0,
    currentCoefficient: 0, treatmentVariable: '', treatmentCoefficient: 0 };
  return { ...graph, study: { times, clusterSizes, observation }, nodes: graph.nodes.map((n) => {
    const scope = scopeOf(n);
    const next = { ...n, scope, repeated: { ...zero, ...(n.repeated ?? {}) }, lags: [...(n.lags ?? [])] };
    if (!clusters) next.repeated.clusterSD = 0;
    if (!clusters && scope === 'cluster') next.scope = 'subject';
    if (!visits && (scope === 'visit' || scope === 'time')) { next.scope = 'subject'; next.repeated = { ...zero }; next.lags = []; }
    return next;
  }) };
}

// An arrow between two boxes: from the side facing the other box, curved, optionally bowed sideways by `bow`
// (lagged links drawn beside a same-time arrow).  Returns the path and its midpoint.
export function edgePath(source, target, bow = 0) {
  const W = NODE_W, H = NODE_H;
  const dx = target.x - source.x, dy = target.y - source.y;
  const vertical = Math.abs(dy) > Math.abs(dx) * 0.6;
  const x1 = source.x + (vertical ? W / 2 : dx >= 0 ? W : 0);
  const y1 = source.y + (vertical ? dy >= 0 ? H : 0 : H / 2);
  const x2 = target.x + (vertical ? W / 2 : dx >= 0 ? 0 : W);
  const y2 = target.y + (vertical ? dy >= 0 ? 0 : H : H / 2);
  const bend = Math.max(30, Math.abs(vertical ? y2 - y1 : x2 - x1) * 0.4);
  const [c1, c2] = vertical ? [[x1 + bow, y1 + Math.sign(dy || 1) * bend], [x2 + bow, y2 - Math.sign(dy || 1) * bend]]
    : [[x1 + Math.sign(dx || 1) * bend, y1 + bow], [x2 - Math.sign(dx || 1) * bend, y2 + bow]];
  const d = `M ${x1} ${y1} C ${c1[0]} ${c1[1]}, ${c2[0]} ${c2[1]}, ${x2} ${y2}`;
  const mid = [0, 1].map((k) => ([x1, y1][k] + 3 * c1[k] + 3 * c2[k] + [x2, y2][k]) / 8);  // the curve at t = 1/2
  return { d, mid };
}
// A variable's own earlier visit: a loop over the box's top right corner.
export function loopPath(node) {
  const x = node.x + NODE_W - 36, y = node.y;
  return { d: `M ${x} ${y} C ${x - 6} ${y - 52}, ${x + 62} ${y - 52}, ${node.x + NODE_W} ${y + 24}`, mid: [x + 22, y - 40] };
}

// Columns by observation level, then generation depth; rows ordered to put each box near its parents and children
// (a few barycentre sweeps).  Structured designs leave room for the plates around the levels.
export function layout(graph) {
  if (!graph.nodes.length) return graph;
  const parents = new Map(graph.nodes.map((n) => [n.id, []])), children = new Map(graph.nodes.map((n) => [n.id, []]));
  for (const e of graph.edges) { parents.get(e.target).push(e.source); children.get(e.source).push(e.target); }
  const depth = new Map(), pending = new Set(graph.nodes.map((n) => n.id));
  while (pending.size) for (const id of pending) {
    if (parents.get(id).some((p) => !depth.has(p))) continue;
    depth.set(id, Math.max(0, ...parents.get(id).map((p) => depth.get(p) + 1))); pending.delete(id);
  }
  const level = new Map(graph.nodes.map((n) => [n.id, levelOf(graph, n)]));
  const keys = [...new Set(graph.nodes.map((n) => `${level.get(n.id)}:${String(depth.get(n.id)).padStart(3, '0')}`))].sort();
  const columnOf = new Map(graph.nodes.map((n) => [n.id, keys.indexOf(`${level.get(n.id)}:${String(depth.get(n.id)).padStart(3, '0')}`)]));
  const columns = keys.map((_, c) => graph.nodes.filter((n) => columnOf.get(n.id) === c).map((n) => n.id));
  // start in the current vertical order, then sweep
  const original = new Map(graph.nodes.map((n) => [n.id, n.y]));
  for (const col of columns) col.sort((a, b) => original.get(a) - original.get(b));
  const rank = new Map();
  const index = () => columns.forEach((col) => col.forEach((id, i) => rank.set(id, i - (col.length - 1) / 2)));
  index();
  const centre = (ids) => ids.length ? ids.reduce((s, id) => s + rank.get(id), 0) / ids.length : null;
  for (let sweep = 0; sweep < 4; sweep++) {
    const order = sweep % 2 ? [...columns].reverse() : columns;
    for (const col of order) {
      const at = new Map(col.map((id) => [id, centre(sweep % 2 ? children.get(id) : parents.get(id)) ?? rank.get(id)]));
      col.sort((a, b) => at.get(a) - at.get(b) || rank.get(a) - rank.get(b));
      index();
    }
  }
  const structured = Boolean(graph.study), gapX = NODE_W + 90, gapY = NODE_H + 52;
  const levelGap = structured ? 2 * PAD + 30 : 0;
  const top = structured ? 40 + (graph.study.clusterSizes?.length ? 2 : 1) * (PAD + LABEL) + 2 * EFFECT_GAP : 40;
  const tallest = Math.max(...columns.map((c) => c.length));
  const x = [], levels = keys.map((k) => Number(k.split(':')[0]));
  for (let c = 0; c < columns.length; c++)
    x.push(c ? x[c - 1] + gapX + (levels[c] !== levels[c - 1] ? levelGap : 0) : structured ? 60 : 40);
  const placed = new Map();
  columns.forEach((col, c) => col.forEach((id, i) => placed.set(id, {
    x: Math.min(CANVAS_LIMIT.x, Math.round(x[c])),
    y: Math.min(CANVAS_LIMIT.y, Math.round(top + ((tallest - col.length) / 2 + i) * gapY)) })));
  return { ...graph, nodes: graph.nodes.map((n) => ({ ...n, ...placed.get(n.id) })) };
}

// A free spot for a new box: the first position near `near` (spiralling out on a grid) that overlaps no box.
export function freeSpot(graph, near = { x: 40, y: 40 }) {
  const clash = (x, y) => graph.nodes.some((n) => Math.abs(n.x - x) < NODE_W + 20 && Math.abs(n.y - y) < NODE_H + 20);
  const clamp = (p) => ({ x: Math.max(0, Math.min(CANVAS_LIMIT.x, Math.round(p.x))), y: Math.max(0, Math.min(CANVAS_LIMIT.y, Math.round(p.y))) });
  const start = clamp(near);
  for (let ring = 0; ring < 30; ring++)
    for (let i = -ring; i <= ring; i++) for (let j = -ring; j <= ring; j++) {
      if (Math.max(Math.abs(i), Math.abs(j)) !== ring) continue;
      const p = clamp({ x: start.x + i * (NODE_W + 40), y: start.y + j * (NODE_H + 40) });
      if (!clash(p.x, p.y)) return p;
    }
  return start;
}

function reach(graph, from, forward, skip = null) {
  const found = new Set();
  const stack = [from];
  while (stack.length) {
    const v = stack.pop();
    for (const e of graph.edges) {
      const [a, b] = forward ? [e.source, e.target] : [e.target, e.source];
      if (a === v && b !== skip && !found.has(b)) { found.add(b); stack.push(b); }
    }
  }
  return found;
}
// Whether each variable can take a new arrow out of the variable named source, and why not: itself, already its child, an ancestor
// (a cycle), visit time (set by the schedule) or a coarser level than the source.
export function connectTargets(graph, source) {
  const from = graph.nodes.find((n) => n.id === source), above = reach(graph, source, false);
  return new Map(graph.nodes.map((n) => {
    let reason = '';
    if (n.id === source) reason = 'An arrow needs two different variables.';
    else if (graph.edges.some((e) => e.source === source && e.target === n.id)) reason = `${source} already points to ${n.id}.`;
    else if (above.has(n.id)) reason = `That arrow creates a cycle: ${n.id} already affects ${source}.`;
    else if (graph.study && scopeOf(n) === 'time') reason = 'Visit time comes from the schedule and has no parents.';
    else if (from && levelOf(graph, from) > levelOf(graph, n)) reason = `${source} varies within ${LEVELS[levelOf(graph, n)]}s, so it cannot affect ${n.id}, which is measured once per ${LEVELS[levelOf(graph, n)]}.`;
    return [n.id, { ok: !reason, reason }];
  }));
}

// What each other variable is for the comparison of `exposure` on `outcome`, from the arrows alone: a mediator
// (on a path from one to the other), a confounder (affects both, not only through the exposure), a collider
// (affected by both), an instrument (affects the outcome only through the exposure) or prognostic (affects the
// outcome only).  Variables with none of these roles are left out.
export function inferRoles(graph, exposure, outcome) {
  const roles = new Map(), ids = new Set(graph.nodes.map((n) => n.id));
  if (!ids.has(exposure) || !ids.has(outcome) || exposure === outcome) return roles;
  const belowA = reach(graph, exposure, true), aboveA = reach(graph, exposure, false), belowY = reach(graph, outcome, true);
  const aboveY = reach(graph, outcome, false), aboveYWithoutA = reach(graph, outcome, false, exposure);
  for (const id of ids) {
    if (id === exposure || id === outcome) continue;
    if (belowA.has(id) && aboveY.has(id)) roles.set(id, 'mediator');
    else if (aboveA.has(id) && aboveYWithoutA.has(id)) roles.set(id, 'confounder');
    else if (belowA.has(id) && belowY.has(id)) roles.set(id, 'collider');
    else if (aboveA.has(id) && aboveY.has(id)) roles.set(id, 'instrument');
    else if (aboveYWithoutA.has(id) && !belowA.has(id)) roles.set(id, 'prognostic');
  }
  return roles;
}

// The region the view must show: every box, plate and random-effect circle, with a margin.
export function contentBounds(graph, plates = [], effects = []) {
  const box = bounds([...graph.nodes.map(nodeBox), ...plates, ...effects.map((e) => ({ x: e.cx - 16, y: e.cy - 16, w: 32, h: 32 }))]);
  return box ? { x: box.x - 40, y: box.y - 40, w: box.w + 80, h: box.h + 80 } : { x: 0, y: 0, w: 1100, h: 1100 * ASPECT };
}
// The view (in drawing units) that shows `box` whole at the drawing's aspect, no closer than `minWidth` across.
export function fitView(box, minWidth = 1100) {
  const w = Math.max(box.w, box.h / ASPECT, minWidth), h = w * ASPECT;
  return { x: box.x + box.w / 2 - w / 2, y: box.y + box.h / 2 - h / 2, w };
}
export const ZOOM = { min: 360, max: 9000 };
// The view zoomed by `factor` (> 1 closer) keeping the drawing point `at` where it is on screen.
export function zoomView(view, factor, at = { x: view.x + view.w / 2, y: view.y + view.w * ASPECT / 2 }) {
  const w = Math.min(ZOOM.max, Math.max(ZOOM.min, view.w / factor)), k = w / view.w;
  return { x: at.x - (at.x - view.x) * k, y: at.y - (at.y - view.y) * k, w };
}
export const viewBox = (view) => `${view.x.toFixed(1)} ${view.y.toFixed(1)} ${view.w.toFixed(1)} ${(view.w * ASPECT).toFixed(1)}`;
