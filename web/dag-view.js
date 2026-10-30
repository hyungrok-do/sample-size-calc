// What the DAG canvas draws besides boxes and arrows: a variable's simulated distribution in miniature, and the paths
// between the compared variable and the outcome, with what blocks each one.  Plain functions of data, no DOM.

const fmt = (x) => String(Number(x.toPrecision(3))).replace('-', '−');
function quantiles(values) {
  const sorted = Float64Array.from(values).sort();
  return (p) => sorted[Math.min(sorted.length - 1, Math.max(0, Math.round(p * (sorted.length - 1))))];
}
const bar = (x, width, height, h) => height > 0 ? `M${x.toFixed(1)} ${h}V${(h - height).toFixed(1)}H${(x + width).toFixed(1)}V${h}Z` : '';

// A sample drawn small in a w × h box: bars (histogram, two bars for 0/1, one per count) or, for event times, the
// survival curve.  Returns path data, an optional short label and a tooltip.
export function sketch(kind, values, w, h) {
  const n = values.length;
  if (!n) return null;
  const q = quantiles(values);
  let mean = 0;
  for (const v of values) mean += v;
  mean /= n;
  if (kind === 'binary') {
    const p = mean, width = 14;
    return { bars: bar(0, width, (1 - p) * h, h) + bar(width + 4, width, p * h, h), line: '',
      label: `P(1) ≈ ${fmt(p)}`, title: `Simulated share of 1s: ${fmt(p)} (n = ${n})` };
  }
  if (kind === 'survival') {
    const end = q(0.95) || 1, steps = 48;
    const sorted = Float64Array.from(values).sort();
    let k = 0, d = '';
    for (let i = 0; i <= steps; i++) {
      const t = end * i / steps;
      while (k < n && sorted[k] <= t) k++;
      d += `${i ? 'L' : 'M'}${(w * i / steps).toFixed(1)} ${(h - h * (n - k) / n).toFixed(1)}`;
    }
    return { bars: '', line: d, label: `median ${fmt(q(0.5))}`,
      title: `Simulated survival curve to the 95th percentile (${fmt(end)}); median ${fmt(q(0.5))} (n = ${n})` };
  }
  let sd = 0;
  for (const v of values) sd += (v - mean) ** 2;
  sd = Math.sqrt(sd / Math.max(1, n - 1));
  const title = `Simulated: mean ${fmt(mean)}, SD ${fmt(sd)}, 1st–99th percentile ${fmt(q(0.01))} to ${fmt(q(0.99))} (n = ${n})`;
  let edges;  // bin edges
  if (kind === 'count') {
    const top = Math.min(40, Math.max(1, q(0.99)));
    edges = Array.from({ length: top + 2 }, (_, i) => i - 0.5);
  } else {
    let lo = q(0.01), hi = q(0.99);
    if (!(hi > lo)) { lo -= 0.5; hi += 0.5; }
    edges = Array.from({ length: 25 }, (_, i) => lo + (hi - lo) * i / 24);
  }
  const bins = edges.length - 1, counts = new Array(bins).fill(0), span = edges[bins] - edges[0];
  for (const v of values) {
    const i = Math.floor((v - edges[0]) / span * bins);
    if (i >= 0 && i < bins) counts[i]++;
  }
  const top = Math.max(...counts), width = w / bins;
  return { bars: counts.map((c, i) => bar(i * width, Math.max(0.6, width - 1), c / top * h, h)).join(''), line: '', label: '', title };
}

// Every simple path between two variables in the graph with its arrows' directions ignored (at most `limit`), as
// steps { from, to, forward } where forward means the arrow points from → to.
export function findPaths(graph, from, to, limit = 40) {
  const links = new Map(graph.nodes.map((n) => [n.id, []]));
  for (const e of graph.edges) {
    links.get(e.source).push({ next: e.target, forward: true });
    links.get(e.target).push({ next: e.source, forward: false });
  }
  const found = [], seen = new Set([from]), steps = [];
  let truncated = false;
  (function walk(at) {
    if (truncated) return;
    if (at === to) {
      if (found.length === limit) truncated = true;
      else found.push([...steps]);
      return;
    }
    for (const { next, forward } of links.get(at) ?? []) {
      if (seen.has(next)) continue;
      seen.add(next); steps.push({ from: at, to: next, forward });
      walk(next);
      steps.pop(); seen.delete(next);
    }
  })(from);
  return { paths: found, truncated };
}

// What a path is and whether it is open when the analysis adjusts for `adjusted` (a Set of variable ids): causal
// (every arrow points toward the outcome), back-door (it starts with an arrow into the compared variable) or another
// non-causal path; blocked by an adjusted non-collider, or by a collider that neither it nor any descendant of it
// is adjusted for.
export function classifyPath(graph, path, adjusted = new Set()) {
  const kind = path.every((s) => s.forward) ? 'causal' : path[0].forward ? 'other' : 'backdoor';
  const below = (id) => {
    const found = new Set([id]);
    for (const v of found) for (const e of graph.edges) if (e.source === v) found.add(e.target);
    return found;
  };
  for (let i = 1; i < path.length; i++) {
    const v = path[i].from, collider = path[i - 1].forward && !path[i].forward;
    if (collider && ![...below(v)].some((d) => adjusted.has(d))) return { kind, open: false, by: v, collider: true };
    if (!collider && adjusted.has(v)) return { kind, open: false, by: v, collider: false };
  }
  return { kind, open: true, by: null, collider: false };
}
export const pathText = (path) => path.reduce((text, s) => `${text} ${s.forward ? '→' : '←'} ${s.to}`, path[0]?.from ?? '');
