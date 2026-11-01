import { nodeTerms, termKey, termText as oneTermText, parseAnalysis, analysisTerms, nonlinearIn } from './dag-functions.js?v=dev';
// UI compatibility hints only. The C++ engine independently validates and selects the actual analysis.
import { validateGraph, topologicalOrder } from './dag-model.js?v=dev';
import { validateEstimand, MEASURES } from './estimand-model.js?v=dev';

import { DEFAULT_POWER, parsePowerSettings } from './power-settings.js?v=dev';
export { DEFAULT_POWER, parsePowerSettings };
const fail = (message) => { throw new Error(message); };
function ancestors(graph, id) {
  const found = new Set([id]), queue = [id];
  for (let i = 0; i < queue.length; i++) for (const edge of graph.edges)
    if (edge.target === queue[i] && !found.has(edge.source)) { found.add(edge.source); queue.push(edge.source); }
  return found;
}
function descendants(graph, id) {
  const found = new Set([id]), queue = [id];
  for (let i = 0; i < queue.length; i++) for (const edge of graph.edges)
    if (edge.source === queue[i] && !found.has(edge.target)) { found.add(edge.target); queue.push(edge.target); }
  return found;
}
export function calibrationEdges(graph, estimand) {
  if (!estimand || estimand.kind === 'custom') return [];
  if (estimand.kind === 'coefficient') return graph.edges.filter((e) => e.source === estimand.exposure && e.target === estimand.outcome);
  const before = ancestors(graph, estimand.outcome), after = descendants(graph, estimand.exposure);
  return graph.edges.filter((e) => after.has(e.source) && before.has(e.target)).sort((a, b) =>
    Number(b.source === estimand.exposure && b.target === estimand.outcome) - Number(a.source === estimand.exposure && a.target === estimand.outcome));
}
const termText = (terms) => terms.map(oneTermText).join(' + ');
const censored = (node) => node.kind === 'survival' && node.censoring.mode !== 'none';
const latentNote = (id) => `${id} is a censored event time; the analysis uses its latent values as if they were measured.`;
const informativeNote = (node) => `${node.censoring.variable} censors ${node.id} but depends on ${node.id} itself, so no analysis of the observed data is unbiased; the results show the bias.`;
const addNote = (notes, text) => { if (!notes.includes(text)) notes.push(text); };
const reverseNote = (a, y) => `${y} lies upstream of ${a}, so ${a} cannot change ${y}: the effect is zero by construction and no study can detect one. Swap the two roles to study the effect of ${y} on ${a}.`;
const noPathNote = (a, y) => `No arrow path leads from ${a} to ${y}, so ${a} has no effect on ${y} here; the rejection rate is the test's false-positive rate.`;
const measureScale = (estimand, outcome) => {
  const m = estimand.measure;
  return m === 'risk_ratio' ? 'log risk ratio' : m === 'odds_ratio' ? 'log odds ratio'
    : m === 'survival_difference' ? `survival difference at time ${estimand.horizon}`
    : m === 'rmst_difference' ? `RMST difference through time ${estimand.horizon}`
    : m === 'risk_difference' || outcome.kind === 'binary' ? 'risk difference' : 'mean difference';
};
// Mirrors the engine's AIPW plan for a binary comparison variable with parents: the exposure's own terms as
// the propensity model; within each arm, an outcome model on the exposure's parents plus the outcome's
// parents that the exposure does not affect.
function aipwPlan(graph, estimand, exposure, outcome, notes) {
  const node = (id) => graph.nodes.find((n) => n.id === id);
  const propensity = nodeTerms(graph, exposure);
  const affected = descendants(graph, exposure.id);
  const chosen = new Set();
  for (const t of propensity) for (const id of [t.left, t.right]) if (id) {
    if (censored(node(id))) addNote(notes, latentNote(id));
    chosen.add(id);
  }
  for (const t of nodeTerms(graph, outcome)) if (!t.right && !affected.has(t.left) && !censored(node(t.left))) chosen.add(t.left);
  const terms = graph.nodes.filter((n) => chosen.has(n.id)).map((n) => ({ left: n.id }));
  for (const raw of [...nodeTerms(graph, exposure), ...nodeTerms(graph, outcome)]) {
    let t = { ...raw };
    if (t.left === exposure.id && t.right) t = { left: t.right, ...(t.rightTransform ? { transform: t.rightTransform } : {}) };
    else if (t.right === exposure.id) { delete t.right; delete t.rightTransform; }
    if (chosen.has(t.left) && (!t.right || chosen.has(t.right)) && !terms.some((u) => termKey(u) === termKey(t))) terms.push(t);
  }
  const outcomeModel = outcome.family === 'bernoulli' ? 'bernoulli' : outcome.family === 'poisson' ? 'poisson'
    : ['gamma', 'lognormal'].includes(outcome.family) ? 'gamma' : 'gaussian';
  const scale = estimand.measure === 'risk_ratio' ? 'log risk ratio' : estimand.measure === 'odds_ratio' ? 'log odds ratio'
    : estimand.measure === 'risk_difference' || outcome.kind === 'binary' ? 'risk difference' : 'mean difference';
  const model = { gaussian: 'linear regression', bernoulli: 'logistic regression', poisson: 'Poisson log-link regression', gamma: 'Gamma log-link regression' }[outcomeModel];
  return { method: 'aipw', label: `Doubly robust AIPW Wald test of the marginal ${scale}, adjusting for ${termText(terms)}`, terms,
    propensityTerms: propensity, outcomeModel, notes,
    formula: `${outcome.id}: ${exposure.id} = ${estimand.active} versus ${estimand.reference}; propensity ${exposure.id} ~ ${termText(propensity)} (logistic); outcome ${outcome.id} ~ ${termText(terms)} within each arm (${model})` };
}
// Mirrors the engine's parametric g-formula: refit the outcome and every variable the comparison affects on
// the way to it, each on its own parents, and standardize over the study's rows.
function gformulaPlan(graph, estimand, exposure, outcome, notes) {
  const label = `Parametric g-formula Wald test of the marginal ${measureScale(estimand, outcome)}`;
  const levels = `${outcome.id}: ${exposure.id} = ${estimand.active} versus ${estimand.reference}`;
  // the outcome causes the comparison variable: nothing to refit, and the contrast is zero in every study
  if (descendants(graph, outcome.id).has(exposure.id)) return { method: 'gformula', label, terms: [], models: [], notes, formula: `${levels}; zero by construction` };
  const after = descendants(graph, exposure.id), before = ancestors(graph, outcome.id);
  const node = (id) => graph.nodes.find((n) => n.id === id);
  const models = topologicalOrder(graph).filter((id) => id !== exposure.id && after.has(id) && before.has(id))
    .map((id) => ({ node: id, terms: nodeTerms(graph, node(id)) }));
  if (!after.has(outcome.id)) models.push({ node: outcome.id, terms: [{ left: exposure.id }, ...nodeTerms(graph, outcome)] });
  for (const m of models) for (const t of m.terms) for (const id of [t.left, t.right])
    if (id && id !== exposure.id && censored(node(id))) addNote(notes, latentNote(id));
  for (const m of models) {
    const v = node(m.node);
    if (v.kind === 'survival' && v.censoring.mode === 'variable' && descendants(graph, v.id).has(v.censoring.variable)) addNote(notes, informativeNote(v));
  }
  if (outcome.kind === 'survival' && outcome.censoring.mode === 'administrative' && estimand.horizon > outcome.censoring.time)
    addNote(notes, `The horizon ${estimand.horizon} lies past the follow-up cutoff ${outcome.censoring.time}; the Weibull model extrapolates beyond the data.`);
  return { method: 'gformula', label, terms: models.at(-1).terms, models, notes,
    formula: `${levels}; refitted ${models.map((x) => `${x.node} ~ ${termText(x.terms)}`).join(', ')}; other variables as observed` };
}
function defaultPlan(graph, estimand) {
  validateGraph(graph); validateEstimand(graph, estimand);
  if (estimand.kind === 'custom') fail('A custom text definition needs an executable estimator and test. Choose a supported target to run power.');
  const outcome = graph.nodes.find((n) => n.id === estimand.outcome);
  const exposure = graph.nodes.find((n) => n.id === estimand.exposure);
  const edges = calibrationEdges(graph, estimand);
  const common = { targetIndex: -1, outcome: outcome.id, exposure: exposure.id, family: outcome.family, calibrationEdges: edges,
    testScale: MEASURES[estimand.measure].positive ? 'log' : 'original' };
  let method, label, terms = [], targetIndex = -1;
  const notes = [];
  if (estimand.kind === 'marginal') {
    // A two-level comparison of a confounded binary variable uses AIPW; any other comparison variable, or a
    // confounded one with an event-time outcome, the g-formula.
    const upstream = descendants(graph, outcome.id).has(exposure.id);
    if (upstream) notes.push(reverseNote(exposure.id, outcome.id));
    else if (!descendants(graph, exposure.id).has(outcome.id)) notes.push(noPathNote(exposure.id, outcome.id));
    const confounded = graph.edges.some((e) => e.target === exposure.id) || exposure.interactions.length > 0;
    if (upstream || exposure.kind !== 'binary' || (confounded && outcome.kind === 'survival')) return { ...common, ...gformulaPlan(graph, estimand, exposure, outcome, notes) };
    if (confounded) return { ...common, ...aipwPlan(graph, estimand, exposure, outcome, notes) };
    const p = 1 / (1 + Math.exp(-exposure.intercept));
    if (!(p > 0 && p < 1)) fail('The comparison variable gives an empty arm. Use a nondegenerate allocation probability.');
    if (estimand.measure === 'mean_difference' && outcome.kind !== 'binary') { method = 'welch'; label = 'Welch t test of the marginal mean difference'; }
    else if (['risk_difference', 'mean_difference'].includes(estimand.measure)) {
      method = estimand.nullValue === 0 ? 'risk_score' : 'risk_wald';
      label = estimand.nullValue === 0 ? 'Two-arm proportion score test of the risk difference' : 'Wald test of the risk difference';
    } else if (estimand.measure === 'risk_ratio') { method = 'risk_ratio'; label = 'Wald test of the marginal log risk ratio'; }
    else if (estimand.measure === 'odds_ratio') { method = 'odds_ratio'; label = 'Wald test of the marginal log odds ratio'; }
    else if (estimand.measure === 'survival_difference') { method = 'km_survival'; label = `Kaplan–Meier / Greenwood Wald test at time ${estimand.horizon}`; }
    else if (estimand.measure === 'rmst_difference') { method = 'km_rmst'; label = `Kaplan–Meier RMST Wald test through time ${estimand.horizon}`; }
    // Kaplan–Meier needs follow-up to the horizon and censoring independent of the event time within arms;
    // otherwise the g-formula.
    if (outcome.kind === 'survival') {
      const rule = outcome.censoring;
      let shared = false;
      if (rule.mode === 'variable') {
        const t = ancestors(graph, outcome.id), c = ancestors(graph, rule.variable);
        shared = [...t].some((id) => id !== exposure.id && c.has(id));
      }
      if (shared || (rule.mode === 'administrative' && estimand.horizon > rule.time)) return { ...common, ...gformulaPlan(graph, estimand, exposure, outcome, notes) };
    }
  } else {
    terms = nodeTerms(graph, outcome);
    terms.forEach((t, j) => {
      if (t.right || t.transform) return;
      if (t.left === exposure.id) targetIndex = j;
      const node = graph.nodes.find((n) => n.id === t.left);
      if (censored(node)) addNote(notes, latentNote(node.id));
    });
    if (['gaussian', 'uniform', 'lognormal'].includes(outcome.family)) {
      method = outcome.family === 'lognormal' ? 'ols_log' : 'ols';
      label = `${outcome.family === 'lognormal' ? 'Log-response linear' : 'Linear'} regression coefficient t test${outcome.family === 'uniform' ? ' (approximate for Uniform errors)' : ''}`;
    } else if (['bernoulli', 'poisson', 'gamma'].includes(outcome.family)) {
      method = 'glm'; label = `${outcome.family === 'bernoulli' ? 'Logistic' : outcome.family === 'poisson' ? 'Poisson log-link' : 'Gamma log-link'} regression coefficient Wald test`;
    } else if (outcome.kind === 'survival') {
      method = 'cox'; label = 'Cox PH coefficient Wald test (Breslow ties)';
      if (outcome.censoring.mode === 'variable' && ancestors(graph, outcome.censoring.variable).has(outcome.id)) addNote(notes, informativeNote(outcome));
    }
  }
  if (!method) fail('No executable analysis is available for this estimand.');
  return { ...common, method, label, terms, targetIndex, notes,
    formula: estimand.kind === 'coefficient' ? `${method === 'ols_log' ? `log(${outcome.id})` : outcome.id} ~ ${termText(terms)}`
      : `${outcome.id}: ${exposure.id} = ${estimand.active} versus ${estimand.reference}` };
}
export function compilePowerPlan(graph, estimand, rawAnalysis = {}) {
  const spec = parseAnalysis(rawAnalysis);
  for (const id of Object.keys(spec.models ?? {})) if (!graph.nodes.some((n) => n.id === id)) fail('An analysis override refers to a missing variable.');
  let p = defaultPlan(graph, estimand);
  const a = graph.nodes.find((n) => n.id === estimand.exposure), y = graph.nodes.find((n) => n.id === estimand.outcome);
  if (estimand.kind !== 'marginal' && spec.method !== 'auto') fail('Conditional and joint targets use a regression test.');
  if (spec.method === 'gformula') p = { ...p, ...gformulaPlan(graph, estimand, a, y, p.notes) };
  if (['aipw', 'tmle'].includes(spec.method)) {
    if (a.kind !== 'binary' || y.kind === 'survival' || descendants(graph, y.id).has(a.id))
      fail('AIPW/TMLE require a binary exposure, a nonsurvival outcome, and an outcome downstream or unrelated to exposure.');
    if (spec.method === 'tmle' && y.kind !== 'binary') fail('This TMLE implementation supports binary outcomes.');
    p = { ...p, ...aipwPlan(graph, estimand, a, y, p.notes), method: spec.method };
  }
  p.analysis = spec; p.folds = spec.folds;
  if (p.folds > 1 && !['aipw', 'tmle'].includes(p.method)) fail('Cross-fitting is available for AIPW and binary-outcome TMLE. Select one of those methods.');
  p.terms = analysisTerms(graph, p.terms, spec, y.id);
  if (['aipw', 'tmle'].includes(p.method)) {
    p.propensityTerms = analysisTerms(graph, p.propensityTerms, spec, a.id, true);
    const affected = descendants(graph, a.id);
    if ([...p.terms, ...p.propensityTerms].some((t) => affected.has(t.left) || t.right && affected.has(t.right))) fail('AIPW/TMLE adjustment cannot include the exposure or its descendants.');
    p.label = `${p.method === 'tmle' ? 'Binary-outcome TMLE' : 'Doubly robust AIPW'} Wald test of the marginal ${measureScale(estimand, y)}, adjusting for ${termText(p.terms)}` + (p.folds > 1 ? ` with ${p.folds}-fold cross-fitting` : '');
    const model = { gaussian: 'linear regression', bernoulli: 'logistic regression', poisson: 'Poisson log-link regression', gamma: 'Gamma log-link regression' }[p.outcomeModel];
    p.formula = `${y.id}: ${a.id} = ${estimand.active} versus ${estimand.reference}; propensity ${a.id} ~ ${termText(p.propensityTerms)} (logistic); outcome ${y.id} ~ ${termText(p.terms)} within each arm (${model})`;
  }
  if (p.method === 'gformula') {
    p.models = p.models.map((c) => ({ node: c.node, terms: analysisTerms(graph, c.terms, spec, c.node) }));
    for (const c of p.models) {
      const later = descendants(graph, c.node);
      if (c.terms.some((t) => later.has(t.left) || t.right && later.has(t.right))) fail('A g-computation model cannot condition on itself or its descendants.');
    }
    if (p.models.length) {
      p.terms = p.models.at(-1).terms;
      p.formula = `${y.id}: ${a.id} = ${estimand.active} versus ${estimand.reference}; refitted ${p.models.map((c) => `${c.node} ~ ${termText(c.terms)}`).join(', ')}; other variables as observed`;
    }
  }
  if (['coefficient', 'joint'].includes(estimand.kind)) {
    p.targetIndex = p.terms.findIndex((t) => t.left === a.id && !t.right && !t.transform);
    if (estimand.kind === 'coefficient' && p.targetIndex < 0) fail('The planned model omits the selected linear coefficient.');
    if (estimand.kind === 'joint') {
      p.testedIndices = p.terms.flatMap((t, i) => (estimand.measure === 'nonlinearity' ? nonlinearIn(t, a.id) : t.left === a.id || t.right === a.id) ? [i + 1] : []);
      if (!p.testedIndices.length) fail('The planned model has no terms for the selected joint test.');
      p.method = 'joint'; p.label = `Joint Wald test of ${estimand.measure === 'nonlinearity' ? 'nonlinearity' : 'all predictor terms'} (${p.testedIndices.length} df)`;
      addNote(p.notes, y.kind === 'survival' ? 'The joint survival test fits a parametric Weibull proportional-hazards model.' : 'The joint Wald test uses a sandwich covariance and a large-sample chi-square reference distribution.');
    }
    p.formula = `${y.family === 'lognormal' ? `log(${y.id})` : y.id} ~ ${termText(p.terms)}`;
  }
  if (spec.basis !== 'match' || spec.propensity || spec.models) addNote(p.notes, 'Analysis terms are specified separately from the generating model. Inspect bias, interval coverage and null rejection when they differ.');
  if (estimand.kind === 'marginal' && !['gformula', 'aipw', 'tmle'].includes(p.method) && (spec.basis !== 'match' || spec.propensity || spec.models))
    addNote(p.notes, 'The automatic two-arm test does not fit covariate models. Select g-computation, AIPW or TMLE to use the specified regression terms.');
  return p;
}
export function powerCompatibility(graph, estimand, analysis = {}) {
  if (!estimand) return { supported: false, reason: 'Save an estimand above to choose its analysis and calculate power.' };
  try { return { supported: true, ...compilePowerPlan(graph, estimand, analysis) }; }
  catch (error) { return { supported: false, reason: error.message }; }
}

