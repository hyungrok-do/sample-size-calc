// UI compatibility hints only. The C++ engine independently validates and selects the actual analysis.
import { validateGraph } from './dag-model.js?v=dev';
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
const termText = (terms) => terms.map((t) => t.right ? `${t.left} × ${t.right}` : t.left).join(' + ');
// Mirrors the engine's AIPW plan for a binary comparison variable with parents: the exposure's own terms as
// the propensity model; within each arm, an outcome model on the exposure's parents plus the outcome's
// parents that the exposure does not affect.
function aipwPlan(graph, estimand, exposure, outcome) {
  if (outcome.kind === 'survival')
    fail('A survival outcome with a confounded comparison variable needs a confounding- and censoring-adjusted estimator, which is not available yet. Remove the arrows into the comparison variable or use a conditional hazard-ratio target.');
  const node = (id) => graph.nodes.find((n) => n.id === id);
  const censored = (id) => node(id).kind === 'survival' && node(id).censoring.mode !== 'none';
  const intoA = graph.edges.filter((e) => e.target === exposure.id);
  const propensity = [...intoA.map((e) => ({ left: e.source })), ...exposure.interactions.map((t) => ({ left: t.left, right: t.right }))];
  if (intoA.some((e) => censored(e.source)))
    fail('A parent of the comparison variable is a censored event time. Its latent values cannot be used to adjust for confounding.');
  const affected = descendants(graph, exposure.id);
  const chosen = new Set(intoA.map((e) => e.source));
  for (const e of graph.edges) if (e.target === outcome.id && !affected.has(e.source) && !censored(e.source)) chosen.add(e.source);
  const terms = graph.nodes.filter((n) => chosen.has(n.id)).map((n) => ({ left: n.id }));
  for (const t of [...exposure.interactions, ...outcome.interactions])
    if (chosen.has(t.left) && chosen.has(t.right) && !terms.some((u) => u.right && [u.left, u.right].sort().join() === [t.left, t.right].sort().join()))
      terms.push({ left: t.left, right: t.right });
  const outcomeModel = outcome.family === 'bernoulli' ? 'bernoulli' : outcome.family === 'poisson' ? 'poisson'
    : ['gamma', 'lognormal'].includes(outcome.family) ? 'gamma' : 'gaussian';
  const scale = estimand.measure === 'risk_ratio' ? 'log risk ratio' : estimand.measure === 'odds_ratio' ? 'log odds ratio'
    : estimand.measure === 'risk_difference' || outcome.kind === 'binary' ? 'risk difference' : 'mean difference';
  const model = { gaussian: 'linear regression', bernoulli: 'logistic regression', poisson: 'Poisson log-link regression', gamma: 'Gamma log-link regression' }[outcomeModel];
  return { method: 'aipw', label: `Doubly robust AIPW Wald test of the marginal ${scale}, adjusting for ${termText(terms)}`, terms,
    propensityTerms: propensity, outcomeModel,
    formula: `${outcome.id}: ${exposure.id} = ${estimand.active} versus ${estimand.reference}; propensity ${exposure.id} ~ ${termText(propensity)} (logistic); outcome ${outcome.id} ~ ${termText(terms)} within each arm (${model})` };
}
export function compilePowerPlan(graph, estimand) {
  validateGraph(graph); validateEstimand(graph, estimand);
  if (estimand.kind === 'custom') fail('A custom text definition needs an executable estimator and test. Choose a supported target to run power.');
  const outcome = graph.nodes.find((n) => n.id === estimand.outcome);
  const exposure = graph.nodes.find((n) => n.id === estimand.exposure);
  const edges = calibrationEdges(graph, estimand);
  let method, label, terms = [], targetIndex = -1;
  if (estimand.kind === 'marginal') {
    if (exposure.kind !== 'binary')
      fail('Marginal arm-comparison power currently requires a binary comparison variable. A continuous or count intervention needs a different estimator.');
    if (graph.edges.some((e) => e.target === exposure.id))
      return { ...aipwPlan(graph, estimand, exposure, outcome), targetIndex: -1, outcome: outcome.id, exposure: exposure.id,
        family: outcome.family, calibrationEdges: edges, testScale: MEASURES[estimand.measure].positive ? 'log' : 'original' };
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
    if (outcome.kind === 'survival') {
      const rule = outcome.censoring;
      if (rule.mode === 'administrative' && estimand.horizon > rule.time)
        fail('The estimand horizon exceeds the administrative follow-up cutoff. It cannot be estimated by this KM analysis.');
      if (rule.mode === 'variable') {
        const t = ancestors(graph, outcome.id), c = ancestors(graph, rule.variable);
        const shared = [...t].filter((id) => id !== exposure.id && c.has(id));
        if (shared.length) fail(`The event and censoring mechanisms share ${shared.join(', ')} beyond the randomized arm. Unadjusted KM does not identify this marginal target; change the censoring DGP or use a censoring-adjusted estimator.`);
      }
    }
  } else {
    const parents = graph.edges.filter((e) => e.target === outcome.id);
    terms = [...parents.map((e) => ({ left: e.source })), ...outcome.interactions.map((t) => ({ left: t.left, right: t.right }))];
    for (const parent of parents) {
      const node = graph.nodes.find((n) => n.id === parent.source);
      if (node.kind === 'survival' && node.censoring.mode !== 'none')
        fail('The planned regression includes a censored event-time predictor. It cannot use its latent values as observed covariates.');
    }
    targetIndex = parents.findIndex((e) => e.source === exposure.id);
    if (['gaussian', 'uniform', 'lognormal'].includes(outcome.family)) {
      method = outcome.family === 'lognormal' ? 'ols_log' : 'ols';
      label = `${outcome.family === 'lognormal' ? 'Log-response linear' : 'Linear'} regression coefficient t test${outcome.family === 'uniform' ? ' (approximate for Uniform errors)' : ''}`;
    } else if (['bernoulli', 'poisson', 'gamma'].includes(outcome.family)) {
      method = 'glm'; label = `${outcome.family === 'bernoulli' ? 'Logistic' : outcome.family === 'poisson' ? 'Poisson log-link' : 'Gamma log-link'} regression coefficient Wald test`;
    } else if (outcome.kind === 'survival') {
      method = 'cox'; label = 'Cox PH coefficient Wald test (Breslow ties)';
      if (outcome.censoring.mode === 'variable' && ancestors(graph, outcome.censoring.variable).has(outcome.id))
        fail('The censoring variable is a descendant of the event time. The proposed Cox analysis requires independent censoring given its fitted parents.');
    }
  }
  if (!method) fail('No executable analysis is available for this estimand.');
  return { method, label, terms, targetIndex, outcome: outcome.id, exposure: exposure.id, family: outcome.family,
    calibrationEdges: edges, testScale: MEASURES[estimand.measure].positive ? 'log' : 'original',
    formula: estimand.kind === 'coefficient' ? `${method === 'ols_log' ? `log(${outcome.id})` : outcome.id} ~ ${termText(terms)}`
      : `${outcome.id}: ${exposure.id} = ${estimand.active} versus ${estimand.reference}` };
}
export function powerCompatibility(graph, estimand) {
  if (!estimand) return { supported: false, reason: 'Save an estimand above to choose its analysis and calculate power.' };
  try { return { supported: true, ...compilePowerPlan(graph, estimand) }; }
  catch (error) { return { supported: false, reason: error.message }; }
}
