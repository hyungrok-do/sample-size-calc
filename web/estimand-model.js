import { DISTRIBUTIONS, parseGraph } from './dag-model.js?v=dev';
import { parsePowerSettings } from './power-settings.js?v=dev';

export const MEASURES = {
  mean_difference: { label: 'Difference in marginal means', kinds: ['continuous', 'count', 'binary'], nullValue: 0 },
  risk_difference: { label: 'Risk difference', kinds: ['binary'], nullValue: 0, bounded: true },
  risk_ratio: { label: 'Risk ratio', kinds: ['binary'], nullValue: 1, positive: true },
  odds_ratio: { label: 'Marginal odds ratio', kinds: ['binary'], nullValue: 1, positive: true },
  survival_difference: { label: 'Survival probability difference at a time', kinds: ['survival'], nullValue: 0, bounded: true, horizon: true },
  rmst_difference: { label: 'Restricted mean survival time difference', kinds: ['survival'], nullValue: 0, horizon: true },
  coefficient: { label: 'Conditional predictor coefficient', nullValue: 0 },
  exponentiated_coefficient: { label: 'Exponentiated conditional coefficient', nullValue: 1, positive: true },
  custom: { label: 'Custom estimand', nullValue: 0 },
};
const finite = (x) => typeof x === 'number' && Number.isFinite(x);
const fail = (text) => { throw new Error(text); };
const logLinks = ['bernoulli', 'poisson', 'gamma', 'lognormal', 'weibull', 'exponential'];

export function availableMeasures(graph, kind, outcome) {
  const node = graph.nodes.find((n) => n.id === outcome);
  if (kind === 'custom') return ['custom'];
  if (!node) return [];
  if (kind === 'coefficient') return ['coefficient', ...(logLinks.includes(node.family) ? ['exponentiated_coefficient'] : [])];
  return Object.entries(MEASURES).filter(([, spec]) => spec.kinds?.includes(node.kind)).map(([key]) => key);
}
export function validateEstimand(graph, value) {
  if (!value || value.version !== 1 || !['marginal', 'coefficient', 'custom'].includes(value.kind)) fail('Choose an estimand type.');
  if (typeof value.population !== 'string' || !value.population.trim() || value.population.length > 1000)
    fail('Describe the target population (up to 1,000 characters).');
  if (typeof value.definition !== 'string' || value.definition.length > 3000) fail('The estimand definition must be text of up to 3,000 characters.');
  if (typeof value.outcome !== 'string' || typeof value.exposure !== 'string') fail('Estimand variable identifiers must be text.');
  if (!finite(value.reference) || !finite(value.active)) fail('Comparison levels must be finite numbers.');
  if (value.horizon !== null && !finite(value.horizon)) fail('The time horizon must be a finite number or left unset.');
  const outcome = graph.nodes.find((n) => n.id === value.outcome);
  const exposure = graph.nodes.find((n) => n.id === value.exposure);
  if (value.kind === 'custom') {
    if (!value.definition.trim()) fail('Write the custom estimand definition.');
    if (value.measure !== 'custom') fail('Use the custom measure for a custom estimand.');
  } else {
    if (!outcome || !exposure) fail('Select an outcome variable and a comparison variable.');
    if (outcome.id === exposure.id) fail('The outcome and comparison variable must differ.');
    if (!availableMeasures(graph, value.kind, value.outcome).includes(value.measure)) fail('Choose an effect measure compatible with the outcome distribution.');
    if (value.kind === 'coefficient') {
      if (!graph.edges.some((e) => e.source === exposure.id && e.target === outcome.id))
        fail('A conditional coefficient requires an arrow from the selected predictor to the outcome.');
    } else {
      if (!finite(value.reference) || !finite(value.active) || value.reference === value.active) fail('Enter two different, finite comparison levels.');
      for (const level of [value.reference, value.active]) {
        if (exposure.kind === 'binary' && ![0, 1].includes(level)) fail('Binary comparison levels must be 0 or 1.');
        if (exposure.kind === 'count' && (!Number.isInteger(level) || level < 0)) fail('Count comparison levels must be nonnegative integers.');
        if (exposure.kind === 'survival' && level <= 0) fail('Survival-time intervention levels must be positive.');
      }
    }
  }
  const spec = MEASURES[value.measure];
  if (!spec) fail('Choose an effect measure.');
  if (spec.horizon && (!finite(value.horizon) || value.horizon <= 0)) fail('Enter a positive time horizon in the event-time variable’s units.');
  for (const [label, number] of [['Null value', value.nullValue], ['Planning alternative', value.alternative]]) {
    if (number === null && label === 'Planning alternative') continue;
    if (!finite(number)) fail(`${label} must be a finite number.`);
    if (spec.positive && number <= 0) fail(`${label} must be positive on a ratio scale.`);
    if ((spec.bounded || value.kind === 'marginal' && outcome?.kind === 'binary' && value.measure === 'mean_difference') && (number < -1 || number > 1)) fail(`${label} must lie between −1 and 1.`);
    if (value.measure === 'rmst_difference' && Math.abs(number) > value.horizon) fail(`${label} cannot exceed the RMST horizon in magnitude.`);
  }
  if (value.alternative !== null && value.alternative === value.nullValue) fail('The planning alternative must differ from the null value.');
  return value;
}

export function copyEstimand(value) {
  return Object.fromEntries(['version', 'kind', 'population', 'definition', 'outcome', 'exposure', 'measure',
    'reference', 'active', 'horizon', 'nullValue', 'alternative'].map((key) => [key, value[key]]));
}
export function renameEstimand(value, from, to) {
  return value ? { ...value, outcome: value.outcome === from ? to : value.outcome,
    exposure: value.exposure === from ? to : value.exposure } : null;
}

export function describeEstimand(graph, value) {
  validateEstimand(graph, value);
  const y = value.outcome, a = value.exposure, lo = value.reference, hi = value.active, tau = value.horizon;
  let formula, note;
  if (value.kind === 'custom') {
    formula = value.definition;
    note = 'This definition is recorded for planning. Its estimator and simulation evaluation must be specified explicitly.';
  } else if (value.kind === 'coefficient') {
    const node = graph.nodes.find((n) => n.id === y);
    const coefficient = graph.edges.find((e) => e.source === a && e.target === y).coefficient;
    const ratio = value.measure === 'exponentiated_coefficient';
    formula = ratio ? `θ = exp(β_${a}→${y}); current DGP value = ${Math.exp(coefficient)}`
      : `θ = β_${a}→${y} on the ${DISTRIBUTIONS[node.family].link} scale; current DGP value = ${coefficient}`;
    note = 'This is a conditional structural coefficient, rather than a marginal total effect. With interactions, it applies when the interacting parents equal zero. The analysis must target this coefficient.';
  } else {
    const intervention = (term, level) => `${term} under do(${a} = ${level})`;
    if (value.measure === 'mean_difference') formula = `θ = E[${intervention(y, hi)}] − E[${intervention(y, lo)}]`;
    if (value.measure === 'risk_difference') formula = `θ = P(${intervention(`${y} = 1`, hi)}) − P(${intervention(`${y} = 1`, lo)})`;
    if (value.measure === 'risk_ratio') formula = `θ = P(${intervention(`${y} = 1`, hi)}) / P(${intervention(`${y} = 1`, lo)})`;
    if (value.measure === 'odds_ratio') formula = `θ = odds(${y} = 1 under do(${a} = ${hi})) / odds(${y} = 1 under do(${a} = ${lo}))`;
    if (value.measure === 'survival_difference') formula = `θ = P(${intervention(`${y} > ${tau}`, hi)}) − P(${intervention(`${y} > ${tau}`, lo)})`;
    if (value.measure === 'rmst_difference') formula = `θ = E[min(${y}, ${tau}) under do(${a} = ${hi})] − E[min(${y}, ${tau}) under do(${a} = ${lo})]`;
    note = 'This is a marginal total-effect contrast under interventions in the specified DGP, averaged over its population distribution. Mediators follow their generating models in each arm; they are not held fixed or automatically adjusted for.';
    if (MEASURES[value.measure].horizon) note += ' The target concerns latent event times; censoring affects the analysis, not this estimand.';
  }
  return { formula, note, population: value.population, hypothesis: `H₀: θ = ${value.nullValue}`,
    alternative: value.alternative === null ? 'Planning alternative: use the effect implied by the graph as drawn.'
      : `Planning alternative: θ = ${value.alternative}. ${value.kind === 'custom' ? 'An executable estimator and test are required for power.' : 'Power calibrates the run DGP to this target and verifies its implied effect.'}` };
}
export function exportDesign(graph, estimand, power = null) {
  if (estimand) validateEstimand(graph, estimand);
  return { version: 1, graph: parseGraph(graph), estimand: estimand ? copyEstimand(estimand) : null,
    ...(power === null ? {} : { power: parsePowerSettings(power) }) };
}
export function parseDesign(input) {
  const raw = typeof input === 'string' ? JSON.parse(input) : input;
  if (raw?.nodes) return { graph: parseGraph(raw), estimand: null }; // standalone DAG files remain valid
  if (!raw || raw.version !== 1 || !raw.graph) fail('Unsupported simulation design file.');
  const graph = parseGraph(raw.graph), estimand = raw.estimand ? copyEstimand(raw.estimand) : null;
  if (estimand) validateEstimand(graph, estimand);
  return { graph, estimand, ...(raw.power == null ? {} : { power: parsePowerSettings(raw.power) }) };
}
