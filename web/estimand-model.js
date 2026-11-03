import { contrastWeights } from './dag-study.js?v=dev';
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
  overall_effect: { label: 'Overall predictor effect (joint test)', nullValue: 0 },
  nonlinearity: { label: 'Nonlinearity beyond linear terms (joint test)', nullValue: 0 },
  interaction: { label: 'Effect modification: interaction terms (joint test)', nullValue: 0 },
  custom: { label: 'Custom estimand', nullValue: 0 },
};
const finite = (x) => typeof x === 'number' && Number.isFinite(x);
const fail = (text) => { throw new Error(text); };
const logLinks = ['bernoulli', 'beta', 'poisson', 'negbin', 'gamma', 'lognormal', 'weibull', 'exponential'];

export function availableMeasures(graph, kind, outcome) {
  const node = graph.nodes.find((n) => n.id === outcome);
  if (kind === 'custom') return ['custom'];
  if (!node) return [];
  if (graph.study && kind === 'marginal') return node.kind === 'binary' ? ['risk_difference'] : node.kind === 'survival' ? [] : ['mean_difference'];
  if (kind === 'joint') return ['overall_effect', 'nonlinearity', 'interaction'];
  if (kind === 'coefficient') return ['coefficient', ...(logLinks.includes(node.family) ? ['exponentiated_coefficient'] : [])];
  return Object.entries(MEASURES).filter(([, spec]) => spec.kinds?.includes(node.kind)).map(([key]) => key);
}
export function validateEstimand(graph, value) {
  if (!value || value.version !== 1 || !['marginal', 'coefficient', 'joint', 'custom'].includes(value.kind)) fail('Choose an estimand type.');
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
    if (outcome.measured === false || exposure.measured === false) fail('The compared variable and the outcome must be measured.');
    if (!availableMeasures(graph, value.kind, value.outcome).includes(value.measure)) fail('Choose an effect measure compatible with the outcome distribution.');
    if (value.kind === 'coefficient' || value.kind === 'joint') {
      if (value.kind === 'joint' && (value.nullValue !== 0 || value.alternative !== null)) fail('Joint tests use a zero coefficient-block null and the graph as drawn.');
      if (value.kind === 'joint' && value.measure === 'nonlinearity' && exposure.kind === 'binary') fail('A binary predictor has no nonlinear dose-response component.');
      if (!graph.edges.some((e) => e.source === exposure.id && e.target === outcome.id))
        fail('A conditional coefficient requires an arrow from the selected predictor to the outcome.');
      if (outcome.unknown?.inputs?.includes(exposure.id))
        fail(`The effect of ${exposure.id} on ${outcome.id} has an unknown form, so no coefficient describes it. Compare groups or doses instead.`);
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
    if ((spec.bounded || value.kind === 'marginal' && outcome?.kind === 'binary' && value.measure === 'mean_difference') && !graph.study && (number < -1 || number > 1)) fail(`${label} must lie between −1 and 1.`);
    if (value.measure === 'rmst_difference' && Math.abs(number) > value.horizon) fail(`${label} cannot exceed the RMST horizon in magnitude.`);
  }
  if (value.alternative !== null && value.alternative === value.nullValue) fail('The planning alternative must differ from the null value.');
  if (graph.study && value.kind === 'marginal') contrastWeights(graph, value);
  return value;
}

export function copyEstimand(value) {
  const out = Object.fromEntries(['version', 'kind', 'population', 'definition', 'outcome', 'exposure', 'measure',
    'reference', 'active', 'horizon', 'nullValue', 'alternative'].map((key) => [key, value[key]]));
  if (value.repeated) out.repeated = { contrast: value.repeated.contrast ?? 'final', ...(value.repeated.visit != null ? { visit: value.repeated.visit } : {}), ...(value.repeated.weights ? { weights: [...value.repeated.weights] } : {}) };
  return out;
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
  } else if (value.kind === 'joint') {
    formula = `H₀: ${value.measure === 'nonlinearity' ? 'all nonlinear terms in' : value.measure === 'interaction' ? 'all interaction terms of' : 'all terms involving'} ${a} in the model for ${y} equal zero`;
    note = 'This is a conditional regression block test. Its hypothesis depends on the specified analysis basis. It is not a marginal causal contrast.';
  } else if (value.kind === 'coefficient') {
    const node = graph.nodes.find((n) => n.id === y);
    const coefficient = graph.edges.find((e) => e.source === a && e.target === y).coefficient;
    const ratio = value.measure === 'exponentiated_coefficient';
    formula = ratio ? `θ = exp(β_${a}→${y}). Current DGP value = ${Math.exp(coefficient)}`
      : `θ = β_${a}→${y} on the ${DISTRIBUTIONS[node.family].link} scale; current DGP value = ${coefficient}`;
    note = 'This is a conditional structural coefficient, rather than a marginal total effect. With ordinary products, it applies when the other interacting parents equal zero. Additional response terms can also contribute to effects and derivatives; this coefficient alone is not the full effect. The analysis must target this coefficient.';
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
  if (graph.study && value.kind === 'marginal') formula = `θ = Σ wⱼ {E[${y}ⱼ under do(${a} = ${hi})] − E[${y}ⱼ under do(${a} = ${lo})]}, w = [${contrastWeights(graph, value).map((x) => Number(x.toPrecision(4))).join(', ')}]`;
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


// What the power analysis targets, in plain words, for the card above the question: the target's name, what it
// means, over whom it averages, its formula, the hypotheses, and what the planned analysis estimates. `plan` is the
// compiled analysis plan (or null) and `side` the test's direction.
export function explainEstimand(graph, value, plan = null, side = 'two-sided') {
  validateEstimand(graph, value);
  const node = (id) => graph.nodes.find((n) => n.id === id);
  const a = node(value.exposure), y = node(value.outcome);
  const name = (n) => n.label && n.label !== n.id ? `${n.label} (${n.id})` : n.id;
  const A = name(a), Y = name(y);
  // "X1", "X1 and X2", "X1, X2 and X3", with labels
  const list = (ids) => { const names = ids.map((id) => name(node(id) ?? { id })); return names.length < 2 ? names.join('') : `${names.slice(0, -1).join(', ')} and ${names.at(-1)}`; };
  const relation = side === 'greater' ? '>' : side === 'less' ? '<' : '≠';
  const rows = [];
  if (value.kind === 'marginal') {
    const treatment = a.kind === 'binary' && value.reference === 0 && value.active === 1;
    const lead = treatment ? 'Average treatment effect (ATE)' : `Average causal effect of setting ${A} to ${value.active} versus ${value.reference}`;
    let what = `the mean of ${Y}`, words = `how much the average ${Y}`;
    if (graph.study) {
      const r = value.repeated ?? { contrast: 'final' }, times = graph.study.times;
      const at = r.contrast === 'change' ? `the change in ${Y} from time ${times[0]} to ${times.at(-1)}` : r.contrast === 'slope' ? `the time trend of ${Y}`
        : r.contrast === 'average' ? `${Y} averaged over the visits` : r.contrast === 'custom' ? `a weighted combination of ${Y} over the visits`
          : `${Y} at time ${times[r.visit ?? times.length - 1]}`;
      what = `the mean of ${at}`; words = `how much the average of ${at}`;
    }
    const m = value.measure, h = value.horizon;
    const target = {
      mean_difference: y.kind === 'binary' ? `${lead} on the risk of ${Y} (risk difference)` : `${lead} on ${what}`,
      risk_difference: `${lead} on the risk of ${Y} (risk difference)`,
      risk_ratio: `Causal risk ratio of ${Y}`,
      odds_ratio: `Marginal causal odds ratio of ${Y}`,
      survival_difference: `${lead} on the probability of being free of ${Y} at time ${h}`,
      rmst_difference: `${lead} on the restricted mean ${Y} time up to ${h}`,
    }[m];
    const verb = {
      mean_difference: y.kind === 'binary' ? `how much the share of the population with ${Y} = 1 would differ` : `${words} would differ`,
      risk_difference: `how much the share of the population with ${Y} = 1 would differ`,
      risk_ratio: `how many times larger the share with ${Y} = 1 would be`,
      odds_ratio: `how many times larger the odds of ${Y} = 1 would be`,
      survival_difference: `how much the share still free of the event at time ${h} would differ`,
      rmst_difference: `how much the average event-free time up to ${h} would differ`,
    }[m];
    rows.push(['Target', target],
      ['In words', `${verb[0].toUpperCase()}${verb.slice(1)} if everyone had ${A} = ${value.active} rather than everyone ${A} = ${value.reference}. The variables ${A} affects change with it; everything else keeps the distribution the graph gives it.`],
      ['Population', graph.study ? `Every subject the graph describes, each counted once${graph.study.clusterSizes?.length ? ', across clusters' : ''}; not a subgroup.` : 'Everyone the graph describes, at the distributions it generates; not a subgroup.']);
  } else if (value.kind === 'coefficient') {
    const ratio = value.measure === 'exponentiated_coefficient';
    const others = graph.edges.filter((e) => e.target === y.id && e.source !== a.id).map((e) => e.source);
    const scale = DISTRIBUTIONS[y.family].link;
    rows.push(['Target', `${ratio ? 'Exponentiated regression coefficient' : 'Regression coefficient'} of ${A} in the model for ${Y}: a conditional effect, not a population average`],
      ['In words', `The change in the ${scale} of ${Y} per unit of ${A}${ratio ? ', as a ratio' : ''}${others.length ? `, holding ${list(others)} fixed` : ''}. It equals the average effect only when the model for ${Y} is linear in ${A} with no interactions or curved terms involving it.`],
      ['Population', others.length ? `People with the same ${list(others)}: the effect is conditional on them.` : 'Everyone the graph describes.']);
  } else {
    const which = value.measure === 'nonlinearity' ? `the terms that bend the effect of ${A}` : value.measure === 'interaction' ? `the terms in which ${A} multiplies another variable` : `all terms involving ${A}`;
    rows.push(['Target', `No single effect: a test of whether ${which} in the model for ${Y} are all zero`],
      ['In words', value.measure === 'interaction' ? `Whether the effect of ${A} on ${Y} depends on another variable (effect modification).`
        : value.measure === 'nonlinearity' ? `Whether the effect of ${A} on ${Y} is curved rather than a straight line.` : `Whether ${A} is related to ${Y} at all, once the other variables in the model are accounted for.`],
      ['Population', 'Conditional on the other variables in the model.']);
  }
  const d = describeEstimand(graph, value);
  rows.push(['Formula', d.formula]);
  rows.push(['Hypotheses', value.kind === 'joint' ? 'H₀: all of those coefficients are 0; H₁: at least one is not.' : `H₀: θ = ${value.nullValue}; H₁: θ ${relation} ${value.nullValue}.`]);
  if (plan?.supported !== false && plan?.label) {
    const used = [...new Set([...(plan.terms ?? []), ...(plan.propensityTerms ?? [])].flatMap((t) => [t.left, ...(t.right ? [t.right] : [])]))].filter((id) => id !== a.id);
    const how = ['welch', 'risk_score', 'risk_wald', 'risk_ratio', 'odds_ratio', 'km_survival', 'km_rmst'].includes(plan.method)
      ? `It compares the two groups directly, which estimates θ when ${A} is assigned at random (nothing in the graph causes it).`
      : ['aipw', 'tmle'].includes(plan.method) ? `It adjusts for ${used.length ? list(used) : 'nothing measured'} and estimates θ if ${used.length === 1 ? 'it blocks' : 'they block'} every back-door path from ${A} to ${Y}${plan.learner === 'superlearner' ? '; its nuisance models are learned, not assumed' : ' and one of its models is right'}.`
        : plan.method === 'gformula' ? `It refits the model${(plan.models?.length ?? 1) > 1 ? 's' : ''} for ${list((plan.models ?? []).map((x) => x.node)) || Y} and averages its predictions with ${A} set to each level; it estimates θ if those models are right.`
          : ['mmrm', 'lmm', 'trajectory'].includes(plan.method) ? `It compares the randomized groups' outcomes over the visits, which estimates θ.`
            : value.kind === 'marginal' ? '' : `It fits the model ${plan.formula ?? ''}.`;
    const caution = (plan.notes ?? []).some((n) => /unmeasured|unknown form that the planned|lies upstream|No arrow path|censors .* depends on/.test(n))
      ? ' The notes below say what it cannot account for; the bias and coverage in the results show the cost.' : '';
    rows.push(['Analysis', `${plan.label}. ${how}${caution}`.trim()]);
  }
  return rows;
}
