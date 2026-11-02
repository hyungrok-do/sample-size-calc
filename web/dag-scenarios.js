import { exampleGraph, parseGraph } from './dag-model.js?v=dev';
import { defaultTarget } from './dag-power.js?v=dev';
import { DEFAULT_POWER, parsePowerSettings } from './power-settings.js?v=dev';

// Each example is a whole design: loading it never inherits another example's
// estimand, estimator or sample-size axis. Values are illustrative assumptions.
export const SCENARIOS = [
  ['trial', 'Common designs', 'Two-group continuous outcome', 'One measurement per person; randomized treatment.', 'Welch inference. Change the effect and residual variation to fit your study.'],
  ['repeated', 'Repeated and clustered', 'Repeated continuous outcome', 'Four visits with subject intercepts and slopes.', 'MMRM with an unstructured covariance and approximate t inference.'],
  ['clustered-repeated', 'Repeated and clustered', 'Clustered repeated outcome', 'Treatment assigned to groups of unequal size.', 'Sample size counts clusters. Robust inference is approximate with few clusters.'],
  ['repeated-binary', 'Repeated and clustered', 'Repeated binary outcome', 'Correlated binary responses at four visits.', 'Available-visit risk contrasts with subject-level inference.'],
  ['mediation', 'Common designs', 'Continuous mediation', 'A baseline treatment affects an intermediate variable and the outcome.', 'The target is the total effect, not a natural direct or indirect effect.'],
  ['binary', 'Common designs', 'Binary outcome', 'Randomized treatment and a prognostic baseline covariate.', 'A marginal risk difference; rare events may need much larger samples.'],
  ['survival', 'Common designs', 'Survival with censoring', 'Weibull event times and independently generated censoring.', 'Single-event survival and RMST are supported; competing risks are not.'],
  ['nonlinear-dose', 'Model sensitivity', 'Nonlinear dose response', 'A curved response to a continuous exposure.', 'Compare two doses or change the question to an overall/nonlinearity test.'],
  ['nonlinear-confounding', 'Model sensitivity', 'Nonlinear confounding', 'Treatment and outcome depend nonlinearly on a baseline covariate.', 'Change the two nuisance models separately to inspect bias from misspecification.'],
  ['heterogeneous-null', 'Model sensitivity', 'Zero average, heterogeneous effects', 'Positive and negative effects cancel on average.', 'An average-effect null can retain subgroup effects.'],
  ['tmle-continuous', 'TMLE', 'TMLE · continuous outcome', 'Confounded baseline treatment, nonlinear outcome and three-fold cross-fitting.', 'Identity fluctuation for an unbounded mean; quadratic nuisance regressions.'],
  ['tmle-binary', 'TMLE', 'TMLE · binary outcome', 'Adjusted marginal risk difference with three-fold cross-fitting.', 'Logistic targeting preserves probability bounds; Wald inference needs adequate events.'],
  ['tmle-bounded', 'TMLE', 'TMLE · bounded score', 'A score with known 0–100 limits under randomized treatment.', 'Logistic fluctuation on a prespecified scale. Bounds must cover every possible outcome.'],
  ['tmle-repeated', 'TMLE', 'TMLE · confounded repeated outcomes', 'Average of complete visits with baseline confounding.', 'Point-treatment TMLE of a subject summary; no time-varying treatment or dropout adjustment.'],
  ['tmle-clustered', 'TMLE', 'TMLE · unequal clusters', 'Average repeated outcome under cluster-level treatment.', 'Subject-weighted effect. Cross-fitting holds out entire clusters; sample size counts clusters.'],
  ['tmle-repeated-binary', 'TMLE', 'TMLE · average repeated risk', 'Average risk across complete repeated binary responses.', 'The subject summary is bounded between zero and one; dependence is retained in inference.'],
  ['rare-binary', 'Stress tests', 'Rare binary events', 'A low event probability and a treatment effect on the odds scale.', 'TMLE can fail with too few events. Inspect failures, null rejection and coverage, not power alone.'],
  ['weak-overlap', 'Stress tests', 'Weak treatment overlap', 'Strong selection into treatment and an explicit propensity bound.', 'Compare bounds and sample sizes. Truncation trades extreme weights for potential bias.'],
  ['unequal-allocation', 'Common designs', 'Unequal treatment allocation', 'About one participant in five receives active treatment.', 'Sample sizes are totals across both groups; allocation is generated, not forced to fixed counts.'],
  ['skewed-outcome', 'Stress tests', 'Strongly skewed positive outcome', 'Lognormal outcomes with a long right tail.', 'Unbounded continuous TMLE targets a mean. A large sample may be needed for reliable Wald inference.'],
  ['count-confounding', 'Common designs', 'Confounded count outcome', 'Poisson counts and nonlinear baseline confounding.', 'AIPW targets a marginal count difference. This is not a negative-binomial or zero-inflated model.'],
  ['irregular-visits', 'Repeated and clustered', 'Irregular visits and correlation', 'Unequally spaced visits with residual correlation and random slopes.', 'The contrast uses actual visit times. Examine covariance misspecification with LMM versus MMRM.'],
  ['mar-dropout', 'Missing observations', 'History-dependent dropout', 'Missing later visits depend on a previously observed outcome.', 'MMRM assumes an adequate mean/covariance model and ignorable missingness. TMLE with dropout is not implemented.'],
  ['mnar-dropout', 'Missing observations', 'Unobserved-outcome dropout', 'Dropout depends on the current unobserved response.', 'A deliberate bias sensitivity scenario. Neither MMRM nor the available TMLE identifies this effect without extra assumptions.'],
  ['intermittent-visits', 'Missing observations', 'Intermittent missed visits', 'Subjects can miss one visit and return later.', 'LMM uses observed visits. The monotone MMRM and complete-summary TMLE paths reject this design.'],
  ['lagged-response', 'Repeated and clustered', 'Lagged outcome response', 'Each response depends on its previous value.', 'A lagged generating model with a final-visit contrast; change the lag strength and visit schedule.'],
  ['few-clusters', 'Stress tests', 'Few, unequal clusters', 'A small number of independent groups with imbalanced sizes.', 'Small-sample stress test, not a validated correction. Check null rejection and coverage.'],
  ['empty', 'Common designs', 'Empty canvas', 'Build a custom data-generating graph.', 'Add a comparison, outcome and analysis before computing power.'],
].map(([id, category, label, description, assumption]) => ({ id, category, label, description, assumption }));

export function scenarioDesign(id) {
  if (!SCENARIOS.some((s) => s.id === id)) throw new Error('Unknown study scenario.');
  const base = {
    'tmle-continuous': 'nonlinear-confounding', 'tmle-binary': 'binary', 'tmle-bounded': 'trial',
    'tmle-repeated': 'repeated', 'tmle-clustered': 'clustered-repeated', 'tmle-repeated-binary': 'repeated-binary',
    'rare-binary': 'binary', 'weak-overlap': 'nonlinear-confounding', 'unequal-allocation': 'trial',
    'skewed-outcome': 'trial', 'count-confounding': 'nonlinear-confounding', 'irregular-visits': 'repeated',
    'mar-dropout': 'repeated', 'mnar-dropout': 'repeated', 'intermittent-visits': 'repeated', 'lagged-response': 'repeated', 'few-clusters': 'clustered-repeated',
  }[id] ?? id;
  const graph = exampleGraph(base), a = graph.nodes.find((n) => n.id === 'A'), y = graph.nodes.find((n) => n.id === 'Y');
  const power = structuredClone(DEFAULT_POWER);
  power.sampleSizes = graph.study ? graph.study.clusterSizes.length ? [24, 48, 96] : [100, 200, 400] : [100, 200, 400];
  if (id.startsWith('tmle-') || ['rare-binary', 'weak-overlap', 'skewed-outcome'].includes(id)) power.analysis = { method: 'tmle', folds: 3 };
  if (id === 'tmle-continuous') power.analysis.basis = 'quadratic';
  if (id === 'tmle-binary' || id === 'tmle-repeated') graph.edges.push({ source: 'X', target: 'A', coefficient: .5 });
  if (id === 'tmle-bounded') {
    y.family = 'uniform'; y.intercept = 50; y.width = 50; graph.edges[0].coefficient = 5;
    power.analysis.outcomeBounds = [0, 100];
  }
  if (id === 'tmle-clustered') { graph.study.clusterSizes = [4, 12, 24]; power.sampleSizes = [30, 60, 100]; }
  if (id === 'rare-binary') {
    y.intercept = -4; graph.edges.find((e) => e.source === 'A').coefficient = .7;
    power.sampleSizes = [500, 1000, 2000]; power.analysis.folds = 1;
  }
  if (id === 'weak-overlap') {
    a.functions = []; a.intercept = 0; graph.edges.find((e) => e.target === 'A').coefficient = 3.5;
    power.analysis.propensityBound = .025; power.analysis.basis = 'quadratic'; power.sampleSizes = [250, 500, 1000];
  }
  if (id === 'unequal-allocation') a.intercept = Math.log(.2 / .8);
  if (id === 'skewed-outcome') { y.family = 'lognormal'; y.sd = 1.2; graph.edges[0].coefficient = .2; }
  if (id === 'count-confounding') {
    y.kind = 'count'; y.family = 'poisson'; y.intercept = .2; y.functions[0].coefficient = .15;
    graph.edges.find((e) => e.source === 'A').coefficient = .25; power.analysis = { method: 'aipw', folds: 3 };
  }
  if (id === 'irregular-visits') { graph.study.times = [0, 1, 4, 10]; y.repeated.residualCorrelation = .45; }
  if (['mar-dropout', 'mnar-dropout', 'intermittent-visits'].includes(id)) {
    graph.study.observation = { ...graph.study.observation, mode: id === 'intermittent-visits' ? 'intermittent' : 'dropout', intercept: -2,
      historyVariable: 'Y', historyCoefficient: id === 'mar-dropout' ? .5 : 0, currentCoefficient: id === 'mnar-dropout' ? .5 : 0 };
    power.analysis = { method: id === 'intermittent-visits' ? 'lmm' : 'mmrm', timeBasis: 'categorical', randomSlope: true };
  }
  if (id === 'lagged-response') y.lags = [{ source: 'Y', lag: 1, coefficient: .4, initial: 0 }];
  if (id === 'few-clusters') { graph.study.clusterSizes = [3, 8, 25]; power.sampleSizes = [8, 12, 20]; }
  const parsed = parseGraph(graph), estimand = defaultTarget(parsed);
  if (id.startsWith('tmle-') && graph.study) estimand.repeated = { contrast: 'average' };
  if (id === 'heterogeneous-null') power.nullMode = 'average_effect';
  return { graph: parsed, estimand, power: parsePowerSettings(power) };
}
