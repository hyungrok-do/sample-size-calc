import { exampleGraph, parseGraph, makeNode, DAG_VERSION } from './dag-model.js?v=dev';
import { MEASURES, validateEstimand } from './estimand-model.js?v=dev';
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
  ['ancova', 'Common designs', 'Baseline-adjusted outcome (pre–post)', 'The outcome is measured before and after randomized treatment.', 'ANCOVA: the treatment coefficient adjusted for baseline. Change the question to a mean difference to see the unadjusted test.'],
  ['factorial', 'Common designs', '2 × 2 factorial trial', 'Two randomized treatments, given alone or together.', 'The main effect of A averages over B. Change the question to an interaction test to plan for synergy.'],
  ['noninferiority', 'Common designs', 'Non-inferiority', 'A new treatment expected to be as good as control, with a margin of 0.3.', 'One-sided test against the margin −0.3; power is the chance of showing non-inferiority when the true difference is zero.'],
  ['effect-modification', 'Common designs', 'Effect modification', 'The treatment effect differs by a measured biomarker.', 'Tests the treatment × biomarker interaction. Interactions typically need about four times the sample of a main effect.'],
  ['cluster-trial', 'Repeated and clustered', 'Cluster randomized trial', 'Clinics are randomized; each patient is measured once.', 'Sample size counts clinics of 20 patients; the intracluster correlation is about 0.08. Robust inference is approximate with few clusters.'],
  ['staggered-entry', 'Outcome types', 'Staggered entry, common end date', 'Participants enrol over two years and are followed until the study closes.', 'Follow-up is uniform from 2 to 4 years and independent of events; a horizon of 2 keeps Kaplan–Meier unbiased.'],
  ['overdispersed-count', 'Outcome types', 'Overdispersed counts', 'Event counts more variable than Poisson (negative binomial).', 'Rate ratio by quasi-Poisson regression with a sandwich standard error, which stays valid under overdispersion.'],
  ['rate-followup', 'Outcome types', 'Event rates with unequal follow-up', 'Counts accumulate over follow-up times from 6 to 18 months.', 'log(follow-up) enters the rate with coefficient 1 (an offset); the analysis estimates its coefficient alongside the rate ratio.'],
  ['heavy-tails', 'Outcome types', 'Heavy-tailed outcome', 'Student t errors with 3 degrees of freedom: occasional extreme values.', 'Welch’s test stays valid, but outliers cost power; compare with the two-group continuous design.'],
  ['proportion-outcome', 'Outcome types', 'Proportion outcome (0 to 1)', 'A share such as the fraction of days adherent, with a prognostic covariate.', 'Beta-distributed values; the effect is a difference in mean proportion on the 0–1 scale.'],
  ['unmeasured-confounding', 'Unmeasured variables', 'Unmeasured confounding', 'An unrecorded U affects both treatment and outcome; X is recorded.', 'AIPW adjusts for X only. The Checks column shows the bias, coverage and false positives that U leaves.'],
  ['measurement-error', 'Unmeasured variables', 'Confounder measured with error', 'Only an error-prone version W of the true severity is recorded.', 'Adjusting for W removes part of the confounding; the rest shows as bias. Lower W’s error SD to see it shrink.'],
  ['zero-inflated', 'Unmeasured variables', 'Zero-inflated counts', 'An unrecorded 30% are never at risk and always count zero.', 'The difference in mean counts stays valid; a Poisson model that ignores the extra zeros would understate its standard error.'],
  ['responders', 'Unmeasured variables', 'Responders and non-responders', 'Only an unrecorded 40% respond to treatment.', 'Responders improve by 0.8, so the average effect is 0.32. The trial estimates the average, not who responds.'],
  ['empty', 'Common designs', 'Empty canvas', 'Build a custom data-generating graph.', 'Add a comparison, outcome and analysis before computing power.'],
].map(([id, category, label, description, assumption]) => ({ id, category, label, description, assumption }));

// A variable for the designs below; `extra` sets its distribution and parameters.
const v = (id, label, kind, x, y, extra = {}) => ({ ...makeNode(id, kind, x, y), label, ...extra });
const logit = (p) => Math.log(p / (1 - p));
// Designs drawn here rather than derived from an example graph: the graph, and the question and power settings that
// differ from the defaults.
const DRAWN = {
  ancova: () => ({ graph: { version: DAG_VERSION, nodes: [v('B', 'Baseline score', 'continuous', 60, 60), v('A', 'Treatment', 'binary', 60, 280), v('Y', 'Follow-up score', 'continuous', 520, 170, { sd: .8 })],
    edges: [{ source: 'B', target: 'Y', coefficient: .6 }, { source: 'A', target: 'Y', coefficient: .3 }] },
  estimand: { kind: 'coefficient', measure: 'coefficient' } }),
  factorial: () => ({ graph: { version: DAG_VERSION, nodes: [v('A', 'Drug', 'binary', 60, 60), v('B', 'Exercise', 'binary', 60, 280),
    v('Y', 'Outcome', 'continuous', 520, 170, { interactions: [{ left: 'A', right: 'B', coefficient: .2 }] })],
    edges: [{ source: 'A', target: 'Y', coefficient: .4 }, { source: 'B', target: 'Y', coefficient: .3 }] } }),
  noninferiority: () => ({ graph: { version: DAG_VERSION, nodes: [v('A', 'New treatment', 'binary', 240, 180), v('Y', 'Outcome', 'continuous', 680, 180)],
    edges: [{ source: 'A', target: 'Y', coefficient: 0 }] },
  estimand: { nullValue: -.3 }, power: { side: 'greater', sampleSizes: [200, 400, 800] } }),
  'effect-modification': () => ({ graph: { version: DAG_VERSION, nodes: [v('X', 'Biomarker', 'binary', 60, 60), v('A', 'Treatment', 'binary', 60, 280),
    v('Y', 'Outcome', 'continuous', 520, 170, { interactions: [{ left: 'A', right: 'X', coefficient: .5 }] })],
    edges: [{ source: 'X', target: 'Y', coefficient: .2 }, { source: 'A', target: 'Y', coefficient: .3 }] },
  estimand: { kind: 'joint', measure: 'interaction' }, power: { sampleSizes: [300, 600, 1200] } }),
  'cluster-trial': () => ({ graph: { version: DAG_VERSION, study: { times: [0], clusterSizes: [20] },
    nodes: [v('A', 'Clinic arm', 'binary', 60, 180, { scope: 'cluster' }), v('X', 'Age', 'continuous', 340, 60, { scope: 'subject' }),
      v('Y', 'Outcome', 'continuous', 600, 180, { scope: 'visit', repeated: { subjectSD: 0, slopeSD: 0, slopeCorrelation: 0, clusterSD: .3, residualCorrelation: 0 } })],
    edges: [{ source: 'A', target: 'Y', coefficient: .3 }, { source: 'X', target: 'Y', coefficient: .3 }] },
  power: { sampleSizes: [10, 20, 40] } }),
  'staggered-entry': () => ({ graph: { version: DAG_VERSION, nodes: [v('A', 'Treatment', 'binary', 60, 60),
    v('T', 'Time to event', 'survival', 520, 60, { shape: 1.3, scale: 5, censoring: { mode: 'variable', time: 20, variable: 'F' } }),
    v('F', 'Follow-up', 'continuous', 520, 280, { family: 'uniform', intercept: 3, width: 2 })],
    edges: [{ source: 'A', target: 'T', coefficient: Math.log(.65) }] },
  estimand: { measure: 'survival_difference', horizon: 2 }, power: { sampleSizes: [300, 600, 1200] } }),
  'overdispersed-count': () => ({ graph: { version: DAG_VERSION, nodes: [v('A', 'Treatment', 'binary', 240, 180),
    v('Y', 'Exacerbations', 'count', 680, 180, { family: 'negbin', intercept: Math.log(2), overdispersion: .8 })],
    edges: [{ source: 'A', target: 'Y', coefficient: Math.log(.75) }] },
  estimand: { kind: 'coefficient', measure: 'exponentiated_coefficient', nullValue: 1 }, power: { sampleSizes: [200, 400, 800] } }),
  'rate-followup': () => ({ graph: { version: DAG_VERSION, nodes: [v('A', 'Treatment', 'binary', 60, 60),
    v('F', 'Follow-up', 'continuous', 60, 280, { family: 'uniform', intercept: 1, width: 1 }),
    v('Y', 'Events', 'count', 520, 170, { intercept: Math.log(1.5), functions: [{ left: 'F', transform: { kind: 'log' }, coefficient: 1 }] })],
    edges: [{ source: 'A', target: 'Y', coefficient: Math.log(.7) }, { source: 'F', target: 'Y', coefficient: 0 }] },
  estimand: { kind: 'coefficient', measure: 'exponentiated_coefficient', nullValue: 1 }, power: { sampleSizes: [100, 200, 400] } }),
  'heavy-tails': () => ({ graph: { version: DAG_VERSION, nodes: [v('A', 'Treatment', 'binary', 240, 180),
    v('Y', 'Outcome', 'continuous', 680, 180, { family: 'student', sd: .8, df: 3 })],
    edges: [{ source: 'A', target: 'Y', coefficient: .5 }] } }),
  'proportion-outcome': () => ({ graph: { version: DAG_VERSION, nodes: [v('X', 'Prior adherence', 'continuous', 60, 60), v('A', 'Reminder app', 'binary', 60, 280),
    v('Y', 'Days adherent', 'continuous', 520, 170, { family: 'beta', intercept: logit(.6), precision: 8 })],
    edges: [{ source: 'X', target: 'Y', coefficient: .4 }, { source: 'A', target: 'Y', coefficient: .35 }] } }),
  'unmeasured-confounding': () => ({ graph: { version: DAG_VERSION, nodes: [v('U', 'Frailty', 'continuous', 300, 30, { measured: false }),
    v('X', 'Age', 'continuous', 60, 150), v('A', 'Treatment', 'binary', 60, 330, { intercept: 0 }), v('Y', 'Outcome', 'continuous', 600, 230)],
    edges: [{ source: 'U', target: 'A', coefficient: .8 }, { source: 'U', target: 'Y', coefficient: .5 }, { source: 'X', target: 'A', coefficient: .5 },
      { source: 'X', target: 'Y', coefficient: .4 }, { source: 'A', target: 'Y', coefficient: .3 }] } }),
  'measurement-error': () => ({ graph: { version: DAG_VERSION, nodes: [v('S', 'True severity', 'continuous', 60, 30, { measured: false }),
    v('W', 'Recorded severity', 'continuous', 340, 30, { sd: .7 }), v('A', 'Treatment', 'binary', 60, 300), v('Y', 'Outcome', 'continuous', 600, 220)],
    edges: [{ source: 'S', target: 'W', coefficient: 1 }, { source: 'S', target: 'A', coefficient: .8 }, { source: 'S', target: 'Y', coefficient: .6 },
      { source: 'A', target: 'Y', coefficient: .3 }] },
  // the analyst adjusts for the recorded proxy; the automatic analysis adjusts only for measured causes of treatment
  power: { analysis: { method: 'aipw', basis: 'custom', terms: [{ left: 'W' }], propensity: { basis: 'custom', terms: [{ left: 'W' }] } } } }),
  'zero-inflated': () => ({ graph: { version: DAG_VERSION, nodes: [v('R', 'At risk', 'binary', 60, 60, { intercept: logit(.7), measured: false }),
    v('A', 'Treatment', 'binary', 60, 280), v('Y', 'Events', 'count', 520, 170, { intercept: Math.log(2) - 6 })],
    edges: [{ source: 'R', target: 'Y', coefficient: 6 }, { source: 'A', target: 'Y', coefficient: Math.log(.7) }] },
  power: { sampleSizes: [200, 400, 800] } }),
  responders: () => ({ graph: { version: DAG_VERSION, nodes: [v('R', 'Responder', 'binary', 60, 60, { intercept: logit(.4), measured: false }),
    v('A', 'Treatment', 'binary', 60, 280), v('Y', 'Outcome', 'continuous', 520, 170, { interactions: [{ left: 'A', right: 'R', coefficient: .8 }] })],
    edges: [{ source: 'R', target: 'Y', coefficient: .2 }, { source: 'A', target: 'Y', coefficient: 0 }] } }),
};

export function scenarioDesign(id) {
  if (!SCENARIOS.some((s) => s.id === id)) throw new Error('Unknown study scenario.');
  if (DRAWN[id]) {
    const drawn = DRAWN[id](), graph = parseGraph(drawn.graph), power = { ...structuredClone(DEFAULT_POWER), ...drawn.power };
    const start = defaultTarget(graph), estimand = { ...start, ...drawn.estimand };
    if (drawn.estimand?.kind && drawn.estimand.kind !== 'marginal') Object.assign(estimand, { reference: 0, active: 1, horizon: null, alternative: null,
      nullValue: drawn.estimand.nullValue ?? MEASURES[estimand.measure].nullValue });
    return { graph, estimand: validateEstimand(graph, estimand), power: parsePowerSettings(power) };
  }
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
