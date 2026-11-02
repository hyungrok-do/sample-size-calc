// Declarative study structure and UI validation. All sampling and fitting run in C++.
const fail = (message) => { throw new Error(message); };
const finite = Number.isFinite;
export const scopeOf = (node) => node.scope ?? 'subject';
export const sizeUnit = (graph) => graph.study?.clusterSizes?.length ? 'clusters' : 'subjects';
export const repeatedKeys = ['subjectSD', 'slopeSD', 'slopeCorrelation', 'clusterSD', 'residualCorrelation'];
export function parseStudy(raw) {
  if (!raw || typeof raw !== 'object' || Array.isArray(raw)) fail('Study structure must be an object.');
  const times = raw.times;
  if (!Array.isArray(times) || times.length < 1 || times.length > 12 || times.some((t, i) => !finite(t) || i && t <= times[i - 1]))
    fail('Specify 1–12 finite, increasing visit times.');
  const clusterSizes = raw.clusterSizes ?? [];
  if (!Array.isArray(clusterSizes) || clusterSizes.length > 100 || clusterSizes.some((x) => !Number.isInteger(x) || x < 1 || x > 1000)) fail('Cluster sizes must be integers from 1 to 1,000.');
  const o = raw.observation ?? {};
  if (!o || typeof o !== 'object' || Array.isArray(o)) fail('Observation settings must be an object.');
  const observation = { mode: o.mode ?? 'none', intercept: o.intercept ?? -3, historyVariable: o.historyVariable ?? '',
    historyCoefficient: o.historyCoefficient ?? 0, currentCoefficient: o.currentCoefficient ?? 0,
    treatmentVariable: o.treatmentVariable ?? '', treatmentCoefficient: o.treatmentCoefficient ?? 0 };
  if (!['none', 'dropout', 'intermittent'].includes(observation.mode)) fail('Choose no missing visits, dropout or intermittent missing visits.');
  if (['intercept', 'historyCoefficient', 'currentCoefficient', 'treatmentCoefficient'].some((k) => !finite(observation[k]))) fail('Missing-visit coefficients must be finite.');
  return { times: [...times], clusterSizes: [...clusterSizes], observation };
}
export function validateStudy(graph) {
  if (!graph.study) {
    if (graph.nodes.some((n) => n.scope && n.scope !== 'subject' || n.lags?.length || repeatedKeys.some((k) => n.repeated?.[k]))) fail('Enable study structure before setting repeated or cluster variables.');
    return;
  }
  const s = parseStudy(graph.study), byId = new Map(graph.nodes.map((n) => [n.id, n]));
  const level = (n) => scopeOf(n) === 'cluster' ? 0 : scopeOf(n) === 'subject' ? 1 : 2;
  const reserved = new Set(['subject_id', 'cluster_id', 'visit_index', 'visit_time']);
  for (const n of graph.nodes) { if (reserved.has(n.id)) fail('A variable ID conflicts with a study identifier column.'); reserved.add(n.id); }
  for (const n of graph.nodes) {
    const scope = scopeOf(n), r = n.repeated ?? {};
    if (!['subject', 'cluster', 'visit', 'time'].includes(scope)) fail('Choose cluster, subject, visit or time scope.');
    if (scope === 'cluster' && !s.clusterSizes.length) fail('Cluster variables require cluster sizes.');
    if (scope === 'time' && (n.family !== 'gaussian' || graph.edges.some((e) => e.target === n.id))) fail('Visit time must be a continuous Normal root; its values come from the schedule.');
    if (scope === 'visit' && n.kind === 'survival') fail('Repeated event histories need a separate model. Use a subject-level survival variable.');
    if (repeatedKeys.some((k) => !finite(r[k] ?? 0))) fail('Repeated-measure parameters must be finite.');
    if (['subjectSD', 'slopeSD', 'clusterSD'].some((k) => (r[k] ?? 0) < 0) || Math.abs(r.slopeCorrelation ?? 0) > 1 || (r.residualCorrelation ?? 0) < 0 || (r.residualCorrelation ?? 0) >= 1)
      fail('Random-effect SDs must be nonnegative; slope correlation in [-1,1]; residual correlation in [0,1).');
    if (repeatedKeys.some((k) => r[k]) && scope !== 'visit') fail('Repeated random effects require a visit-level variable.');
    if (r.clusterSD && !s.clusterSizes.length) fail('A cluster random effect requires cluster sizes.');
    if (r.residualCorrelation && !['gaussian', 'lognormal'].includes(n.family)) fail('Residual correlation is available for Normal and log-normal errors.');
    if (n.lags != null && (!Array.isArray(n.lags) || n.lags.length > 40)) fail('Use at most 40 lagged relationships per variable.');
    const used = new Set();
    for (const l of n.lags ?? []) {
      const key = `${l.source}:${l.lag}`;
      if (scope !== 'visit' || !byId.has(l.source) || !Number.isInteger(l.lag) || l.lag < 1 || l.lag > 12 || !finite(l.coefficient) || !finite(l.initial ?? 0) || used.has(key)) fail('Lagged relationships need a visit-level child, valid source, lag 1–12 and finite coefficients.');
      used.add(key);
    }
    if (scope === 'visit') for (const suffix of ['_observed', '_recorded']) {
      if (reserved.has(n.id + suffix)) fail('A repeated observation column conflicts with a variable ID.'); reserved.add(n.id + suffix);
    }
  }
  for (const e of graph.edges) if (level(byId.get(e.source)) > level(byId.get(e.target))) fail('A variable cannot depend on a finer observation level.');
  const o = s.observation;
  if (o.historyVariable && scopeOf(byId.get(o.historyVariable) ?? {}) !== 'visit') fail('Missing visits must reference a visit-level history variable.');
  if (o.treatmentVariable && (!byId.has(o.treatmentVariable) || level(byId.get(o.treatmentVariable)) > 1)) fail('Choose a subject- or cluster-level treatment variable.');
  if ((o.historyCoefficient || o.currentCoefficient) && !o.historyVariable) fail('Choose the variable determining missing visits.');
  if (o.treatmentCoefficient && !o.treatmentVariable) fail('Choose the treatment variable determining missing visits.');
}
export function contrastWeights(graph, estimand) {
  const t = graph.study.times, r = estimand.repeated ?? {}, kind = r.contrast ?? 'final', w = t.map(() => 0);
  if (kind === 'final') { const v = r.visit ?? t.length - 1; if (!Number.isInteger(v) || v < 0 || v >= t.length) fail('Choose a visit in the schedule.'); w[v] = 1; }
  else if (kind === 'average') w.fill(1 / t.length);
  else if (kind === 'change') { if (t.length < 2) fail('Change needs at least two visits.'); w[0] = -1; w[t.length - 1] = 1; }
  else if (kind === 'slope') {
    if (t.length < 2) fail('A trend needs at least two visits.');
    const mean = t.reduce((a, b) => a + b, 0) / t.length, sum = t.reduce((a, b) => a + (b - mean) ** 2, 0);
    t.forEach((x, i) => { w[i] = (x - mean) / sum; });
  } else if (kind === 'custom') {
    if (!Array.isArray(r.weights) || r.weights.length !== t.length || r.weights.some((x) => !finite(x))) fail('Custom contrast weights must match the visit schedule.');
    w.splice(0, w.length, ...r.weights);
  } else fail('Choose final visit, average, change, linear trend or custom contrast.');
  if (w.every((x) => x === 0)) fail('At least one contrast weight must be nonzero.');
  return w;
}
export function studyPlan(graph, e, spec) {
  const a = graph.nodes.find((n) => n.id === e.exposure), y = graph.nodes.find((n) => n.id === e.outcome);
  if (e.kind !== 'marginal' || !['mean_difference', 'risk_difference'].includes(e.measure)) fail('Structured studies currently estimate marginal mean/risk contrasts across visits.');
  if (a.kind !== 'binary' || !['subject', 'cluster'].includes(scopeOf(a))) fail('Structured power requires a binary baseline subject- or cluster-level exposure. Time-varying treatment needs a separate longitudinal estimator.');
  if (y.kind === 'survival' || scopeOf(y) === 'time') fail('Select a nonsurvival measured outcome for a repeated contrast.');
  const grouped = Boolean(graph.study.clusterSizes?.length), T = graph.study.times.length, continuous = y.kind === 'continuous';
  const method = spec.method === 'auto' ? !grouped && continuous && scopeOf(y) === 'visit' && T > 1 && graph.study.observation?.mode !== 'intermittent' ? 'mmrm' : 'trajectory' : spec.method;
  if (method !== 'tmle' && graph.edges.some((edge) => edge.target === a.id)) fail('For baseline confounding, select TMLE with complete outcomes. Other structured analyses require randomized exposure.');
  if (method === 'tmle' && graph.study.observation?.mode && graph.study.observation.mode !== 'none') fail('Trajectory TMLE requires complete outcomes. Dropout and intermittent missingness need a longitudinal missing-data estimator; this TMLE does not provide one.');
  if (!['mmrm', 'lmm', 'trajectory', 'tmle'].includes(method)) fail('Choose MMRM, linear mixed model or trajectory means for a structured study.');
  if (['mmrm', 'lmm'].includes(method) && (grouped || !continuous || scopeOf(y) !== 'visit')) fail('Mixed-model trajectory analyses require continuous outcomes and independent subjects. Use trajectory means for clustered or discrete outcomes.');
  if (method === 'mmrm' && graph.study.observation?.mode === 'intermittent') fail('This MMRM fitter supports monotone dropout. Use LMM or trajectory means for intermittent missing visits.');
  if (method === 'lmm' && (T < 2 || (spec.timeBasis ?? 'categorical') === 'categorical' && T > 6)) fail('The LMM supports 2–6 categorical visits or 2–12 visits with linear time.');
  if (method !== 'tmle' && (spec.folds !== 1 || spec.basis !== 'match' || spec.propensity || spec.models)) fail('Structured analyses use their trajectory settings, not cross-sectional nuisance-model terms.');
  const notes = ['Power is conditional on the generating distributions, response functions, dependence and observation process. Monte Carlo intervals measure simulation error, not uncertainty in those assumptions.',
    grouped ? 'Sample size counts independent clusters. Effects average over subjects; unequal cluster sizes follow the saved repeating size pattern.' : 'Sample size counts independent subjects; visits are repeated observations of those subjects.'];
  if (method === 'tmle') notes.push('TMLE targets the subject-weighted mean of the complete outcome summary under baseline treatment. It adjusts only for baseline covariates. Folds and influence-function inference preserve independent subjects or entire clusters. It is not longitudinal TMLE for time-varying treatment or dropout.');
  if (method === 'mmrm') notes.push('MMRM fits arm-by-visit means and a common unstructured covariance; t inference uses last-visit sample size minus two, not Kenward–Roger degrees of freedom.');
  if (method === 'lmm') notes.push('LMM estimates covariance by REML with Satterthwaite inference. The planned time model and random effects are separate from the generator.');
  if (method === 'trajectory') notes.push('Trajectory means use a CR1 sandwich at the highest independent unit with units minus two degrees of freedom; inference is approximate with few units.');
  if (graph.study.observation?.mode && graph.study.observation.mode !== 'none') notes.push(method === 'trajectory' ? 'Available-visit means require noninformative missingness within arm and visit. History-dependent dropout can bias this analysis.' : 'Likelihood analyses rely on the planned mean/covariance model and ignorable missingness. Selection on an unobserved current outcome can bias them.');
  return { method, label: { mmrm: 'MMRM contrast of arm-by-visit means', lmm: 'Linear mixed-model trajectory contrast', trajectory: 'Trajectory means with cluster-robust inference' }[method],
    formula: `${y.id} ~ arm × ${method === 'lmm' && spec.timeBasis === 'linear' ? 'time' : 'visit'}`, notes, weights: contrastWeights(graph, e), terms: [], analysis: spec, folds: 1, testScale: 'identity', family: y.family, exposure: a.id, outcome: y.id };
}
export function updateStudyReferences(study, from, to) {
  if (!study) return undefined;
  const s = parseStudy(study), o = s.observation;
  if (o.historyVariable === from) { o.historyVariable = to ?? ''; if (!to) o.historyCoefficient = o.currentCoefficient = 0; }
  if (o.treatmentVariable === from) { o.treatmentVariable = to ?? ''; if (!to) o.treatmentCoefficient = 0; }
  return s;
}
