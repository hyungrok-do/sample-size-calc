// Plain-language summaries of the design and of the question, for the canvas heading and the power panel.

export function studySummary(graph) {
  if (!graph.nodes.length) return 'No variables yet · choose a starting point or add a variable';
  const s = graph.study, plural = (k, word) => `${k} ${word}${k === 1 ? '' : 's'}`;
  return `${plural(graph.nodes.length, 'variable')} · ${plural(graph.edges.length, 'arrow')} · ${s ? `${s.times.length} scheduled visit${s.times.length === 1 ? '' : 's'} · sample size counts ${s.clusterSizes.length ? 'clusters' : 'subjects'}` : 'sample size counts observations'}`;
}

export function questionSummary(graph, target) {
  if (!target) return 'Choose a compared variable and an outcome to define your question.';
  const name = (id) => graph.nodes.find((n) => n.id === id)?.label || id;
  const outcome = name(target.outcome), exposure = name(target.exposure);
  if (target.kind === 'coefficient') return `Estimate the ${target.measure === 'exponentiated_coefficient' ? 'exponentiated ' : ''}coefficient of ${exposure} in the model for ${outcome}.`;
  if (target.kind === 'joint') return target.measure === 'interaction' ? `Test whether the effect of ${exposure} on ${outcome} depends on another variable (its interaction terms).`
    : `Test ${target.measure === 'nonlinearity' ? 'nonlinearity' : 'the overall effect'} of ${exposure} in the model for ${outcome}.`;
  const measures = { mean_difference: 'mean difference', risk_difference: 'risk difference', risk_ratio: 'risk ratio', odds_ratio: 'odds ratio', survival_difference: 'survival probability difference', rmst_difference: 'restricted mean survival time difference' };
  let response = outcome, when = '';
  if (graph.study) {
    const r = target.repeated ?? { contrast: 'final' }, times = graph.study.times;
    if (r.contrast === 'change') response = `the change in ${outcome} from time ${times[0]} to ${times.at(-1)}`;
    else if (r.contrast === 'slope') response = `the linear trend in ${outcome} over time`;
    else when = r.contrast === 'average' ? ' averaged across scheduled visits'
      : r.contrast === 'custom' ? ' using your custom visit weights' : ` at time ${times[r.visit ?? times.length - 1]}`;
  } else if (target.horizon != null) when = ` ${target.measure === 'rmst_difference' ? 'through' : 'at'} time ${target.horizon}`;
  return `Compare ${response}${when} if ${exposure} were set to ${target.active} versus ${target.reference}, using the ${measures[target.measure] ?? target.measure.replaceAll('_', ' ')}.`;
}
