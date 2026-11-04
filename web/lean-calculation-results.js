// Calculation output is executable output, never a proof certificate. Keep the
// actual source and the independent kernel declaration report alongside it.
const prefix = 'SSC_CALCULATION ';
const probability = x => typeof x === 'number' && Number.isFinite(x) && x >= 0 && x <= 1;
const estimate = e => e && ['exact', 'numerical', 'simulation'].includes(e.method) && probability(e.value)
  && (e.method !== 'simulation' || (probability(e.lower) && probability(e.upper)
    && e.lower <= e.value && e.value <= e.upper && Number.isSafeInteger(e.replicates) && e.replicates > 0
    && probability(e.confidence) && e.confidence > 0
    && Number.isSafeInteger(e.successes) && e.successes >= 0 && e.successes <= e.replicates
    && typeof e.seed === 'string'));
const evaluation = e => e && Object.hasOwn(e, 'design') && estimate(e.typeI) && estimate(e.power);

export function extractCalculations(report) {
  const calculations = [];
  const strip = text => text.split('\n').filter(line => {
    if (!line.startsWith(prefix)) return true;
    try {
      const r = JSON.parse(line.slice(prefix.length));
      if (r.schema !== 'ssc.calculation.v1') return true;
      if (r.kind === 'evaluation') { if (!evaluation(r.evaluation)) return true; }
      else if (r.kind === 'search') {
        if (!Array.isArray(r.rows) || r.rows.length > 10000 || !r.rows.every(row => evaluation(row.evaluation)
          && ['meets', 'fails', 'unresolved'].includes(row.decision))) return true;
        if (r.selected !== null && (!Number.isSafeInteger(r.selected) || r.selected < 0
          || r.selected >= r.rows.length || r.rows[r.selected].decision !== 'meets')) return true;
        if (!Array.isArray(r.earlierUnresolved) || !r.earlierUnresolved.every(i => Number.isSafeInteger(i)
          && i >= 0 && i < r.rows.length && r.rows[i].decision === 'unresolved')) return true;
      } else return true;
      calculations.push(r);
      return false;
    } catch { return true; }
  }).join('\n');
  const output = strip(report.output || '');
  const diagnostics = report.diagnostics.map(d => d.severity === 'information' ? { ...d, data: strip(d.data) } : d)
    .filter(d => d.data.trim());
  return { ...report, output, diagnostics, calculations };
}

const element = (tag, text, className) => {
  const node = document.createElement(tag);
  if (text != null) node.textContent = text;
  if (className) node.className = className;
  return node;
};
const designText = value => typeof value === 'object' ? JSON.stringify(value) : String(value);
const probText = e => e.method === 'simulation'
  ? `${e.value.toFixed(4)}\n[${e.lower.toFixed(4)}, ${e.upper.toFixed(4)}]`
  : e.value.toFixed(6);

export function renderCalculations(container, reports) {
  container.replaceChildren(); container.hidden = !reports.length;
  for (const r of reports) {
    const card = element('section', null, 'calculation');
    card.append(element('h4', r.kind === 'search' ? 'Design search' : 'Power evaluation'));
    const rows = r.kind === 'search' ? r.rows : [{ evaluation: r.evaluation, decision: 'evaluated' }];
    if (r.kind === 'search') {
      card.append(element('p', `Target power ≥ ${r.targetPower}${r.maxTypeI ? ` · Type-I error ≤ ${r.maxTypeI}` : ' · Type-I error reported separately'}`, 'small'));
      if (r.selected !== null) {
        card.append(element('p', `Selected design: ${designText(rows[r.selected].evaluation.design)}`, 'calculation-answer'));
        card.append(element('p', r.earlierUnresolved.length
          ? `${r.earlierUnresolved.length} earlier candidates remain unresolved. A minimum has not been established.`
          : 'First candidate meeting the criteria in the supplied order, at the displayed evaluation accuracy.', 'small'));
      } else card.append(element('p', r.earlierUnresolved.length
        ? 'No selection: some candidates remain unresolved. Increase simulation precision or refine the supplied candidate set.'
        : 'No supplied candidate met the criteria.', 'calculation-answer'));
      card.append(element('p', 'Search scope is the supplied candidate set. Numerical decisions and simulation intervals are not kernel-checked minimality proofs.', 'small muted'));
    }
    const details = element('details'); details.open = r.kind !== 'search';
    details.append(element('summary', `${rows.length} candidate evaluation${rows.length === 1 ? '' : 's'}`));
    const wrap = element('div', null, 'trace-scroll'), table = element('table');
    const head = element('tr');
    const columns = ['Design', 'Type-I error', 'Power', 'Method'];
    if (r.kind === 'search') columns.push('Decision');
    for (const title of columns) head.append(element('th', title));
    const thead = element('thead'); thead.append(head); table.append(thead);
    const body = element('tbody');
    for (const { evaluation: e, decision } of rows) {
      const tr = element('tr');
      const method = e.typeI.method === e.power.method ? e.power.method : `${e.typeI.method} / ${e.power.method}`;
      const cells = [designText(e.design), probText(e.typeI), probText(e.power), method];
      if (r.kind === 'search') cells.push(decision);
      for (const value of cells) tr.append(element('td', value));
      body.append(tr);
    }
    table.append(body); wrap.append(table); details.append(wrap); card.append(details);
    if (rows.some(row => row.evaluation.typeI.method !== row.evaluation.power.method))
      card.append(element('p', 'When methods differ, they are shown as null / alternative.', 'small muted'));
    const simulated = rows.flatMap(row => [row.evaluation.typeI, row.evaluation.power]).filter(e => e.method === 'simulation');
    if (simulated.length) {
      const confidence = r.kind === 'search' ? r.simultaneousConfidence : 1 - 2 * (1 - simulated[0].confidence);
      card.append(element('p', `Simulation: ${simulated[0].replicates.toLocaleString()} replicates per probability. Brackets show Hoeffding intervals with ${(100 * confidence).toFixed(1)}% simultaneous confidence across the planned Monte Carlo evaluations. Seeds and success counts are in the downloaded report.`, 'small'));
    }
    const notes = [...new Set(rows.flatMap(row => [row.evaluation.typeI.reason, row.evaluation.power.reason]).filter(Boolean))];
    const scope = element('details'); scope.append(element('summary', 'Evaluation method and scope'));
    for (const note of notes) scope.append(element('p', note));
    scope.append(element('p', r.scope || 'Executable output; theorem proofs are reported separately.'));
    card.append(scope); container.append(card);
  }
}
