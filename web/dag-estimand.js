import { MEASURES, availableMeasures, validateEstimand, describeEstimand, copyEstimand } from './estimand-model.js?v=dev';

function el(tag, attrs = {}, ...children) {
  const node = document.createElement(tag);
  for (const [key, value] of Object.entries(attrs)) {
    if (value == null || value === false) continue;
    if (key.startsWith('on')) node.addEventListener(key.slice(2), value);
    else node.setAttribute(key, value === true ? '' : value);
  }
  node.append(...children.flat()); return node;
}
function field(title, input, help = '') {
  return el('div', { class: 'dag-field' }, el('label', { for: input.id }, title), input,
    ...(help ? [el('p', { class: 'small muted' }, help)] : []));
}
function select(id, options, value) {
  const input = el('select', { id }, options.map(([key, label]) => el('option', { value: key }, label)));
  input.value = value; return input;
}
const numeric = (id, value) => el('input', { id, type: 'number', step: 'any', value: value ?? '' });
const number = (input, optional = false) => {
  if (input.value.trim() === '') return optional ? null : NaN;
  return Number(input.value);
};

export function mountEstimand(root, { getGraph, onSave, onDraftChange = () => {} }) {
  function refresh(saved) {
    const graph = getGraph();
    const defaultValue = { version: 1, kind: 'marginal', population: 'Population represented by the data generating process',
      definition: '', exposure: '', outcome: '', measure: '', reference: 0, active: 1,
      horizon: null, nullValue: 0, alternative: null };
    const value = saved ?? defaultValue;
    const badge = el('span', { class: 'dag-badge' }, saved ? 'Saved target' : 'Define the target');
    const kind = select('estimand-kind', [['marginal', 'Marginal total-effect contrast'], ['coefficient', 'Conditional model coefficient'], ['custom', 'Custom definition']], value.kind);
    const choices = [['', 'Choose a variable'], ...graph.nodes.map((n) => [n.id, `${n.id} · ${n.label}`])];
    const outcome = select('estimand-outcome', choices, value.outcome);
    const exposure = select('estimand-exposure', choices, value.exposure);
    const measure = select('estimand-measure', [], value.measure);
    const population = el('textarea', { id: 'estimand-population', rows: 2, maxlength: 1000 }, value.population);
    const definition = el('textarea', { id: 'estimand-definition', rows: 3, maxlength: 3000 }, value.definition);
    const reference = numeric('estimand-reference', value.reference), active = numeric('estimand-active', value.active);
    const horizon = numeric('estimand-horizon', value.horizon), nullValue = numeric('estimand-null', value.nullValue);
    const alternative = numeric('estimand-alternative', value.alternative);
    const variableFields = el('div', { class: 'dag-estimand-grid' }, field('Comparison / intervention variable', exposure), field('Outcome variable', outcome));
    const levels = el('div', { class: 'dag-estimand-grid' }, field('Reference level', reference), field('Active level', active));
    const horizonField = field('Time horizon', horizon, 'Use the same time units as the survival variable.');
    const customField = field('Estimand definition', definition, 'Specify the population quantity, comparison, and any time horizon.');
    const preview = el('div', { id: 'estimand-preview', class: 'dag-estimand-preview', 'aria-live': 'polite' });
    const form = el('form', { id: 'estimand-form', novalidate: true },
      el('div', { class: 'dag-estimand-grid' }, field('Target population', population,
        'The executable population is the graph’s distribution. Describe that population here; text does not subset or reweight generated data.'), field('Estimand type', kind)),
      variableFields, field('Effect measure', measure), levels, horizonField, customField,
      el('div', { class: 'dag-estimand-grid' }, field('Null value', nullValue), field('Planning alternative (optional)', alternative,
        'Sets the effect used in power simulations. Leave blank to use the graph’s implied effect.')),
      el('button', { type: 'submit', class: 'primary' }, 'Save estimand'));
    function candidate() {
      return { version: 1, kind: kind.value, population: population.value.trim(), definition: definition.value.trim(),
        exposure: kind.value === 'custom' ? '' : exposure.value, outcome: kind.value === 'custom' ? '' : outcome.value,
        measure: measure.value, reference: number(reference), active: number(active),
        horizon: MEASURES[measure.value]?.horizon ? number(horizon) : null,
        nullValue: number(nullValue), alternative: number(alternative, true) };
    }
    function showPreview() {
      try {
        const target = candidate(), text = describeEstimand(graph, target);
        preview.replaceChildren(el('p', { class: 'eyebrow' }, 'Target definition'), el('p', { class: 'dag-estimand-formula' }, text.formula),
          el('p', {}, `Population: ${text.population}`), el('p', {}, text.hypothesis), el('p', {}, text.alternative), el('p', { class: 'small muted' }, text.note),
          el('p', { class: 'small muted' }, 'Save this target to select its compatible analysis and run repeated simulations in the power section below.'));
      } catch (error) { preview.replaceChildren(el('p', { class: 'muted' }, error.message)); }
    }
    function updateMeasures(resetNull = true) {
      const wanted = measure.value || value.measure;
      const options = availableMeasures(graph, kind.value, outcome.value);
      measure.replaceChildren(...options.map((key) => el('option', { value: key }, MEASURES[key].label)));
      measure.value = options.includes(wanted) ? wanted : options[0] ?? '';
      if (resetNull) { nullValue.value = MEASURES[measure.value]?.nullValue ?? 0; alternative.value = ''; }
      visibility();
    }
    function visibility() {
      variableFields.hidden = kind.value === 'custom';
      levels.hidden = kind.value !== 'marginal';
      customField.hidden = kind.value !== 'custom';
      horizonField.hidden = !MEASURES[measure.value]?.horizon;
      showPreview();
    }
    kind.addEventListener('change', () => updateMeasures());
    outcome.addEventListener('change', () => updateMeasures());
    measure.addEventListener('change', () => { nullValue.value = MEASURES[measure.value]?.nullValue ?? 0; alternative.value = ''; visibility(); });
    const edited = () => { badge.textContent = 'Unsaved changes'; onDraftChange(true); showPreview(); };
    form.addEventListener('input', edited); form.addEventListener('change', edited);
    form.addEventListener('submit', (event) => {
      event.preventDefault();
      try { const target = candidate(); validateEstimand(graph, target); onSave(copyEstimand(target)); }
      catch (error) { preview.replaceChildren(el('p', { class: 'error', role: 'alert' }, error.message)); }
    });
    updateMeasures(false);
    root.replaceChildren(el('div', { class: 'dag-estimand-heading' }, el('h3', {}, 'What is the estimand?'),
      badge),
      el('p', {}, 'Specify the population quantity your power analysis should detect. Variable roles are chosen for this target; the DAG variables remain generic.'),
      el('div', { class: 'dag-estimand-layout' }, form, preview));
  }
  return { refresh };
}
