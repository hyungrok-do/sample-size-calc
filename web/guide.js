import { guideData } from './guide-data.js?v=dev';
import { GuideState } from './guide-state.js?v=dev';
import { prepareSelection } from './selection.js?v=dev';

const el = (tag, className = '', text) => {
  const node = document.createElement(tag);
  node.className = className;
  if (text !== undefined) node.textContent = text;
  return node;
};

export function mountGuide(root, { schemas, onChange, onApply }) {
  const state = new GuideState(guideData);
  const byId = (id) => root.querySelector(`#${id}`);
  const host = byId('guide-questions');
  const requirements = byId('guide-requirements');
  const result = byId('guide-recommendation-content');
  const groups = new Map();
  let revision = 0;
  let applied = null;
  let applyButton = null;
  let status = null;

  function chip(name, value, text, type) {
    const label = el('label', 'guide-chip');
    const input = el('input');
    input.type = type; input.name = name; input.value = value;
    label.append(input, el('span', '', text));
    return label;
  }
  for (const feature of [...guideData.features, { id: 'none', label: 'None beyond choices above' }, { id: 'unsure', label: 'I’m unsure' }]) {
    byId('guide-requirement-chips').append(chip('guide-extra', feature.id, feature.label, 'checkbox'));
  }

  function group(node) {
    if (!groups.has(node.node_id)) {
      const fieldset = el('fieldset');
      const legend = el('legend', '', node.prompt);
      legend.append(el('span', 'guide-hint', 'Choose one'));
      const chips = el('div', 'guide-chips');
      node.choices.forEach((choice, i) => chips.append(chip(node.node_id, String(i), choice.label, 'radio')));
      chips.append(chip(node.node_id, 'unsure', 'I’m unsure', 'radio'));
      fieldset.append(legend, chips);
      groups.set(node.node_id, fieldset);
    }
    const fieldset = groups.get(node.node_id);
    fieldset.querySelectorAll('input').forEach((input) => { input.checked = state.answers.get(node.node_id) === input.value; });
    return fieldset;
  }

  function setting(proc, key, value) {
    const p = proc.params.find((p) => p.id === key);
    return [p?.label ?? key, p?.options?.find((o) => o.value === value)?.label ?? value];
  }

  function message(title, detail) {
    result.append(el('p', 'guide-message', title));
    if (detail) result.append(el('p', 'muted', detail));
  }

  function analysisName(id, opts, fallback) {
    if (id === 'prop.two') return {
      fisher: 'Fisher exact test for two proportions',
      zpool: 'Pooled score test for two proportions',
      zunpool: 'Wald test for two proportions',
      zpoolcc: 'Continuity-corrected test for two proportions',
    }[opts.test] ?? fallback;
    if (id === 'means.one.t' && opts.test === 'z') return 'One-sample z test with known variance';
    if (id === 'means.two.t' && opts.test === 'z') return 'Two-sample z test with known variance';
    if (id === 'sim.surv.nph') return {
      logrank: 'Logrank under nonproportional hazards',
      fh: 'Fleming–Harrington weighted logrank',
      maxcombo: 'MaxCombo under nonproportional hazards',
      rmst: 'Restricted mean survival time (RMST)',
    }[opts.test] ?? fallback;
    return fallback;
  }

  function present(snapshot) {
    const { node, path, selection } = snapshot;
    result.replaceChildren();
    applyButton = null;
    status = null;
    if (node.type === 'question') {
      message(snapshot.unsure ? 'The analysis is not yet determined.' : 'Complete the study choices to see a recommendation.',
        snapshot.unsure ? 'This choice changes the analysis. Clarify the study objective or design before continuing.' : node.prompt);
      return;
    }
    if (node.type === 'no_matching_calculator') {
      message('No matching calculator for this selection.', node.reason);
      return;
    }
    const entry = guideData.catalog[node.procedure_id];
    const proc = schemas.get(node.procedure_id);
    if (!entry || !proc) { message('This procedure is unavailable in the current engine.'); return; }
    if (!state.extra.size || state.extra.has('unsure')) {
      message('Confirm the additional study requirements.', 'Choose all that apply, or “None beyond choices above”.');
      return;
    }
    if (snapshot.blockers.length) {
      message('No recommendation for the current combination.', `${entry.label} does not cover these additional requirements: ${snapshot.blockers.map((f) => f.label.toLowerCase()).join(', ')}. Revise the study choices to look for a matching design.`);
      return;
    }
    try { prepareSelection(proc, selection); }
    catch (error) { message('The recommendation needs review.', error.message); return; }
    const opts = selection.inputs;
    result.append(el('h3', 'guide-analysis-name', analysisName(proc.id, opts, entry.label)));
    if (Object.keys(opts).length) {
      const settings = el('dl', 'guide-settings');
      for (const [key, value] of Object.entries(opts)) {
        const [label, name] = setting(proc, key, value);
        const row = el('div');
        row.append(el('dt', '', label), el('dd', '', name));
        settings.append(row);
      }
      result.append(settings);
    }
    result.append(el('h4', '', 'Why this fits'), el('p', '', path.join(' · ')));
    result.append(el('h4', '', 'Scope and assumptions'), el('p', '', entry.scope_and_checks));
    const details = el('details');
    details.append(el('summary', '', 'Inputs needed for this calculation'), el('p', '', entry.inputs_to_request));
    result.append(details);
    applyButton = el('button', 'primary', 'Use this analysis');
    applyButton.type = 'button';
    applyButton.dataset.procedure = proc.id;
    status = el('p', 'muted small', 'Your analysis choices will be applied to the calculator below. Review its example values before calculating.');
    status.setAttribute('role', 'status');
    result.append(applyButton, status);
    applyButton.addEventListener('click', async () => {
      if (applyButton.getAttribute('aria-disabled') === 'true') return;
      const started = revision;
      applyButton.setAttribute('aria-disabled', 'true');
      try {
        const opened = await onApply(selection);
        if (started !== revision) return;
        if (opened) {
          applied = selection;
          applyButton.textContent = 'Reapply this analysis';
          status.textContent = 'Analysis applied below. Review the study inputs and choose what to calculate.';
        }
      } catch (error) {
        if (started === revision) status.textContent = error.message;
      } finally {
        if (started === revision) applyButton.setAttribute('aria-disabled', 'false');
      }
    });
  }

  function render(announce) {
    const snapshot = state.snapshot();
    const desired = snapshot.questions.map(group);
    // Preserve the currently focused chip; replacing the whole form loses keyboard focus.
    for (const child of Array.from(host.children)) if (!desired.includes(child)) child.remove();
    desired.forEach((child, i) => { if (host.children[i] !== child) host.insertBefore(child, host.children[i] ?? null); });
    requirements.hidden = snapshot.node.type !== 'candidate';
    root.querySelectorAll('input[name="guide-extra"]').forEach((input) => { input.checked = state.extra.has(input.value); });
    present(snapshot);
    if (announce) byId('guide-announcement').textContent = `Recommended analysis updated. ${result.firstElementChild?.textContent ?? ''}`;
  }

  function changed() {
    revision++;
    applied = null;
    onChange();
    render(true);
  }
  root.addEventListener('change', (event) => {
    const input = event.target;
    if (!(input instanceof HTMLInputElement)) return;
    if (input.name === 'guide-extra') state.toggle(input.value, input.checked);
    else state.choose(input.name, input.value);
    changed();
  });
  byId('guide-reset').addEventListener('click', () => {
    state.reset(); changed();
    host.querySelector('input')?.focus();
  });
  render(false);
  return {
    clearApplied() {
      applied = null;
      if (status) status.textContent = 'The calculator was changed outside the guide. Reapply this analysis to restore these settings.';
    },
    matches(proc, lookup) {
      return applied?.procedure === proc?.id && Object.entries(applied.inputs).every(([key, value]) => lookup(key) === value);
    },
  };
}
