// Shared by the schema-driven form and the guide. Never silently replace an explicit guide choice.
export function visible(rule, lookup) {
  if (!rule || !rule.length) return true;
  return rule.some((clauses) => clauses.every((c) => (c.values ?? []).includes(lookup(c.key)) !== Boolean(c.negate)));
}

export function prepareSelection(proc, selection) {
  if (!proc || proc.id !== selection.procedure) throw new Error('This calculator is unavailable.');
  const solve = selection.solve_for ?? proc.solve_for[0].value;
  if (!proc.solve_for.some((s) => s.value === solve)) throw new Error('This calculation goal is unavailable.');
  const values = Object.fromEntries(proc.params.map((p) => [p.id, p.default ?? '']));
  // Guide presets are analysis choices, never guessed numerical study assumptions.
  for (const [key, value] of Object.entries(selection.inputs ?? {})) {
    const p = proc.params.find((p) => p.id === key);
    if (p?.type !== 'select' || !p.options.some((o) => o.value === value)) {
      // the name shown comes from the calculator, never from the address, which anyone can write
      throw new Error(p?.type === 'select'
        ? `The recommended setting “${p.label}” is unavailable. Please review the study choices.`
        : 'The address names a setting this calculator does not have. Please review the study choices.');
    }
    values[key] = value;
  }
  const lookup = (key) => key === '$solve' ? solve : values[key];
  for (const [key, value] of Object.entries(selection.inputs ?? {})) {
    const p = proc.params.find((p) => p.id === key);
    if (!visible(p.visible, lookup) || !visible(p.options.find((o) => o.value === value).visible, lookup)) {
      throw new Error(`The recommended setting “${p.label}” does not apply to this combination.`);
    }
  }
  return { solve, values };
}

export function selectionHash(selection) {
  const options = new URLSearchParams(selection.inputs ?? {}).toString();
  return `#${encodeURIComponent(selection.procedure)}${options ? `?${options}` : ''}`;
}

export function parseSelectionHash(hash) {
  const [id, options = ''] = hash.replace(/^#/, '').split('?');
  try {
    return { procedure: decodeURIComponent(id), inputs: Object.fromEntries(new URLSearchParams(options)) };
  } catch { return { procedure: '', inputs: {} }; }
}
