// A single choice at each branch; independent additional requirements can be selected together.
export class GuideState {
  constructor(data) {
    this.data = data;
    this.answers = new Map();
    this.extra = new Set();
    this.nodes = new Map();
    const register = (node) => {
      this.nodes.set(node.node_id, node);
      if (node.type === 'question') node.choices.forEach((c) => register(c.next));
    };
    register(data.root);
  }

  choose(id, answer) {
    const node = this.nodes.get(id);
    if (node?.type !== 'question' || (answer !== 'unsure' && !node.choices[Number(answer)])) return;
    const clear = (n) => {
      this.answers.delete(n.node_id);
      if (n.type === 'question') n.choices.forEach((c) => clear(c.next));
    };
    node.choices.forEach((c) => clear(c.next));
    this.answers.set(id, String(answer));
  }

  toggle(id, checked) {
    if (!['none', 'unsure', ...this.data.features.map((f) => f.id)].includes(id)) return;
    if (checked) {
      if (id === 'none' || id === 'unsure') this.extra.clear();
      else { this.extra.delete('none'); this.extra.delete('unsure'); }
      this.extra.add(id);
    } else this.extra.delete(id);
  }

  reset() { this.answers.clear(); this.extra.clear(); }

  // The additional requirements chosen; "none" and "unsure" ask for nothing beyond the questions.
  wanted() { return this.data.features.filter((f) => this.extra.has(f.id)); }

  // The analyses below a node, in the guide's order: one per procedure, with the presets every route to it
  // agrees on (a choice the routes differ on is left to the calculator, for the next questions to settle).
  candidatesUnder(node) {
    const byProcedure = new Map();
    const walk = (n) => {
      if (n.type === 'question') n.choices.forEach((c) => walk(c.next));
      if (n.type !== 'candidate') return;
      const seen = byProcedure.get(n.procedure_id);
      if (!seen) byProcedure.set(n.procedure_id, { procedure: n.procedure_id, inputs: { ...n.proposed_option_values } });
      else for (const [key, value] of Object.entries(seen.inputs)) if (n.proposed_option_values[key] !== value) delete seen.inputs[key];
    };
    walk(node);
    return [...byProcedure.values()];
  }

  snapshot() {
    const questions = [], path = [];
    const wanted = this.wanted();
    const fits = (id) => wanted.every((f) => f.allow.includes(id));
    let node = this.data.root;
    while (node.type === 'question') {
      questions.push(node);
      const answer = this.answers.get(node.node_id);
      if (answer === undefined || answer === 'unsure') {
        // Before the path ends, every analysis still possible: shown once anything narrows the whole catalog.
        const candidates = path.length || wanted.length ? this.candidatesUnder(node).filter((c) => fits(c.procedure)) : [];
        return { node, questions, path, wanted, candidates, unsure: answer === 'unsure', selection: null };
      }
      const choice = node.choices[Number(answer)];
      path.push(choice.label);
      node = choice.next;
    }
    const blockers = node.type === 'candidate' ? wanted.filter((f) => !fits(node.procedure_id)) : [];
    const ready = node.type === 'candidate' && !blockers.length;
    return { node, questions, path, wanted, candidates: [], blockers, selection: ready ? {
      procedure: node.procedure_id, inputs: { ...node.proposed_option_values },
    } : null };
  }
}
