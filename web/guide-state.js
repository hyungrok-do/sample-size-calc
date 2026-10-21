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

  snapshot() {
    const questions = [], path = [];
    let node = this.data.root;
    while (node.type === 'question') {
      questions.push(node);
      const answer = this.answers.get(node.node_id);
      if (answer === undefined || answer === 'unsure') return { node, questions, path, unsure: answer === 'unsure' };
      const choice = node.choices[Number(answer)];
      path.push(choice.label);
      node = choice.next;
    }
    const blockers = this.data.features.filter((f) => this.extra.has(f.id) && !f.allow.includes(node.procedure_id));
    const ready = node.type === 'candidate' && this.extra.size > 0 && !this.extra.has('unsure') && !blockers.length;
    return { node, questions, path, blockers, selection: ready ? {
      procedure: node.procedure_id, inputs: { ...node.proposed_option_values },
    } : null };
  }
}
