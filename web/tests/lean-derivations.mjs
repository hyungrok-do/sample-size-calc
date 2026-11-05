import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { examples } from '../lean-derivation-examples.js';

const close = (x, y) => assert.ok(Math.abs(x - y) < 2e-13, `${x} differs from ${y}`);
const contains = (estimate, p) => {
  assert.equal(estimate.method, 'simulation');
  assert.ok(estimate.lower <= p && p <= estimate.upper,
    `${p} is outside [${estimate.lower}, ${estimate.upper}]`);
};
const checkedTheorem = (report, name) => {
  const theorem = report.results.find(result => result.name === name);
  assert.ok(theorem, `Missing theorem ${name}`);
  assert.equal(theorem.kind, 'theorem');
  assert.ok(!theorem.axioms.some(a => /sorry|trustCompiler|ofReduce/.test(a)));
  return theorem;
};

const reference = JSON.parse(await readFile(new URL('./fixtures/lean-derivations.json', import.meta.url), 'utf8'));
const compareExactRows = (search, name) => {
  for (const { evaluation } of search.rows) {
    const expected = reference[name].find(x => x.n === evaluation.design);
    assert.ok(expected);
    assert.equal(evaluation.typeI.exact, expected.typeI);
    assert.equal(evaluation.power.exact, expected.power);
  }
};

export async function checkDerivations({ compile, report, page, reject }) {
  await compile(examples['bernoulli-observations']);
  let result = await report('bernoulli-observations.json');
  let [search, simulation] = result.calculations;
  let answer = search.rows[search.selected].evaluation;
  assert.equal(answer.design, 25);
  assert.equal(search.earlierUnresolved.length, 0);
  assert.ok(search.rows.every(row => row.evaluation.power.method === 'exact'));
  compareExactRows(search, 'bernoulli-observations');
  close(answer.typeI.value, .0333999446114952529);
  close(answer.power.value, .8065115579119348040);
  contains(simulation.evaluation.typeI, answer.typeI.value);
  contains(simulation.evaluation.power, answer.power.value);
  checkedTheorem(result, 'ObservedTrial.distribution_from_observations');
  checkedTheorem(result, 'ObservedTrial.count_recurrence');
  checkedTheorem(result, 'ObservedTrial.target_from_participant_model');
  const minimum = checkedTheorem(result, 'ObservedTrial.minimum_sample_size');
  assert.match(minimum.statement, /meetsTarget 25/);
  assert.ok(minimum.dependencies.includes('ObservedTrial.compiled_iff_observation_target'));
  console.log('Observation model -> statistic distribution -> exact power -> minimum n: all kernel checked.');

  await compile(examples['paired-categorical']);
  result = await report('paired-categorical.json');
  [search, simulation] = result.calculations;
  answer = search.rows[search.selected].evaluation;
  assert.equal(answer.design, 23);
  compareExactRows(search, 'paired-categorical');
  close(answer.typeI.value, .0424170238910503507);
  close(answer.power.value, .8165689646545179111);
  contains(simulation.evaluation.typeI, answer.typeI.value);
  contains(simulation.evaluation.power, answer.power.value);
  assert.match(simulation.evaluation.power.reason, /derivation exceeded.*original observation program/);
  checkedTheorem(result, 'PairedCategorical.score_distribution_from_pairs');
  checkedTheorem(result, 'PairedCategorical.calibration_from_certificate');
  checkedTheorem(result, 'PairedCategorical.selected_cutoff');
  checkedTheorem(result, 'PairedCategorical.target_from_joint_observations');
  for (const name of ['checkedNull', 'checkedAlternative']) {
    const certificate = result.results.find(x => x.name === `PairedCategorical.${name}`);
    assert.equal(certificate.kind, 'definition');
    assert.ok(certificate.dependencies.some(x => x.includes('_finiteCheckpoint')));
    assert.ok(!certificate.axioms.some(x => /sorry|trustCompiler|ofReduce/.test(x)));
  }
  console.log('Paired observations: exact cutoff and power bounds kernel checked.');

  await compile(examples['dependent-observations']);
  result = await report('dependent-observations.json');
  const [direct, simulated] = result.calculations.map(x => x.evaluation);
  const expectedMarkov = reference['dependent-observations'].find(x => x.n === 12);
  assert.equal(direct.typeI.exact, expectedMarkov.typeI);
  assert.equal(direct.power.exact, expectedMarkov.power);
  close(direct.typeI.value, .031284633835);
  close(direct.power.value, .21798912);
  contains(simulated.typeI, direct.typeI.value);
  contains(simulated.power, direct.power.value);
  checkedTheorem(result, 'DependentObservations.derived_power');
  checkedTheorem(result, 'DependentObservations.dependent_probability');
  console.log('Dependent paired and Markov observations: exact derivations and original-model simulations agree.');

  const generic = `import SSC
open SSC.Finite SSC.Derivation
open SSC.Calculation (emitEvaluation)
set_option maxRecDepth 10000
def coin : Program Bool := .draw (.bernoulli (1/3) (by decide +kernel))
def step (s : Nat) (x : Bool) := s + if x then 1 else 0
def count := coin.iidFold step 0 25
theorem compact_support : (compileWithin count).map (fun c => c.distribution.atoms.length) = some 26 := by
  decide +kernel
theorem no_partial_certificate : (compileWithin count 1).isNone = true := by decide +kernel
def rule : Rejection Nat := .both (.test (fun n => n ≥ 5)) (.test (fun n => n ≤ 7))
def spec : Specification Unit Nat where
  null := fun _ => count
  alternative := fun _ => count
  reject := fun _ => rule
#eval emitEvaluation (spec.evaluate () { replicates := 2000, seed := 41004 } (budget := 1))
#eval emitEvaluation (spec.evaluate () { route := .simulation, replicates := 2000, seed := 41004 })
def small := coin.iidFold step 0 3
def staged : Certificate small := by finite_derive
theorem staged_observation_probability : small.probability (fun k => k == 3) = 1/27 := by
  rw [← staged.probability_eq]
  decide +kernel
theorem exact_observation_probability : small.probability (fun k => k == 3) = 1/27 := by
  apply derivedProbability_sound small (fun k => k == 3) 1000000
  decide +kernel
`;
  await compile(generic);
  result = await report('derivation-generic.json');
  const [fallback, forced] = result.calculations.map(x => x.evaluation);
  for (const which of ['typeI', 'power']) {
    assert.equal(fallback[which].successes, forced[which].successes);
    assert.equal(fallback[which].seed, forced[which].seed);
    assert.equal(fallback[which].lower, forced[which].lower);
  }
  checkedTheorem(result, 'compact_support');
  checkedTheorem(result, 'exact_observation_probability');
  checkedTheorem(result, 'staged_observation_probability');
  console.log('State compaction, bounded certificate soundness, joint rules and identical fallback draws checked.');

  await compile(await readFile(new URL('./fixtures/lean-performance.lean', import.meta.url), 'utf8'));
  result = await report('performance-semantics.json');
  checkedTheorem(result, 'family_correct');
  checkedTheorem(result, 'memoized_law_is_original');
  checkedTheorem(result, 'quotientCertificate_correct');
  checkedTheorem(result, 'quotientBindCertificate_correct');
  checkedTheorem(result, 'arbitrary_order_preserves_distinct_outcomes');
  assert.equal(result.results.find(d => d.name === 'contextualOrdering').parameters.length, 1);
  assert.ok(result.results.find(d => d.name === 'quotientCertificate').axioms.includes('Quot.sound'));
  const messages = [result.output, ...result.diagnostics.map(d => d.data)].join('\n');
  assert.match(messages, /SSC_FINITE_CHECKPOINTS 9(?:\s|$)/);
  assert.match(messages, /SSC_ORDERED_CHECKPOINTS 8(?:\s|$)/);
  assert.match(messages, /SSC_PERFORMANCE_SEMANTICS_OK/);
  assert.ok(result.performance.totalMs >= 0);
  for (const key of ['headerMs', 'parseMs', 'kernelReplayMs', 'axiomAuditMs',
    'declarationRenderMs', 'traceRenderMs', 'commandMs']) assert.ok(result.performance[key] >= 0, key);
  assert.ok(result.performance.commands.every(command => command.line > 0 && command.ms >= 0));
  console.log('Shared finite families and cached sampling preserve certificates, draws, budgets and failures.');

  await reject(`import SSC
open SSC.Derivation
def tooWide : (n : Fin 257) → Certificate (Program.pure n.val) := by finite_derive`, /256 finite-family cases/);

  // Change an observation assumption; the previous mathematical claim must fail.
  const changed = examples['bernoulli-observations'].replaceAll('(3/10)', '(2/10)');
  await reject(changed.slice(0, changed.indexOf('#eval emitSearch')) + '\nend ObservedTrial',
    /false|kernel|unsolved|failed/i);
  await reject(`import SSC
open SSC.Finite SSC.Derivation
def falseCertificate : Certificate (Program.pure (0 : Nat)) where
  distribution := Distribution.pure 1
  correct := by intro f; rfl`, /failed|definitionally|unsolved/i);
  await reject(`import SSC
open SSC.Derivation
local instance (priority := 2000) : Lean.ToExpr Nat where
  toTypeExpr := Lean.mkConst \`\`Nat
  toExpr _ := Lean.mkNatLit 0
def invalidProposal : Certificate (Program.pure (1 : Nat)) := by finite_derive`,
    /kernel|type mismatch|false|failed/i);
  await reject(`import SSC
open SSC.Derivation
local instance (priority := 2000) : Lean.ToExpr Nat where
  toTypeExpr := Lean.mkConst \`\`Nat
  toExpr _ := Lean.mkNatLit 0
def invalidBranchProposal : Certificate
  ((Program.pure (0 : Nat)).bind (fun n => Program.pure (n + 1))) := by finite_derive`,
    /kernel|type mismatch|false|failed/i);
  await reject(`import SSC
open SSC.Finite
def invalid : Distribution Nat where
  atoms := [(0, 2)]
  nonnegative := by intro x hx; simp only [List.mem_singleton] at hx; subst x; decide +kernel
  normalized := by decide +kernel`, /false|kernel|failed/i);
  await reject(`import SSC
open SSC.Derivation
open SSC.Calculation (emitEvaluation)
def p : Program Nat := .pure 1
def spec : Specification Unit Nat where
  null := fun _ => p
  alternative := fun _ => p
  reject := fun _ => .test (fun _ => true)
#eval emitEvaluation (spec.evaluate () { route := .exact } (budget := 0))`,
    /derivation exceeded.*budget/i);
  assert.equal(await page.locator('#lean-error').isVisible(), true);
  console.log('Changed observation models, false certificates, invalid mass and forced unavailable exact routes rejected.');
}
