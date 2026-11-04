import assert from 'node:assert/strict';

export async function checkVerification({ compile, reject, report }) {
  const tactic = `import SSC
open Lean Elab Tactic
elab "bogusProof" : tactic => do
  let goal ← getMainGoal
  goal.assign (mkConst \`\`True.intro)
  replaceMainGoal []
`;
  // Global and scoped options must both be caught by independent kernel replay.
  for (const source of [
    tactic + 'set_option debug.skipKernelTC true\ntheorem impossible : False := by bogusProof',
    tactic + 'set_option debug.skipKernelTC true in\ntheorem impossible : False := by bogusProof',
  ]) await reject(source, /replaying|kernel|type mismatch/i);
  await reject(`import SSC
open SSC.Derivation
set_option debug.skipKernelTC true
local instance (priority := 2000) : Lean.ToExpr Nat where
  toTypeExpr := Lean.mkConst \`\`Nat
  toExpr _ := Lean.mkNatLit 0
def invalidProposal : Certificate (Program.pure (1 : Nat)) := by finite_derive`,
  /kernel|type mismatch|false|failed/i);
  // Correct source remains usable after each rejected module, including
  // inductives whose constructors/recursors are reconstructed during replay.
  await compile(`import SSC
set_option debug.skipKernelTC true
inductive Bit where | zero | one deriving DecidableEq
theorem distinct : Bit.zero ≠ Bit.one := by decide +kernel`);
  const checked = await report('independent-kernel-replay.json');
  assert.ok(checked.results.some(d => d.name === 'distinct' && d.kind === 'theorem'));
  console.log('Independent kernel replay rejects unchecked false proofs and preserves valid inductives.');

  await compile(`import SSC
open SSC.Derivation
open SSC.Calculation (emitEvaluation)
def coin : Program Bool := .draw (.bernoulli (1/2) (by decide +kernel))
def count (n : Nat) := coin.iidFold (fun (s : Nat) x => s + if x then 1 else 0) 0 n
def spec : Specification Nat Nat where
  null := count
  alternative := count
  reject := fun n => .test (fun k => k ≥ n/2)
#eval emitEvaluation (spec.evaluate 1000 { replicates := 100, seed := 1 } (budget := 0))
#eval emitEvaluation (spec.evaluate 1000 { replicates := 100, seed := 1 })
#eval emitEvaluation (spec.evaluate 1000 { route := .simulation, replicates := 100, seed := 1 })
`);
  const results = (await report('large-observation-simulation.json')).calculations;
  assert.equal(results.length, 3);
  // Independent integer SplitMix64 reference: 100 replicates of 1,000 draws
  // per law, rejecting when at least 500 are true. Check every draw is consumed.
  assert.equal(results[0].evaluation.typeI.successes, 55);
  assert.equal(results[0].evaluation.typeI.seed, '1');
  assert.equal(results[0].evaluation.power.successes, 61);
  assert.equal(results[0].evaluation.power.seed, '7357944848429976353');
  for (const key of ['typeI', 'power']) {
    const baseline = results[0].evaluation[key];
    assert.equal(baseline.method, 'simulation');
    for (const result of results.slice(1)) {
      const actual = result.evaluation[key];
      assert.equal(actual.method, 'simulation');
      assert.equal(actual.successes, baseline.successes);
      assert.equal(actual.seed, baseline.seed);
      assert.equal(actual.lower, baseline.lower);
      assert.equal(actual.upper, baseline.upper);
    }
  }
  console.log('1,000 original observations: automatic, exhausted-budget and forced simulation agree.');
  await compile(`import SSC
open SSC.Derivation
def coin : Program Bool := .draw (.bernoulli (1/2) (by decide +kernel))
def large := coin.iidFold (fun (s : Nat) x => s + if x then 1 else 0) 0 10000
def dependent := Program.iterate (fun (s : Nat) => Program.pure (s+1)) 0 10000
#eval (show IO Unit from do
  match (SSC.Probability.sample large.toLaw).run ⟨1, 100000⟩ with
  | .ok (count, _) =>
    if count != 4836 then throw (IO.userError "Wrong 10,000-observation reference count")
  | .error message => throw (IO.userError message)
  match (SSC.Probability.sample dependent.toLaw).run ⟨1, 100000⟩ with
  | .ok (count, _) =>
    if count != 10000 then throw (IO.userError "Dependent iteration skipped observations")
  | .error message => throw (IO.userError message))
`);
  console.log('10,000-step IID and dependent programs execute without browser recursion overflow.');
}
