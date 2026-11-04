import assert from 'node:assert/strict';

export async function checkEnumeration({ compile, reject, report }) {
  // Independent hand-enumerated branch masses, with changing intermediate
  // types and a zero-mass atom. Exactly 20 units of enumeration work are needed.
  await compile(`import SSC
open SSC.Probability
theorem exact_computation_remains_reducible :
    (match exactProbability (Law.pure true) (Event.ofBool id) 1 with
     | .ok p => p == 1 | .error _ => false) = true := by decide +kernel
def mixed : Law Nat := do
  let flag ← Law.finite #[(false, 2/3), (true, 1/3), (false, 0)]
  if flag then
    let x ← Law.finite #[(10, 1/4), (20, 3/4)]
    pure (x+1)
  else
    let text ← Law.pure "zero"
    pure text.length
#eval (show IO Unit from do
  match (enumerate mixed).run 20 with
  | .ok (xs, remaining) =>
    if xs != #[(4, 2/3), (11, 1/12), (21, 1/4)] || remaining != 0 then
      throw (IO.userError "Wrong branch weights, order or work consumption")
  | .error _ => throw (IO.userError "Exact enumeration rejected a sufficient budget")
  for budget in [:20] do
    match (enumerate mixed).run budget with
    | .error (.unavailable _) => pure ()
    | _ => throw (IO.userError "Enumeration returned a partial table"))
`);
  console.log('Enumeration kernel reduction, rational branch weights and exact work boundary passed.');

  const source = `import SSC
open SSC.Probability SSC.Calculation
def spec : Specification Nat (List Nat) where
  null := fun n => (Law.pure 1).iid n
  alternative := fun n => (Law.pure 1).iid n
  reject := fun n => Event.ofBool (fun xs => xs.length == n)
`;
  await compile(source + `
#eval emitEvaluation (evaluate spec 1000 { replicates := 100, seed := 1 })
#eval emitEvaluation (evaluate spec 10000 { route := .exact, exactBudget := 100000 })
`);
  const rows = (await report('deep-law-enumeration.json')).calculations.map(x => x.evaluation);
  console.log('1,000- and 10,000-level laws enumerate exactly without recursion overflow.');
  for (const route of ['automatic', 'simulation']) {
    await compile(source + `\n#eval emitEvaluation (evaluate spec 10000 { route := .${route}, replicates := 100, seed := 1 })`);
    rows.push((await report(`deep-law-${route}.json`)).calculations[0].evaluation);
    console.log(`10,000-level law: ${route} evaluation completed.`);
  }
  assert.equal(rows.length, 4);
  for (const which of ['typeI', 'power']) {
    for (const row of rows.slice(0, 2)) {
      assert.equal(row[which].method, 'exact');
      assert.equal(row[which].exact, '1');
    }
    for (const row of rows.slice(2)) {
      assert.equal(row[which].method, 'simulation');
      assert.equal(row[which].successes, 100);
    }
    for (const key of ['value', 'seed', 'lower', 'upper'])
      assert.equal(rows[2][which][key], rows[3][which][key]);
  }
  await reject(source + '#eval emitEvaluation (evaluate spec 10000 { route := .exact })',
    /enumeration exceeded its work budget/i);
  // A clean budget failure must leave the same worker usable.
  await compile('import SSC\ntheorem stillUsable : 1 = 1 := rfl');
  console.log('Exact enumeration: branch weights, work limits, 10,000-level laws and unchanged simulation fallback passed.');
}
