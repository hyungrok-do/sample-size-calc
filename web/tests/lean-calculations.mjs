import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { examples } from '../lean-calculation-examples.js';

const header = 'import SSC\nopen Lean SSC.Numeric SSC.Probability SSC.Calculation\nset_option maxRecDepth 10000\nset_option maxHeartbeats 8000000\n';

export async function checkCalculations({ compile, report, page, reject }) {
  const fixture = JSON.parse(await readFile(new URL('./fixtures/lean-numerics.json', import.meta.url), 'utf8'));
  let source = header + `def printNumber (r : Result Float) : IO Unit := match r with
  | .ok x => IO.println (floatJson x).compress
  | .error e => throw (IO.userError e)\n`;
  for (const c of fixture.cases) {
    const args = c.arguments.map(x => `(${x})`);
    const expression = ['chiSquare', 'fisherF'].includes(c.function)
      ? `(ScalarLaw.${c.function} ${args.slice(0, -1).join(' ')}).cdf ${args.at(-1)}`
      : c.function === 'normalCDF' ? `finite (normalCDF ${args.join(' ')})` : `${c.function} ${args.join(' ')}`;
    source += `#eval printNumber (${expression})\n`;
  }
  await compile(source);
  const values = (await page.locator('#lean-evaluation').textContent()).trim().split(/\s+/).map(Number);
  assert.equal(values.length, fixture.cases.length);
  fixture.cases.forEach((c, i) => {
    const tolerance = /Quantile/.test(c.function) ? 2e-7 : 2e-8;
    assert.ok(Math.abs(values[i] - c.expected) < tolerance,
      `${c.function}(${c.arguments}): ${values[i]} vs ${c.expected}`);
  });
  console.log(`${values.length} distribution/quantile cases agree with ${fixture.reference}.`);

  const results = {};
  for (const name of Object.keys(examples)) {
    await compile(examples[name]);
    const checked = await report(`calculation-${name}.json`);
    assert.ok(checked.calculations.length > 0, name);
    results[name] = checked.calculations;
    assert.ok(await page.locator('#lean-calculations').isVisible());
    if (name === 'normal-means') {
      if (process.env.SCREENSHOTS) await page.screenshot({ path: resolve(process.env.SCREENSHOTS, 'lean-calculations.png'), fullPage: true });
      await page.setViewportSize({ width: 390, height: 844 });
      await page.locator('#lean-calculations details').first().evaluate(d => { d.open = true; });
      assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth), false,
        'calculation tables scroll within the mobile viewport');
      await page.setViewportSize({ width: 1400, height: 1000 });
    }
    console.log('Shared-specification calculation passed:', name);
  }
  const selection = r => r.rows[r.selected].evaluation.design;
  assert.equal(selection(results['normal-means'][0]), 63);
  assert.equal(selection(results['student-means'][0]), 64);
  assert.equal(selection(results['student-means'][1]), 0.5);
  assert.equal(selection(results['exact-binomial'][0]), 25);
  // Matched manual-tab baseline inputs: means.two.t (equal allocation) and
  // prop.one (upper, exact), alpha=.05 and target power=.80. This is deliberately
  // two baseline comparisons, not a claim of full procedure/option coverage.
  assert.ok(Math.abs(results['student-means'][0].rows.at(-1).evaluation.power.value - 0.8014595579222547) < 2e-8);
  assert.ok(Math.abs(results['exact-binomial'][0].rows.at(-1).evaluation.power.value - 0.8065115579119346) < 2e-12);
  assert.equal(results['exact-binomial'][0].rows.at(-1).evaluation.power.method, 'exact');
  for (const name of ['normal-means', 'exact-binomial']) {
    const direct = results[name][0].rows.at(-1).evaluation;
    const mc = results[name][1].evaluation;
    for (const key of ['power', 'typeI']) {
      assert.equal(mc[key].method, 'simulation');
      assert.ok(mc[key].lower <= direct[key].value && direct[key].value <= mc[key].upper,
        `${name}: same-specification ${key} direct result inside MC interval`);
      assert.equal(mc[key].successes / mc[key].replicates, mc[key].value);
    }
  }
  for (const name of ['paired-lognormal', 'regression-simulation']) {
    assert.ok(results[name][0].rows.every(row => row.evaluation.power.method === 'simulation'));
    assert.match(results[name][0].rows[0].evaluation.power.reason, /composed law/);
  }

  const samplerSource = header + `
def checkDistribution (law : Law Float) (event : Event Float) : IO Unit :=
  let spec : Specification Nat Float := ⟨fun _ => law, fun _ => law, fun _ => event⟩
  emitEvaluation (evaluate spec 1 { route := .simulation, replicates := 8000, seed := 90210 })
#eval checkDistribution (Law.normal 0 1) (.outside (-1.959963984540054) 1.959963984540054)
#eval checkDistribution (Law.gamma 2 1) (.lower 3)
#eval checkDistribution (Law.beta 2 3) (.lower 0.4)
#eval checkDistribution (Law.studentT 5 2) (.upper 2)
#eval checkDistribution (Law.chiSquare 4 5) (.lower 9)
#eval checkDistribution (Law.fisherF 3 40 5) (.lower 2)
#eval checkDistribution (Law.exponential 2) (.upper 1)
#eval checkDistribution (Law.uniform 0 1) (.inside 0.25 0.75)
`;
  await compile(samplerSource);
  let checked = await report('samplers.json');
  // Independent SciPy values for the eight events above.
  const references = [0.05, 0.8008517265285442, 0.5248, 0.5362932775854147,
    0.5692367217539747, 0.4066443072619022, 0.1353352832366127, 0.5];
  assert.equal(checked.calculations.length, references.length);
  checked.calculations.forEach((r, i) => {
    for (const key of ['power', 'typeI']) {
      const e = r.evaluation[key];
      assert.ok(e.lower < references[i] && references[i] < e.upper, `sampler ${i}, ${key}`);
    }
  });
  console.log('Eight continuous samplers agree with independent event probabilities.');

  await compile(header + `
def expensiveIntegral : Specification Nat Float where
  null := fun _ => Law.fisherF 3 40
  alternative := fun _ => Law.fisherF 3 40 700
  reject := fun _ => .upper 200
#eval emitEvaluation (evaluate expensiveIntegral 1 { replicates := 2000, seed := 71 })
def countModel : Specification Nat Nat where
  null := fun _ => Law.poisson 20
  alternative := fun _ => Law.poisson 30
  reject := fun _ => .ofBool (fun k => k >= 35)
#eval emitEvaluation (evaluate countModel 1 { replicates := 4000, seed := 72 })
`);
  checked = await report('numerical-fallback.json');
  const integrated = checked.calculations[0].evaluation;
  assert.equal(integrated.typeI.method, 'numerical');
  assert.equal(integrated.power.method, 'simulation');
  assert.match(integrated.power.reason, /Noncentral mixture exceeds/);
  assert.ok(integrated.power.lower < 0.7724280622765395 && 0.7724280622765395 < integrated.power.upper);
  const counted = checked.calculations[1].evaluation;
  assert.equal(counted.power.method, 'simulation');
  assert.match(counted.power.reason, /unbounded support/);
  assert.ok(counted.power.lower < 0.20269167451688266 && 0.20269167451688266 < counted.power.upper);
  console.log('Numerical-budget and unbounded-support fallbacks preserve the specified distributions.');

  const same = header + `
def composed : Specification Nat (List Bool) where
  null := fun _ => (Law.bernoulli (1/2)).iid 8
  alternative := fun _ => (Law.bernoulli (3/4)).iid 8
  reject := fun _ => .ofBool (fun xs => (xs.filter id).length >= 5)
#eval emitEvaluation (evaluate composed 1 { exactBudget := 20, replicates := 4000, seed := 77 })
#eval emitEvaluation (evaluate composed 1 { exactBudget := 20, replicates := 4000, seed := 77 })
#eval emitEvaluation (evaluate composed 1 { route := .exact, exactBudget := 10000 })
`;
  await compile(same);
  checked = await report('enumeration-fallback.json');
  const [a, b, exact] = checked.calculations.map(c => c.evaluation);
  assert.deepEqual(a, b, 'same seed and specification reproduce identical estimates');
  assert.match(a.power.reason, /work budget/);
  for (const key of ['power', 'typeI']) {
    assert.ok(a[key].lower <= exact[key].value && exact[key].value <= a[key].upper);
  }

  await compile(header + `
def nonmonotone : Specification Nat Bool where
  null := fun _ => pure false
  alternative := fun n => pure (n == 7 || n == 12)
  reject := fun _ => .ofBool id
#eval emitSearch (search nonmonotone #[5, 6, 7, 8, 12])
def uncertain : Specification Nat Bool where
  null := fun _ => pure false
  alternative := fun n => if n == 1 then Law.bernoulli (4/5) else pure true
  reject := fun _ => .ofBool id
#eval emitSearch (search uncertain #[1, 2] {} { route := .simulation, replicates := 100, seed := 7 })
#eval emitSearch (search uncertain #[1] {} { route := .simulation, replicates := 100, seed := 7 })
`);
  checked = await report('search-uncertainty.json');
  assert.equal(selection(checked.calculations[0]), 7);
  assert.equal(selection(checked.calculations[1]), 2);
  assert.deepEqual(checked.calculations[1].earlierUnresolved, [0]);
  assert.equal(checked.calculations[2].selected, null);
  assert.match(await page.locator('#lean-calculations').textContent(), /minimum has not been established/);

  const invalid = (law, event = '.ofBool (fun _ => true)', options = '{}') => header + `
def invalidSpec : Specification Nat Float := ⟨fun _ => ${law}, fun _ => ${law}, fun _ => ${event}⟩
#eval emitEvaluation (evaluate invalidSpec 1 ${options})`;
  await reject(invalid('Law.normal 0 (-1)'), /positive SD/);
  await reject(invalid('Law.normal 0 1', '.invalid "invalid rejection rule"'), /invalid rejection/);
  await reject(invalid('Law.normal 0 1', '.test (fun _ => .error "estimator did not converge")'), /estimator did not converge/);
  await reject(invalid('Law.normal 0 1', '.upper 2', '{ route := .simulation, maxDraws := 5 }'), /budget exhausted/);
  await reject(invalid('Law.normal 0 1', '.upper 2', '{ replicates := 0 }'), /replicates/);
  await reject(invalid('Law.finite #[(1.0, 2)]'), /sum exactly to one/);
  await reject(invalid('Law.normal 0 1', '.ofBool (fun x => x > 2)', '{ route := .numerical }'), /No direct integral/);
  console.log('Invalid models, invalid rules, failed estimators and exhausted budgets fail explicitly.');
}
