# Native Lean sample-size compiler

The third tab is for **statistical study planning**. A Lean study specifies a
model, an estimand, an analysis, allocation/counting units, and a power criterion.
The compiler derives the effect and variance coefficient, constructs the
rejection probability, and executes generated C++ to find an admissible integer
sample size. The result shows the statistical reasoning before the optional
C++ source. Low-level algebra commands remain available as implementation tools.

## Run the native compiler

For a hosted application that compiles directly from the browser, use the
[Docker deployment](../packaging/lean/README.md). It bundles all three tabs and
their WASM engines, HTTPS/authentication, and disposable native compiler workers.
GitHub Pages alone remains an editor because it cannot execute the Lean backend.
The following setup is the local development option.

With the pinned Lean 4.19.0 toolchain, a C++17 compiler, and Python 3.9+ on POSIX:

```sh
cd lean
lake exe cache get
lake build
cd ..
python scripts/serve_lean.py
```

Open `http://127.0.0.1:8765/lean.html`. The GitHub Pages version is an editor;
it does not run the native Lean or C++ toolchains. This local service is for
trusted source: Lean modules can perform IO. It checks the loopback origin and
host, accepts one job at a time, and bounds source size, output, and total elapsed
time across Lean checking, C++ compilation, and execution. Stop cancels the
active process group, including the numerical phase.

## Specify a statistical design

```lean
open DesignCompiler.Planning

def study : Study := {
  name := "Treatment-by-visit contrast"
  estimand := "Between-arm difference in mean change"
  analysis := .knownCovarianceZ
  sizeUnit := "participants per arm"
  groups := [{ label := "Control" }, { label := "Treatment" }]
  model := .gaussianContrast {
    labels := ["Control baseline", "Control follow-up", "Treatment baseline", "Treatment follow-up"]
    means := [0, 0, 0, 1/2]
    coefficients := [1, -1, -1, 1]
    covariance := [[1, 3/5, 0, 0], [3/5, 1, 0, 0],
                   [0, 0, 1, 3/5], [0, 0, 3/5, 1]]
  }
  criterion := { alpha := 1/20, targetPower := 4/5, minN := 2, maxN := 10000 }
}

solve sample_size repeated_measures using study
solve power at_one_hundred using study at 100
solve effect detectable_change using study at 100
```

The contrast and matrix are arbitrary numerical Lean expressions, not a lookup
by example or test name. Ordinary user functions can construct them. The example
derives the effect `1/2` and variance coefficient `8/5`, hence estimator variance
`(8/5)/n`. It evaluates both tails of the shifted-normal rejection region and
selects **51 participants per arm (102 total)**: power is 0.80591384, versus
0.79817620 at 50 per arm. The dominant-tail closed-form approximation is displayed
for comparison, but the final answer uses the full two-sided power function.

The covariance describes the vector in **one independent sampling block**.
In this example, a block contains one participant from each arm, with two
measurements per participant. Repeated observations are correlated within a
participant; participants and arms are independent. Allocation entries map one
increment of `n` into arm counts and total participants. They are part of the
user's design specification, not inferred from labels.

A covariance matrix must be dimensionally conformable, symmetric and positive
semidefinite. Validation uses exact rational Schur complements and permits
singular matrices when the requested contrast still has positive variance.
Changing the covariance changes the derived variance and sample size. An
incompatible analysis method is rejected rather than silently relabeled.

## Supported statistical models

| Model | What is derived | Statistical basis |
| --- | --- | --- |
| Arbitrary Gaussian linear contrast | `psi = c^T mu`, `V = c^T Sigma c`, `Var(psi_hat) = V/n` | Exact z-test power under independent jointly Gaussian blocks with known population covariance |
| Finite-stratum ATE | Weighted ATE, efficient influence function, outcome-noise and treatment-heterogeneity variance terms | Asymptotic Wald power for AIPW/TMLE attaining the specified EIF |
| User-defined influence model | Effect and variance computed by arbitrary closed Lean expressions | Conditional on the supplied asymptotically linear representation and variance source |

The examples cover repeated measures, unequal allocation, a cluster-randomized
trial with fixed cluster size and exchangeable ICC, targeted ATE planning, and a
custom pilot influence distribution. Cluster counts are reported separately
from participants. **The Gaussian cluster example assumes known covariance**;
it is not a small-cluster t approximation. The generic covariance representation
can describe many contrasts, but choosing valid sampling blocks and an appropriate
analysis remains part of specifying the statistical model.

For an ATE, the compiler calculates

```text
psi = E[mu1(W) - mu0(W)]
V = E[v1(W)/g(W) + v0(W)/(1-g(W))] + Var[mu1(W)-mu0(W)]
```

The source provides the baseline-stratum probabilities, treatment propensities,
conditional means and variances. The generated explanation shows the EIF and
both variance components. Positivity, consistency/exchangeability for a causal
interpretation, nuisance-estimation conditions, asymptotic linearity and standard
error consistency are stated explicitly. Double robustness alone does not imply
efficiency. Solving for an ATE effect while holding this variance fixed is
rejected; change the planning strata and recompile instead.

All statistical queries currently use two-sided Gaussian/Wald inference.
Estimated-covariance t/F inference, survival likelihoods, and automatic lowering
of arbitrary Lean stochastic models to Monte Carlo simulation are not implemented
in this route. The separate browser study sandbox retains its simulation engine.
An unproved CLT or a suitable estimator is not inferred from arbitrary source code.

## Derivation, numerical execution, and verification

The displayed steps are: estimand/estimator, covariance or EIF variance,
null hypothesis/rejection rule, distribution under the alternative, both-tail
power, and the requested numerical inversion. The generated C++ runs as part of
compilation, so the main output is the requested sample size, power, or detectable
effect, rather than merely a downloadable program. Every statistical result
includes its assumptions and counting units.

```sh
python scripts/lean_equations.py examples/lean-studies/targeted-ate.lean \
  --output result.json --cpp study.cpp
c++ -std=c++17 -O2 study.cpp -o study
./study
```

Lean produces reusable declarations such as `repeated_measures_variance` and
`repeated_measures_effect` proving the exact rational arithmetic in the actual
source module. General checked results include finite-contrast second-moment
propagation, independent-average `V/n` scaling, exchangeable-cluster variance,
and the shifted-normal rejection identity. These theorems are conditional on
their mathematical hypotheses. They do not certify that an actual study satisfies
the sampling assumptions, prove the EIF/CLT for arbitrary nuisance fitting, or
certify the generated floating-point program.

Integer search evaluates every admissible count in ascending order and reports
the selected count's predecessor power. A comparison within `1e-12` of the target
returns `numerically_unresolved`. This is an ambiguity guard, not a certified
error bound. Minimality is numerical and restricted to the requested domain.
Detectable-effect inversion holds the supplied variance coefficient fixed.

Implementation: `lean/DesignCompiler/Planning.lean`, `Planning/Theory.lean`,
`Planning/Command.lean`, and `scripts/lean_study_cpp.py`. Validation runs the native
Lean kernel and emitted C++, compares sample sizes/power/effect inversion with
SciPy, checks a multivariate Gaussian simulation independently, and tests
invalid covariance matrices, analysis mismatches, ATE positivity, theorem reuse,
multiple commands, counting units, and the real browser/service pipeline.

```sh
python scripts/check_lean_studies.py
python scripts/check_lean_equations.py
cd web && npm run test:lean-native-browser
```

Related sources: the pinned Mathlib probability/variance and Gaussian modules;
[Designing efficient randomized trials: power and sample size calculation when
using semiparametric efficient estimators](https://pubmed.ncbi.nlm.nih.gov/34364314/);
and [Estimating the efficiency gain of covariate-adjusted analyses in future
clinical trials using external data](https://pmc.ncbi.nlm.nih.gov/articles/PMC10433828/).

## Underlying equation engine

The service supplies `import DesignCompiler.Planning.Command` (which includes the equation engine), so no header is required in
the editor. Additional imports may refer to installed Lean modules.

```lean
noncomputable def response (x : ℝ) : ℝ := Real.exp (2*x + 1)

solve response_solution for x : ℝ where response x = 5
#check response_solution
```

The generated declaration has type:

```lean
response_solution (x : ℝ) :
  response x = 5 ↔ x = (Real.log 5 - 1) / 2
```

Definitions, namespaces, structures, notation, and ordinary Lean theorems are
elaborated by Lean. Function application uses `response x`, not the custom study
language's `response(x)` convention. Mathematical real functions often require
`noncomputable def`. `solve for x : ℝ where ...` also works, with an automatically
assigned theorem name. A module can contain multiple solve commands.

Parameters and assumptions belong in ordinary sections:

```lean
section
variable (a b y : ℝ) (ha : a ≠ 0) (hy : 0 < y)
solve formula for x : ℝ where Real.exp (a*x+b) = y
-- formula includes a, b, y, ha, hy and proves x = (Real.log y - b) / a.
end
```

Without those hypotheses, the compiler can retain conditions and degenerate
branches instead of dividing by a potentially zero coefficient. Square equations
return both roots, or `False` when no real root exists. A result of `True` means
all values satisfy the equation under the stated assumptions.

## Teach the compiler a new function

`@[equation_rule]` registers a **proved Lean equality or equivalence** for search.
For example, `examples/lean-equations/custom-rule.lean` defines
`cubeLink x = x^3 + x`, proves `cubeLink x = 0 ↔ x = 0`, and solves
`cubeLink (2*x-3) = 0` as `x = 3/2`. This requires no C++ parser or solver change.
The generated theorem can be used by subsequent proofs in the source module.

The built-in rules cover arithmetic rearrangement, exponential/logarithm
inversion with hypotheses, squares, square roots, absolute values, products,
reciprocals, zero-coefficient cases, and completing the square. These are mathematical
primitives, not complete-example printouts. Search composes their proofs on the
actual elaborated expression; Lean simplification, exact arithmetic, and ring
normalization provide further transformations. Each displayed step is drawn
from the proof-producing transformation that ran.

Search explores at most 96 states, favors expressions closer to isolating the
unknown (preferring fewer branches at equal depth), and has heartbeat/time limits. It is not a complete computer algebra
system or an arbitrary theorem prover. It currently rewrites the first matching
occurrence per rule, so rule ordering and expression form can affect discovery.
An **unresolved** result is a checked partial equivalence, not a proof that an
analytic solution is impossible. Unsupported functions can be handled by adding
proved rules or by extending numerical lowering.

The final theorem is kernel checked with kernel checking enabled. Its transitive
axiom dependencies must be a subset of `propext`, `Classical.choice`, and
`Quot.sound`; derivations relying on `sorryAx` or user axioms are rejected. As with
ordinary Lean use, this trusts the installed toolchain, imported modules, and
local metaprogramming environment. Pretty-printed proof terms are explanatory;
the downloadable module reproduces the actual checked declaration.

## Command line and C++

```sh
python scripts/lean_equations.py examples/lean-equations/composition.lean \
  --output result.json --olean Checked.olean --cpp equation.cpp --run
```

Parameterized formulas compile once and accept named real inputs at runtime:

```sh
python scripts/lean_equations.py examples/lean-equations/parameters.lean \
  --cpp equation.cpp --run --param a=2 --param b=1 --param y=5
c++ -std=c++17 -O2 equation.cpp -o equation
./equation a=2 b=1 y=5
./equation y=10 a=3 b=0
```

Every real section variable becomes an input, and every section hypothesis in
the generated theorem becomes a runtime guard. The first example checks `a ≠ 0`
and `0 < y` before evaluating `(log y - b) / a`. Missing, duplicated, unknown,
non-finite, and out-of-range inputs are rejected. Parameter names retain their
Lean spelling, including Unicode; generated C++ identifiers are separate.
An unsupported parameter type or hypothesis disables C++ export with an
explanation. It does not silently drop a theorem assumption.

The result preserves conditional solution sets. With `a*x = y` and no nonzero
assumption, the same executable returns a finite root for `a=2 y=6`, an empty set
for `a=0 y=6`, and all reals for `a=0 y=0`. JSON reports `solutionSet` as
`finite`, `empty`, or `all_reals`, with `roots` populated only for finite sets.
The exported report carries branch predicates separately from their roots.
Conjunctions intersect solution sets, and disjunctions take their union.

For a precision-planning example:

```lean
variable (variance targetVariance : ℝ) (hv : 0 < variance) (ht : 0 < targetVariance)
solve precision_size for n : ℝ where variance / n = targetVariance
-- n = variance / targetVariance
```

```sh
python scripts/lean_equations.py examples/lean-equations/precision.lean \
  --run --param variance=2 --param targetVariance=0.01
```

This returns `n = 200`. It solves the continuous variance equation conditional
on the stipulated variance model; it does not establish integer minimality or
derive a test's power from a data-generating model. Reciprocal rules preserve
Lean's totalized division: without positivity restrictions, `1/x = 0` has the
Lean solution `x = 0`. Add appropriate hypotheses or domain conditions when
that is outside the intended statistical model.

The Python adapter supplies the prelude, invokes Lean, and transports its report.
It does not parse equations. The Lean backend lowers elaborated real expressions
to numeric C++ expressions; Python packages these into a standalone executable.
The generated code checks for finite intermediate values. Standard real arithmetic
and order instances are checked before lowering, and helpers preserve Lean's totalized
`Real.log`, `Real.sqrt`, and division semantics (including exceptional arguments).

For an unresolved equation that has numerical lowering:

```sh
python scripts/lean_equations.py examples/lean-equations/numerical.lean \
  --cpp equation.cpp --run --lower 0 --upper 1
```

This attempts `cos x = x` with bracketed bisection. Parameterized unresolved
equations also accept `--param NAME=VALUE`; their standalone C++ takes the
positional bounds first, followed by the named inputs. It requires a continuous,
sign-changing bracket, checks the residual, and reports failure if the bracket
or precision is inadequate. It does not establish uniqueness or prove the
absence of other roots. C++17 and a `c++` command are required for `--run`.

The C++ backend supports real parameters, arithmetic, constant natural
powers, exp/log/sqrt, sin/cos, and absolute value. Guards support real equality,
inequality, order comparisons, and Boolean combinations. Runtime guard
comparisons use binary64, so rounding can affect branch or assumption decisions
near a boundary. They are not formal proofs of the supplied hypotheses.
Binary64 evaluation, numerical roots, the
transport layer, and the C++ compiler are **not certified by the Lean theorem**.
The theorem certifies exact symbolic equivalence over the stated Lean types.

Monte Carlo models remain in the separate [browser study sandbox](study-language.md),
accessible from this tab at `web/study.html`. The native equation compiler does
not yet lower Lean probability distributions to simulation. The sandbox retains
its previous offline study language and C++/WASM runtime.

## Implementation and validation

- `lean/DesignCompiler/Solve.lean`: command elaboration, proof search, checked
  theorem construction, and report export.
- `lean/DesignCompiler/Solve/Attribute.lean`: extensible theorem registry.
- `lean/DesignCompiler/Solve/Rules.lean`: independently proved transformations.
- `lean/DesignCompiler/Solve/Lower.lean`: lowering from typed Lean expressions.
- `scripts/lean_equations.py`, `lean_cpp.py`, `serve_lean.py`: native transport,
  C++ packaging, and local editor service.

```sh
cd lean && lake build && cd ..
python scripts/check_lean_equations.py
cd web && npm ci && npx playwright install chromium
npm run test:lean-native-browser
npm test
npm run test:lean-browser
```

Tests invoke the actual Lean kernel and compile emitted C++ independently. They
cover user-defined functions/rules, namespaces/structures/notation, theorem
reuse, executable parameters, checked runtime guards, degenerate cases, both
roots, empty and all-real solution sets, intersections, reciprocals, Unicode
parameter names, unsupported hypotheses, invalid numeric input, overflow, invalid proofs and
axioms, unresolved equations, numerical evaluation, and process cancellation.
The native browser test checks the actual service, mobile layout, CSP, multiple
queries, stale-result invalidation, and the static-host unavailable state.

## Earlier Lean declaration templates


The earlier `NormalInput` and `InfluenceInput` examples remain available:

1. Write a `NormalInput` or `InfluenceInput` declaration, or load an example.
2. Choose **sample size**, **power**, or **detectable effect**. Power/effect require
   a fixed count: per arm for the normal model, independent units for the
   influence-function model.
3. **Compile** validates the syntax and parameters and generates a versioned
   plan and standalone C++17 source.
4. **Solve** executes the plan in the C++/WebAssembly runtime.
5. Save the `.lean`, download the `.cpp`, or export the result. Editing source
   or solve settings invalidates previous results. The tab works offline after caching.

```lean
import DesignCompiler.Input
open DesignCompiler

def study : NormalInput := {
  name := "Two-arm normal comparison"
  alternative := 0.5
  nullValue := 0
  sd := 1
  alpha := 0.05
  targetPower := 0.80
  minN := 1
  maxN := 10000
}
```

`sd` is the **known** common population SD, and allocation is equal and fixed.
This is a z-test, not a t-test with estimated variance. Two-sided power is
`Phi(-z-lambda) + Phi(lambda-z)`, where
`lambda = abs(alternative-nullValue) * sqrt(n / (2*sd^2))` and
`z = inversePhi(1-alpha/2)`. `n` counts participants **per arm**.

The example selects 63 per arm (126 total), with power 0.8013023941055797.
At 62 per arm, power is 0.7950080284018122.

For `InfluenceInput`, `variance` is the per-independent-unit coefficient `V` in
`Var(estimator) ≈ V/N`. Specify `.participant` or `.cluster`, the estimand,
estimator, and variance source. This is an asymptotic approximation conditional
on the correct influence representation for the planned estimator. Solving for
effect holds `V` fixed and reports an absolute alternative-minus-null difference.

## Legacy declaration language and execution

The source is valid Lean 4 using `lean/DesignCompiler/Input.lean`. The **browser**
uses a C++ compiler for a closed subset:

- One `def study` structure, optionally preceded by the shown import and `open`.
- `NormalInput` or `InfluenceInput` with their declared fields and defaults.
- Integer/decimal literals, signed values, integer fractions, optional parentheses
  around a number, quoted strings, line comments and nested block comments.
- Fields separated by commas or new lines. Duplicate or incompatible fields fail.
- Other imports, arithmetic expressions, tactics, additional definitions and
  executable commands produce a line/column error.

This older declaration route is **syntax and parameter compilation**, not a Lean kernel invocation in the
browser. Compilation explicitly returns `leanKernelChecked: false`. The static
site contains neither a general Lean compiler nor a C++ JIT. **Solve** interprets
the compiled plan in the precompiled C++/WebAssembly engine. The downloadable C++
implements the same calculation independently using the standard library.

Native Lean checking of the same input also compares the two exported plans:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
cd lean
lake exe cache get Mathlib.Probability.Distributions.Gaussian Mathlib.Probability.CDF Mathlib.Tactic
lake build
cd ..
python scripts/compile_lean.py examples/design/normal.lean --cpp study.cpp
c++ -std=c++17 -O2 study.cpp -o study
./study
```

C++ is emitted only if Lean and the C++ compiler agree on the plan. Use
`--solve power --n 63` or `--solve effect --n 63` for other targets.
Lean 4.19.0 and an exact Mathlib commit are pinned. Library authors may also use
the proof-carrying `NormalStudy`/`InfluenceStudy` types and
`lake exe designc normal` to emit plans for `build/ssc_cli design <file|->`.

## Legacy Lean library verification scope

| Component | What is checked |
| --- | --- |
| Typed study library | Probability, positivity and finite search-domain constraints |
| `normal_rejection_formula` | Exact rejection probability for a shifted standard-normal statistic through its CDF |
| `normal_ncp_lowering` | Equality of source and variance-coefficient noncentrality formulas over exact reals |
| `scan_found` | Every smaller candidate was inadequate before a found candidate was returned |
| `scan_exhausted` | Every candidate in an exhausted finite domain was inadequate |
| Lower/upper bound lemmas | Conditional conclusions when genuine power bounds are supplied |
| ATE variance-term lemma | Nonnegativity given positivity and nonnegative outcome variances |
| JSON, C++ and WASM execution | Tested implementations, without a formal refinement proof |

The raw iid-data-to-normal-statistic derivation, quantile specification,
normal-power monotonicity, floating-point interval bounds, and C++ refinement
are not formalized here. Normal-theory power is exact under its model and
evaluated in ordinary binary64. Input JSON is not proof-authenticated, and a
`formalVerified` flag is rejected rather than trusted.

Sample-size search checks integers in ascending order; it needs no monotonicity
theorem. A comparison within `1e-12` of the target stops as
`numerically_unresolved`. This guard is **not a certified error bound**. The
search never skips an unresolved count. Minimality refers only to the supplied
domain and floating-point evaluations. Searches are bounded at ten million
units and support cancellation. Detectable effect is a numerical inversion
with fixed variance, returning a bracket and endpoint powers.

The versioned `ssc.design-plan` API also supports `ate_eif` through Lean/JSON:
finite-stratum variance is
`sum_w p(w) [var(Y|1,w)/g(w) + var(Y|0,w)/(1-g(w)) + (tau(w)-ATE)^2]`.
It requires the planned AIPW/TMLE to attain this influence function under its
nuisance and remainder conditions. There is no finite-sample power guarantee.
It is not yet an editor template. DAG simulation compilation is deferred.

The C API is `ssc_design_json`, exposed as `Engine.design()` and `ssc_cli design`.
Compilation requests use `op: "compile"`, `source`, `solveFor`, and (for
power/effect) `n`. Plans accept an optional
`solve: {for: "sample_size" | "power" | "effect", n?: integer}`.

## Validation

```sh
build/ssc_tests
python scripts/check_design_compiler.py
python scripts/check_study_compiler.py
cd lean && lake env lean Audit.lean
cd ../web && npm test && npm run test:lean-browser
```

The axiom audit permits only `propext`, `Classical.choice`, and `Quot.sound`,
rejecting placeholders and additional axioms. Checks cover search boundaries,
unit distinctions, invalid source/parameters, effect inversion, ATE overlap and
heterogeneity, independent SciPy references and independently compiled C++.
Expression tests also cover shared versus independent draws, function argument
evaluation, algebraic and numerical roots, simulation bounds, invalid domains,
plan tampering, and generated-C++ agreement with independent SciPy references.
The browser test uses the actual WASM engine for code-based queries and legacy
targets, method explanations, simulation traces, stale-result invalidation,
error handling, persistence, mobile layout, CSP and offline use.
