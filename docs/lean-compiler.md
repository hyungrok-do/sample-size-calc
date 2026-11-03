# Lean study compiler

The **Lean compiler** is a separate third tab at `web/lean.html`:
**Lean source → Compile → generated C++ and execution plan → Solve**.
The Calculator and DAG simulator retain their existing workflows.

## Use the tab

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

## Supported language and execution

The source is valid Lean 4 using `lean/DesignCompiler/Input.lean`. The **browser**
uses a C++ compiler for a closed subset:

- Exactly the shown import and `open`, followed by one `def study` structure.
- `NormalInput` or `InfluenceInput` with their declared fields and defaults.
- Integer/decimal literals, signed values, integer fractions, optional parentheses
  around a number, quoted strings, line comments and nested block comments.
- Fields separated by commas or new lines. Duplicate or incompatible fields fail.
- Other imports, arithmetic expressions, tactics, additional definitions and
  executable commands produce a line/column error.

This is **syntax and parameter compilation**, not a Lean kernel invocation in the
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

## Verification scope

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
cd lean && lake env lean Audit.lean
cd ../web && npm test && npm run test:lean-browser
```

The axiom audit permits only `propext`, `Classical.choice`, and `Quot.sound`,
rejecting placeholders and additional axioms. Checks cover search boundaries,
unit distinctions, invalid source/parameters, effect inversion, ATE overlap and
heterogeneity, independent SciPy references and independently compiled C++.
The browser test uses the actual WASM engine for three solve targets, stale-result
invalidation, error handling, persistence, mobile layout, CSP and offline use.
