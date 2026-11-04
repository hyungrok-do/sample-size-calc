> This is the earlier custom language at `web/study.html`. The third tab now
> opens the [Lean compiler in WebAssembly](lean-compiler.md), which accepts
> actual Lean source and checks generated theorems in Lean's kernel.

# Browser study language sandbox

The **browser study sandbox** at `web/study.html` accepts user-defined formulas and random
models. Imports and headers are optional: the mathematical runtime is built in.
The syntax borrows `def`, typed parameters, `let`, `if`, and `fun` from Lean.
`solve for`, `evaluate`, `probability`, and `simulation` are this compiler's own
syntax. This is not an unrestricted Lean 4 implementation or a proof assistant.

## Write how to solve in the source

```lean
def objective (x : Real) : Real := x^2 + 2*x

solve for x : Real in [0, 10]
where objective(x) = 8
```

The answer is `x = 2`, obtained algebraically and checked in the original
expression. **Solution derivation** shows normalization, completing the square,
both roots, and restriction to the domain. C++ is an optional collapsed output.
A query specifies the variable, domain, and condition:

| Query | Meaning |
| --- | --- |
| `evaluate objective(2)` | Calculate one value |
| `solve for x : Real in [a,b] where f(x) = c` | Find a real equation root |
| `solve for n : Nat in [a,b] where f(n) >= c` | First qualifying integer |
| `solve for x : Grid in [a,b] step h where f(x) < c` | First qualifying point on a finite grid |

Bounds and settings must be constant expressions. `Nat` permits a positive
integer `step` too; its default is 1. The solve variable must be a fresh name.
`Real` accepts `=` only. Both sides may contain the unknown. Optional
`tolerance 0.000001` follows the bounds/step and precedes `where`. Each file has
definitions, optionally one simulation configuration, and one final query.

## Lean proofs of the derivation

For real equations built from rational arithmetic with degree at most two,
the compiler emits exact algebraic steps and a `.lean` proof. For example:

```text
x² + 2x = 8
x² + 2x − 8 = 0
(2x + 2)² = 36
x = 2 or x = −4
x ∈ [0,10] ⇒ x = 2
```

The proof establishes normalization, the complete set of roots (including
irrational square-root expressions), and, for rational roots, the filtered
domain solution. Negative discriminants and constant equations have explicit
proofs too. Source decimals are interpreted as exact rationals here: `0.1`
means exactly `1/10`. A separate bounded rational layer preserves source
constants, so a binary64 approximation is never presented as exact algebra.
Exact rational numerators/denominators are currently limited to 10^9 and source
numeric mantissas to nine digits; unsupported expressions retain the numerical
route and are labelled accordingly. The native Lean proof is the final check
of each generated symbolic claim.

```sh
python scripts/compile_study.py examples/programs/equation.study \
  --lean derivation.lean --check-lean
```

Use the [pinned Lean/Mathlib setup](lean-compiler.md) first. The browser exports
the proof and says **not kernel checked**; the native command actually invokes
Lean and reports `leanDerivationKernelChecked: true` only after success.
These theorems concern exact real arithmetic, not floating-point execution,
tolerance-based root acceptance, quadrature accuracy, or a stochastic model's
validity. Gaussian cases show their analytic probability derivation and identify
numerical inversion separately. Simulation shows the estimator and its bounds;
it is never labelled an analytic solution or a Lean proof.

## Probability and automatic simulation

```lean
def power (n : Real) : Real :=
  let statistic := normal(0.5 * sqrt(n / 2), 1)
  probability(abs(statistic) > normalQuantile(0.975))

solve for n : Nat in [2, 10000]
where power(n) >= 0.80
```

For this equal-allocation, known-variance normal test, `n` counts participants
per group. The compiler derives the two normal tails and finds 63 per group
(126 total), with power approximately 0.8013023941. No power formula is hardcoded
to this function name: changing the supplied event changes the derived expression.

For a nonlinear event, use the same query mechanism:

```lean
def reject (n : Nat) : Bool :=
  let control := sample(n, fun i => exp(normal(0, 0.6)))
  let treated := sample(n, fun i => 0.8 + exp(normal(0, 0.6)))
  let se := sqrt(variance(control) / n + variance(treated) / n)
  abs(mean(treated) - mean(control)) / se > normalQuantile(0.975)

simulation { replicates := 5000, seed := 2026, confidence := 0.95 }

solve for n : Nat in [10, 60] step 10
where probability(reject(n)) >= 0.80
```

Here the user defines a rejection rule using estimated standard errors. Its
nominal threshold does not itself establish finite-sample type-I error control.
The compiler simulates the model and estimates the probability of that rule.
To examine size under a null model, write the null-generating model and evaluate
the same rejection event.

`method := auto` is the default. `analytic` rejects events outside the available
analytic rules; `simulation` forces Monte Carlo. The automatic rules handle
comparisons of affine Gaussian expressions and absolute-value thresholds.
Shared draws are tracked: `let x := normal(0,1); x+x` has variance 4; two separate
calls have variance 2. Function arguments are evaluated once, left to right.
Each sample lambda is evaluated separately for each observation.

Simulation currently requires `probability(event)` as the whole query quantity,
possibly wrapped by user functions and `let` bindings. It does not support
nested probabilities, arbitrary arithmetic combining nonanalytic probabilities,
exact stochastic equality, or continuous stochastic root finding. Unsupported
forms produce an error. Use `Nat` or `Grid` with `<`, `<=`, `>`, or `>=`.

## Expressions and built-ins

Define functions before use; recursion and forward references are unsupported.
Supported types are `Real` (`Float` and `ℝ` aliases), `Nat` (`ℕ` alias), `Bool`,
and `Array Real`. Numeric values use binary64; `Nat` annotations add integer
checks. `f(a,b)` and curried `f a b` calls are supported. Names may use Unicode.
Use `--` for line comments and `/- … -/` for nested block comments.

| Category | Operations |
| --- | --- |
| Arithmetic | `+ - * / ^`, unary signs, `pi`, parentheses |
| Boolean | `< <= > >= = !=`, `&&`, `||`, `!`, `true`, `false`; `≤ ≥ ≠ ∧ ∨ ¬` aliases |
| Local control | `let x := value; body` (newline or `in` also accepted), `if condition then a else b` |
| Math | `sqrt log exp abs sin cos floor ceil min max pow logistic` |
| Normal functions | `normalCDF(x)`, `normalSF(x)`, `normalQuantile(p)` |
| Scalar draws | `normal(mean,sd)`, `uniform(lo,hi)`, `bernoulli(p)`, `exponential(rate)`, `gamma(shape,scale)`, `binomial(n,p)` |
| Arrays | `[a,b,c]`, zero-based `xs[i]`, `sample(n, fun i => expression)`, `map(xs, fun x => expression)` |
| Sample helpers | `normalSample(n,mean,sd)`, `lognormalSample(n,logMean,logSd)`, `bernoulliSample(n,p)`, `exponentialSample(n,rate)` |
| Summaries | `sum mean variance sd length dot quantile`; variance uses denominator `n-1`, quantile uses linear interpolation |
| Integration | `integrate(fun x => expression, lo, hi, evenPanels)`; composite Simpson rule |

Top-level values must be deterministic. Put draws in parameterized functions
or local bindings inside a probability expression. Domain errors, unknown names,
wrong argument types, out-of-range indices, and non-finite arithmetic fail.
Scalar/array random output must be converted to a Boolean event inside
`probability`; it is not silently treated as a deterministic objective.

## Solution methods and uncertainty

Linear and quadratic polynomial objectives use algebraic candidates, validated
against the original expression. Other real equations use bisection with an
absolute residual tolerance (default `1e-8`). An endpoint sign change is needed
unless an endpoint already meets the tolerance. Continuity is a model
assumption, uniqueness is not claimed, and a narrowed interval without a small
residual reports `numerically_unresolved`. Failure to bracket is not proof that
no root exists. Quadrature error is separate from root residual: choose enough
panels and check convergence for integration-based objectives.

Finite-grid searches scan in ascending order, without assuming monotonicity.
Gaussian probability comparisons within `1e-12 * max(1, abs(target))` report
numerical uncertainty instead of skipping the candidate. This heuristic is not
a certified floating-point bound. Deterministic results otherwise refer to the
evaluated floating-point expressions and supplied search domain.

Monte Carlo uses fixed replicates `R` per candidate and deterministic random
streams indexed by seed, grid position and replicate. With `K` supplied grid
points, the two-sided Hoeffding radius is
`sqrt(log(2*K/(1-confidence)) / (2*R))`. Clipped intervals cover all candidate
probabilities simultaneously under independent, identically distributed
replicates. The full grid size is used even if the search stops early.
These conservative bounds cover simulation error, not model error, numerical
error, or the validity of a statistical decision rule.

A minimum is reported only when the candidate is confirmed adequate and all
preceding grid points are confirmed inadequate. If earlier candidates remain
uncertain, the result gives a range for the minimum (`bounded`). If there is no
confirmed adequate upper point, the result is `unresolved` or
`not_found_in_domain`, depending on whether any candidates remain possible.
Increase replicates or refine the grid as needed; point-estimate minima are
labelled separately. Saved results include the evaluated candidates and bounds.

## Build and execute

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
python scripts/compile_study.py examples/programs/normal-power.study --cpp study.cpp --run
```

The compiler performs lexing, parsing, type checks, function expansion,
exact algebraic derivation and Lean-proof emission where supported, probability
reduction, solver selection, and C++ emission. Generated C++ contains
the emitted operations and shared numerical runtime, with no source parser or
AST interpreter. Compile it independently with `c++ -std=c++17 -O2 study.cpp`.
The browser uses the same runtime through a typed-expression interpreter in
precompiled WASM; it does not run a C++ JIT. Both paths are checked against each
other and independent numerical references by `scripts/check_study_compiler.py`.
`scripts/check_study_derivations.py` kernel-checks generated proofs and audits
their axioms across linear, quadratic, rational, irrational and degenerate cases.

The existing `ssc_design_json` / `ssc_cli design` entry point accepts
`{"op":"compile","source":"…","solveFor":"code"}` (`solveFor` may be omitted).
The result contains generated C++, method explanations, and an
`ssc.expression-plan` v1. Executing a saved plan reparses its source and checks
that the stored typed IR matches; edited or unrecognized plans fail.

Limits: 128 KiB source, 30,000 expanded expression nodes, nesting depth 120,
100,000 observations/panels per array/integration, ten million deterministic
grid points, 2,000 simulation grid points, 100–1,000,000 replicates per point,
ten million simulated trials per query, and a runtime work budget. Stop cancels
the browser worker. Seeds are integers from 0 through 2^53−1. New study programs
use `.study`; legacy `.lean` declaration templates remain supported with their
separate [Lean checking workflow](lean-compiler.md).
