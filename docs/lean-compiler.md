# Lean mathematical compiler in WebAssembly

The independent third tab accepts ordinary Lean modules: definitions, hypotheses,
structures, recursive functions, theorem statements, and proof scripts. There is
no `Study`, `Model`, `Analysis`, or `solve` dispatcher in the browser compiler.
Statistical equations live in the source, not in compiler cases.

Lean elaborates that source, constructs proof terms through the specified tactics,
and checks those terms in its kernel. `#eval` runs computable definitions using
Lean's interpreter inside the same WebAssembly runtime. There is no separate
C++ statistical backend and no compiler-server request.

## Write a module

Start with `import SSC`, then write ordinary Lean. The bundled library contains
`Init`, `Std`, `Lean`, and `SSC`. SSC imports reusable numerical, probability,
and calculation modules, and provides reporting and a generic proof-carrying
bounded search. There is no additional source grammar to learn and no dispatch
on statistical procedure names.

Use a theorem to state a goal. Supply hypotheses as theorem parameters, and use
Lean proof constructs such as `by`, `have`, `calc`, `exact`, `simp`, `omega`,
`grind`, and `decide +kernel`. Proof automation is bounded and may require lemmas
or intermediate arguments. An unresolved goal is a compilation failure with
Lean's actual goal and context; no approximate answer substitutes for a proof.

The compiler preserves symbolic parameters. It does not require equations to be
closed numeric terms. User-defined types, namespaces, notation and functions
are handled by Lean's parser and elaborator.

## One specification, direct calculation and simulation

`SSC.Calculation.Specification Design Outcome` contains three typed functions:

```lean
null        : Design → Law Outcome
alternative : Design → Law Outcome
reject      : Design → Event Outcome
```

`Design` can be a sample count, a proposed effect, or a record of allocation,
cluster count, cluster size and other parameters. `Outcome` can be a test
statistic or a complete dataset. `Law` has ordinary Lean `pure`, `bind` and
`iid` operations, so latent variables, dependence, transformations and custom
estimators can be expressed in ordinary source. There is no separate simulation
specification to keep in sync.

`evaluate spec design options` calculates both null rejection probability and
power. `search spec candidates requirements options` searches an explicitly
ordered candidate array. The default route uses:

1. Exact rational enumeration of finite probability programs within a work budget.
2. Numerical CDF evaluation for supported scalar laws and interval/tail events.
3. Simulation of the same law, applying the same rejection event.

Unavailable integration and exhausted enumeration budgets can trigger fallback.
Invalid probabilities, invalid distribution parameters, invalid rejection rules,
failed estimators, Lean errors and failed proofs do not become successful
simulation results. Estimator failures abort; replicates are never silently
discarded. Simulation has explicit replicate, random-draw and browser-time limits.

Set `route := .simulation` to compare direct evaluation with simulation using the
identical specification. `.exact` and `.numerical` require that route and report
unavailability instead of falling back. An arbitrary rejection function uses
`Event.test : (Outcome → Except String Bool) → Event Outcome`; this lets an
estimator report nonconvergence. `Event.ofBool` is the simpler total-function form.

The numerical library includes normal CDFs and quantiles, regularized incomplete
gamma and beta, central t quantiles, noncentral t integration, and central and
noncentral chi-square/F CDFs. Composable samplers include these distribution
families plus uniform, exponential, gamma, beta, binomial, Bernoulli, Poisson and
arbitrary finite rational laws. A new test requires a source definition, not a
new compiler branch. A distribution of a statistic supplied as an asymptotic
approximation remains an approximation; the library does not infer its validity.

Simulation reports success counts, fixed replicate counts, the initial stream
state and Hoeffding intervals. The random generator is SplitMix64. A search
allocates its error probability across both null/alternative evaluations at
every planned candidate. Thus the Monte Carlo intervals are simultaneous under
the usual independent-draw assumptions, including when the search stops at its
first acceptable candidate. They do not account for model misspecification,
finite-precision sampling or deterministic numerical approximation error.

Search uses the lower power endpoint and, when requested, the upper type-I-error
endpoint. A candidate near a constraint can remain **unresolved**. A later
acceptable candidate is reported with all earlier unresolved indices; no minimum
is claimed. No monotonicity assumption or bisection shortcut is used. Search
makes no claim about candidates outside the supplied set. Numerical comparisons
use a configurable decision guard, not a certified error bound.

`Requirements.maxTypeI` is optional. When absent, type-I error is reported but
does not constrain the search. A Monte Carlo interval around a nominal 5% null
rejection probability will generally not establish a strict ≤5% bound; increase
precision or supply a justified calibration theorem rather than ignoring this
uncertainty. Exact rational search comparisons do not round the constraints.

Use standard Lean `#eval emitEvaluation (...)` or `#eval emitSearch (...)` to
produce a structured result table and downloadable JSON. Binary64 probability
values and interval endpoints retain all bits in the report. Numerical results
are execution output, separate from theorem proofs.

## Calculation examples and coverage

`examples/lean-calculations/` contains editable complete modules:

| Module | Calculation | Routes exercised |
| --- | --- | --- |
| `normal-means.lean` | Two normal means with known variance, per-group sample size | Numerical CDF and same-specification simulation |
| `student-means.lean` | Equal-allocation pooled t test, sample size and detectable-effect grid | Noncentral t integration |
| `exact-binomial.lean` | Null-calibrated one-sided binomial test, power and sample size | Exact rational evaluation and same-specification simulation; a small probability theorem is kernel checked |
| `paired-lognormal.lean` | Correlated lognormal pairs, studentized paired differences | Automatic simulation of actual observations; null size is measured |
| `regression-simulation.lean` | Random covariates, heteroscedastic errors, HC1 slope test | Automatic simulation with user-defined estimation and rejection |

These examples establish an extensible computation path. They are not a claim
of parity with all **71 manual procedures / 163 procedure–target combinations**.
Full parity across every procedure option, and formal derivations of every
statistical approximation, remain separate work. The library currently supplies
neither a general mixed-model optimizer nor a complete formal probability theory.

`web/tests/fixtures/lean-numerics.json` records 53 independent SciPy reference
values for the numerical primitives. Browser tests execute the actual Lean WASM
interpreter, compare direct and simulated probabilities, check eight continuous
samplers, exercise arbitrary composed laws and nonmonotone designs, and test
unresolved candidates, reproducibility and explicit errors. These are statistical
computation tests, distinct from the existing frontend regression suite.

## Proof and evaluation reports

The report reads declarations from the kernel-checked environment. It shows:

- The actual statement and pretty-printed proof or defining expression.
- Parameters and formal hypotheses in the declaration's type.
- Source tactics and their real before/after goal states from Lean's elaborator.
- Direct declaration dependencies and transitive foundational axiom dependencies.
- Separate output from `#eval`, `#check` and other source commands.

The report is an observation of elaboration; it does not generate mathematical
claims or fill explanation templates. A tactic trace is not itself the proof;
the completed proof term is kernel checked. A direct term proof can have no
source tactic trace. Traces are limited to 256 entries per source command and
explicitly marked when truncated. Pretty printing may abbreviate large terms; the downloaded
source is the reproducible proof script.

`sorry`, unproved user axioms, and `native_decide`/compiler-trust axioms cannot
produce a successful report. The accepted foundational axioms are `propext`,
`Classical.choice`, and `Quot.sound`, and their actual use is shown. Safe
definitions are distinguished from unsafe/partial executable definitions.

An evaluated number is not automatically a certified statistical answer. A
separate theorem must establish its intended property. In particular, source
hypotheses remain conditions, a supplied formula is not automatically derived
from a probability model, and floating-point error bounds are not inferred.

## Generic finite search

`SSC.leastBelow P bound`, for any decidable `P : Nat → Prop`, produces a value of
`SSC.SearchResult P bound`:

- `found` carries a candidate, a proof it is below the bound, a proof of `P`,
  and a proof that every smaller natural number fails `P`.
- `exhausted` carries a proof that every candidate below the bound fails `P`.

`result.value?` exposes the computable witness; `result.sound` exposes the
corresponding proposition. The implementation is proved by recursion, with no
monotonicity, continuity, model-family or statistical-method assumptions. A
search that exhausts its bound makes no claim about candidates beyond it.

This is one optional library construction, not a required way to derive answers.
Parametric proofs and other user-defined algorithms compile through the same
path. A theorem about a computed witness can use `decide +kernel`, so the
reported numeral is checked through kernel reduction rather than trusting native
evaluation.

## Statistical examples

The proof-oriented source modules remain in `examples/lean-proofs/`:

| Module | What is established | Scope |
| --- | --- | --- |
| `exact-trial.lean` | Minimum sample size 25, excluding every smaller size over all specified threshold rejection rules; the cutoff 6 meets the exact error/power inequalities | Independent Bernoulli null/alternative weights and rejection equations are supplied in the source; the finite model's interpretation is not a separately formalized probability-space theorem |
| `certified-bounds.lean` | Feasibility and minimality from arbitrary lower/upper power bounds | All power functions, domains and bound hypotheses remain symbolic; numerical bounds must be supplied and justified separately |
| `variance-propagation.lean` | Symmetric covariance expansion and its application to supplied ratio-linearization coefficients | Does not prove differentiability, an asymptotic normal limit, or a delta-method theorem |

Changing these equations needs no compiler change. The tests also compile fresh
non-monotone predicates, exhausted search domains, namespaces and structures;
mutating a trial equation invalidates its old conclusion.

## Deployment and limits

The actual Lean 4.33.0-pre WASM runtime is pinned to commit
`62b6a2291302d4bbeace37642a066b7510d0145c`. The matching module artifacts are
bundled, and all modules under `lean/browser/` are compiled by that same WASM binary.
Native 64-bit `.olean` files cannot substitute for wasm32 files.

The first load is approximately 175 MiB, cached for later and offline use.
A current desktop browser with shared-memory WebAssembly is required. On static
hosting, the service worker establishes cross-origin isolation and may reload
once. Stop terminates the worker; restart reuses authenticated cached assets.
Compilation has a 120-second browser limit and a 128-KiB source limit.

Mathlib is not bundled. General real analysis, probability theory, automatic
statistical model derivation, and a universal automated mathematician are not
implemented. Real-valued noncomputable definitions require computable bridges
before numerical execution. These are explicit limits of the available library
and automation, not a hidden fallback to predefined statistical routines.

Build with Node 24+ (or the Node runtime configured by Emscripten):

```sh
python scripts/build_lean_wasm.py
```

The complete calculator/DAG site still uses its existing C++ engines. To assemble
and verify the entire site with Emscripten 6.0.13:

```sh
web/build.sh
python scripts/build_lean_wasm.py
web/assemble.sh _site development
cd web && npm ci && npm run test:lean-wasm-browser -- ../_site
```

The browser regression uses the actual WASM compiler on a static subpath without
server isolation headers. It checks proof traces, symbolic hypotheses, changed
equations, failed proofs, axiom auditing, numerical-output separation, cancellation
and offline restart. API requests and the removed C++ statistical module are
explicitly excluded from this compiler path.

The earlier native Lean 4.19/Mathlib routines and probability audits remain
separate [legacy command-line tools](lean-native-tools.md). The old
`?runtime=native` UI adapter is retired.
