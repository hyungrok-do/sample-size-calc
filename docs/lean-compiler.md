# Lean mathematical compiler in WebAssembly

The independent third tab accepts ordinary Lean modules: definitions, hypotheses,
structures, recursive functions, theorem statements, and proof scripts. There is
no `Study`, `Model`, `Analysis`, or `solve` dispatcher in the browser compiler.
Statistical equations live in the source, not in compiler cases.

Lean elaborates that source, constructs proof terms through the specified tactics,
and checks those terms in its kernel. Before reporting success, all local
declarations are independently replayed from the saved import environment with
unconditional kernel checking. This also checks declarations produced while
`debug.skipKernelTC` was enabled, including scoped options inside tactics.
`#eval` runs computable definitions using
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

## Derive a finite distribution from observations

`SSC.Finite` formalizes finite probability distributions over rational weights.
Every `Distribution α` includes proofs that its weights are nonnegative and sum
to one. Expectation, probability bounds, composition, and the collection of equal
outcomes are proved in Lean. No user axioms or floating-point results are used in
these proofs.

Finite event probabilities use balanced rational summation to keep kernel
reduction shallow when weights have large numerators and denominators. Its
equivalence to the defining expectation is proved for every table and event;
the optimization changes neither the model nor the exact probability.
Use `Distribution.probability_eq_expect` to connect a probability expression
to an expectation proof.

`SSC.Derivation.Program α` describes observations and their transformations using
ordinary Lean definitions:

| Operation | Meaning |
| --- | --- |
| `Program.draw distribution` | Draw from a proved finite probability distribution. |
| `Program.pure value` | A deterministic value. |
| `program.bind next` | Continue with a distribution that may depend on the previous outcome. |
| `program.map statistic` | Apply a deterministic statistic. |
| `program.iid n` | Retain a complete independent dataset. |
| `program.iidFold update initial n` | Process independent observations with an arbitrary state update. |
| `Program.iterate transition initial n` | Run dependent or adaptive state transitions. |

Bind needs decidable equality of its intermediate type so equal states can be
combined. These are typed library functions, not additional source grammar.
Programs can use arbitrary user-defined finite outcomes, states and predicates.

`compileWithin program budget` returns either no result or a `Certificate program`:

```lean
distribution : SSC.Finite.Distribution α
correct : ∀ f, distribution.expect f = program.expect f
```

The proof applies to **every** rational-valued function of the output, including
every rejection indicator. It follows from finite additivity and composition,
not from identifying a statistical procedure. Compilation combines equal
intermediate states, so a count based on 25 binary observations has 26 output
states without enumerating all 33,554,432 datasets. The compiler constructs its
certificate alongside that computation. `derive` is the unbounded reference
version intended for mathematical definitions and proofs.

For a larger closed program, `def evidence : Certificate program := by
finite_derive` checks intermediate distributions separately. The tactic
evaluates each proposed table, then requires a kernel-checked equality before
using it in the next certificate. Balanced merging reduces the depth of each
check. Evaluation and quotation are untrusted proposal mechanisms: neither can
introduce an axiom or bypass the equality proof. This also avoids asking the
browser kernel to reduce an entire long computation in a single proof.

Intermediate types need `Lean.ToExpr` instances to propose literal values.
The tactic stages a sequential bind chain, with a limit of 256 checkpoints;
individual continuations must still fit the browser's resource limits.
Parameterized programs use the symbolic correctness theorems above. Staging
changes proof construction, not the observation program or simulation route.

`Program.iidFold_eq_iid` proves that the state update gives the same expectation
as applying the fold to the full dataset. Thus the reduction from observations
to a statistic is explicit and checked. A statistic must have the stated fold
representation to use this reduction; the compiler does not infer an efficient
sufficient statistic from an arbitrary function on datasets.

The corresponding `SSC.Derivation.Specification Design Outcome` defines
`null`, `alternative` and `reject`. The two laws are observation programs, and
the rejection rule is a `Rejection Outcome`. Construct it with `.test predicate`,
or compose rules with `.both`, `.either`, `.all` and `.any`. Preparing a rejection
rule at a design is separated from applying it to each replicate, so calibration
can be computed once.

Its `evaluate` and `search` methods compile these same programs. If the compiler
exhausts its budget, they simulate the original observations and statistic. A
forced `.simulation` route also uses the original program, not the compiled
distribution. Invalid proofs do not trigger fallback. Budget exhaustion does
not truncate, renormalize, or discard outcomes. The default compilation budget
is 1,000,000 units; it accounts for branch operations and bounds on table-merging
work. Exact compilation also limits structural recursion to 64 levels and
returns no certificate when that limit is reached. Program builders use Lean's
proved tail-recursive iteration, and simulation uses suspended conversion and
a typed explicit continuation stack to process long observation sequences.
The general probability-law enumerator also uses an explicit continuation stack,
preserving branch weights, order and work accounting for deeply nested laws.
If its work budget is insufficient, automatic evaluation simulates the same law;
an explicitly requested exact route returns an error instead of a partial table.
These budgets do not bound arbitrary source computations or elapsed time. The browser execution
limit still applies. Calibration functions written using unbounded `derive`
must themselves be computationally feasible.

The three complete modules in `examples/lean-derivations/` are available in the
editor:

| Module | Checked statistical connection |
| --- | --- |
| `bernoulli-observations.lean` | Individual Bernoulli observations imply the count recurrence and derived distribution. A kernel-checked theorem proves minimum n = 25 for the specified calibrated upper-tail rules, and the cutoff 6 satisfies type-I error ≤ 0.05 and power ≥ 0.80. |
| `paired-categorical.lean` | An arbitrary joint distribution of binary pairs implies the distribution of their summed differences. Under the supplied dependent joint laws, n = 23 and cutoff 5 satisfy the stated bounds. |
| `dependent-observations.lean` | State-dependent transitions imply the distribution of a repeated-observation count. A concrete probability is kernel checked without assuming independent observations. |

The observation-to-distribution theorems are symbolic and reusable. The examples
also check concrete probability or sample-size claims with `decide +kernel`,
using staged certificates for the larger paired calculation.
The structured calculation tables remain execution output: only the separate
theorem declarations certify their stated mathematical claims. Monte Carlo and
floating-point samplers are tested implementations, not formally verified
samplers.

The browser regression compares every exact design row in these examples against
independent rational references. `scripts/generate_lean_derivation_references.py`
uses combinatorial binomial and multinomial formulas, and exhaustive Markov
trajectories, rather than the compiler's state-merging algorithm. It also tests
false certificates, changed observation assumptions, invalid probability mass,
bounded compilation, and equality of seeded original-program simulations across
automatic and explicitly requested routes.

This layer establishes finite probability semantics and proof-producing
distribution derivation. Continuous observation-to-statistic theorems,
noncentral-t derivations, asymptotic normality, and certified numerical
integration are not yet included. The numerical examples below continue to
assume their supplied sampling distributions or approximations.

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
neither a general mixed-model optimizer nor a complete formal theory for continuous probability.

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

The browser prepares a trusted Lean frontend function once and calls it through
the pinned runtime ABI, passing source bytes and the report token separately.
The driver parses the input in a separate environment and returns each native
frontend action for the host to execute directly. This keeps command elaboration
off an interpreted driver's call stack. After the last command it calls the audit
directly. The audit is not appended to user source, so `#exit` and user macros
cannot skip or replace it. The report identifier is absent from the input's
source and declaration environment. Printed report-shaped text is ordinary
output; tactic steps come directly from the elaborator's information trees.
Clean import environments are cached by import list and module mode, so repeated
requests reuse initialized libraries without accumulating imported regions.

The report is an observation of elaboration; it does not generate mathematical
claims or fill explanation templates. A tactic trace is not itself the proof;
the completed proof term is kernel checked. A direct term proof can have no
source tactic trace. Traces are limited to 256 entries per declaration and
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
Source and output limits count UTF-8 bytes; output is limited to 8 MiB.
Large valid proofs can reach the time limit on slower or freshly initialized
browsers. A timeout is a resource failure, not a mathematical counterexample;
the worker is discarded and a subsequent request starts a fresh runtime.

This deployment runs Lean's parser, elaborator, dependent type checker, kernel
and executable interpreter. The statistical constructions are libraries in Lean,
not additional grammar or a model-name dispatcher. Its proof-report policy is
deliberately stricter than the language: legal Lean modules using additional
axioms or native-evaluation proof axioms are rejected. The browser does not
provide the complete Lake/project/native-code toolchain or an independent
implementation of Lean's kernel. Replaying declarations is an extra kernel
check in the same runtime, not an external proof-checker certification.

Mathlib is not bundled. General real analysis, continuous probability theory,
unrestricted statistical derivation, and a universal automated mathematician are
not implemented. Real-valued noncomputable definitions require computable bridges
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

Run the differential and resource stress suite against the assembled site:

```sh
python scripts/generate_lean_stress_references.py
cd web
SSC_STRESS_DEADLINE=1 npm run test:lean-stress -- ../_site
```

The seeded reference generator exhaustively enumerates observation trajectories
with Python rational arithmetic. The browser compares 32 generated dependent
models against 160 exact results and invariants, with eight kernel-checked
certificates using intermediate checkpoints. A separate large evaluation expression
may succeed or exhaust the browser stack; the test requires a clean rejection
and restart in the latter case. Raising Lean's recursion budget does not enlarge
the browser's host stack. Additional tests cover dependent types, type classes, macros,
symbolic proofs, invalid syntax and types, forged report/trace text, an early
exit after bypassed kernel checking, UTF-8 boundaries, 1,000 declarations,
100 alternating successful/failed requests across `SSC`, `Std` and implicit
`Init` imports with a heap-growth check, a 100,000-level deterministic exact
program, concurrent requests, cancellation, runtime exit, output flooding and
recovery. `SSC_STRESS_DEADLINE=1`
also waits for the actual 120-second deadline and verifies restart. These are
regression and stress checks of a pinned implementation, not a claim of full
language conformance or exhaustive statistical coverage.
CI runs the complete stress suite including the real deadline; publication also
runs the stress suite before deploying, with the deadline test left to CI.

The earlier native Lean 4.19/Mathlib routines and probability audits remain
separate [legacy command-line tools](lean-native-tools.md). The old
`?runtime=native` UI adapter is retired.
