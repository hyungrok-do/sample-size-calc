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
`Init`, `Std`, `Lean`, and `SSC`. The last module provides reporting and a generic
proof-carrying bounded search. It has no statistical formulas or custom tactics.
There is no additional source grammar to learn.

Use a theorem to state a goal. Supply hypotheses as theorem parameters, and use
Lean proof constructs such as `by`, `have`, `calc`, `exact`, `simp`, `omega`,
`grind`, and `decide +kernel`. Proof automation is bounded and may require lemmas
or intermediate arguments. An unresolved goal is a compilation failure with
Lean's actual goal and context; no approximate answer substitutes for a proof.

The compiler preserves symbolic parameters. It does not require equations to be
closed numeric terms. User-defined types, namespaces, notation and functions
are handled by Lean's parser and elaborator.

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

The complete source modules are in `examples/lean-proofs/`:

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
bundled, and `lean/browser/SSC.lean` is compiled by that same WASM binary.
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
