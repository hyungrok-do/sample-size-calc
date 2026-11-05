# Compiler performance and verification

The browser executes the pinned Lean parser, elaborator, interpreter and kernel
in WebAssembly. Performance changes preserve the observation programs, theorem
statements, probability calculations and simulation random stream.

## What a request measures

| Report field | Work included |
| --- | --- |
| `headerMs` | Native header parsing and import-environment lookup/loading |
| `parseMs` | Native command parsing, including user syntax extensions |
| `commands` | Source line, actual parser kind and elapsed time for every command |
| `commandMs` | Host measurement of command elaboration, compilation, initial kernel checks and any evaluation |
| `kernelReplayMs` | Unconditional replay of **all** local declarations against the clean imported environment |
| `axiomAuditMs` | Transitive dependency audit, including private/generated declarations |
| `declarationRenderMs` | Statements, parameters, proof bodies and dependencies |
| `traceRenderMs` | Actual elaborator goal traces |
| `totalMs` | Worker request time, excluding initial download/startup and client rendering |

Command time is not labelled kernel-only time. A `#eval` can itself compile code,
search candidate designs and simulate datasets. Header, parser and reporting
measurements do not include that work. Millisecond rounding and host/diagnostic
overhead mean the smaller counters need not sum exactly to `totalMs`.

The UI exposes these measurements in **Compilation timing**, and the downloaded
JSON includes them. Progress messages identify the current command number and
elapsed request time. The 120-second production deadline is unchanged.

## Shared finite derivations

`finite_derive` accepts a closed `Certificate program` or a family such as:

```lean
def checked : (n : Fin 26) → Certificate (observations n.val) := by
  finite_derive
```

The tactic enumerates the finite domain with Lean's `Fin.cases` and `Fin.elim0`.
For each closed program it reduces the outer constructor, recursively obtains
the left-hand certificate of a bind, and applies the generic composition theorem.
It proposes a literal table through evaluation, proves equality with the
computed table using ordinary kernel reduction, and registers the checked
checkpoint. Direct kernel conversion checks reflexivity first, avoiding an
extra interpreted `DecidableEq` comparison over the table. If a lawful quotation
is only propositionally equal (for example another quotient representative),
the tactic restores its saved state and uses the original checked decision
proof. Both paths are tested. No evaluator-trust axiom is introduced.

Rational table entries are stored as normalized data constants. The ordinary
`ToExpr Rat` quotation represents fractions with division expressions, which
causes later reductions to normalize them again. `checkedRational` instead
takes proposed integer numerator/denominator data, proves that the denominator
is nonzero and the components are coprime, and registers a safe `Rat.mk'`
definition. Equal rational values share that definition within the invocation.
The table-equality proof still checks these values against the actual derived
table: valid rational representation alone does not establish a probability.
All generated rational definitions and their proofs participate in final kernel
replay and the dependency audit.

`deriveAtoms` computes the same table as the reference compiler without building
recursive certificate records along the computational path. Its theorem
`deriveAtoms_eq` proves equality with `(derive p).distribution.atoms` for every
program, by structural recursion. Staged bind composition uses this data-only
path; positivity, normalization and all-expectations correctness follow from
that equality and the existing generic composition theorems. This separates
data computation from proof construction without weakening either check.

For result types with an `Ord` instance, composition sorts each branch table
and merges ordered streams. With a lawful ordering, merging two sorted tables
takes at most their combined number of stream steps, avoiding repeated scans
of the accumulated support. Balanced merging combines the branches. Types
without `Ord` retain the equality-only composition path.
An ordering supplied as a local parameter also uses that path, keeping the
evaluated certificate closed even when its surrounding theorem has parameters.

The optimization does not assume that an arbitrary user comparator is lawful.
It checks outcome equality independently before combining weights. Generic
theorems prove conservation of every expectation and of nonnegative weights
for any comparator, even one that always returns `equal`. Poor ordering can
reduce compaction efficiency but cannot identify distinct outcomes or alter
their probabilities. Tests cover that comparator, the missing-instance path,
quotient-valued outcomes and malicious literal quotation.

For detailed diagnosis, `set_option ssc.profileFinite true` reports each
checkpoint's support size and evaluation, quotation, kernel-checking and
registration times. It is off by default and does not affect proof acceptance.

Two generic algebraic reductions also remove unnecessary rational arithmetic.
`multiplyWeight` handles multiplication by zero or one using a proved equality
to ordinary multiplication. `eventMassBalanced` sums only accepted atom masses,
with an inductive proof that it equals the original expectation of the Boolean
indicator. Neither optimization recognizes a statistical model or changes an
exact probability. Balanced summation retains its existing depth bound and
exact zero-fuel fallback.

A memo local to that one invocation maps reduced program expressions to already
checked checkpoint constants. Expression equality resolves hash collisions.
Repeated prefixes therefore reuse proofs; changed expressions do not reuse
their old answers. The memo disappears after the invocation. It cannot retain
declarations from another request or a different imported environment.

There are at most 256 finite-family leaves and 256 recursive checkpoints on a
program path. Exhaustion fails the tactic. It never returns partial evidence.
Other tactic goals are preserved. Empty and nested finite families are tested.

The Bernoulli observation example uses two such families. Its calibration lemma
equates the original rule with the same rule calculated from a certificate.
Another lemma equates the original target with the target computed from checked
tables. The minimum-sample-size theorem reuses the selected design's established
feasibility theorem, and proves that every smaller design fails the power
requirement. It still depends on the
observation-model equivalence theorem, and changing the observation probability
invalidates the original numerical claim.

## Sampling preparation

Sampling repeatedly validates finite weights and converts their cumulative
exact rational masses to binary64 thresholds. `countSamples` keeps a cache of
this deterministic preparation across replicates. Its entries contain the
weight array, threshold array and a proof that the preparation function produced
those thresholds. The cache holds at most eight arrays of at most 1,024 weights.
Larger supports remain valid and bypass storage.

The key contains the ordered exact weights. Outcome values and their types are
not needed to reuse thresholds. Adaptive models with different weights miss the
cache; eviction changes only execution time. Sampling still follows every bind
and observation in the original program, consumes the same uniform variates,
uses the same strict threshold comparisons and charges the same work budget.
Estimator errors, invalid laws and exhausted budgets still return failure with
no partial estimate. `RandomState` and the existing `sample` interface retain
their original two-field state contract.

The regression fixture compares cached sampling against a separate copy of the
original interpreter, including seeds, remaining budget, errors, zero weights,
adaptive laws and more distinct tables than the cache can hold.

`Program.toLaw` also uses Lean's `Thunk` to memoize deferred construction.
Repeated replicates reuse the constructed prefix. `Program.toLaw_bind` is an
ordinary reflexivity proof of the original defining equation: the same deferred
nodes, bind order and interpreter work charges remain. The regression fixture
proves equality to an independent copy of the original translation for arbitrary
programs, and compares seeded execution including exhausted budgets.

## Loading and builds

Browser assets are cached by their authenticated content hashes. An unchanged
runtime or library pack can be reused after a reporter or website update. Cache
reads are rehashed before use. Writes overlap decompression and imports, and
finish before readiness is reported so immediate offline reload still works.
Packs have stable names and boundaries for Init, Lean, SSC and Std. Growing an
SSC module therefore cannot shift unchanged standard-library bytes into new
packs. A packaging regression test crosses an SSC pack boundary, checks that
the other pack hashes remain unchanged, and reconstructs every file from its
manifest offsets.

Build artifacts are keyed by actual runtime/standard-library inputs, the Node
driver, module source, module name and prior dependency outputs. Cached outputs
are hash-verified before reuse. A missing/corrupt entry or changed input rebuilds
the module. A failed invocation cannot reuse an old output file, and reported
Lean errors fail the build even if the process exits with status zero.

Trace assembly indexes steps once by declaration instead of repeatedly scanning
the full trace for every declaration. This changes report assembly from a
declaration-by-trace product to one indexing pass plus per-declaration sorting.

## Reproduce

```sh
python scripts/build_lean_wasm.py
python scripts/test_build_lean_wasm.py
npm test --prefix web
npm run test:lean-wasm-browser --prefix web
npm run test:lean-stress --prefix web
SSC_PERFORMANCE_REPORT=/tmp/lean-performance.json npm run test:lean-performance --prefix web
```

`PLAYWRIGHT_CHROMIUM_EXECUTABLE_PATH` selects a fixed browser binary. The
benchmark accepts a static root followed by optional Lean files, separates
startup from warm requests, and records 30 small requests plus the selected
workloads. By default it includes 100/500/1,000 declarations and all three
observation examples. For a baseline that already times out, explicitly set
`SSC_PERFORMANCE_DEADLINE_MS=600000`. Such a result is marked diagnostic and does
not establish that the workload passes the production deadline.

The report records both elapsed time and renderer-process CPU seconds. They
measure different things: elapsed time includes scheduling delays, while the
CPU counter also includes renderer background and compilation threads. CPU
time is `null` when the sandbox does not expose usable counters. Source hashes
identify the exact inputs used for each measurement. CPU
profiling is opt-in and marked in the report because it can perturb timings.
CI runs the standard benchmark at the production deadline and uploads the
`lean-performance` artifact, including measurements completed before a failure.

Kernel replay, axiom auditing, source/output limits, cancellation and failure
recovery remain part of the browser and stress suites. Performance claims must
name the browser, workload and deadline, and must compare the complete workload
without removing proofs, search rows or simulation replicates.
