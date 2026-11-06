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

New downloads use lossless Brotli packs, with gzip retained for compatibility.
An existing authenticated gzip cache is reused without downloading a second
encoding. Compressed bytes are checked before decoding, and decoded bytes are
checked against their original size and SHA-256 before installation. The
package preserves all 6,930 library files, including interpreted IR bodies:
having a native module initializer does not mean every function has a native
implementation. A browser regression dynamically evaluates an imported
function that needs that IR fallback.

Browsers without native Brotli decoding use a pinned, small decoder in a
temporary worker. Its JavaScript has a content-addressed filename, its WASM
bytes are verified before transfer, and its retained heap is released before
Lean imports begin. The ordinary same-origin static-host trust model for
JavaScript remains unchanged. Runtime/glue loading and three bounded library
loaders overlap. Libraries can be installed as soon as the pinned Emscripten
filesystem exists, while WebAssembly initialization continues. No Lean code
runs until both are ready. Immutable pack views are installed directly into
MEMFS, avoiding a second copy of every library file.

The ready event includes an `initialization` report for runtime assets,
WebAssembly initialization, library loading, file installation, imports and
cache completion. These stages overlap, and file installation is a substage
of library loading, so their durations must not be added together. Readiness
still waits for the complete runtime, libraries, bootstrap and cache writes.

Build artifacts are keyed by actual runtime/standard-library inputs, the Node
driver, module source, module name and prior dependency outputs. Cached outputs
are hash-verified before reuse. A missing/corrupt entry or changed input rebuilds
the module. A failed invocation cannot reuse an old output file, and reported
Lean errors fail the build even if the process exits with status zero.

Compression outputs have a separate verified cache keyed by raw content,
codec, options and compressor version. Unchanged modules therefore avoid both
recompilation and recompression. Corrupt or missing cache entries regenerate
their output. Compression caches stay outside published assets.

Trace assembly indexes steps once by declaration instead of repeatedly scanning
the full trace for every declaration. This changes report assembly from a
declaration-by-trace product to one indexing pass plus per-declaration sorting.

### Measured startup

The [paired Chromium 157 report](benchmarks/lean-startup-chromium157.json)
compares the previous loader with this package on one local HTTP server in the
same browser context. Each package receives a warmup, followed by three cached
starts in alternating order. All three measured pairs improved. The baseline
uses commit `1d48e4cf299681e0335c3413db9da6274b540a48` plus the included
[timing-only instrumentation](benchmarks/lean-startup-baseline.patch).

| Measurement | Before | After |
| --- | ---: | ---: |
| Selected compressed download | 176.00 MiB | 142.36 MiB |
| Cached startup median | 14.72 s | 12.88 s |
| Library files | 6,930 | 6,930 |
| Uncompressed library bytes | 453,465,980 | 453,465,980 |

This is a 19.11% reduction in selected transfer and a 12.49% reduction in
observed cached startup time. It is not an Internet-speed benchmark or a
cross-device guarantee. Gzip and Brotli are both hosted, so total hosting
storage increases even though each new browser downloads less. The Lean WASM
binary and library contents are unchanged. Startup still takes seconds:
restarting a worker reconstructs the runtime and imported environment even
when compressed assets are cached.

The same run passed a checked theorem and dynamic imported-IR evaluation after
every startup, offline startup with zero requests, corrupt-cache rejection,
cancellation while the decoder was active followed by child-worker termination
and restart, and reuse of 30 existing gzip assets with no replacement pack
downloads. The [decoder probe](benchmarks/lean-brotli-decoders-chromium157.json)
records why the smaller decoder-only package was selected. Full startup
measurements, rather than isolated decompression times, determine the result.

## Reproduce

```sh
python scripts/build_lean_wasm.py
python scripts/test_build_lean_wasm.py
npm test --prefix web
npm run test:lean-wasm-browser --prefix web
npm run test:lean-stress --prefix web
npm run test:lean-decompress --prefix web
npm run test:lean-startup-page --prefix web
SSC_STARTUP_REPORT=/tmp/lean-startup.json npm run test:lean-startup --prefix web
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
profiling is opt-in and marked in the report because it can perturb timings.
CI runs the standard benchmark at the production deadline and uploads the
`lean-performance` artifact, including measurements completed before a failure.

Kernel replay, axiom auditing, source/output limits, cancellation and failure
recovery remain part of the browser and stress suites. Performance claims must
name the browser, workload and deadline, and must compare the complete workload
without removing proofs, search rows or simulation replicates.

## Prepared parameterized execution

The browser now exposes an opt-in reusable execution interface. In the Lean tab,
load **Repeated calculations: parameterized normal sample-size search**, expand
**Repeated calculations**, and prepare the model. Parameters are JSON data passed
to an ordinary Lean definition:

```lean
def sscEvaluate (input : Lean.Json) : Except String Lean.Json := do
  -- Validate inputs, construct the specification, and evaluate it here.
  return input
```

The entry's type must be definitionally equal to
`Lean.Json → Except String Lean.Json`; it must be a safe, monomorphic definition.
The interface has no statistical-model dispatcher and introduces no grammar.
The example uses the existing general probability and candidate-search library.
Other definitions can use exact arithmetic or simulate the same specification.

```javascript
const report = await compiler.prepare(source, 'sscEvaluate');
const answer = await compiler.evaluate(report.prepared.id, {
  effectMilli: 500, alphaMillion: 50000, targetMilli: 800, maxN: 200
});
```

Preparation runs the complete source, then independently replays **all** local
Lean declarations and audits their axiom dependencies. Only after those checks
and entry-type validation does the driver retain a runtime function. Every warm
call executes that function with new data; it does not parse or elaborate new
Lean source. Changing parameters does not retrieve a previously cached answer.

The warm reply is explicitly `verification: "compiled-execution"` and
`concreteResultKernelChecked: false`. The module's checked theorems retain their
original statements and hypotheses. Runtime execution, even of a function with
a proved specification, is not a new kernel certificate for its returned
number. A user requiring a fresh theorem about each concrete result must still
use the ordinary compilation/proof path. This feature does not claim to speed
up that path or arbitrary proof search.

One worker retains at most one prepared function. The handle identifies the
exact source SHA-256, package version and entry name. Every compile or prepare
attempt invalidates the previous handle, including invalid source and failed
preparation. Source/entry edits invalidate it in the UI. Cancellation, fatal
runtime failure and worker restart discard the entire session. Parameter errors
leave a valid prepared function available. Requests are serialized across async
hashing so a handle cannot refer to another source's function. Runtime result
JSON is nested under `value`, separate from host-owned verification metadata.

The existing 120-second deadline and 8 MiB output limit remain in effect.
Parameter JSON is limited to 128 KiB. Input and output reject nonfinite numbers
and integers outside JavaScript's safe integer range, rather than silently
rounding exact values. Encode exact integers outside that range and exact
rationals as strings. Other finite decimal JSON numbers use binary64.
A safe Lean definition is not a sandbox against compiler overrides
or expensive computation; the worker/deadline provide execution isolation.

`npm run test:lean-prepared` runs real browser WASM checks for parameter changes,
reference execution, seeded simulation, replay/audit failures, entry signatures,
handle replacement and restart, and data isolation. It reports startup,
preparation, and 30 full warm sample-size searches separately. Set
`SSC_PREPARED_REPORT=/tmp/lean-prepared.json` to retain the report. A controlled
runner can enforce an explicit latency target with `SSC_PREPARED_P95_MS`.
Ordinary CI reports measurements without a hardware-dependent timing assertion.

### Measured prepared execution

The checked-in [Chromium 157 report](benchmarks/lean-prepared-compact-chromium157.json)
records the actual browser WebAssembly run on a local HTTP server. Thirty
queries varied the standardized effect from 0.500 to 0.645, with two-sided
alpha 0.05, target power 0.80, and an ordered search of 2 through 200 participants
per group, stopping at the first qualifying candidate. These are full numerical
searches on each call, not cached answers.

| Measurement | Observed time |
| --- | --- |
| Browser runtime initialization | 13.37 s |
| Model preparation, including replay and audit | 2.56 s |
| Warm API round trip, median | 18.20 ms |
| Warm API round trip, 95th percentile | 28.18 ms |
| Warm API round trip, maximum | 31.55 ms |
| One simulation call, 100 replicates for each law and 3 observations per replicate | 35.48 ms |

The normal-search results matched ordinary `#eval` for the reference parameter
sets. Simulation matched ordinary execution, including successes, seeds and
intervals. The same run passed source/type/axiom rejection, original independent
replay regression checks, 20 model replacements, active/queued cancellation and
restart, and JSON precision/output-limit recovery.

These timings exclude DOM rendering, are not a cross-device guarantee, and do
not benchmark the expensive paired categorical certificate workload. The warm
calls do not create new kernel proofs of their returned numbers. Large
simulations and fresh exact certificates still require separate optimization.
The 13.37 s initialization was measured locally and is not an Internet-download
benchmark. The earlier [prepared-execution report](benchmarks/lean-prepared-chromium157.json)
is retained for reference. Those separate runs are not a controlled latency
comparison; use the paired startup experiment above for the loader comparison.
