# Engine interface

The C++ engine is a shared library (`libssc_core.so`, `.dylib` or `.dll`) with a
small C interface declared in [`core/include/ssc/ssc.h`](../core/include/ssc/ssc.h).
Requests and results are UTF-8 JSON, so any language with a C foreign-function
interface can use the engine. The Python package (`ssc.run`, `ssc-calc`,
`ssc-mcp`) is one such client; the web page calls the same functions in the
engine compiled to WebAssembly (see [`web/README.md`](../web/README.md)).

## C functions

| Function | Purpose |
|---|---|
| `const char* ssc_version(void)` | Version string (static, do not free). |
| `int ssc_hardware_threads(void)` | Number of hardware threads. |
| `char* ssc_catalog_json(void)` | Categories and the list of procedures. |
| `char* ssc_procedure_json(const char* id)` | Schema of one procedure, or `NULL` for an unknown id. |
| `char* ssc_search_json(const char* query, int limit)` | Procedures ranked for a free-text query (see [Search](#search)); `limit` ≤ 0 returns every match. |
| `char* ssc_run_json(const char* request, ssc_token* token)` | Run a calculation; `token` may be `NULL`. |
| `void ssc_free(char* p)` | Release any `char*` returned by the library. |
| `ssc_token* ssc_token_new(void)` | Create a cancellation/progress token. |
| `void ssc_token_cancel(ssc_token*)` | Ask a running calculation to stop (thread-safe). |
| `double ssc_token_progress(ssc_token*)` | Fraction of the work done, 0 to 1 (thread-safe). |
| `void ssc_token_free(ssc_token*)` | Release a token. |
| `int ssc_eval(const char* name, const double* args, int nargs, double* out)` | Evaluate a distribution function, e.g. `pnt`, `qchisq`, `pcorr`; returns 0 on success. |

All functions are thread-safe. A calculation runs on the calling thread and,
for its scenarios, power-curve points and simulated studies, on an internal
pool of worker threads that calculations running at the same time share.
Another thread can watch progress and cancel it through the token.

## Calculation request

```json
{
  "procedure": "means.two.t",
  "solve_for": "n",
  "inputs": {"power": "0.8 0.9", "mu1": "1", "mu2": "0", "sigma": "1"},
  "options": {"threads": 0, "power_curve": true, "curve_points": 41, "max_scenarios": 20000,
              "approximate": false}
}
```

- `solve_for` is one of the procedure's `solve_for` values (usually `n`,
  `power` or an effect-size parameter).
- `inputs` maps parameter ids to text. Missing parameters take their defaults,
  and parameters hidden by the current choices are ignored. A numeric input can
  be a single value or a list separated by spaces (`0.8 0.9`; a comma before the
  space is fine, `0.8, 0.9`, but `0.8,0.9` is refused, since `2,500` could be a
  thousands separator). It can also be a range: `10 to 100 by 10`, or `10 to 100`
  with step 1; `10-20` is refused rather than read as 10 and −20. The minus sign
  `−` and the spaces of pasted text (no-break, thin) are read like `-` and a
  space. A list has at most 10,000 values and a text at most 1,000,000
  characters. Counts (whole-number inputs) go up to 1,000,000,000. Every other
  number (each value of a list too) is 0 or between 1e-30 and 1e30 in size (of
  either sign where the input's own limits allow negative values): no study's
  quantity lies beyond in any sensible unit, and a value beyond is refused with
  a message in the terms of the input's limits (`values must be between 1e-30
  and 1e30`, `values must be 0 or at least 1e-30` ...). A vector input
  (for example group means `0 0.5 1`) can hold several alternatives separated
  by `;`.
- A request nested deeper than 64 levels of objects and arrays is refused.
- Every combination of list values is one **scenario**. The last parameter
  varies fastest.
- `options` (all optional):
  - `threads`: 0 means all cores.
  - `power_curve`: compute power as a function of the sample size for up to
    eight distinct designs.
  - `curve_points`: number of points per curve.
  - `max_scenarios`: a lower limit on the number of scenarios than the engine's
    own (20,000); it cannot raise it. Group-sequential requests count by the
    time they take: a sample-size search counts 1 + K/2 scenarios for K
    analyses, and each distinct design 1 + K/8 more (1 + K/2 with β-spending
    futility).
  - `approximate`: for a simulation procedure, return its quick analytic
    approximation instead of simulating (live previews use it); other
    procedures ignore it.

## Result

```json
{
  "ok": true,
  "procedure": "means.two.t", "name": "Two-Sample T-Tests Assuming Equal Variance",
  "solve_for": "n", "solve_label": "Sample Size",
  "scenario_count": 2, "threads": 4, "elapsed_ms": 0.8,
  "columns": [{"key": "power", "label": "Actual Power", "fmt": "prob", "tip": "…"}, …],
  "rows": [
    {"n1": 17, "n2": 17, "n": 34, "target_power": 0.8, "power": 0.807036…,
     "_statement": "Group sample sizes of 17 and 17 …"},
    …
  ],
  "dropout": {"columns": […], "rows": […]},
  "plots": [
    {"kind": "result", "title": "…", "x_label": "…", "y_label": "…", "x_key": "mu1", "y_key": "n",
     "series": [{"name": "Target Power = 0.8", "x": […], "y": […]}]},
    {"kind": "curve", "title": "Power vs N1 (Group 1 Sample Size)", …,
     "series": [{"name": "…", "x": […], "y": […], "markers": [[17, 0.807]]}]}
  ],
  "definitions": [["N1, N2", "Number of subjects in groups 1 and 2; N = N1 + N2."], …],
  "references": ["…"],
  "notes": ["…"]
}
```

- **`columns`** lists the result table's columns in order. `fmt` is one of:

  | `fmt` | Meaning |
  |---|---|
  | `prob` | power and other probabilities |
  | `int` | sample sizes |
  | `num` | general numbers |
  | `alpha` | significance levels |
  | `pct` | a fraction shown as a percentage |
  | `ratio` | ratios |
  | `text` | text |

- **`rows`** holds one object per scenario. Keys are column keys, and values are
  numbers or text. Some keys are reserved:

  | Key | Meaning |
  |---|---|
  | `_error` | why this scenario could not be computed; the others are unaffected |
  | `_statement` | a one-paragraph summary of the scenario |
  | `_notes` | remarks specific to this scenario |
  | `_tables` | extra tables, for example group-sequential boundaries: `{"title", "columns": [{key, label, fmt}], "rows": [[…], …]}`; the first 1,000 rows carry them |

- **Simulation procedures** (`"simulation": true` in the result and the
  schema) estimate power from simulated studies. Their rows add `method`
  (`"simulation"` or `"approximation"`), and for a simulation `power_lo` and
  `power_hi` (95% Wilson interval), `power_se` (Monte Carlo standard error),
  `approx_power` (the analytic approximation), `reps`, `seed`, `failures`
  (studies whose analysis failed, as when its model could not be fitted,
  counted as not rejecting; when none could be analyzed, the row has an
  `_error` instead of a power, and so it has when more than half could not
  be in the simulations whose analyses fail only when the arithmetic breaks
  down: the longitudinal, cluster and stepped-wedge mixed models and the
  two-sample tests on skewed or heavy-tailed data. The MMRM, GLMM, ordinal
  and win-ratio simulations count failures that the design brings, such as a
  group with no patients left at a visit, apart from those of the
  arithmetic, which alone are held to that rule, and say why in
  `failures_why`, as their statements do) and `mean_df`. Those of them that
  run other tests on the same studies give each one's power (`power_mmrm`,
  `power_wald`, …; `null` when that test could not be computed in any study,
  or its arithmetic broke down in more than half), and for a test other than
  the one chosen whose analysis failed in some studies, their number and
  reasons in `<key>_failures` and `<key>_why` (`power_wald_failures`). For
  the ordinal Wald test an approximation gives `separated`, the expected
  share of separated studies, when it is 1% or more. Columns
  without a value in any row are left out, so an approximation has no Monte
  Carlo columns. The result reports
  `"approximate": true` for an approximation. Results depend only on the seed:
  every simulated cluster or subject draws from its own random stream, so the
  thread count never changes them, and a larger study contains the units of a
  smaller one. The token's progress advances smoothly while a simulation runs.
- **`dropout`** appears for procedures with a dropout-rate input. It lists sample
  sizes inflated by 1/(1 − rate) and the expected number of dropouts.
- **`plots`** has two kinds of entry:

  | `kind` | What it plots |
  |---|---|
  | `result` | the solved quantity against a varied input, one series per combination of the other inputs |
  | `curve` | power against the sample size; markers show the solutions |

  `x_key` names the input along the x axis (for a `curve`, the sample-size
  input such as `n1` or `k1`; `n2` when N1 is given and N2 solved, `n1` when
  N2 is given, `ntotal` with a percentage allocation), and a `result` plot's
  `y_key` the solved quantity. A `curve` also has `x_with`, the choices to
  make along with that input to reproduce a point (`{"alloc": "sep"}` when one
  group's size was given alone), so that a click on the curve can set the form.

Inputs that the calculation did not use are reported rather than dropped:
`unknown_inputs` lists names that are no input of the procedure,
`ignored_inputs` those that the chosen options leave out (`n2` with equal
groups, `power` when solving for power), and `replaced_choices` maps a choice
not available in this mode to `{"given": ..., "used": ...}` (`alloc=sep` when
solving for the sample size). Each is present only when it has entries.

If the request cannot be run, `ok` is `false` and `error` describes the
problem. Invalid inputs also produce `input_errors`, a map from parameter id to
message. A cancelled calculation returns `"cancelled": true`.

## Procedure schema

`ssc_procedure_json` describes everything a client needs to build an input
form:

```json
{
  "id": "means.two.t", "name": "…", "category": "Means", "subcategory": "Two Independent Means",
  "summary": "…", "keywords": ["…"], "simulation": false,
  "solve_for": [{"value": "n", "label": "Sample Size", "ykey": "n"}, {"value": "power", "label": "Power", "ykey": "power"}, …],
  "groups": ["Test", "Power and Alpha", "Sample Size", "Effect Size"],
  "params": [
    {"id": "alpha", "label": "Alpha (Significance Level)", "symbol": "Alpha", "group": "Power and Alpha",
     "help": "…", "type": "number", "default": "0.05", "multi": true, "integer": false,
     "min": 0, "max": 1, "min_open": true, "max_open": true,
     "options": [], "suggest": ["0.05", "0.01", …], "visible": []},
    …
  ],
  "example": {"solve_for": "n", "inputs": {…}},
  "definitions": […], "references": […], "notes": […]
}
```

- `ykey` of a `solve_for` option is the result row's key of the solved
  quantity (for example `k1` when solving for the number of clusters).
- `type` is one of:

  | `type` | Input |
  |---|---|
  | `number` | one value, or a list or range if `multi` is true |
  | `select` | one of `options` |
  | `bool` | true or false |
  | `vector` | a list of numbers, such as group means |

- `visible` (on a parameter or a select option) says when it applies. It is a
  list of alternatives; each alternative is a list of clauses that must all
  hold. A clause `{"key": "alternative", "values": ["upper", "lower"],
  "negate": false}` holds when the select parameter `alternative` has one of
  those values. The key `$solve` refers to the solve-for choice. An empty list
  means always visible.

## Search

`ssc_search_json` (`ssc.search` in Python, `ssc-calc search` on the command
line) finds procedures from a few words, the way the web page's search box
does:

```json
{"query": "equivalance two means", "corrected": "equivalence two means", "partial": false,
 "hits": [{"id": "means.two.equiv", "name": "Equivalence Tests for Two Means (TOST)",
           "category": "Means", "score": 8.68, "terms": ["equivalence", "two", "means"]}, ...]}
```

- **What is indexed.** Each procedure's name, keywords, topic, id, parameter
  and option labels, summary, definitions and references. A word counts more
  in the name and keywords than in the summary or references.
- **Normalization.** Text is lowercased; accents are removed, Greek letters
  spelled out (κ → kappa, χ² → chi2), plurals reduced (analyses → analysis)
  and common words dropped. Hyphenated words are indexed apart and joined
  (log-rank → log, rank, logrank), so "log rank", "log-rank" and "logrank"
  all match.
- **Synonyms and abbreviations.** About 45 groups of equivalent phrases
  (non-inferiority, NI; confidence interval, CI, margin of error; Mann–Whitney,
  rank sum, U test; ...) match each other.
- **Inexact words.** A word that is not in the index matches longer words it
  begins (still being typed: "bioequ") and, from four letters on, words
  within one edit (two from eight letters on; a swap of neighboring letters
  counts as one): "mcnemmar" finds McNemar and sets `corrected`. Candidates
  come from a trigram index over the vocabulary, so a typo never scans every
  word, and each candidate's edit distance is computed with an early exit.
- **Ranking.** BM25 over the weighted fields, plus a bonus when the query reads
  as a phrase of a procedure's name; a procedure id typed exactly comes first.
  Every distinctive word must match; if no procedure has them all, the
  procedures matching most of them are returned with `partial` set.
- **Long text.** Only the first 2,000 bytes of a query are read (cut at a
  character boundary), and only its first 64 words; `query` gives the text
  as read.

The index is built once, on first use, in about 20 ms. A query then takes a
few hundredths of a millisecond, including the call from Python. The relevance benchmark in
`core/tests/fixtures/search.json` (about 150 queries a user might type, with
the procedure expected first or among the first few) is part of the C++ tests.

## Python

```python
import ssc

ssc.catalog()                                   # dict: categories, procedures
ssc.schema("prop.two")                          # dict: the procedure schema
res = ssc.run("prop.two", "power", {"n1": "20 to 200 by 20", "p1": 0.6, "p2": 0.4},
              threads=0, power_curve=True, curve_points=41, approximate=False, token=None)
ssc.evaluate("pbinom", 3, 10, 0.5, 1)           # distribution functions
ssc.search("logrnk survival", limit=5)          # ranked procedures (see Search)

with ssc.Token() as token:                      # cancel from another thread:
    ...                                         #   token.cancel(); token.progress()
```

Lists and numbers in `inputs` are converted to text automatically.
`ssc.run` releases the GIL, so several calculations can run in parallel
threads.
