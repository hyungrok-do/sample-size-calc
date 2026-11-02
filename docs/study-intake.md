# Study intake schema (draft)

Status: design draft for discussion. Nothing here is implemented yet.

The DAG simulator asks people to think like the simulator: draw nodes, connect
arrows, then type coefficients on predictor scales. Clustered and repeated
designs are hard to express that way, because a DAG says which variables cause
which but says nothing about how units are nested, how time is indexed or where
correlation comes from.

This document proposes a different front door. The user first describes the
study in their own terms (what is randomized, how many visits, what the
outcome looks like in the control arm, what difference matters) and an agent
fills a structured **study intake** from that description. The intake is
checked for sufficiency: the agent says exactly what is still missing and does
not run until every required field is either given by the user or proposed by
the agent and confirmed. The intake then compiles into the existing engine
inputs (the `{version, graph, estimand, power}` design envelope, or a
calculator procedure), and the generated DAG is shown as a readable view of
what was assumed rather than the thing the user has to draw.

```
user's description + data summaries
        │  agent asks only for what is missing
        ▼
study intake (this schema)  ──  sufficiency check: missing / assumed / given
        │  compile
        ├──► calculator procedure (closed form or specialized simulation)
        └──► DAG design envelope {graph + study, estimand, power} ──► ssc_dag_json power
```

## Principles

1. **Study terms, not engine terms.** The user gives ICCs, control-arm
   prevalences, SDs, dropout by the last visit and effects on the outcome
   scale. Converting these to random-effect SDs, intercepts and coefficients is
   the compiler's job, and every conversion is shown.
2. **Every field has a status.** `missing` (blocks the run), `assumed` (the
   agent proposed a default; listed for confirmation before the run) or
   `given` (the user stated it or it came from their data). Assumed values are
   recorded in `assumed` as JSON pointers.
3. **Requirements depend on the design.** A field can be required for one
   design type and meaningless for another. The table in
   [Sufficiency](#sufficiency) is the single source of truth for that.
4. **Sources are kept.** Pilot data, registries and papers the user cites are
   recorded in `sources`, and a value can point to the source it came from.
   The engine never reads raw data; the intake stores the summaries.
5. **The intake is the reproducible artifact.** It is exported with the
   result next to the compiled envelope, so a reviewer can read the
   assumptions without reading a graph.

## Layers

The intake is one JSON object with eight layers. The machine-readable shape is
in [`study-intake.schema.json`](study-intake.schema.json); conditional
requirements are listed under [Sufficiency](#sufficiency) because they are
easier to check, and to explain to a user, in code than in JSON Schema.

### 1. `study` and `design`: the question and the design type

| Field | Type | Notes |
| --- | --- | --- |
| `study.title` | text ≤ 200 | |
| `study.question` | text ≤ 3000 | The scientific question in one or two sentences. Becomes `estimand.definition`. |
| `study.population` | text ≤ 500 | Who the effect averages over. |
| `sources[]` | `{id, kind, description, n?}` | `kind`: `pilot`, `prior_trial`, `registry`, `ehr`, `literature`, `expert`. |
| `design.type` | enum | `parallel_rct`, `cluster_rct`, `longitudinal_rct`, `cluster_longitudinal_rct`, `stepped_wedge`, `crossover`, `observational`. |
| `design.solveFor` | enum | `sample_size` (default), `power`, `detectable_effect`. |
| `design.randomizationLevel` | enum | `subject` or `cluster`. Implied by the type for every RCT design. |
| `design.allocation` | `{ratio: [a, b], scheme}` | `ratio` default `[1, 1]`; `scheme` `simple` (independent coin flips) or `fixed` (exact split). |

The design type fixes the unit levels: `parallel_rct` and `observational` are
subject only (observational may add clusters or visits), `cluster_rct` is
cluster > subject, `longitudinal_rct` is subject > visit,
`cluster_longitudinal_rct` is cluster > subject > visit, `stepped_wedge` is
cluster > period > subject, and `crossover` is subject > period.

### 2. `units`: nesting, sizes and attrition

| Field | Type | Notes |
| --- | --- | --- |
| `units.cluster.label` | text | "practice", "school", "ward". |
| `units.cluster.size` | one of `{fixed}`, `{pattern}`, `{mean, cv}` | Subjects per cluster. `{mean, cv}` compiles to a deterministic repeating pattern of gamma quantiles. |
| `units.cluster.count` | integer list | Candidate numbers of clusters (only when solving for power, or to fix the grid). |
| `units.subject.label` | text | "patient", "student". |
| `units.subject.count` | integer list | Candidate numbers of subjects for non-clustered designs. |
| `units.attrition.mode` | enum | `none`, `dropout` (monotone), `intermittent`, `missing_outcome` (single measurement). |
| `units.attrition.rate` | number in [0, 0.9) | Proportion without the primary outcome at the last scheduled visit. |
| `units.attrition.dependsOn` | enum | `nothing` (default), `treatment`, `observed_history`, `current_value` (MNAR sensitivity). |
| `units.attrition.strength` | odds ratio > 0 | Required unless `dependsOn` is `nothing`: how much the named factor changes the odds of missing a visit. |

### 3. `schedule`: timing

| Field | Type | Notes |
| --- | --- | --- |
| `schedule.unit` | enum | `days`, `weeks`, `months`, `years`. |
| `schedule.times` | 1–12 increasing numbers | Measurement times. A single time is a cross-sectional endpoint. |
| `schedule.baselineBeforeTreatment` | boolean | Whether the first time is measured before treatment starts. Decides whether a baseline difference of zero is imposed. |
| `schedule.periods` | integer | Stepped-wedge and crossover only: number of periods. |
| `schedule.followUp` | number | Time-to-event outcomes: planned follow-up (administrative censoring). |

### 4. `correlation`: where dependence comes from

| Field | Type | Notes |
| --- | --- | --- |
| `correlation.icc` | number in [0, 1) | Intra-cluster correlation of the primary outcome. |
| `correlation.iccScale` | enum | `observed` or `latent`. Binary outcomes default to `latent` (logit scale with residual variance π²/3); continuous outcomes are on the outcome scale. |
| `correlation.withinSubject.structure` | enum | `random_intercept` (default), `random_intercept_slope`, `exponential` (correlation decays with the time gap), `intercept_plus_exponential`. |
| `correlation.withinSubject.correlation` | number in [0, 1) | Correlation between two measurements of the same subject (for `exponential`, between measurements one time unit apart). |
| `correlation.withinSubject.slopeSD` | number ≥ 0 | Between-subject SD of the per-unit-time slope, on the outcome scale. |
| `correlation.withinSubject.interceptSlopeCorrelation` | number in (−1, 1) | |
| `correlation.clusterAutocorrelation` | number in (0, 1] | Stepped-wedge only: correlation of cluster effects between periods. |

### 5. `variables`: what is measured

One entry per variable the user cares about. Roles exist only in the intake;
the engine still treats every node generically. Roles drive the compiler's
defaults (a treatment in an RCT becomes an unconfounded root, a confounder gets
arrows into both exposure and outcome).

| Field | Type | Notes |
| --- | --- | --- |
| `id` | identifier | Same rules as the engine: ≤ 32 letters, digits or underscores, starting with a letter. |
| `label` | text ≤ 80 | |
| `role` | enum | `treatment`, `exposure`, `outcome`, `covariate`, `confounder`, `mediator`, `effect_modifier`, `censoring`. |
| `type` | enum | `continuous`, `binary`, `count`, `time_to_event`. |
| `level` | enum | `cluster`, `subject`, `visit`. Baseline covariates are `subject` (or `cluster`); a repeated outcome is `visit`. |
| `distribution` | by type | See below. For the outcome these are **reference-arm** (control or unexposed) values at the primary time. |
| `source` | source id | Optional link to `sources[]`. |

Distribution by type:

- `continuous`: `{mean, sd, shape}`; `shape` is `symmetric` (Normal, default), `right_skewed` (log-normal) or `positive_skewed` (Gamma).
- `binary`: `{prevalence}`; for treatment or exposure this is the share in the active group, which also follows from `design.allocation`.
- `count`: `{mean}` (Poisson).
- `time_to_event`: `{median}` or `{survivalAt: {time, probability}}`, plus `hazardShape` (`constant` default, `increasing`, `decreasing`) and `loss` (a proportion lost to follow-up by `schedule.followUp`).

### 6. `relationships` and `effect`: the within-slice DAG and the primary effect

`relationships[]` lists the arrows other than the primary effect, each as
`{from, to, strength}` on a scale the user can reason about:

| `to` type | Accepted `strength` |
| --- | --- |
| binary | `{oddsRatio, per}` (per one unit, or per `per` units of a continuous `from`) |
| continuous | `{difference, per}` or `{standardized}` (outcome SDs per `from` SD) |
| count | `{rateRatio, per}` |
| time-to-event | `{hazardRatio, per}` |
| any | `weak`, `moderate`, `strong`: allowed only as an agent proposal, compiled to explicit numbers and always listed in `assumed` |

Optional extras: `lags[]` (`{variable, lag, strength}` for a previous visit's
value) and `modifiers[]` (`{by, ratioOfEffects | differenceOfEffects}` for
treatment-by-covariate or treatment-by-time heterogeneity).

`effect` is the primary comparison:

| Field | Type | Notes |
| --- | --- | --- |
| `effect.compare` | variable id | The treatment or exposure. |
| `effect.outcome` | variable id | |
| `effect.measure` | enum | `mean_difference`, `risk_difference`, `risk_ratio`, `odds_ratio`, `hazard_ratio`, `survival_difference`, `rmst_difference`. |
| `effect.value` | number (or list for a grid) | Required unless solving for the detectable effect. |
| `effect.direction` | enum | `increase`, `decrease`, `either`. Sets the sign check and, with `power.sides`, the rejection region. |
| `effect.timing` | object | Repeated outcomes: `contrast` (`final`, `visit`, `average`, `change`, `slope`, `custom`), `visit` when `contrast` is `visit`, `weights` when it is `custom`, and `shape`: how the effect builds up over visits (`immediate`, `linear`, `final_only`). Time-to-event: `horizon` for survival or RMST differences. |
| `effect.justification` | text | Minimal clinically important difference, pilot estimate, prior trial. Optional but prompted. |

### 7. `analysis` and `power`

| Field | Type | Notes |
| --- | --- | --- |
| `analysis.estimator` | enum | `auto` (default) or a named estimator: `welch`, `two_proportion`, `logrank`, `cox`, `regression`, `mmrm`, `lmm`, `trajectory_means`, `glmm`, `aipw`, `tmle`, `gformula`. |
| `analysis.adjustFor` | variable ids | Baseline covariates in the planned analysis. |
| `power.alpha` | (0, 0.5) | Default 0.05. |
| `power.sides` | enum | `two-sided` (default), `one-sided`. |
| `power.target` | (alpha, 1) | Default 0.8. |
| `power.sampleSizes` | integer list | Candidate N (clusters when clustered). The agent proposes a grid around a closed-form approximation. |
| `power.repetitions` | 10–10,000 | Default 1,000 for a final answer, 200 for a quick look. |
| `power.seed` | uint32 | Default 2026. |

### 8. Bookkeeping

`schemaVersion` (1), `assumed[]` (JSON pointers the agent filled and the user
has not confirmed) and `notes` (free text the compiler copies into the report).

## Sufficiency

Legend: **R** required, **C** required under the stated condition, **D** the
agent may propose a default (recorded in `assumed`), **–** not used.

| Field | parallel RCT | longitudinal RCT | cluster RCT | cluster longitudinal | stepped wedge | crossover | observational |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `study.question` | R | R | R | R | R | R | R |
| `design.solveFor` | D | D | D | D | D | D | D |
| `design.allocation` | D | D | D | D | – | D | – |
| `units.cluster.size` | – | – | R | R | R | – | C: if clustered |
| `units.attrition.mode`, `rate` | R | R | R | R | R | R | R |
| `schedule.times` | D: one time | R, ≥ 2 | D: one time | R, ≥ 2 | – | – | C: if repeated |
| `schedule.baselineBeforeTreatment` | – | R | – | R | – | – | C: if repeated |
| `schedule.periods` | – | – | – | – | R | R | – |
| `schedule.followUp` | C: time-to-event | – | C: time-to-event | – | – | – | C: time-to-event |
| `correlation.icc` | – | – | R | R | R | – | C: if clustered |
| `correlation.withinSubject.correlation` | – | R | – | R | – | R | C: if repeated |
| `correlation.withinSubject.structure` | – | D | – | D | – | – | C: if repeated |
| `correlation.clusterAutocorrelation` | – | – | – | – | D: 1 | – | – |
| outcome `distribution` | R | R | R | R | R | R | R |
| exposure `distribution.prevalence` | – | – | – | – | – | – | R |
| confounders with both arrows | – | – | – | – | – | – | R: at least one, or an explicit "none measured" |
| `effect.measure` | R | R | R | R | R | R | R |
| `effect.value` | C: unless detectable effect | C | C | C | C | C | C |
| `effect.timing.contrast`, `shape` | – | R | – | R | – | – | C: if repeated |
| `effect.timing.horizon` | C: survival/RMST measure | – | C | – | – | – | C |
| `analysis.estimator` | D | D | D | D | D | D | D |
| `power.alpha`, `sides`, `target` | D | D | D | D | D | D | D |
| `power.sampleSizes` | C: if solving for power | C | C | C | C | C | C |

Validation rules, checked before the sufficiency report:

- Probabilities in (0, 1); prevalences in (0.001, 0.999); SDs, medians and
  follow-up > 0; ratios > 0; ICC < the within-subject correlation when both are
  given (a subject's two visits cannot be less alike than two members of the
  same cluster).
- A risk difference keeps both arms' risks in (0, 1) given the reference
  prevalence; a mean difference larger than 3 outcome SDs is flagged as
  implausible but allowed.
- `effect.direction` agrees with the sign of `effect.value` (or the ratio side
  of 1).
- `effect.compare` is binary for every RCT design, and baseline
  (subject or cluster level) for every repeated design. A visit-level exposure
  is a time-varying treatment and is rejected with that explanation.
- Every `relationships[]` endpoint exists; arrows run from coarser to finer
  levels or within a level; the graph within a slice is acyclic.
- In an RCT, no arrow enters the treatment. In an observational design, at
  least one confounder has arrows into both the exposure and the outcome, or the
  user has stated that none is measured (then the run shows the unadjusted bias
  rather than refusing).
- Mediators are never put in `analysis.adjustFor` for a total effect.

### The missing-information report

The agent turns the check into at most three questions at a time, ordered by
how much each answer moves power: the effect size first, then outcome
variability (SD or reference prevalence), then ICC or within-subject
correlation, then attrition, then everything else. Each question names the
default it would use, so the user can answer with "fine":

```json
{
  "ready": false,
  "missing": [
    {"path": "/correlation/icc", "ask": "How similar are patients within the same practice? An ICC of 0.01–0.05 is typical for primary-care process outcomes.", "proposal": 0.03},
    {"path": "/units/attrition/rate", "ask": "What share of patients will have no 6-month outcome?", "proposal": 0.1}
  ],
  "assumed": ["/power/alpha", "/power/sides", "/analysis/estimator"],
  "invalid": []
}
```

The run starts only when `missing` and `invalid` are empty and the user has
seen the `assumed` list once.

## Intake template (what the user sees)

The agent opens with an invitation to describe the study and any data in their
own words, then works through the parts below, skipping questions the
description already answered and parts the design type does not use.

1. **Your question.** What are you comparing, in whom, and on what outcome?
   Do you have pilot data, a prior trial or published estimates to base the
   numbers on?
2. **Design.** Is the comparison randomized? If so, what is randomized: people,
   or groups such as practices, schools or wards? Does everyone get both
   treatments in sequence (crossover), or do groups switch from control to
   intervention at staggered times (stepped wedge)? Equal allocation?
3. **Groups** (clustered designs). What is a group? About how many people per
   group, and how much does that vary? How similar are outcomes within a group
   (ICC)?
4. **Timing.** When is the outcome measured? Is the first measurement before
   treatment starts? For time to an event: how long is follow-up?
5. **Repeated measurements.** How correlated are two measurements of the same
   person? Do people differ in how fast they change?
6. **The outcome in the control group.** Continuous: mean, SD and whether it is
   skewed. Yes/no: the proportion with the outcome. Count: the mean count. Time
   to event: the median time, or the proportion event-free at a time.
7. **The effect you want to detect.** On the outcome's scale (a difference in
   means or risks, a ratio). Which direction counts as benefit? For repeated
   outcomes: at which visit, or averaged over visits, and does the effect
   appear at once or build up? Why this size?
8. **Other variables** (observational, or adjusted RCT analyses). Which factors
   affect both who is exposed and the outcome? For each: how common or how
   variable, and how strongly it relates to exposure and to outcome?
9. **Missing data.** What share will have no outcome at the final visit? Could
   that depend on treatment or on how the person is doing?
10. **Analysis and error rates.** A planned analysis if you have one; otherwise
    I will propose one. Alpha, one- or two-sided, target power.

Before running, the agent shows the completed intake as a short table, the
`assumed` list, the compiled DAG with its parameter table, and the analysis it
will fit.

## Worked examples

### A. Cluster RCT with a binary outcome

"Primary-care practices randomized 1:1 to a BP-control program, about 20
patients per practice (varies a lot), outcome is controlled BP at 6 months: 40%
in usual care, we hope for 50%. ICC around 0.05. Maybe 10% have no 6-month
reading."

```json
{
  "schemaVersion": 1,
  "study": {"title": "BP control program", "question": "Does the program raise the proportion of patients with controlled BP at 6 months?"},
  "design": {"type": "cluster_rct", "solveFor": "sample_size", "randomizationLevel": "cluster", "allocation": {"ratio": [1, 1], "scheme": "simple"}},
  "units": {
    "cluster": {"label": "practice", "size": {"mean": 20, "cv": 0.4}},
    "subject": {"label": "patient"},
    "attrition": {"mode": "missing_outcome", "rate": 0.1, "dependsOn": "nothing"}
  },
  "schedule": {"unit": "months", "times": [6]},
  "correlation": {"icc": 0.05, "iccScale": "latent"},
  "variables": [
    {"id": "Program", "role": "treatment", "type": "binary", "level": "cluster", "distribution": {"prevalence": 0.5}},
    {"id": "BPcontrol", "role": "outcome", "type": "binary", "level": "visit", "distribution": {"prevalence": 0.4}}
  ],
  "effect": {"compare": "Program", "outcome": "BPcontrol", "measure": "risk_difference", "value": 0.1, "direction": "increase"},
  "analysis": {"estimator": "auto"},
  "power": {"alpha": 0.05, "sides": "two-sided", "target": 0.8},
  "assumed": ["/units/cluster/size/cv", "/correlation/iccScale", "/analysis/estimator", "/power/alpha", "/power/sides", "/power/target"]
}
```

Compiles to: `study.times = [6]`; `study.clusterSizes` = a repeating pattern
with mean 20 and CV 0.4 (for example `[10, 14, 18, 24, 34]`), shrunk by 10% for
the missing outcomes; `Program` at cluster scope; `BPcontrol` at visit scope
with `repeated.clusterSD = sqrt(0.05 / 0.95 × π²/3) = 0.416`; the
`BPcontrol` intercept calibrated so the marginal usual-care risk is 0.40; the
`Program → BPcontrol` arrow calibrated to a marginal risk difference of 0.10
(`estimand.alternative`); trajectory means with a cluster-level CR1 sandwich.
The agent also offers `prop.two.cluster` as a closed-form cross-check.

### B. Longitudinal RCT with a continuous outcome

"Depression score at baseline, 3 and 6 months, SD about 6. Treatment should
lower the 6-month score by 2.5 points, building up over time. Visits correlate
around 0.6. 20% dropout by month 6."

```json
{
  "schemaVersion": 1,
  "study": {"title": "Depression program", "question": "Does the program lower the depression score at 6 months?"},
  "design": {"type": "longitudinal_rct", "solveFor": "sample_size", "randomizationLevel": "subject", "allocation": {"ratio": [1, 1], "scheme": "simple"}},
  "units": {"subject": {"label": "participant"}, "attrition": {"mode": "dropout", "rate": 0.2, "dependsOn": "nothing"}},
  "schedule": {"unit": "months", "times": [0, 3, 6], "baselineBeforeTreatment": true},
  "correlation": {"withinSubject": {"structure": "random_intercept", "correlation": 0.6}},
  "variables": [
    {"id": "Treat", "role": "treatment", "type": "binary", "level": "subject", "distribution": {"prevalence": 0.5}},
    {"id": "Score", "role": "outcome", "type": "continuous", "level": "visit", "distribution": {"mean": 15, "sd": 6, "shape": "symmetric"}}
  ],
  "effect": {"compare": "Treat", "outcome": "Score", "measure": "mean_difference", "value": -2.5, "direction": "decrease",
             "timing": {"contrast": "final", "shape": "linear"}, "justification": "Half the minimal important difference reported for the scale."},
  "analysis": {"estimator": "mmrm"},
  "power": {"alpha": 0.05, "sides": "two-sided", "target": 0.9},
  "assumed": ["/variables/1/distribution/mean", "/correlation/withinSubject/structure"]
}
```

Compiles to: `study.times = [0, 3, 6]`; `Score` at visit scope, Normal with
`repeated.subjectSD = sqrt(0.6 × 36) = 4.65` and residual `sd = sqrt(0.4 × 36) =
3.79`; a time node and a `Treat × time` function with coefficient
`−2.5 / 6 = −0.417` per month (zero difference at baseline, linear build-up);
`study.observation = {mode: "dropout", intercept: logit(1 − sqrt(0.8)) = −2.14}`
(two post-baseline visits); `estimand.repeated.contrast = "final"`; MMRM.

### C. Observational study with a confounded binary exposure

"EHR cohort. About 30% on the drug. One-year MACE is 8% in non-users. Older
patients and diabetics are more likely to get the drug (OR 1.5 per 10 years,
OR 2.0) and to have MACE (OR 1.6 per 10 years, OR 1.8). Age mean 62, SD 10;
25% diabetic. We want to detect RR 0.75."

```json
{
  "schemaVersion": 1,
  "study": {"title": "Drug and one-year MACE", "question": "What is the effect of the drug on one-year MACE risk, as a risk ratio?"},
  "sources": [{"id": "ehr", "kind": "ehr", "description": "Health-system EHR extract, 2018–2024", "n": 12000}],
  "design": {"type": "observational", "solveFor": "sample_size"},
  "units": {"subject": {"label": "patient"}, "attrition": {"mode": "missing_outcome", "rate": 0, "dependsOn": "nothing"}},
  "schedule": {"unit": "years", "times": [1]},
  "variables": [
    {"id": "Age", "role": "confounder", "type": "continuous", "level": "subject", "distribution": {"mean": 62, "sd": 10}, "source": "ehr"},
    {"id": "Diabetes", "role": "confounder", "type": "binary", "level": "subject", "distribution": {"prevalence": 0.25}, "source": "ehr"},
    {"id": "Drug", "role": "exposure", "type": "binary", "level": "subject", "distribution": {"prevalence": 0.3}, "source": "ehr"},
    {"id": "MACE", "role": "outcome", "type": "binary", "level": "subject", "distribution": {"prevalence": 0.08}, "source": "ehr"}
  ],
  "relationships": [
    {"from": "Age", "to": "Drug", "strength": {"oddsRatio": 1.5, "per": 10}},
    {"from": "Diabetes", "to": "Drug", "strength": {"oddsRatio": 2.0}},
    {"from": "Age", "to": "MACE", "strength": {"oddsRatio": 1.6, "per": 10}},
    {"from": "Diabetes", "to": "MACE", "strength": {"oddsRatio": 1.8}}
  ],
  "effect": {"compare": "Drug", "outcome": "MACE", "measure": "risk_ratio", "value": 0.75, "direction": "decrease"},
  "analysis": {"estimator": "aipw", "adjustFor": ["Age", "Diabetes"]},
  "power": {"alpha": 0.05, "sides": "two-sided", "target": 0.8},
  "assumed": ["/power/alpha", "/power/sides", "/power/target"]
}
```

Compiles to: a cross-sectional graph (no `study` block) with `Age` Normal(62,
10), `Diabetes` Bernoulli, arrows `Age → Drug` = log(1.5)/10 = 0.0405 per year,
`Diabetes → Drug` = log 2 = 0.693, `Age → MACE` = 0.0470, `Diabetes → MACE` =
0.588; `Drug` and `MACE` intercepts calibrated to marginal prevalences 0.30 and
0.08 (unexposed); `Drug → MACE` calibrated to a marginal risk ratio of 0.75;
AIPW.

## What the engine supports today

| Layer | Current engine field or procedure | Status |
| --- | --- | --- |
| Design type: parallel, cluster, longitudinal, cluster longitudinal, observational | `graph`, `graph.study`, node `scope` | Supported |
| Design type: stepped wedge, crossover | `cluster.sw`, `sim.sw`, `means.crossover` (calculator) | Supported as calculator procedures only; the DAG engine has no time-varying treatment |
| Allocation ratio, simple | Bernoulli treatment root, intercept `logit(share)` | Supported |
| Allocation, fixed split | none | **New**: fixed or blocked assignment in `study.cpp` |
| Cluster sizes, fixed or pattern | `study.clusterSizes` (≤ 100 entries, 1–1,000) | Supported |
| Cluster sizes, `{mean, cv}` | compiled to a pattern | **New** (compiler only) |
| Dropout, intermittent missingness, MNAR | `study.observation` | Supported; the compiler converts "rate by last visit" to the per-visit intercept |
| Missing outcome at a single measurement | none (first visit always observed) | **New**: either a first-visit missingness rule in the engine, or the compiler shrinks sizes (MCAR only) |
| Visit schedule | `study.times` (1–12) | Supported |
| Baseline before treatment | treatment × time function with zero at baseline | Supported via compiled functions |
| ICC | `repeated.clusterSD` on the predictor scale | **New** conversion: observed-scale ICC for binary outcomes needs a numerical solve; latent scale is closed form |
| Within-subject correlation, slopes, exponential decay | `repeated.subjectSD`, `slopeSD`, `slopeCorrelation`, `residualCorrelation` | Supported; conversion from a correlation is new compiler code |
| Variables, types, levels | nodes with `kind`, `family`, `scope` | Supported; roles are intake-only |
| Marginal (reference-arm) prevalence, mean or median | node `intercept`, `scale` | **New**: intercept calibration to a marginal target (the engine only calibrates arrows); can reuse the `population` operation |
| Relationship strengths (OR, HR, standardized) | arrow `coefficient` | Supported; conversion is compiler code |
| Primary effect on the marginal scale | `estimand.alternative` + `calibrationEdge` | Supported for mean, risk difference, risk ratio, odds ratio, survival difference, RMST |
| Hazard ratio as the primary effect | conditional coefficient target, Cox | Supported for a randomized comparison; marginal HR is not a target |
| Repeated contrasts | `estimand.repeated` (`final` with an optional `visit`, `average`, `change`, `slope`, `custom`) | Supported; the intake's `visit` compiles to `final` with that visit |
| Estimators | `power.analysis.method` (auto, AIPW, TMLE, g-formula, MMRM, LMM, trajectory means) | Supported; GLMM for clustered binary outcomes is calculator-only (`sim.glmm.cluster`) |
| Uncertain inputs as ranges | one envelope per scenario | **New**: the compiler expands a grid and runs it |

## Choices to settle

1. **Routing.** Send stepped-wedge, crossover and other designs the DAG engine
   cannot express to the calculator procedures (proposed), or extend the
   engine with time-varying treatment.
2. **Effect scale.** Elicit the effect on the marginal outcome scale and
   calibrate the arrow (proposed, and what the engine already does), or let
   advanced users give coefficients directly in the intake.
3. **Fixed allocation.** Add exact or blocked assignment to the engine, or
   keep Bernoulli assignment and say so in the report.
4. **ICC for binary outcomes.** Ask on the latent scale (closed form,
   proposed default) or the observed scale (more familiar to applied users,
   needs a numerical solve).
5. **Ranges.** Accept ranges for uncertain inputs (ICC 0.02–0.08) and run a
   small sensitivity grid by default, or keep single values.
6. **Intercept calibration.** Add a `calibrateIntercept` operation to
   `ssc_dag_json`, or let the agent iterate on `population`.
