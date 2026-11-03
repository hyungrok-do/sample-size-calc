# Calculator reference coverage audit

Inventory date: 2026-11-03. Source commit: `a24694320f69b1cb0caabed5ffbbd7f82cace209`.

**Status: bibliography import is blocked. This is a coverage inventory, not a completed citation verification or numerical-method audit.**

## Findings

- All 71 registered catalog procedures expose nonempty reference lists. These lists contain 125 distinct citation strings, copied verbatim below from the engine schema.
- The tracked tree at the source commit contains no `.bib` file, including `references.bib`. `CITATION.cff` cites the calculator software itself; it is not a bibliography of statistical methods.
- Catalog references are embedded as strings in C++ procedure definitions, serialized in the schema and result JSON, and rendered in the standard calculator and Word reports. There are no BibTeX keys or recorded citation-export provenance.
- The DAG and Lean engines are outside the 71-procedure registry. Their result objects do not expose corresponding bibliographic reference lists. The DAG documentation has some methodological links, but these do not establish citation coverage in exported results.
- Many catalog procedures contain selectable test, power-approximation, or interval variants. A nonempty procedure-level list does not establish that every variant is supported by the cited sources.
- A direct Google Scholar citation-export attempt was blocked by an unusual-traffic response. No Scholar BibTeX was retrieved. No bibliography entries have been synthesized or imported from another source.

## Integration points

| Surface | Current source | Required follow-up after the supplied bibliography is located |
| --- | --- | --- |
| Catalog schema and results | `core/src/engine/runner.cpp`, `core/src/engine/procedure.hpp`, `core/src/procedures/*.cpp` | Associate each procedure and relevant method variant with verified citation keys. |
| Standard calculator | `web/app.js` | Render records resolved from those keys; retain plain-text references for backward compatibility. |
| Word report | `python/ssc/docx.py` | Use the same resolved records as the engine and web result. |
| DAG analysis | `core/src/dag/power.cpp`, `core/src/dag/study.cpp` | Attach citations for the actual selected estimator, test, variance calculation, and nuisance learner. |
| Lean design plans | `core/src/design/plan.cpp` | Cite known-variance normal power separately from influence-function Wald power and the ATE EIF. Statistical references do not certify the Lean-to-C++ implementation. |

## Procedure-to-reference inventory

The `R` labels below are local inventory numbers, not BibTeX keys. Existing bibliographic text has not been treated as a retrieved record. Full references follow the table.

| Procedure | Implemented catalog name | Existing reference inventory | Source |
| --- | --- | --- | --- |
| `means.crossover` | Two-Treatment 2×2 Cross-Over Design (Difference of Means) | R001, R002, R003 | [means_designs.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_designs.cpp) |
| `means.be.ratio` | Bioequivalence: Ratio of Two Means (TOST, Log Scale) | R004, R005, R006, R001 | [means_designs.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_designs.cpp) |
| `means.anova` | One-Way Analysis of Variance (ANOVA) F-Test | R007, R008, R009, R010 | [means_designs.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_designs.cpp) |
| `means.anova2` | Two-Way Factorial ANOVA (Main Effects and Interaction) | R007, R011 | [means_more.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_more.cpp) |
| `means.one.t` | Tests for One Mean (One-Sample T-Test) | R001, R002, R012, R013, R010 | [means_t.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_t.cpp) |
| `means.one.ni` | Non-Inferiority Tests for One Mean (or Paired Means) | R001, R002, R014 | [means_ni_equiv.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_ni_equiv.cpp) |
| `means.one.equiv` | Equivalence Tests for One Mean (or Paired Means) | R001, R002, R004, R005 | [means_ni_equiv.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_ni_equiv.cpp) |
| `means.paired.t` | Tests for Paired Means (Paired T-Test) | R001, R002, R012 | [means_t.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_t.cpp) |
| `means.two.repeated` | Two Means in a Repeated-Measures Design (Time-Averaged Response) | R015, R016 | [means_designs.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_designs.cpp) |
| `means.rm.anova` | Repeated-Measures ANOVA (One Within-Subject Factor) | R007, R017 | [means_more.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_more.cpp) |
| `means.two.slopes` | Two Mean Slopes in a Longitudinal Design (Random Slopes) | R018, R015, R019 | [means_more.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_more.cpp) |
| `means.two.t` | Two-Sample T-Tests Assuming Equal Variance | R001, R002, R012, R013, R010 | [means_t.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_t.cpp) |
| `means.two.welch` | Two-Sample T-Tests Allowing Unequal Variance (Welch) | R001, R002, R012, R013, R010, R020 | [means_t.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_t.cpp) |
| `means.two.ni` | Non-Inferiority Tests for Two Means | R001, R002, R014 | [means_ni_equiv.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_ni_equiv.cpp) |
| `means.two.ancova` | Two Means Adjusted for Covariates (ANCOVA) | R021, R001 | [means_more.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_more.cpp) |
| `means.two.equiv` | Equivalence Tests for Two Means (TOST) | R001, R002, R004 | [means_ni_equiv.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_ni_equiv.cpp) |
| `means.coprimary` | Two Means with Co-Primary Endpoints (All Must Be Significant) | R022, R023, R024 | [means_coprimary.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_coprimary.cpp) |
| `prop.cc` | Unmatched Case-Control Study (Odds Ratio) | R025, R026, R027 | [proportions_two.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_two.cpp) |
| `prop.cc.matched` | Matched Case-Control Study (1:M Matching) | R028, R025, R029 | [proportions_two.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_two.cpp) |
| `prop.chisq` | Chi-Square Tests (Goodness-of-Fit and Contingency Tables) | R007, R010, R030 | [proportions_more.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_more.cpp) |
| `prop.one` | Tests for One Proportion | R001, R026, R012, R031 | [proportions_one.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_one.cpp) |
| `prop.paired.mcnemar` | Tests for Two Correlated Proportions (McNemar Test) | R032, R033, R026 | [proportions_more.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_more.cpp) |
| `prop.trend` | Cochran–Armitage Test for a Trend in Proportions | R034, R035 | [trend_stepped_wedge.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/trend_stepped_wedge.cpp) |
| `prop.one.simon` | Simon's Two-Stage Designs (Single-Arm Phase II) | R036, R037 | [phase2_simon.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/phase2_simon.cpp) |
| `prop.two` | Tests for Two Proportions | R026, R001, R012, R038, R039 | [proportions_two.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_two.cpp) |
| `prop.two.ni` | Non-Inferiority Tests for Two Proportions | R040, R041, R001, R002 | [proportions_more.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_more.cpp) |
| `prop.two.equiv` | Equivalence Tests for Two Proportions | R040, R041, R001 | [proportions_more.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_more.cpp) |
| `corr.one` | Tests for One Correlation (Pearson) | R042, R043, R007, R012 | [correlation_regression.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/correlation_regression.cpp) |
| `corr.two` | Tests for Two Independent Correlations | R007, R043 | [correlation_regression.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/correlation_regression.cpp) |
| `reg.multiple` | Multiple Regression (R² Increase F-Test) | R007, R044, R010 | [correlation_regression.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/correlation_regression.cpp) |
| `reg.slope` | Simple Linear Regression (Test of the Slope) | R045, R012 | [correlation_regression.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/correlation_regression.cpp) |
| `reg.logistic` | Logistic Regression (One Covariate of Interest) | R046, R047 | [correlation_regression.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/correlation_regression.cpp) |
| `reg.prediction` | Developing a Clinical Prediction Model (Minimum Sample Size) | R048, R049, R050, R051, R052, R053 | [prediction_models.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/prediction_models.cpp) |
| `reg.cox` | Cox Regression (One Covariate of Interest) | R054, R055 | [correlation_regression.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/correlation_regression.cpp) |
| `surv.logrank` | Logrank Tests (Exponential Survival, Accrual and Follow-Up) | R056, R057, R058, R012 | [survival_rates_variances.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/survival_rates_variances.cpp) |
| `surv.logrank.ni` | Non-Inferiority Logrank Tests (Hazard Ratio Margin) | R056, R059, R012 | [survival_rates_variances.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/survival_rates_variances.cpp) |
| `rates.one.poisson` | Tests for One Poisson Rate | R060, R061, R012 | [survival_rates_variances.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/survival_rates_variances.cpp) |
| `rates.two.poisson` | Tests for the Ratio of Two Poisson Rates | R061, R012 | [survival_rates_variances.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/survival_rates_variances.cpp) |
| `rates.two.negbin` | Tests for the Ratio of Two Negative Binomial Rates | R062, R063 | [survival_rates_variances.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/survival_rates_variances.cpp) |
| `var.one` | Tests for One Variance (Chi-Square Test) | R064, R012 | [survival_rates_variances.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/survival_rates_variances.cpp) |
| `var.two` | Tests for Two Variances (F-Test) | R008, R012 | [survival_rates_variances.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/survival_rates_variances.cpp) |
| `nonpar.signedrank` | Wilcoxon Signed-Rank Test (One Sample or Paired) | R065, R066 | [nonparametric_ci.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/nonparametric_ci.cpp) |
| `nonpar.sign` | Sign Test (Exact Binomial) | R067, R012 | [nonparametric_ci.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/nonparametric_ci.cpp) |
| `nonpar.kruskal` | Kruskal–Wallis Test | R065, R068 | [nonparametric_ci.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/nonparametric_ci.cpp) |
| `nonpar.mannwhitney` | Mann–Whitney U / Wilcoxon Rank-Sum Test | R067, R065 | [nonparametric_ci.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/nonparametric_ci.cpp) |
| `gs.logrank` | Group-Sequential Logrank Tests (Hazard Ratio) | R069, R070, R071, R072, R073, R074, R075, R076, R056 | [groupseq_procedures.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/groupseq_procedures.cpp) |
| `gs.means` | Group-Sequential Tests for Two Means | R069, R070, R071, R072, R073, R074, R075, R076 | [groupseq_procedures.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/groupseq_procedures.cpp) |
| `gs.props` | Group-Sequential Tests for Two Proportions | R069, R070, R071, R072, R073, R074, R075, R076 | [groupseq_procedures.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/groupseq_procedures.cpp) |
| `agree.icc` | Tests for an Intraclass Correlation (Reliability) | R077, R078 | [agreement_diagnostic.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/agreement_diagnostic.cpp) |
| `agree.kappa` | Tests for Cohen's Kappa (Two Raters, Binary Ratings) | R079, R080 | [agreement_diagnostic.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/agreement_diagnostic.cpp) |
| `diag.auc` | Tests for the Area Under a ROC Curve | R081, R082 | [agreement_diagnostic.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/agreement_diagnostic.cpp) |
| `diag.sens` | Tests for Sensitivity and Specificity (One Diagnostic Test) | R083, R084 | [agreement_diagnostic.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/agreement_diagnostic.cpp) |
| `diag.paired` | Comparing Two Diagnostic Tests on the Same Subjects (Paired Sensitivity and Specificity) | R032, R085, R084 | [agreement_diagnostic.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/agreement_diagnostic.cpp) |
| `corr.ci` | Confidence Interval for a Correlation | R012 | [correlation_regression.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/correlation_regression.cpp) |
| `means.one.ci` | Confidence Interval for One Mean | R086, R012 | [nonparametric_ci.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/nonparametric_ci.cpp) |
| `means.two.ci` | Confidence Interval for the Difference of Two Means | R086, R012 | [nonparametric_ci.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/nonparametric_ci.cpp) |
| `prop.one.ci` | Confidence Interval for One Proportion | R087, R088, R089, R026 | [proportions_one.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_one.cpp) |
| `prop.two.ci` | Confidence Interval for the Difference of Two Proportions | R090, R091, R026 | [proportions_two.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_two.cpp) |
| `means.two.cluster` | Two Means in a Cluster-Randomized Design | R092, R093, R012 | [means_designs.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/means_designs.cpp) |
| `prop.two.cluster` | Two Proportions in a Cluster-Randomized Design | R092, R094 | [proportions_more.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/proportions_more.cpp) |
| `cluster.sw` | Stepped-Wedge Cluster-Randomized Designs | R095, R096 | [trend_stepped_wedge.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/trend_stepped_wedge.cpp) |
| `sim.lmm.cluster` | Two Means in a Cluster-Randomized Design with a Mixed Model, by Simulation | R092, R097, R098, R099 | [sim_lmm.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_lmm.cpp) |
| `sim.glmm.cluster` | Two Proportions in a Cluster-Randomized Design with a Logistic Mixed Model, by Simulation | R100, R097, R101, R099 | [sim_glmm.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_glmm.cpp) |
| `sim.sw` | Stepped-Wedge Trial with a Mixed Model (Cluster-Period Effects), by Simulation | R095, R102, R103, R099 | [sim_stepped_wedge.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_stepped_wedge.cpp) |
| `sim.winratio` | Win Ratio for a Hierarchical Composite Endpoint (Death, Then a Non-Fatal Event, Then a Score), by Simulation | R104, R105, R106, R107, R099 | [sim_winratio.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_winratio.cpp) |
| `sim.lmm.longitudinal` | Two Mean Slopes with a Mixed Model, by Simulation | R108, R109, R098, R099, R110 | [sim_lmm.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_lmm.cpp) |
| `sim.mmrm` | Two Means at the Last Visit with Dropout: Mixed Model for Repeated Measures (MMRM), by Simulation | R111, R112, R113, R099 | [sim_mmrm.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_mmrm.cpp) |
| `sim.means.robust` | Two Means with Skewed or Heavy-Tailed Data (t, Welch, Mann–Whitney), by Simulation | R067, R114, R020, R099 | [sim_robust.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_robust.cpp) |
| `sim.ordinal` | Two Groups with an Ordinal Outcome (Proportional Odds Model or Mann–Whitney), by Simulation | R115, R116, R117, R099 | [sim_ordinal.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_ordinal.cpp) |
| `sim.multiarm` | Several Treatments against a Control with a Multiplicity Adjustment (Dunnett, Holm, Hochberg), by Simulation | R118, R119, R120, R121, R099 | [sim_multiarm.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_multiarm.cpp) |
| `sim.surv.nph` | Two Survival Curves with Non-Proportional Hazards, by Simulation | R122, R123, R124, R056, R125, R099 | [sim_survival.cpp](https://github.com/hyungrok-do/sample-size-calc/blob/a24694320f69b1cb0caabed5ffbbd7f82cace209/core/src/procedures/sim_survival.cpp) |

## Selectable test and approximation variants

These labels are copied from the engine schema. They identify submethods that need source-level mapping; they do not assert that any particular citation validates all choices.

| Procedure | Setting | Available choices |
| --- | --- | --- |
| `means.crossover` | Test Statistic | T-test / Z-test (σ known) |
| `means.one.t` | Test Statistic | T-test (σ estimated from the data) / Z-test (σ known) |
| `means.one.ni` | Test Statistic | T-test (σ estimated) / Z-test (σ known) |
| `means.one.equiv` | Test Statistic | T-tests (σ estimated) / Z-tests (σ known) |
| `means.paired.t` | Test Statistic | T-test (σ estimated from the data) / Z-test (σ known) |
| `means.two.repeated` | Test Statistic | Z-test (large sample) / T-test (df = N − 2) |
| `means.two.slopes` | Test Statistic | T-test (variances estimated) / Z-test (variances known) |
| `means.two.t` | Test Statistic | T-test (σ estimated from the data) / Z-test (σ known) |
| `means.two.welch` | Test Statistic | T-test (σ estimated from the data) / Z-test (σ known) |
| `means.two.ni` | Test Statistic | T-test (σ estimated) / Z-test (σ known) |
| `means.two.equiv` | Test Statistic | T-tests (σ estimated) / Z-tests (σ known) |
| `means.coprimary` | Test | T-tests (variances estimated) / Z-tests (known variances) |
| `prop.cc` | Test Type | Chi-square test (pooled z-test) / Chi-square test with continuity correction (Yates) / Z-test, unpooled variance / Fisher's exact test |
| `prop.cc` | Power Calculation | Exact (binomial enumeration) / Normal approximation |
| `prop.one` | Test Type | Exact binomial test / Z-test using P0 (score) / Z-test using P̂ (Wald) / Z-test using P0 with continuity correction / Z-test using P̂ with continuity correction |
| `prop.one` | Power Calculation | Exact (binomial enumeration) / Normal approximation |
| `prop.paired.mcnemar` | Test / Power Calculation | Exact conditional test, exact power / Normal approximation (Connor) |
| `prop.two` | Test Type | Z-test, pooled variance (= chi-square) / Z-test, unpooled variance / Z-test with continuity correction (Yates) / Fisher's exact test |
| `prop.two` | Power Calculation | Exact (binomial enumeration) / Normal approximation |
| `prop.two.ni` | Power Calculation | Normal approximation / Exact (binomial enumeration) |
| `prop.two.equiv` | Power Calculation | Normal approximation / Exact (binomial enumeration) |
| `corr.one` | Calculation Method | Exact distribution of r / Fisher's z transformation |
| `surv.logrank` | Events Formula | Schoenfeld (1981) / Freedman (1982) |
| `rates.one.poisson` | Test Type | Exact Poisson test / Z-test using λ0 |
| `rates.one.poisson` | Power Calculation | Exact (Poisson tail probabilities) / Normal approximation |
| `rates.two.poisson` | Variance | Pooled rate under H0 (score-type) / Rates under H1 (Wald) |
| `rates.two.negbin` | Variance under H0 | Method 3: pooled (H0 MLE) rate / Method 2: true rates / Method 1: control rate |
| `nonpar.mannwhitney` | Method | Location shift with ARE adjustment / Noether: P(X1 > X2) |
| `gs.logrank` | Efficacy Boundary | O'Brien–Fleming type spending (Lan–DeMets) / Pocock type spending (Lan–DeMets) / Power family spending α·t^ρ (Kim–DeMets) / Hwang–Shih–DeCani spending (γ) / O'Brien–Fleming (classical) / Pocock (classical) / Wang–Tsiatis (Δ) |
| `gs.logrank` | Futility Boundary | None / β-spending (Lan–DeMets type) / Fixed z-values |
| `gs.means` | Efficacy Boundary | O'Brien–Fleming type spending (Lan–DeMets) / Pocock type spending (Lan–DeMets) / Power family spending α·t^ρ (Kim–DeMets) / Hwang–Shih–DeCani spending (γ) / O'Brien–Fleming (classical) / Pocock (classical) / Wang–Tsiatis (Δ) |
| `gs.means` | Futility Boundary | None / β-spending (Lan–DeMets type) / Fixed z-values |
| `gs.props` | Efficacy Boundary | O'Brien–Fleming type spending (Lan–DeMets) / Pocock type spending (Lan–DeMets) / Power family spending α·t^ρ (Kim–DeMets) / Hwang–Shih–DeCani spending (γ) / O'Brien–Fleming (classical) / Pocock (classical) / Wang–Tsiatis (Δ) |
| `gs.props` | Futility Boundary | None / β-spending (Lan–DeMets type) / Fixed z-values |
| `diag.paired` | Test / Power Calculation | Exact conditional test, exact power / Normal approximation (Connor) |
| `prop.one.ci` | Interval Method | Exact (Clopper–Pearson) / Wilson score / Wilson score with continuity correction / Agresti–Coull / Jeffreys / Simple asymptotic (Wald) |
| `prop.two.ci` | Interval Method | Newcombe hybrid score (Wilson) / Newcombe hybrid score with continuity correction / Agresti–Caffo / Simple asymptotic (Wald) / Wald with continuity correction |
| `means.two.cluster` | Test Statistic | T-test on cluster means (df = k1 + k2 − 2) / Z-test |
| `sim.glmm.cluster` | Test | Wald t-test (k1 + k2 − 2 df) / Wald z-test (normal approximation) |
| `sim.sw` | Analysis Model | Random cluster intercept (Hussey & Hughes) / Random cluster and cluster-period intercepts (Hooper et al.) |
| `sim.winratio` | Test | Finkelstein–Schoenfeld (net benefit, permutation variance) / Win ratio (log win ratio, U-statistic variance) |
| `sim.mmrm` | Analysis | MMRM: all visits, unstructured covariance (REML) / t-test on the patients who complete the study |
| `sim.means.robust` | Test | Welch's t-test (unequal variances) / Pooled t-test (equal variances) / Mann–Whitney (Wilcoxon rank-sum) test |
| `sim.ordinal` | Test | Proportional odds model, Wald test / Proportional odds model, likelihood-ratio test / Mann–Whitney (Wilcoxon rank-sum) test with ties |
| `sim.multiarm` | Multiplicity Adjustment | Dunnett's test (single-step) / Holm's step-down procedure / Hochberg's step-up procedure / Bonferroni's correction / None (each test at level α) |
| `sim.surv.nph` | Test | Logrank / Fleming–Harrington weighted logrank FH(ρ, γ) / MaxCombo: largest of FH(0, 0), FH(0, 1), FH(1, 0), FH(1, 1) / Difference in restricted mean survival time (RMST) |

## Sources located for DAG and Lean follow-up

The following primary-source pages were retrieved to establish candidate relevance. These are search findings, **not Google Scholar BibTeX imports**. Import the actual Scholar-exported record before assigning a production bibliography key. Matching a paper to a method is a review judgment and does not itself validate this implementation.

| Component | Retrieved source | Intended role |
| --- | --- | --- |
| AIPW | [Robins, Rotnitzky and Zhao, Estimation of Regression Coefficients When Some Regressors are not Always Observed](https://www.tandfonline.com/doi/full/10.1080/01621459.1994.10476818) | Foundational augmented estimating equations and semiparametric efficiency; already named in the DAG documentation. |
| AIPW | [Bang and Robins, Doubly Robust Estimation in Missing Data and Causal Inference Models](https://onlinelibrary.wiley.com/doi/10.1111/j.1541-0420.2005.00377.x) | Doubly robust causal estimation; already named in the DAG documentation. |
| Parametric g-formula | [Robins, A new approach to causal inference in mortality studies with a sustained exposure period](https://www.sciencedirect.com/science/article/pii/0270025586900886) | G-formula foundation; this is not a source validating every numerical approximation in this calculator. |
| TMLE | [van der Laan and Rubin, Targeted Maximum Likelihood Learning](https://biostats.bepress.com/ucbbiostat/paper213/) | Foundational targeted estimation and efficient influence equations. |
| TMLE implementation | [Gruber and van der Laan, tmle: An R Package for Targeted Maximum Likelihood Estimation](https://www.jstatsoft.org/article/view/v051i13) | Point-treatment effects and variance estimation. Does not imply the calculator implements every option of the R package. |
| Bounded-outcome TMLE | [Gruber and van der Laan, A Targeted Maximum Likelihood Estimator of a Causal Effect on a Bounded Continuous Outcome](https://pmc.ncbi.nlm.nih.gov/articles/PMC3126669/) | Bounded logistic fluctuation and influence-function information. Keep distinct from the unbounded identity-link fluctuation. |
| Super Learner | [van der Laan, Polley and Hubbard, Super learner](https://pubmed.ncbi.nlm.nih.gov/17910531/) | Cross-validated ensemble nuisance estimation. |
| Cross-validated TMLE | [Zheng and van der Laan, Asymptotic Theory for Cross-validated Targeted Maximum Likelihood Estimation](https://biostats.bepress.com/ucbbiostat/paper273/) | Sample splitting and asymptotic theory. Choose the working-paper or published-chapter record deliberately and check that the implementation matches its estimator. |

Additional mapping is needed for DAG Kaplan–Meier/Greenwood inference, RMST, GLMs and joint Wald tests, Cox/Breslow inference, structured repeated/clustered analyses, and Lean influence-function planning. Existing catalog sources may cover some of these after checking the actual formula and assumptions. No comprehensive source verification is claimed here.

## Existing reference text

These entries are **unverified legacy citation strings copied from the code**, not newly composed references or retrieved BibTeX records. The labels identify the procedure table above.

**R001.** Chow, S.-C., Shao, J., Wang, H. & Lokhnygina, Y. (2018). Sample Size Calculations in Clinical Research, 3rd ed. Chapman & Hall/CRC.

**R002.** Julious, S. A. (2010). Sample Sizes for Clinical Trials. Chapman & Hall/CRC.

**R003.** Senn, S. (2002). Cross-over Trials in Clinical Research, 2nd ed. Wiley.

**R004.** Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and the power approach for assessing the equivalence of average bioavailability. Journal of Pharmacokinetics and Biopharmaceutics 15, 657–680.

**R005.** Owen, D. B. (1965). A special case of a bivariate non-central t-distribution. Biometrika 52, 437–446.

**R006.** Diletti, E., Hauschke, D. & Steinijans, V. W. (1991). Sample size determination for bioequivalence assessment by means of confidence intervals. Int. J. Clin. Pharmacol. Ther. Toxicol. 29, 1–8.

**R007.** Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences, 2nd ed. Lawrence Erlbaum.

**R008.** Desu, M. M. & Raghavarao, D. (1990). Sample Size Methodology. Academic Press.

**R009.** Fleiss, J. L. (1986). The Design and Analysis of Clinical Experiments. Wiley.

**R010.** Benton, D. & Krishnamoorthy, K. (2003). Computing discrete mixtures of continuous distributions: noncentral chisquare, noncentral t and the distribution of the square of the sample multiple correlation coefficient. Computational Statistics & Data Analysis 43, 249–267.

**R011.** Montgomery, D. C. (2017). Design and Analysis of Experiments, 9th ed. Wiley.

**R012.** Machin, D., Campbell, M. J., Tan, S. B. & Tan, S. H. (2018). Sample Sizes for Clinical, Laboratory and Epidemiology Studies, 4th ed. Wiley-Blackwell.

**R013.** Lenth, R. V. (1989). Algorithm AS 243: cumulative distribution function of the non-central t distribution. Applied Statistics 38, 185–189.

**R014.** Blackwelder, W. C. (1982). Proving the null hypothesis in clinical trials. Controlled Clinical Trials 3, 345–353.

**R015.** Diggle, P. J., Heagerty, P., Liang, K.-Y. & Zeger, S. L. (2002). Analysis of Longitudinal Data, 2nd ed. Oxford University Press.

**R016.** Liu, G. & Liang, K.-Y. (1997). Sample size calculations for studies with correlated observations. Biometrics 53, 937–947.

**R017.** Muller, K. E. & Barton, C. N. (1989). Approximate power for repeated-measures ANOVA lacking sphericity. Journal of the American Statistical Association 84, 549–555.

**R018.** Fitzmaurice, G. M., Laird, N. M. & Ware, J. H. (2011). Applied Longitudinal Analysis, 2nd ed. Wiley.

**R019.** Ahn, C., Heo, M. & Zhang, S. (2015). Sample Size Calculations for Clustered and Longitudinal Outcomes in Clinical Research. Chapman & Hall/CRC.

**R020.** Welch, B. L. (1947). The generalization of 'Student's' problem when several different population variances are involved. Biometrika 34, 28–35.

**R021.** Borm, G. F., Fransen, J. & Lemmens, W. A. J. G. (2007). A simple sample size formula for analysis of covariance in randomized clinical trials. Journal of Clinical Epidemiology 60, 1234–1238.

**R022.** Sozu, T., Sugimoto, T. & Hamasaki, T. (2011). Sample size determination in superiority clinical trials with multiple co-primary correlated endpoints. Journal of Biopharmaceutical Statistics 21, 650–668.

**R023.** Sugimoto, T., Sozu, T. & Hamasaki, T. (2012). A convenient formula for sample size calculations in clinical trials with multiple co-primary continuous endpoints. Pharmaceutical Statistics 11, 118–128.

**R024.** Offen, W. et al. (2007). Multiple co-primary endpoints: medical and statistical solutions. Drug Information Journal 41, 31–46.

**R025.** Schlesselman, J. J. (1982). Case-Control Studies: Design, Conduct, Analysis. Oxford University Press.

**R026.** Fleiss, J. L., Levin, B. & Paik, M. C. (2003). Statistical Methods for Rates and Proportions, 3rd ed. Wiley.

**R027.** Kelsey, J. L., Whittemore, A. S., Evans, A. S. & Thompson, W. D. (1996). Methods in Observational Epidemiology, 2nd ed. Oxford University Press.

**R028.** Dupont, W. D. (1988). Power calculations for matched case-control studies. Biometrics 44, 1157–1168.

**R029.** Breslow, N. E. & Day, N. E. (1980). Statistical Methods in Cancer Research, Vol. I: The Analysis of Case-Control Studies. IARC.

**R030.** Guenther, W. C. (1977). Power and sample size for approximate chi-square tests. The American Statistician 31, 83–85.

**R031.** Chernick, M. R. & Liu, C. Y. (2002). The saw-toothed behavior of power versus sample size and software solutions. The American Statistician 56, 149–155.

**R032.** Connor, R. J. (1987). Sample size for testing differences in proportions for the paired-sample design. Biometrics 43, 207–211.

**R033.** Lachin, J. M. (1992). Power and sample size evaluation for the McNemar test with application to matched case-control studies. Statistics in Medicine 11, 1239–1251.

**R034.** Nam, J. (1987). A simple approximation for calculating sample sizes for detecting linear trend in proportions. Biometrics 43, 701–705.

**R035.** Armitage, P. (1955). Tests for linear trends in proportions and frequencies. Biometrics 11, 375–386.

**R036.** Simon, R. (1989). Optimal two-stage designs for phase II clinical trials. Controlled Clinical Trials 10, 1–10.

**R037.** Jung, S.-H., Lee, T., Kim, K. & George, S. L. (2004). Admissible two-stage designs for phase II cancer clinical trials. Statistics in Medicine 23, 561–569.

**R038.** Casagrande, J. T., Pike, M. C. & Smith, P. G. (1978). An improved approximate formula for calculating sample sizes for comparing two binomial distributions. Biometrics 34, 483–486.

**R039.** Fisher, R. A. (1935). The logic of inductive inference. Journal of the Royal Statistical Society 98, 39–54.

**R040.** Farrington, C. P. & Manning, G. (1990). Test statistics and sample size formulae for comparative binomial trials with null hypothesis of non-zero risk difference or non-unity relative risk. Statistics in Medicine 9, 1447–1454.

**R041.** Miettinen, O. & Nurminen, M. (1985). Comparative analysis of two rates. Statistics in Medicine 4, 213–226.

**R042.** Guenther, W. C. (1977). Desk calculation of probabilities for the distribution of the sample correlation coefficient. The American Statistician 31, 45–48.

**R043.** Fisher, R. A. (1921). On the 'probable error' of a coefficient of correlation deduced from a small sample. Metron 1, 3–32.

**R044.** Gatsonis, C. & Sampson, A. R. (1989). Multiple correlation: exact power and sample size calculations. Psychological Bulletin 106, 516–524.

**R045.** Dupont, W. D. & Plummer, W. D. (1998). Power and sample size calculations for studies involving linear regression. Controlled Clinical Trials 19, 589–601.

**R046.** Hsieh, F. Y., Bloch, D. A. & Larsen, M. D. (1998). A simple method of sample size calculation for linear and logistic regression. Statistics in Medicine 17, 1623–1634.

**R047.** Hsieh, F. Y. (1989). Sample size tables for logistic regression. Statistics in Medicine 8, 795–802.

**R048.** Riley, R. D., Ensor, J., Snell, K. I. E., Harrell, F. E., Martin, G. P., Reitsma, J. B., Moons, K. G. M., Collins, G. & van Smeden, M. (2020). Calculating the sample size required for developing a clinical prediction model. BMJ 368, m441.

**R049.** Riley, R. D., Snell, K. I. E., Ensor, J., Burke, D. L., Harrell, F. E., Moons, K. G. M. & Collins, G. S. (2019). Minimum sample size for developing a multivariable prediction model: Part I – Continuous outcomes. Statistics in Medicine 38, 1262–1275.

**R050.** Riley, R. D., Snell, K. I. E., Ensor, J., Burke, D. L., Harrell, F. E., Moons, K. G. M. & Collins, G. S. (2019). Minimum sample size for developing a multivariable prediction model: Part II – Binary and time-to-event outcomes. Statistics in Medicine 38, 1276–1296.

**R051.** Riley, R. D., Van Calster, B. & Collins, G. S. (2021). A note on estimating the Cox-Snell R² from a reported C statistic (AUROC) to inform sample size calculations for developing a prediction model with a binary outcome. Statistics in Medicine 40, 859–864.

**R052.** Ensor, J., Martin, E. C. & Riley, R. D. pmsampsize: Calculates the Minimum Sample Size Required for Developing a Multivariable Prediction Model. R package and Stata module.

**R053.** van Houwelingen, J. C. & Le Cessie, S. (1990). Predictive value of statistical models. Statistics in Medicine 9, 1303–1325.

**R054.** Hsieh, F. Y. & Lavori, P. W. (2000). Sample-size calculations for the Cox proportional hazards regression model with nonbinary covariates. Controlled Clinical Trials 21, 552–560.

**R055.** Schoenfeld, D. A. (1983). Sample-size formula for the proportional-hazards regression model. Biometrics 39, 499–503.

**R056.** Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival distributions. Biometrika 68, 316–319.

**R057.** Freedman, L. S. (1982). Tables of the number of patients required in clinical trials using the logrank test. Statistics in Medicine 1, 121–129.

**R058.** Lachin, J. M. & Foulkes, M. A. (1986). Evaluation of sample size and power for analyses of survival with allowance for nonuniform patient entry, losses to follow-up, noncompliance, and stratification. Biometrics 42, 507–519.

**R059.** Jung, S.-H., Kang, S. J., McCall, L. M. & Blumenstein, B. (2005). Sample size computation for two-sample noninferiority log-rank test. Journal of Biopharmaceutical Statistics 15, 969–979.

**R060.** Garwood, F. (1936). Fiducial limits for the Poisson distribution. Biometrika 28, 437–442.

**R061.** Gu, K., Ng, H. K. T., Tang, M. L. & Schucany, W. R. (2008). Testing the ratio of two Poisson rates. Biometrical Journal 50, 283–298.

**R062.** Zhu, H. & Lakkis, H. (2014). Sample size calculation for comparing two negative binomial rates. Statistics in Medicine 33, 376–387.

**R063.** Keene, O. N., Jones, M. R. K., Lane, P. W. & Anderson, J. (2007). Analysis of exacerbation rates in asthma and chronic obstructive pulmonary disease. Pharmaceutical Statistics 6, 89–97.

**R064.** Ostle, B. & Malone, L. C. (1988). Statistics in Research, 4th ed. Iowa State University Press.

**R065.** Lehmann, E. L. (1975). Nonparametrics: Statistical Methods Based on Ranks. Holden-Day.

**R066.** Al-Sunduqchi, M. S. (1990). Determining the Appropriate Sample Size for Inferences Based on the Wilcoxon Statistics. PhD thesis, University of Wyoming.

**R067.** Noether, G. E. (1987). Sample size determination for some common nonparametric tests. Journal of the American Statistical Association 82, 645–647.

**R068.** Mahoney, M. & Magel, R. (1996). Estimation of the power of the Kruskal-Wallis test. Biometrical Journal 38, 613–630.

**R069.** Jennison, C. & Turnbull, B. W. (2000). Group Sequential Methods with Applications to Clinical Trials. Chapman & Hall/CRC.

**R070.** Lan, K. K. G. & DeMets, D. L. (1983). Discrete sequential boundaries for clinical trials. Biometrika 70, 659–663.

**R071.** O'Brien, P. C. & Fleming, T. R. (1979). A multiple testing procedure for clinical trials. Biometrics 35, 549–556.

**R072.** Pocock, S. J. (1977). Group sequential methods in the design and analysis of clinical trials. Biometrika 64, 191–199.

**R073.** Kim, K. & DeMets, D. L. (1987). Design and analysis of group sequential tests based on the type I error spending rate function. Biometrika 74, 149–154.

**R074.** Hwang, I. K., Shih, W. J. & DeCani, J. S. (1990). Group sequential designs using a family of type I error probability spending functions. Statistics in Medicine 9, 1439–1445.

**R075.** Wang, S. K. & Tsiatis, A. A. (1987). Approximately optimal one-parameter boundaries for group sequential trials. Biometrics 43, 193–199.

**R076.** Pampallona, S., Tsiatis, A. A. & Kim, K. (2001). Interim monitoring of group sequential trials using spending functions for the type I and type II error probabilities. Drug Information Journal 35, 1113–1121.

**R077.** Walter, S. D., Eliasziw, M. & Donner, A. (1998). Sample size and optimal designs for reliability studies. Statistics in Medicine 17, 101–110.

**R078.** Shrout, P. E. & Fleiss, J. L. (1979). Intraclass correlations: uses in assessing rater reliability. Psychological Bulletin 86, 420–428.

**R079.** Fleiss, J. L., Cohen, J. & Everitt, B. S. (1969). Large sample standard errors of kappa and weighted kappa. Psychological Bulletin 72, 323–327.

**R080.** Cantor, A. B. (1996). Sample-size calculations for Cohen's kappa. Psychological Methods 1, 150–153.

**R081.** Hanley, J. A. & McNeil, B. J. (1982). The meaning and use of the area under a receiver operating characteristic (ROC) curve. Radiology 143, 29–36.

**R082.** Obuchowski, N. A. (1998). Sample size calculations in studies of test accuracy. Statistical Methods in Medical Research 7, 371–392.

**R083.** Flahault, A., Cadilhac, M. & Thomas, G. (2005). Sample size calculation should be performed for design accuracy in diagnostic test studies. Journal of Clinical Epidemiology 58, 859–862.

**R084.** Pepe, M. S. (2003). The Statistical Evaluation of Medical Tests for Classification and Prediction. Oxford University Press.

**R085.** Alonzo, T. A., Pepe, M. S. & Moskowitz, C. S. (2002). Sample size calculations for comparative studies of medical tests for detecting presence of disease. Statistics in Medicine 21, 835–852.

**R086.** Kupper, L. L. & Hafner, K. B. (1989). How appropriate are popular sample size formulas? The American Statistician 43, 101–105.

**R087.** Newcombe, R. G. (1998). Two-sided confidence intervals for the single proportion: comparison of seven methods. Statistics in Medicine 17, 857–872.

**R088.** Clopper, C. J. & Pearson, E. S. (1934). The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika 26, 404–413.

**R089.** Agresti, A. & Coull, B. A. (1998). Approximate is better than 'exact' for interval estimation of binomial proportions. The American Statistician 52, 119–126.

**R090.** Newcombe, R. G. (1998). Interval estimation for the difference between independent proportions: comparison of eleven methods. Statistics in Medicine 17, 873–890.

**R091.** Agresti, A. & Caffo, B. (2000). Simple and effective confidence intervals for proportions and differences of proportions result from adding two successes and two failures. The American Statistician 54, 280–288.

**R092.** Donner, A. & Klar, N. (2000). Design and Analysis of Cluster Randomization Trials in Health Research. Arnold.

**R093.** Campbell, M. J. & Walters, S. J. (2014). How to Design, Analyse and Report Cluster Randomised Trials in Medicine and Health Related Research. Wiley.

**R094.** Hayes, R. J. & Moulton, L. H. (2017). Cluster Randomised Trials, 2nd ed. CRC Press.

**R095.** Hussey, M. A. & Hughes, J. P. (2007). Design and analysis of stepped wedge cluster randomized trials. Contemporary Clinical Trials 28, 182–191.

**R096.** Hemming, K., Haines, T. P., Chilton, P. J., Girling, A. J. & Lilford, R. J. (2015). The stepped wedge cluster randomised trial: rationale, design, analysis, and reporting. BMJ 350, h391.

**R097.** Eldridge, S. M., Ashby, D. & Kerry, S. (2006). Sample size for cluster randomized trials: effect of coefficient of variation of cluster size and analysis method. International Journal of Epidemiology 35, 1292–1300.

**R098.** Kuznetsova, A., Brockhoff, P. B. & Christensen, R. H. B. (2017). lmerTest package: tests in linear mixed effects models. Journal of Statistical Software 82(13).

**R099.** Landau, S. & Stahl, D. (2013). Sample size and power calculations for medical studies by simulation when closed form expressions are not available. Statistical Methods in Medical Research 22, 324–345.

**R100.** Li, P. & Redden, D. T. (2015). Comparing denominator degrees of freedom approximations for the generalized linear mixed model in analyzing binary outcome in small sample cluster-randomized trials. BMC Medical Research Methodology 15, 38.

**R101.** Pinheiro, J. C. & Bates, D. M. (1995). Approximations to the log-likelihood function in the nonlinear mixed-effects model. Journal of Computational and Graphical Statistics 4, 12–35.

**R102.** Hooper, R., Teerenstra, S., de Hoop, E. & Eldridge, S. (2016). Sample size calculation for stepped wedge and other longitudinal cluster randomised trials. Statistics in Medicine 35, 4718–4728.

**R103.** Baio, G. et al. (2015). Sample size calculation for a stepped wedge trial. Trials 16, 354.

**R104.** Pocock, S. J., Ariti, C. A., Collier, T. J. & Wang, D. (2012). The win ratio: a new approach to the analysis of composite endpoints in clinical trials based on clinical priorities. European Heart Journal 33, 176–182.

**R105.** Finkelstein, D. M. & Schoenfeld, D. A. (1999). Combining mortality and longitudinal measures in clinical trials. Statistics in Medicine 18, 1341–1354.

**R106.** Bebu, I. & Lachin, J. M. (2016). Large sample inference for a win ratio analysis of a composite outcome based on prioritized components. Biostatistics 17, 178–187.

**R107.** Yu, R. X. & Ganju, J. (2022). Sample size formula for a win ratio endpoint. Statistics in Medicine 41, 950–963.

**R108.** Laird, N. M. & Ware, J. H. (1982). Random-effects models for longitudinal data. Biometrics 38, 963–974.

**R109.** Satterthwaite, F. E. (1946). An approximate distribution of estimates of variance components. Biometrics Bulletin 2, 110–114.

**R110.** Arnold, B. F., Hogan, D. R., Colford, J. M. & Hubbard, A. E. (2011). Simulation methods to estimate design power: an overview for applied research. BMC Medical Research Methodology 11, 94.

**R111.** Mallinckrodt, C. H., Clark, W. S. & David, S. R. (2001). Accounting for dropout bias using mixed-effects models. Journal of Biopharmaceutical Statistics 11, 9–21.

**R112.** Lu, K., Luo, X. & Chen, P.-Y. (2008). Sample size estimation for repeated measures analysis in randomized clinical trials with missing data. International Journal of Biostatistics 4(1), 9.

**R113.** Anderson, T. W. (1957). Maximum likelihood estimates for a multivariate normal distribution when some observations are missing. Journal of the American Statistical Association 52, 200–203.

**R114.** Fagerland, M. W. & Sandvik, L. (2009). The Wilcoxon–Mann–Whitney test under scrutiny. Statistics in Medicine 28, 1487–1497.

**R115.** Whitehead, J. (1993). Sample size calculations for ordered categorical data. Statistics in Medicine 12, 2257–2271.

**R116.** McCullagh, P. (1980). Regression models for ordinal data. Journal of the Royal Statistical Society B 42, 109–142.

**R117.** Roozenbeek, B. et al. (2011). The added value of ordinal analysis in clinical trials: an example in traumatic brain injury. Critical Care 15, R127.

**R118.** Dunnett, C. W. (1955). A multiple comparison procedure for comparing several treatments with a control. Journal of the American Statistical Association 50, 1096–1121.

**R119.** Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics 6, 65–70.

**R120.** Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika 75, 800–802.

**R121.** Horn, M. & Vollandt, R. (1995). Multiple Tests und Auswahlverfahren. Gustav Fischer (sample sizes for many-to-one comparisons).

**R122.** Fleming, T. R. & Harrington, D. P. (1991). Counting Processes and Survival Analysis. Wiley.

**R123.** Lin, R. S., Lin, J., Roychoudhury, S. et al. (2020). Alternative analysis methods for time to event endpoints under nonproportional hazards: a comparative analysis. Statistics in Biopharmaceutical Research 12, 187–198.

**R124.** Royston, P. & Parmar, M. K. B. (2013). Restricted mean survival time: an alternative to the hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. BMC Medical Research Methodology 13, 152.

**R125.** Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of Computational and Graphical Statistics 1, 141–149.

## Completion requirements

1. Locate the user’s existing `references.bib` by repository, branch and path, or receive the file. Preserve existing citation keys and content.
2. Retrieve each required Google Scholar Cite → BibTeX export. Retain the exported text, source link and retrieval date; do not reconstruct records from remembered metadata or search snippets.
3. Resolve duplicates and record method-to-key mappings, distinguishing the test, power formula, numerical algorithm and general background sources.
4. Extend references to DAG and Lean results, and validate that every referenced key exists, that procedure/variant coverage is complete, and that web and document exports agree.
5. Run the existing checks appropriate to the actual code changes. This inventory changes no numerical behavior and is not evidence of statistical validation.
