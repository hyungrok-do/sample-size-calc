// Curated study-design routes; calculation fields and allowed options come from the engine's schemas.
// Keep every procedure reachable. web/tests/guide.test.mjs checks coverage and option compatibility.
const question = (node_id, prompt, choices) => ({
  type: 'question', node_id, prompt, choices: choices.map(([label, next]) => ({ label, next })),
});
const candidate = (node_id, procedure_id, proposed_option_values) => ({
  type: 'candidate', node_id, procedure_id, proposed_option_values,
});
const unavailable = (node_id, reason) => ({ type: 'no_matching_calculator', node_id, reason });

export const guideData = {
  reviewed_commit: "2db563072e23e671b01a453740b9666f33adfcb6",
  root: question("q1", "What is your primary planning objective?", [
    ["Compare outcomes", question("q2", "What kind of comparison?", [
      ["Continuous measurements", question("q3", "How are the continuous measurements organized?", [
        ["One sample versus a reference", question("q4", "What quantity is primary?", [
          ["Mean", question("q5", "How will the population variability be treated?", [
            ["Estimated from the study: t test", candidate("q6", "means.one.t", {"test": "t"})],
            ["Known population variance: z test", candidate("q7", "means.one.t", {"test": "z"})]
          ])],
          ["Median / probability above the reference", candidate("q8", "nonpar.sign", {})],
          ["Symmetric location shift", candidate("q9", "nonpar.signedrank", {})]
        ])],
        ["Two paired measurements", question("q10", "What feature of the paired differences is primary?", [
          ["Mean difference", candidate("q11", "means.paired.t", {})],
          ["Sign / median", candidate("q12", "nonpar.sign", {})],
          ["Symmetric location shift", candidate("q13", "nonpar.signedrank", {})]
        ])],
        ["Two independent groups, one measurement", question("q14", "What comparison is primary?", [
          ["Mean difference adjusted for baseline covariates", candidate("q15", "means.two.ancova", {})],
          ["Unadjusted mean difference", question("q16", "Which variance model will the analysis use?", [
            ["Common variance", question("q17", "How will the population variability be treated?", [
              ["Estimated from the study: t test", candidate("q18", "means.two.t", {"test": "t"})],
              ["Known population variance: z test", candidate("q19", "means.two.t", {"test": "z"})]
            ])],
            ["Unequal variances", candidate("q20", "means.two.welch", {})],
            ["Evaluate the chosen test under a supported non-normal distribution", question("q21", "Which mean-comparison analysis is planned?", [
              ["Welch t test", question("q22", "Which outcome distribution will be simulated?", [
                ["Lognormal", candidate("q23", "sim.means.robust", {"test": "welch", "dist": "lognormal"})],
                ["Gamma", candidate("q24", "sim.means.robust", {"test": "welch", "dist": "gamma"})],
                ["Student t", candidate("q25", "sim.means.robust", {"test": "welch", "dist": "t"})],
                ["Contaminated normal", candidate("q26", "sim.means.robust", {"test": "welch", "dist": "contaminated"})],
                ["Uniform", candidate("q27", "sim.means.robust", {"test": "welch", "dist": "uniform"})],
                ["Normal", candidate("q28", "sim.means.robust", {"test": "welch", "dist": "normal"})]
              ])],
              ["Pooled t test", question("q29", "Which outcome distribution will be simulated?", [
                ["Lognormal", candidate("q30", "sim.means.robust", {"test": "t", "dist": "lognormal"})],
                ["Gamma", candidate("q31", "sim.means.robust", {"test": "t", "dist": "gamma"})],
                ["Student t", candidate("q32", "sim.means.robust", {"test": "t", "dist": "t"})],
                ["Contaminated normal", candidate("q33", "sim.means.robust", {"test": "t", "dist": "contaminated"})],
                ["Uniform", candidate("q34", "sim.means.robust", {"test": "t", "dist": "uniform"})],
                ["Normal", candidate("q35", "sim.means.robust", {"test": "t", "dist": "normal"})]
              ])]
            ])]
          ])],
          ["Several correlated co-primary endpoints; all must succeed", question("q196", "How will each endpoint's variability be treated?", [
            ["Estimated from the study: t tests", candidate("q197", "means.coprimary", {"test": "t"})],
            ["Known population variances: z tests", candidate("q198", "means.coprimary", {"test": "z"})]
          ])],
          ["Rank ordering / location shift", candidate("q36", "nonpar.mannwhitney", {})]
        ])],
        ["Balanced AB/BA crossover, mean difference", question("q37", "How will the population variability be treated?", [
          ["Estimated from the study: t test", candidate("q38", "means.crossover", {"test": "t"})],
          ["Known population variance: z test", candidate("q39", "means.crossover", {"test": "z"})]
        ])],
        ["Three or more independent groups", question("q40", "Which claim must succeed?", [
          ["An omnibus group difference", candidate("q41", "means.anova", {})],
          ["Factorial main effect or interaction", question("q42", "Which factorial effect determines the required power?", [
            ["Factor A", candidate("q43", "means.anova2", {"goal": "A"})],
            ["Factor B", candidate("q44", "means.anova2", {"goal": "B"})],
            ["A × B interaction", candidate("q45", "means.anova2", {"goal": "AB"})],
            ["Minimum marginal power across effects", candidate("q46", "means.anova2", {"goal": "all"})]
          ])],
          ["One or more treatment-control comparisons", question("q47", "Which comparisons must succeed?", [
            ["At least one effective treatment", question("q48", "How will multiplicity be handled?", [
              ["Dunnett", candidate("q49", "sim.multiarm", {"goal": "any", "adjust": "dunnett"})],
              ["Holm", candidate("q50", "sim.multiarm", {"goal": "any", "adjust": "holm"})],
              ["Hochberg", candidate("q51", "sim.multiarm", {"goal": "any", "adjust": "hochberg"})],
              ["Bonferroni", candidate("q52", "sim.multiarm", {"goal": "any", "adjust": "bonferroni"})],
              ["No adjustment", candidate("q53", "sim.multiarm", {"goal": "any", "adjust": "none"})]
            ])],
            ["All effective treatments", question("q54", "How will multiplicity be handled?", [
              ["Dunnett", candidate("q55", "sim.multiarm", {"goal": "all", "adjust": "dunnett"})],
              ["Holm", candidate("q56", "sim.multiarm", {"goal": "all", "adjust": "holm"})],
              ["Hochberg", candidate("q57", "sim.multiarm", {"goal": "all", "adjust": "hochberg"})],
              ["Bonferroni", candidate("q58", "sim.multiarm", {"goal": "all", "adjust": "bonferroni"})],
              ["No adjustment", candidate("q59", "sim.multiarm", {"goal": "all", "adjust": "none"})]
            ])],
            ["Power for each treatment", question("q60", "How will multiplicity be handled?", [
              ["Dunnett", candidate("q61", "sim.multiarm", {"goal": "each", "adjust": "dunnett"})],
              ["Holm", candidate("q62", "sim.multiarm", {"goal": "each", "adjust": "holm"})],
              ["Hochberg", candidate("q63", "sim.multiarm", {"goal": "each", "adjust": "hochberg"})],
              ["Bonferroni", candidate("q64", "sim.multiarm", {"goal": "each", "adjust": "bonferroni"})],
              ["No adjustment", candidate("q65", "sim.multiarm", {"goal": "each", "adjust": "none"})]
            ])]
          ])]
        ])],
        ["Repeated conditions in one group", candidate("q66", "means.rm.anova", {})],
        ["Two groups followed over time", question("q67", "Which time contrast is primary?", [
          ["Average across visits", candidate("q68", "means.two.repeated", {})],
          ["Linear slope difference", question("q69", "How are trajectories analyzed?", [
            ["Complete schedules; analytic participant-slope calculation", candidate("q70", "means.two.slopes", {})],
            ["Random-effect mixed model with supported dropout", candidate("q71", "sim.lmm.longitudinal", {})]
          ])],
          ["Difference at the final visit", candidate("q72", "sim.mmrm", {"analysis": "mmrm"})],
          ["Another treatment-by-time contrast", unavailable("q73", "No general longitudinal contrast engine.")]
        ])]
      ])],
      ["Binary or categorical outcomes", question("q74", "What type of outcome is recorded?", [
        ["Yes/no", question("q75", "How are observations compared?", [
          ["One sample versus a fixed benchmark", question("q76", "Which prespecified one-proportion test?", [
            ["Exact binomial", candidate("q77", "prop.one", {"test": "exact"})],
            ["Score test", candidate("q78", "prop.one", {"test": "zp0"})],
            ["Wald test", candidate("q79", "prop.one", {"test": "zphat"})],
            ["Score with continuity correction", candidate("q80", "prop.one", {"test": "zp0cc"})],
            ["Wald with continuity correction", candidate("q81", "prop.one", {"test": "zphatcc"})]
          ])],
          ["Two independent groups", question("q82", "Which binary analysis will be used?", [
            ["Pooled score / chi-square", candidate("q83", "prop.two", {"test": "zpool"})],
            ["Unpooled Wald", candidate("q84", "prop.two", {"test": "zunpool"})],
            ["Continuity-corrected chi-square", candidate("q85", "prop.two", {"test": "zpoolcc"})],
            ["Fisher exact", candidate("q86", "prop.two", {"test": "fisher", "method": "exact"})]
          ])],
          ["Two paired responses", candidate("q87", "prop.paired.mcnemar", {})],
          ["Cases versus controls (case-control study of an exposure)", question("q199", "How are the controls chosen?", [
            ["Independently of the cases (unmatched)", candidate("q200", "prop.cc", {})],
            ["Individually matched to each case (1 case : M controls)", candidate("q201", "prop.cc.matched", {})]
          ])],
          ["Ordered independent dose/exposure groups; scored trend", candidate("q88", "prop.trend", {})],
          ["Several independent groups; omnibus association", candidate("q89", "prop.chisq", {"es": "w"})],
          ["Three or more repeated responses", unavailable("q90", "No Cochran-Q or general repeated-binary route.")]
        ])],
        ["Unordered categories", question("q91", "Which null hypothesis?", [
          ["Match specified category probabilities", candidate("q92", "prop.chisq", {"es": "gof"})],
          ["Independent contingency-table association", candidate("q93", "prop.chisq", {"es": "w"})]
        ])],
        ["Ordered outcome categories in two independent groups", question("q94", "Which analysis target?", [
          ["Common cumulative odds, Wald test", candidate("q95", "sim.ordinal", {"test": "wald"})],
          ["Common cumulative odds, likelihood-ratio test", candidate("q96", "sim.ordinal", {"test": "lr"})],
          ["Rank ordering with ties", candidate("q97", "sim.ordinal", {"test": "mw"})]
        ])]
      ])],
      ["Survival, rates or event hierarchy", question("q98", "What is the endpoint?", [
        ["Time to a first event", question("q99", "Which effect and analysis are planned?", [
          ["PH logrank with supported accrual and baseline model", candidate("q100", "surv.logrank", {})],
          ["Supported non-PH model", question("q101", "Which prespecified analysis?", [
            ["Logrank", candidate("q102", "sim.surv.nph", {"test": "logrank"})],
            ["Fleming-Harrington weighted logrank", candidate("q103", "sim.surv.nph", {"test": "fh"})],
            ["MaxCombo", candidate("q104", "sim.surv.nph", {"test": "maxcombo"})],
            ["RMST at a fixed horizon", candidate("q105", "sim.surv.nph", {"test": "rmst"})]
          ])],
          ["Competing risks, recurrent survival events, or another model", unavailable("q106", "No dedicated matching survival route.")]
        ])],
        ["Events per exposure time", question("q107", "Which independent count model?", [
          ["One Poisson rate versus reference", candidate("q108", "rates.one.poisson", {})],
          ["Two Poisson rates", candidate("q109", "rates.two.poisson", {})],
          ["Two overdispersed negative-binomial rates", candidate("q110", "rates.two.negbin", {})]
        ])],
        ["Prioritized death / event / score hierarchy", question("q111", "Which analysis?", [
          ["Finkelstein-Schoenfeld", candidate("q112", "sim.winratio", {"test": "fs"})],
          ["Win ratio", candidate("q113", "sim.winratio", {"test": "wr"})]
        ])]
      ])],
      ["Rank-based comparisons", question("q114", "How are the samples related?", [
        ["One sample or paired differences", question("q115", "Which target fits?", [
          ["Sign probability / median", candidate("q116", "nonpar.sign", {})],
          ["Symmetric location shift", candidate("q117", "nonpar.signedrank", {})]
        ])],
        ["Two independent samples", candidate("q118", "nonpar.mannwhitney", {})],
        ["Several independent samples", candidate("q119", "nonpar.kruskal", {})],
        ["Repeated conditions beyond a pair", unavailable("q120", "No Friedman power procedure.")]
      ])],
      ["Clustered or stepped-wedge study", question("q121", "Which randomized-cluster design?", [
        ["Parallel groups", question("q122", "Which endpoint?", [
          ["Continuous", question("q123", "Which analysis and size model?", [
            ["Equal-size cluster-mean analytic test", candidate("q124", "means.two.cluster", {})],
            ["Random-intercept LMM, possibly unequal sizes", candidate("q125", "sim.lmm.cluster", {})]
          ])],
          ["Binary", question("q126", "Which analysis?", [
            ["Equal-size design-effect z approximation", candidate("q127", "prop.two.cluster", {})],
            ["Random-intercept logistic mixed model", candidate("q128", "sim.glmm.cluster", {})]
          ])],
          ["Other endpoint or combined person-level longitudinal structure", unavailable("q129", "No matching general clustered-outcome template.")]
        ])],
        ["Cross-sectional stepped-wedge", question("q130", "Which endpoint and model?", [
          ["Continuous, analytic Hussey-Hughes", candidate("q131", "cluster.sw", {"outcome": "continuous"})],
          ["Binary, analytic approximation", candidate("q132", "cluster.sw", {"outcome": "binary"})],
          ["Continuous LMM simulation", question("q133", "Which stepped-wedge covariance model?", [
            ["Hussey–Hughes", candidate("q134", "sim.sw", {"analysis": "hh"})],
            ["Hooper cluster-period effects", candidate("q135", "sim.sw", {"analysis": "hooper"})]
          ])]
        ])],
        ["Closed-cohort stepped-wedge or partial nesting", unavailable("q136", "No dedicated matching procedure.")]
      ])],
      ["Variability", question("q137", "One or two normal populations?", [
        ["One variance versus reference", candidate("q138", "var.one", {})],
        ["Two independent variances", candidate("q139", "var.two", {})]
      ])]
    ])],
    ["Association or regression", question("q140", "What feature is tested?", [
      ["Pearson correlation", question("q141", "One sample or independent samples?", [
        ["One correlation versus a reference", candidate("q142", "corr.one", {})],
        ["Two independent correlations", candidate("q143", "corr.two", {})],
        ["Dependent / overlapping correlations", unavailable("q144", "No matching dependent-correlation method.")]
      ])],
      ["Predictor contribution", question("q145", "What is the outcome and target?", [
        ["Continuous, simple-regression slope", candidate("q146", "reg.slope", {})],
        ["Continuous, overall or incremental predictor block", candidate("q147", "reg.multiple", {})],
        ["Binary, one eligible logistic covariate", candidate("q148", "reg.logistic", {})],
        ["Survival, one PH covariate", candidate("q149", "reg.cox", {})]
      ])],
      ["Develop a clinical prediction model", question("q202", "What outcome will the model predict?", [
        ["Binary (logistic regression)", candidate("q203", "reg.prediction", {"outcome": "binary"})],
        ["Time to event", candidate("q204", "reg.prediction", {"outcome": "survival"})],
        ["Continuous (linear regression)", candidate("q205", "reg.prediction", {"outcome": "continuous"})]
      ])]
    ])],
    ["Noninferiority or equivalence", question("q150", "Which effect scale and design describe the margin claim?", [
      ["One continuous mean or paired mean difference", question("q151", "One boundary or two?", [
        ["Noninferiority / superiority margin", candidate("q152", "means.one.ni", {})],
        ["Equivalence limits", question("q153", "How will the population variability be treated?", [
          ["Estimated from the study: t test", candidate("q154", "means.one.equiv", {"test": "t"})],
          ["Known population variance: z test", candidate("q155", "means.one.equiv", {"test": "z"})]
        ])]
      ])],
      ["Two independent continuous means", question("q156", "One boundary or two?", [
        ["Noninferiority / superiority margin", candidate("q157", "means.two.ni", {})],
        ["Equivalence limits", question("q158", "How will the population variability be treated?", [
          ["Estimated from the study: t test", candidate("q159", "means.two.equiv", {"test": "t"})],
          ["Known population variance: z test", candidate("q160", "means.two.equiv", {"test": "z"})]
        ])]
      ])],
      ["Ratio of geometric means for bioequivalence", question("q161", "Which bioequivalence design?", [
        ["Balanced AB/BA crossover", candidate("q162", "means.be.ratio", {"design": "2x2"})],
        ["Parallel groups", candidate("q163", "means.be.ratio", {"design": "parallel"})]
      ])],
      ["Two independent response probabilities", question("q164", "Which claim and scale?", [
        ["Noninferiority on risk difference", candidate("q165", "prop.two.ni", {"scale": "diff"})],
        ["Noninferiority on risk ratio", candidate("q166", "prop.two.ni", {"scale": "ratio"})],
        ["Noninferiority on odds ratio", candidate("q206", "prop.two.ni", {"scale": "or"})],
        ["Equivalence on risk difference", candidate("q167", "prop.two.equiv", {"scale": "diff"})],
        ["Equivalence on risk ratio", candidate("q207", "prop.two.equiv", {"scale": "ratio"})],
        ["Equivalence on odds ratio", candidate("q208", "prop.two.equiv", {"scale": "or"})]
      ])],
      ["Proportional-hazards noninferiority", candidate("q169", "surv.logrank.ni", {})]
    ])],
    ["Confidence-interval precision", question("q170", "Which quantity needs a narrow confidence interval?", [
      ["One mean or paired mean difference", candidate("q171", "means.one.ci", {})],
      ["Two independent means, common SD", candidate("q172", "means.two.ci", {})],
      ["One independent proportion", candidate("q173", "prop.one.ci", {})],
      ["Two independent proportions, risk difference", candidate("q174", "prop.two.ci", {})],
      ["One Pearson correlation", candidate("q175", "corr.ci", {})],
      ["Another quantity or dependence structure", unavailable("q176", "No dedicated matching precision procedure.")]
    ])],
    ["Accuracy or agreement", question("q177", "Accuracy against disease status, or agreement between measurements?", [
      ["Diagnostic accuracy", question("q178", "Which claim?", [
        ["One AUC versus benchmark", candidate("q179", "diag.auc", {})],
        ["Sensitivity threshold", candidate("q180", "diag.sens", {"goal": "se"})],
        ["Specificity threshold", candidate("q181", "diag.sens", {"goal": "sp"})],
        ["Both thresholds must succeed", candidate("q182", "diag.sens", {"goal": "both"})],
        ["Compare two tests on the same patients", question("q209", "Which accuracy must differ between the tests?", [
          ["Sensitivity", candidate("q210", "diag.paired", {"goal": "se"})],
          ["Specificity", candidate("q211", "diag.paired", {"goal": "sp"})],
          ["Both must differ", candidate("q212", "diag.paired", {"goal": "both"})]
        ])]
      ])],
      ["Agreement / reliability", question("q183", "Which measurement model?", [
        ["Quantitative ratings, one-way random-effects ICC", candidate("q184", "agree.icc", {})],
        ["Two binary raters with equal marginal positivity", candidate("q185", "agree.kappa", {})],
        ["Another agreement model", unavailable("q186", "No matching kappa/ICC/limits-of-agreement method.")]
      ])]
    ])],
    ["Interim or two-stage design", question("q187", "What stopping design is planned?", [
      ["Single-arm, binary-response phase II, two stages", question("q188", "Which two-stage design criterion?", [
        ["Minimum expected size under the null", candidate("q189", "prop.one.simon", {"criterion": "optimal"})],
        ["Minimum maximum sample size", candidate("q190", "prop.one.simon", {"criterion": "minimax"})]
      ])],
      ["Two independent groups, planned interim analyses", question("q191", "Which endpoint and test?", [
        ["Continuous, known-SD z", candidate("q192", "gs.means", {})],
        ["Binary, difference z", candidate("q193", "gs.props", {})],
        ["PH logrank", candidate("q194", "gs.logrank", {})],
        ["Other design or combined features", unavailable("q195", "No matching sequential combination.")]
      ])]
    ])]
  ]),
  features: [
    {
      "id": "cluster",
      "label": "Clustering",
      "allow": [
        "means.two.cluster",
        "prop.two.cluster",
        "sim.lmm.cluster",
        "sim.glmm.cluster",
        "cluster.sw",
        "sim.sw"
      ]
    },
    {
      "id": "repeated",
      "label": "3+ visits per participant",
      "allow": [
        "means.rm.anova",
        "means.two.repeated",
        "means.two.slopes",
        "sim.lmm.longitudinal",
        "sim.mmrm"
      ]
    },
    {
      "id": "adjust",
      "label": "Covariate adjustment",
      "allow": [
        "means.two.ancova",
        "reg.multiple",
        "reg.logistic",
        "reg.cox"
      ]
    },
    {
      "id": "interim",
      "label": "Interim stopping",
      "allow": [
        "gs.means",
        "gs.props",
        "gs.logrank",
        "prop.one.simon"
      ]
    },
    {
      "id": "ni",
      "label": "Noninferiority margin",
      "allow": [
        "means.one.ni",
        "means.two.ni",
        "prop.two.ni",
        "surv.logrank.ni"
      ]
    },
    {
      "id": "equiv",
      "label": "Equivalence limits",
      "allow": [
        "means.one.equiv",
        "means.two.equiv",
        "means.be.ratio",
        "prop.two.equiv"
      ]
    }
  ],
  catalog: {
    "means.one.t": {
      "label": "One-sample mean test",
      "inputs_to_request": "Reference and expected mean; outcome SD; one- or two-sided alternative.",
      "scope_and_checks": "Student t uses estimated SD. Choose z only when population SD is treated as known."
    },
    "means.paired.t": {
      "label": "Paired mean test",
      "inputs_to_request": "Expected paired difference; SD of differences, or component SDs and correlation.",
      "scope_and_checks": "N counts pairs. Pairing does not turn a two-arm baseline-adjusted comparison into this test."
    },
    "means.two.t": {
      "label": "Pooled two-sample mean test",
      "inputs_to_request": "Two means; common SD; group allocation.",
      "scope_and_checks": "Use only when the planned analysis pools variances. Known-variance z is a separate test option."
    },
    "means.two.welch": {
      "label": "Welch two-sample mean test",
      "inputs_to_request": "Two means; SD in each group; allocation.",
      "scope_and_checks": "Satterthwaite approximation. Serious skew or outliers may warrant the matching simulation."
    },
    "means.two.ancova": {
      "label": "Baseline-adjusted mean comparison",
      "inputs_to_request": "Means; outcome SD; covariate R-squared; number of covariates; allocation.",
      "scope_and_checks": "Analytic ANCOVA approximation. Many covariates and small samples need particular care."
    },
    "means.crossover": {
      "label": "Two-treatment crossover",
      "inputs_to_request": "Treatment difference; within-person residual SD or period-difference SD.",
      "scope_and_checks": "Specific 2-by-2 design. Assume carryover is negligible; not arbitrary repeated observations."
    },
    "sim.means.robust": {
      "label": "Simulation for non-normal outcomes",
      "inputs_to_request": "Group means and SDs; distribution and shape parameters; allocation; simulation repetitions.",
      "scope_and_checks": "Select the test before sizing. Mann-Whitney targets ordering, not the mean difference."
    },
    "means.anova": {
      "label": "One-way ANOVA",
      "inputs_to_request": "Group means and common SD, or Cohen f; group allocation.",
      "scope_and_checks": "Omnibus power is not power for every pairwise comparison."
    },
    "means.anova2": {
      "label": "Two-way factorial ANOVA",
      "inputs_to_request": "Levels of each factor; cell means; common SD; per-cell sample size.",
      "scope_and_checks": "Minimum marginal power refers to the least-powered nonzero effect. Detecting every effect together is a different success criterion."
    },
    "means.rm.anova": {
      "label": "One-factor repeated-measures ANOVA",
      "inputs_to_request": "Condition means or f; SD; number of measurements; correlation; epsilon.",
      "scope_and_checks": "One within-subject factor only. Not a general treatment-by-time interaction calculator."
    },
    "means.two.repeated": {
      "label": "Between-group time-average comparison",
      "inputs_to_request": "Group average means; SD; visit count; CS or AR(1) correlation.",
      "scope_and_checks": "Averages measurements; does not target a final-visit effect or a slope."
    },
    "means.two.slopes": {
      "label": "Between-group slope comparison",
      "inputs_to_request": "Slope difference; random-slope and residual SDs; visit times.",
      "scope_and_checks": "Analytic test on participant slopes. Incomplete trajectories need a matching longitudinal simulation."
    },
    "sim.lmm.longitudinal": {
      "label": "Longitudinal mixed-model simulation",
      "inputs_to_request": "Slope effect; random-effect SDs/correlation; residual SD; schedule; dropout.",
      "scope_and_checks": "Specific linear random-intercept/slope model. Does not cover every treatment-by-time contrast."
    },
    "sim.mmrm": {
      "label": "Final-visit MMRM simulation",
      "inputs_to_request": "Final-visit difference; SD; CS/AR(1) generating correlation; visits; arm-specific dropout.",
      "scope_and_checks": "Fits unstructured covariance. Dropout generation is MCAR; df choices are last-visit count minus 2 or z, not Kenward-Roger."
    },
    "sim.multiarm": {
      "label": "Several treatments versus one control",
      "inputs_to_request": "Treatment means; control mean; common SD; control allocation; success criterion.",
      "scope_and_checks": "any, all, and each define different powers. Not all pairwise comparisons or arbitrary outcome families."
    },
    "means.one.ni": {
      "label": "One-mean or paired noninferiority",
      "inputs_to_request": "Assumed true mean/difference; margin; SD of observations or differences.",
      "scope_and_checks": "For paired data, use differences and count pairs. Higher-is-better and lower-is-better require opposite directions."
    },
    "means.two.ni": {
      "label": "Two-mean noninferiority",
      "inputs_to_request": "Assumed true difference; margin; common SD; allocation.",
      "scope_and_checks": "Common-variance t or known-variance z setup. No dedicated Welch noninferiority route."
    },
    "means.one.equiv": {
      "label": "One-mean or paired equivalence",
      "inputs_to_request": "Assumed mean/difference; lower and upper bounds; appropriate SD.",
      "scope_and_checks": "TOST is not a nonsignificant difference test. Both one-sided tests must succeed."
    },
    "means.two.equiv": {
      "label": "Two-mean equivalence",
      "inputs_to_request": "Assumed mean difference; bounds; common SD; allocation.",
      "scope_and_checks": "Pooled-SD TOST. Unequal-variance equivalence is not a separately supported procedure."
    },
    "means.be.ratio": {
      "label": "Average bioequivalence on log scale",
      "inputs_to_request": "Expected ratio; ratio limits; CV appropriate to design.",
      "scope_and_checks": "Crossover uses within-subject CV; parallel uses total CV. Not replicated or scaled bioequivalence."
    },
    "prop.two.ni": {
      "label": "Two-proportion noninferiority",
      "inputs_to_request": "Both probabilities; margin on the chosen scale (difference, ratio or odds ratio); direction; allocation.",
      "scope_and_checks": "Farrington-Manning score test (Miettinen-Nurminen for the odds ratio). Exact enumeration describes power of the score test, not a new exact test."
    },
    "prop.two.equiv": {
      "label": "Two-proportion equivalence",
      "inputs_to_request": "Both probabilities; lower and upper limits on the chosen scale (difference, ratio or odds ratio); allocation.",
      "scope_and_checks": "Two one-sided score tests. Exact enumeration describes power of those tests, not a new exact test."
    },
    "surv.logrank.ni": {
      "label": "Survival noninferiority",
      "inputs_to_request": "Expected HR; HR margin; baseline survival; accrual; follow-up; loss.",
      "scope_and_checks": "Specific PH/exponential planning setup. Higher hazard may represent harm or benefit."
    },
    "prop.one": {
      "label": "One-proportion test",
      "inputs_to_request": "Null and expected probability; alternative; analysis test.",
      "scope_and_checks": "Exact binomial option is equal-tailed for two-sided tests. A power-enumeration option does not change the analysis test."
    },
    "prop.two": {
      "label": "Two-proportion test, including Fisher",
      "inputs_to_request": "Control probability and group probability/difference/RR/OR; allocation; analysis test.",
      "scope_and_checks": "The analysis test and the method used to calculate its power are separate choices. Fisher uses exact enumeration here; enumeration can also calculate power of z tests."
    },
    "prop.paired.mcnemar": {
      "label": "McNemar paired-binary test",
      "inputs_to_request": "Both discordant-pair probabilities, or total discordance and discordance odds ratio.",
      "scope_and_checks": "Two marginal response rates alone do not identify power. N is the number of pairs."
    },
    "prop.chisq": {
      "label": "Categorical goodness-of-fit or association",
      "inputs_to_request": "Null/alternative category probabilities, or Cohen w and degrees of freedom.",
      "scope_and_checks": "Asymptotic chi-square. Not an exact general-table method; not for dependent repeated categories."
    },
    "prop.trend": {
      "label": "Cochran-Armitage trend",
      "inputs_to_request": "Group response probabilities; group weights; ordered scores.",
      "scope_and_checks": "Order is in the groups/exposure. This is not proportional-odds regression for an ordinal outcome."
    },
    "sim.ordinal": {
      "label": "Two-group ordinal simulation",
      "inputs_to_request": "Control category probabilities; common OR or treatment probabilities; allocation.",
      "scope_and_checks": "Choose proportional-odds Wald, likelihood-ratio, or rank test explicitly. Ties are represented."
    },
    "nonpar.sign": {
      "label": "Sign test",
      "inputs_to_request": "Probability of a positive observation/difference, or normal-shift effect input.",
      "scope_and_checks": "Exact binomial power. Plan how ties/zero differences affect the usable sample; do not assume all enrolled pairs contribute."
    },
    "nonpar.signedrank": {
      "label": "Wilcoxon signed-rank",
      "inputs_to_request": "Location shift; SD of observations/differences; assumed parent distribution.",
      "scope_and_checks": "Power uses an ARE approximation, not exact signed-rank enumeration. Small samples may need separate verification."
    },
    "nonpar.mannwhitney": {
      "label": "Mann-Whitney / rank-sum",
      "inputs_to_request": "Shift and SD/distribution, or P(group 1 > group 2); allocation.",
      "scope_and_checks": "ARE or Noether approximation. Not generally a test of a mean or median difference; many ties need matching planning."
    },
    "nonpar.kruskal": {
      "label": "Kruskal-Wallis",
      "inputs_to_request": "Expected group locations; common scale; parent distribution; group count.",
      "scope_and_checks": "ANOVA/ARE approximation. Not a repeated-measures Friedman test or a pairwise comparison calculator."
    },
    "corr.one": {
      "label": "One Pearson correlation",
      "inputs_to_request": "Null and expected correlations; sample size or power.",
      "scope_and_checks": "Exact distribution assumes bivariate normal data; Fisher-z approximation is also available."
    },
    "corr.two": {
      "label": "Two independent correlations",
      "inputs_to_request": "Expected correlation in each sample; allocation.",
      "scope_and_checks": "Fisher-z approximation. Not overlapping/dependent correlations from the same people."
    },
    "reg.slope": {
      "label": "Simple linear-regression slope",
      "inputs_to_request": "Null/expected slope; predictor SD; residual SD or total outcome SD.",
      "scope_and_checks": "Simple regression. For adjusted Gaussian predictor blocks, consider the R-squared-increase route."
    },
    "reg.multiple": {
      "label": "Multiple-regression contribution",
      "inputs_to_request": "Numbers of baseline and tested parameters; reduced/full R-squared.",
      "scope_and_checks": "Set baseline count to zero for the overall test. This does not size a prediction model for overfitting or calibration."
    },
    "reg.logistic": {
      "label": "Logistic-regression coefficient",
      "inputs_to_request": "Baseline event probability; OR; covariate prevalence/distribution; covariate-correlation R-squared.",
      "scope_and_checks": "Specific Wald approximation. Not a generic logistic interaction, prediction-development, or causal-weighting calculator."
    },
    "reg.cox": {
      "label": "Cox-regression coefficient",
      "inputs_to_request": "HR; covariate prevalence or SD; covariate R-squared; observed-event probability.",
      "scope_and_checks": "Returns events and subjects. Event probability is an input; not a time-varying coefficient or competing-risk route."
    },
    "var.one": {
      "label": "One normal variance",
      "inputs_to_request": "Null and expected SD; alternative.",
      "scope_and_checks": "Exact chi-square under normality. The fields request SD, not squared SD."
    },
    "var.two": {
      "label": "Two normal variances",
      "inputs_to_request": "Expected variance ratio; allocation; alternative.",
      "scope_and_checks": "F test under normality. Levene/Brown-Forsythe and robust spread tests are not implemented."
    },
    "surv.logrank": {
      "label": "Fixed-design logrank",
      "inputs_to_request": "HR; control survival; accrual; additional follow-up; loss hazard; allocation.",
      "scope_and_checks": "Exponential survival with uniform accrual converts event requirement to enrollment."
    },
    "sim.surv.nph": {
      "label": "Nonproportional-hazards simulation",
      "inputs_to_request": "Two time-segment HRs; change time; control median; accrual; follow-up; loss; RMST horizon if relevant.",
      "scope_and_checks": "Specific two-segment HR model. Not arbitrary hazard functions, competing risks, or recurrent-event survival."
    },
    "rates.one.poisson": {
      "label": "One Poisson rate",
      "inputs_to_request": "Null/expected rate; common exposure duration per participant; alternative.",
      "scope_and_checks": "Exact test or z test. Two-sided exact test is equal-tailed; preserve rate/time units."
    },
    "rates.two.poisson": {
      "label": "Two Poisson rates",
      "inputs_to_request": "Two rates; common exposure duration; allocation; variance convention.",
      "scope_and_checks": "No overdispersion or arbitrary unequal follow-up distribution in this template."
    },
    "rates.two.negbin": {
      "label": "Two negative-binomial rates",
      "inputs_to_request": "Control rate; expected rate ratio; dispersion; common exposure; allocation.",
      "scope_and_checks": "Supports negative-binomial count variation, not every recurrent-event or cluster model."
    },
    "sim.winratio": {
      "label": "Hierarchical composite simulation",
      "inputs_to_request": "Tier-specific event effects; score effect/tie margin; accrual; follow-up; loss.",
      "scope_and_checks": "Unmatched treatment-control pairs and a fixed supported hierarchy. Not a general recurrent-event win-ratio design."
    },
    "means.two.cluster": {
      "label": "Parallel cluster trial: means",
      "inputs_to_request": "Clusters per arm; participants per cluster; ICC; means; common SD.",
      "scope_and_checks": "t accounts for cluster degrees of freedom; z is asymptotic. No built-in arbitrary longitudinal structure."
    },
    "prop.two.cluster": {
      "label": "Parallel cluster trial: proportions",
      "inputs_to_request": "Clusters per arm; common cluster size; ICC; both response rates.",
      "scope_and_checks": "Uses a pooled z approximation. Do not present this as a finite-cluster-corrected GEE calculation."
    },
    "sim.lmm.cluster": {
      "label": "Cluster random-intercept LMM simulation",
      "inputs_to_request": "Cluster counts; size distribution; ICC; group means and SD.",
      "scope_and_checks": "Random cluster intercept; REML and chosen df method. Not a general multilevel/repeated-measures simulator."
    },
    "sim.glmm.cluster": {
      "label": "Cluster logistic mixed-model simulation",
      "inputs_to_request": "Cluster sizes/counts; group probabilities; proportion scale; ICC and ICC scale.",
      "scope_and_checks": "Population-averaged and typical-cluster probabilities differ. Match the scale and t/z inference explicitly."
    },
    "cluster.sw": {
      "label": "Analytic stepped-wedge design",
      "inputs_to_request": "Steps; clusters per step; subjects per cluster-period; ICC; effect and scale.",
      "scope_and_checks": "Hussey-Hughes random-cluster/period-effects model. Not a closed-cohort repeated-person design."
    },
    "sim.sw": {
      "label": "Stepped-wedge LMM simulation",
      "inputs_to_request": "Steps; cluster counts; cluster-period size; effect/SD; ICC; cluster autocorrelation.",
      "scope_and_checks": "Continuous LMM template. Choose cluster-period effects when required; this is not binary stepped-wedge GLMM."
    },
    "diag.auc": {
      "label": "One ROC AUC versus a benchmark",
      "inputs_to_request": "Null/expected AUC; number of cases; controls per case.",
      "scope_and_checks": "Hanley-McNeil approximation. Not a paired comparison of two AUCs, calibration assessment, or a development sample-size calculation."
    },
    "diag.sens": {
      "label": "Sensitivity and/or specificity",
      "inputs_to_request": "Disease prevalence; null/expected sensitivity and specificity; success criterion.",
      "scope_and_checks": "Uses fixed rounded subgroup counts from prevalence. If both must succeed, power is their joint success probability under that planning setup."
    },
    "agree.icc": {
      "label": "One-way ICC reliability test",
      "inputs_to_request": "Null/expected ICC; raters or repeat measurements per subject.",
      "scope_and_checks": "N counts subjects. Does not cover all two-way, fixed-rater, agreement-versus-consistency ICC definitions."
    },
    "agree.kappa": {
      "label": "Binary Cohen kappa test",
      "inputs_to_request": "Null/expected kappa; positive-rating prevalence.",
      "scope_and_checks": "Current calculation assumes equal positive-rating marginals. Weighted, multicategory, or multiple-rater kappa needs another method."
    },
    "means.one.ci": {
      "label": "Precision of one mean",
      "inputs_to_request": "Confidence level; SD; target half-width; optional tolerance probability.",
      "scope_and_checks": "For pairs enter the SD of differences and count pairs. Tolerance probability is available; it is not statistical power."
    },
    "means.two.ci": {
      "label": "Precision of two-mean difference",
      "inputs_to_request": "Confidence level; common SD; allocation; target half-width; optional tolerance probability.",
      "scope_and_checks": "Not the precision of an adjusted or clustered contrast; no separate unequal-variance interval route."
    },
    "prop.one.ci": {
      "label": "Precision of one proportion",
      "inputs_to_request": "Confidence level; anticipated probability; full width or bound; interval method.",
      "scope_and_checks": "Width is evaluated at the anticipated proportion. It is not a general assurance over all possible observed counts."
    },
    "prop.two.ci": {
      "label": "Precision of a risk difference",
      "inputs_to_request": "Confidence level; both probabilities; allocation; full width; method.",
      "scope_and_checks": "Not risk-ratio, odds-ratio, or paired-proportion precision."
    },
    "corr.ci": {
      "label": "Precision of one Pearson correlation",
      "inputs_to_request": "Confidence level; expected Pearson correlation; full width.",
      "scope_and_checks": "Not a correlation-difference interval, dependent correlation, or Spearman interval calculator."
    },
    "gs.means": {
      "label": "Group-sequential two-mean test",
      "inputs_to_request": "Effect; common known SD; allocation; information fractions; efficacy/futility boundaries.",
      "scope_and_checks": "Known-SD z framework. Not sequential t, ANCOVA, or mixed-model power by default."
    },
    "gs.props": {
      "label": "Group-sequential proportions",
      "inputs_to_request": "Both rates; allocation; pooled/unpooled variance; looks and spending/boundaries.",
      "scope_and_checks": "Sequential z approximation. Not an exact Fisher sequential design."
    },
    "gs.logrank": {
      "label": "Group-sequential survival",
      "inputs_to_request": "HR; control median; accrual/follow-up; event fractions; stopping boundaries.",
      "scope_and_checks": "Exponential enrollment conversion. Does not inherit the NPH simulator or automatically implement sequential noninferiority."
    },
    "prop.one.simon": {
      "label": "Simon single-arm two-stage design",
      "inputs_to_request": "Unacceptable and desirable response rates; alpha/power; optimal or minimax criterion.",
      "scope_and_checks": "Optimal minimizes expected N under the null; minimax minimizes maximum N. It is not a randomized two-arm design."
    },
    "means.coprimary": {
      "label": "Co-primary continuous endpoints",
      "inputs_to_request": "Effect per endpoint (standardized, or mean difference and SD); correlation between endpoints; one-sided alpha; allocation.",
      "scope_and_checks": "Every endpoint must be significant, without a multiplicity adjustment. One common correlation for all pairs; with t tests the joint power is approximated from exact single-endpoint powers."
    },
    "prop.cc": {
      "label": "Unmatched case-control study",
      "inputs_to_request": "Exposure prevalence among controls; odds ratio (or exposure among cases); controls per case.",
      "scope_and_checks": "Compares the exposure of independently sampled cases and controls. Not for individually matched controls or confounder-adjusted analyses."
    },
    "prop.cc.matched": {
      "label": "Matched case-control study",
      "inputs_to_request": "Exposure prevalence among controls; odds ratio; controls per case; correlation of exposure within matched sets.",
      "scope_and_checks": "Dupont's approximation to the conditional test on individually matched sets. Frequency matching is an unmatched design."
    },
    "reg.prediction": {
      "label": "Prediction model development",
      "inputs_to_request": "Candidate predictor parameters; anticipated R-squared or C-statistic; prevalence, event rate and follow-up, or outcome mean and SD.",
      "scope_and_checks": "Riley et al. minimum sample size for developing a model, as computed by pmsampsize. Not a sample size for validating an existing model."
    },
    "diag.paired": {
      "label": "Two diagnostic tests on the same patients",
      "inputs_to_request": "Prevalence; both tests' sensitivities and/or specificities; their dependence (independent, correlation or discordant share).",
      "scope_and_checks": "McNemar tests among patients with and without the condition. Tests given to different patients are two independent proportions."
    }
  },
};
