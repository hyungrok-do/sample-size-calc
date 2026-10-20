// t-tests (and known-variance z-tests) for one mean, paired means and two
// independent means.  Power is computed exactly from the noncentral t
// distribution: for a t statistic with ν degrees of freedom and noncentrality
// λ = δ/SE, power = P(T'_{ν,λ} > t_{1−α}) (+ the lower tail when two-sided).
#include "procedures/common.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

ParamDef test_type_param() {
  return select("test", "Test Statistic",
                {{"t", "T-test (σ estimated from the data)", {}}, {"z", "Z-test (σ known)", {}}}, "t")
      .grp(G_TEST)
      .hint("The t-test uses the sample standard deviation; the z-test assumes σ is known.");
}

// ------------------------------------------------------------------------- one mean
// Shared by the one-sample and paired designs: the statistic is
// (mean − μ0)/(σ/√N) with N − 1 degrees of freedom.
double one_mean_power(const Scenario& s, double n, double diff, double sigma) {
  const double ncp = diff * std::sqrt(n) / sigma;
  const double alpha = s("alpha");
  return s.is("test", "z") ? z_power(ncp, alpha, alt(s)) : t_power(ncp, n - 1, alpha, alt(s));
}

// Solve for N / power / the detectable difference, given the null value m0
// and alternative m1 (both on the same scale), and SD sigma.
void run_one_mean(const Scenario& s, Row& r, RunContext& ctx, double m0, const std::string& m1_key, double sigma,
                  double min_n) {
  const Alt a = alt(s);
  if (s.solve == "n") {
    const double m1 = s(m1_key);
    const double diff = m1 - m0;
    if (diff == 0) throw CalcError("the alternative mean equals the null mean; power cannot exceed alpha");
    if ((a == Alt::Upper && diff < 0) || (a == Alt::Lower && diff > 0))
      throw CalcError("the alternative mean lies on the side opposite to the one-sided alternative");
    solve::NSearch o;
    o.lo = min_n;
    o.guess = z_factor(s("alpha"), s("power"), a) * sigma * sigma / (diff * diff) + 1;
    const double n = solve::min_n([&](double x) { return one_mean_power(s, x, diff, sigma); }, s("power"), o, &ctx);
    r.set("n", n);
    put_power(s, r, one_mean_power(s, n, diff, sigma));
    r.set(m1_key, m1);
  } else if (s.solve == "power") {
    const double n = s("n");
    const double diff = s(m1_key) - m0;
    r.set("n", n);
    put_power(s, r, one_mean_power(s, n, diff, sigma));
    r.set(m1_key, s(m1_key));
  } else {  // detectable alternative mean
    const double n = s("n");
    const double sign = a == Alt::Lower ? -1 : 1;
    const double target = s("power");
    if (target <= s("alpha")) throw CalcError("the power must exceed alpha");
    const double d = solve::solve_increasing([&](double dd) { return one_mean_power(s, n, sign * dd, sigma); }, target,
                                             0, sigma, 1e9 * sigma, 1e-12, "the mean difference");
    r.set("n", n);
    put_power(s, r, one_mean_power(s, n, sign * d, sigma));
    r.set(m1_key, m0 + sign * d);
  }
}

Procedure one_sample() {
  Procedure p;
  p.id = "means.one.t";
  p.name = "Tests for One Mean (One-Sample T-Test)";
  p.category = "Means";
  p.subcategory = "One Mean";
  p.order = 1;
  p.summary =
      "Sample size and power for testing whether a population mean differs from a hypothesized value, using the "
      "one-sample t-test (σ estimated) or z-test (σ known). Power is computed exactly from the noncentral t "
      "distribution.";
  p.keywords = {"one-sample", "t-test", "z-test", "single mean", "student", "one group", "mean against a reference value", "known variance"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Sample Size)"},
                 {"power", "Power", "power", "Power"},
                 {"mu1", "Effect Size (μ1)", "mu1", "Detectable μ1"}};
  p.params = {
      alternative("Two-sided (H1: μ ≠ μ0)", "One-sided (H1: μ > μ0)", "One-sided (H1: μ < μ0)"),
      test_type_param(),
      power_param(),
      alpha_param(),
      n_param("n", "N (Sample Size)", "20", 2).sym("N").hint("Number of subjects measured in the single group."),
      dropout_param(),
      num("mu0", "μ0 (Mean under H0)", "0")
          .grp(G_EFFECT)
          .sym("μ0")
          .sugg({"0", "100"})
          .hint("The mean assumed by the null hypothesis, such as a known reference or historical value."),
      num("mu1", "μ1 (Mean under H1)", "0.5")
          .grp(G_EFFECT)
          .sym("μ1")
          .sugg({"0.5", "0.2 0.5 0.8", "1 to 2 by 0.25"})
          .hint("Actual population mean at which power is evaluated.")
          .show(when(not_solving({"mu1"}))),
      num("sigma", "σ (Standard Deviation)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .sugg({"1", "1 1.5 2", "10"})
          .hint("Standard deviation of the response in the population; take it from earlier studies or a pilot."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double sigma = s("sigma");
    run_one_mean(s, r, ctx, s("mu0"), "mu1", sigma, 2);
    const double m0 = s("mu0"), m1 = r.num("mu1");
    r.set("mu0", m0);
    r.set("delta", m1 - m0);
    r.set("sigma", sigma);
    r.set("es", (m1 - m0) / sigma);
    r.set("alpha", s("alpha"));
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : std::vector<Column>{col("n", "N", "int", "Sample size"),
                                             col("mu0", "μ0", "num", "Mean under the null hypothesis"),
                                             col("mu1", "μ1", "num", "Mean under the alternative"),
                                             col("delta", "δ = μ1 − μ0", "num", "Difference to detect"),
                                             col("sigma", "σ", "num", "Standard deviation"),
                                             col("es", "δ/σ", "num", "Standardized effect size (Cohen's d)"),
                                             col("alpha", "Alpha", "alpha", "Significance level")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A sample size of " + fint(r.num("n")) + " subjects " + achieves_singular(s, r) +
           " to detect a difference of " + fnum(r.num("delta")) + " between the null mean μ0 = " +
           fnum(r.num("mu0")) + " and the actual mean μ1 = " + fnum(r.num("mu1")) +
           ", assuming a standard deviation of " + fnum(r.num("sigma")) + ", with a " + alt_words(alt(s)) +
           " one-sample " + (s.is("test", "z") ? "z-test" : "t-test") + " at significance level α = " +
           fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"N", "Number of subjects (the test has N − 1 degrees of freedom)."},
                   {"μ0", "Population mean under the null hypothesis."},
                   {"μ1", "Actual population mean under the alternative."},
                   {"δ", "Difference μ1 − μ0 the study should detect."},
                   {"σ", "Population standard deviation."},
                   {"δ/σ", "Standardized effect size (Cohen's d)."},
                   {"Alpha", "Probability of rejecting a true null hypothesis."},
                   {"Power", "Probability of rejecting the null hypothesis when μ = μ1."}};
  p.references = {REF_CHOW, REF_JULIOUS, REF_MACHIN, REF_LENTH, REF_BK};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.curve_min = 2;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"alternative", "two"}, {"power", "0.80 0.90"}, {"alpha", "0.05"}, {"mu0", "0"},
               {"mu1", "0.3 0.5 0.7"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(one_sample);

Procedure paired() {
  Procedure p;
  p.id = "means.paired.t";
  p.name = "Tests for Paired Means (Paired T-Test)";
  p.category = "Means";
  p.subcategory = "Paired Means";
  p.order = 1;
  p.summary =
      "Sample size and power for a paired t-test (or z-test) comparing two measurements on the same subjects. The "
      "variability can be given as the SD of the paired differences or as two SDs and their correlation.";
  p.keywords = {"paired", "pre-post", "before after", "matched pairs", "repeated", "t-test", "dependent samples", "paired differences", "within-subject", "change from baseline"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Number of Pairs)"},
                 {"power", "Power", "power", "Power"},
                 {"mud", "Effect Size (μD)", "mud", "Detectable μD"}};
  p.params = {
      alternative("Two-sided (H1: μD ≠ μD0)", "One-sided (H1: μD > μD0)", "One-sided (H1: μD < μD0)"),
      test_type_param(),
      power_param(),
      alpha_param(),
      n_param("n", "N (Number of Pairs)", "20", 2)
          .sym("N")
          .hint("Number of pairs: subjects measured twice (e.g. before and after), or matched pairs of subjects."),
      dropout_param(),
      num("mud0", "μD0 (Mean Difference under H0)", "0")
          .grp(G_EFFECT)
          .sym("μD0")
          .sugg({"0"})
          .hint("Mean of the within-pair differences assumed by the null hypothesis; usually 0 (no change)."),
      num("mud", "μD (Mean Difference under H1)", "0.5")
          .grp(G_EFFECT)
          .sym("μD")
          .sugg({"0.5", "0.25 0.5 0.75"})
          .hint("Mean of the within-pair differences that you expect to be true; the study is sized to detect its "
                "distance from μD0.")
          .show(when(not_solving({"mud"}))),
      select("sdtype", "Variability Input",
             {{"sdd", "SD of the paired differences (σD)", {}},
              {"corr", "SD of each measurement and their correlation", {}}},
             "sdd")
          .grp(G_EFFECT)
          .hint("Give the SD of the paired differences directly, or the SD of each of the two measurements and "
                "their correlation, from which it is derived."),
      num("sdd", "σD (SD of Differences)", "1")
          .grp(G_EFFECT)
          .sym("σD")
          .above(0)
          .sugg({"1", "0.5 1 1.5"})
          .hint("Standard deviation of the within-pair differences; take it from earlier studies or a pilot.")
          .show(when(is("sdtype", {"sdd"}))),
      num("sd1", "σ1 (SD of First Measurement)", "1")
          .grp(G_EFFECT)
          .sym("σ1")
          .above(0)
          .hint("Standard deviation of the first measurement of a pair (e.g. before treatment) among subjects.")
          .show(when(is("sdtype", {"corr"}))),
      num("sd2", "σ2 (SD of Second Measurement)", "1")
          .grp(G_EFFECT)
          .sym("σ2")
          .above(0)
          .hint("Standard deviation of the second measurement of a pair (e.g. after treatment) among subjects.")
          .show(when(is("sdtype", {"corr"}))),
      num("rho", "ρ (Correlation between Measurements)", "0.5")
          .grp(G_EFFECT)
          .sym("ρ")
          .range(-1, 1, true, true)
          .sugg({"0.5", "0.3 0.5 0.7"})
          .hint("Correlation between the two measurements of a pair; a higher correlation makes the differences "
                "less variable and the study smaller.")
          .show(when(is("sdtype", {"corr"}))),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    double sdd;
    if (s.is("sdtype", "corr")) {
      const double a = s("sd1"), b = s("sd2"), rho = s("rho");
      const double v = a * a + b * b - 2 * rho * a * b;
      if (!(v > 0)) throw CalcError("the implied SD of the differences is zero");
      sdd = std::sqrt(v);
    } else {
      sdd = s("sdd");
    }
    run_one_mean(s, r, ctx, s("mud0"), "mud", sdd, 2);
    r.set("mud0", s("mud0"));
    r.set("sdd", sdd);
    r.set("es", (r.num("mud") - s("mud0")) / sdd);
    r.set("alpha", s("alpha"));
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : std::vector<Column>{col("n", "N", "int", "Number of pairs"),
                                             col("mud0", "μD0", "num", "Mean difference under H0"),
                                             col("mud", "μD", "num", "Mean difference under H1"),
                                             col("sdd", "σD", "num", "SD of paired differences"),
                                             col("es", "(μD−μD0)/σD", "num", "Standardized effect size"),
                                             col("alpha", "Alpha", "alpha", "Significance level")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return fint(r.num("n")) + " pairs " + achieves(s, r) + " to detect a mean paired difference of μD = " +
           fnum(r.num("mud")) + " (null value " + fnum(r.num("mud0")) + ") when the SD of the differences is " +
           fnum(r.num("sdd")) + ", using a " + alt_words(alt(s)) + " paired " +
           (s.is("test", "z") ? "z-test" : "t-test") + " with α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"N", "Number of pairs (subjects measured twice)."},
                   {"μD0", "Mean of the paired differences under H0 (usually 0)."},
                   {"μD", "Mean of the paired differences under H1."},
                   {"σD", "SD of the paired differences: √(σ1² + σ2² − 2ρσ1σ2)."},
                   {"Alpha", "Significance level."}};
  p.references = {REF_CHOW, REF_JULIOUS, REF_MACHIN};
  p.curve_param = "n";
  p.curve_label = "N (Number of Pairs)";
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.9"}, {"mud", "0.4 0.6"}, {"sdtype", "corr"}, {"sd1", "1"}, {"sd2", "1"},
               {"rho", "0.3 0.6"}};
  return p;
}
SSC_REGISTER(paired);

// ------------------------------------------------------------------------- two means
struct TwoMeanCfg {
  bool welch;
};

double two_mean_power(const Scenario& s, double n1, double n2, double diff, bool welch) {
  const double alpha = s("alpha");
  const Alt a = alt(s);
  if (welch) {
    const double v1 = s("sd1") * s("sd1") / n1, v2 = s("sd2") * s("sd2") / n2;
    const double se = std::sqrt(v1 + v2);
    if (s.is("test", "z")) return z_power(diff / se, alpha, a);
    const double df = (v1 + v2) * (v1 + v2) / (v1 * v1 / (n1 - 1) + v2 * v2 / (n2 - 1));
    return t_power(diff / se, df, alpha, a);
  }
  const double sigma = s("sigma");
  const double se = sigma * std::sqrt(1 / n1 + 1 / n2);
  return s.is("test", "z") ? z_power(diff / se, alpha, a) : t_power(diff / se, n1 + n2 - 2, alpha, a);
}

void run_two_means(const Scenario& s, Row& r, RunContext& ctx, bool welch) {
  const Alt a = alt(s);
  const double mu2 = s("mu2");
  const double var_sum = welch ? s("sd1") * s("sd1") + s("sd2") * s("sd2") : 2 * s("sigma") * s("sigma");
  if (s.solve == "n") {
    const double diff = s("mu1") - mu2;
    if (diff == 0) throw CalcError("μ1 equals μ2; power cannot exceed alpha");
    if ((a == Alt::Upper && diff < 0) || (a == Alt::Lower && diff > 0))
      throw CalcError("the difference μ1 − μ2 has the sign opposite to the one-sided alternative");
    const double guess = z_factor(s("alpha"), s("power"), a) * var_sum / (diff * diff) + 1;
    const TwoN n = solve_two_group_n(
        s, [&](double n1, double n2) { return two_mean_power(s, n1, n2, diff, welch); }, s("power"), 2, guess, true,
        ctx);
    put_two_group(r, n);
    put_power(s, r, two_mean_power(s, n.n1, n.n2, diff, welch));
    r.set("mu1", s("mu1"));
  } else if (s.solve == "power") {
    const TwoN n = two_group_n(s);
    if (n.n1 < 2 || n.n2 < 2) throw CalcError("each group needs at least 2 subjects");
    put_two_group(r, n);
    put_power(s, r, two_mean_power(s, n.n1, n.n2, s("mu1") - mu2, welch));
    r.set("mu1", s("mu1"));
  } else {
    const TwoN n = two_group_n(s);
    if (n.n1 < 2 || n.n2 < 2) throw CalcError("each group needs at least 2 subjects");
    const double sign = a == Alt::Lower ? -1 : 1;
    const double scale = std::sqrt(var_sum);
    const double d = solve::solve_increasing(
        [&](double dd) { return two_mean_power(s, n.n1, n.n2, sign * dd, welch); }, s("power"), 0, scale,
        1e9 * scale, 1e-12, "the mean difference");
    put_two_group(r, n);
    put_power(s, r, two_mean_power(s, n.n1, n.n2, sign * d, welch));
    r.set("mu1", mu2 + sign * d);
  }
  r.set("mu2", mu2);
  r.set("delta", r.num("mu1") - mu2);
  if (welch) {
    r.set("sd1", s("sd1"));
    r.set("sd2", s("sd2"));
  } else {
    r.set("sigma", s("sigma"));
    r.set("es", r.num("delta") / s("sigma"));
  }
  r.set("alpha", s("alpha"));
}

std::vector<ParamDef> two_mean_params(bool welch) {
  std::vector<ParamDef> v = {
      alternative("Two-sided (H1: μ1 ≠ μ2)", "One-sided (H1: μ1 > μ2)", "One-sided (H1: μ1 < μ2)"),
      test_type_param(), power_param(), alpha_param()};
  for (auto& q : two_group_params(2)) v.push_back(q);
  v.push_back(dropout_param());
  v.push_back(num("mu1", "μ1 (Mean of Group 1)", "1")
                  .grp(G_EFFECT)
                  .sym("μ1")
                  .sugg({"1", "0.5 1 1.5", "10"})
                  .hint("Mean of group 1 (e.g. treatment) under the alternative.")
                  .show(when(not_solving({"mu1"}))));
  v.push_back(num("mu2", "μ2 (Mean of Group 2)", "0")
                  .grp(G_EFFECT)
                  .sym("μ2")
                  .sugg({"0", "10"})
                  .hint("Mean of group 2 (e.g. control)."));
  if (welch) {
    v.push_back(num("sd1", "σ1 (SD of Group 1)", "1")
                    .grp(G_EFFECT)
                    .sym("σ1")
                    .above(0)
                    .sugg({"1", "1 1.5 2"})
                    .hint("Standard deviation of the outcome among subjects in group 1."));
    v.push_back(num("sd2", "σ2 (SD of Group 2)", "1.5")
                    .grp(G_EFFECT)
                    .sym("σ2")
                    .above(0)
                    .sugg({"1.5", "1 2"})
                    .hint("Standard deviation of the outcome among subjects in group 2; it may differ from σ1."));
  } else {
    v.push_back(num("sigma", "σ (Common Standard Deviation)", "1")
                    .grp(G_EFFECT)
                    .sym("σ")
                    .above(0)
                    .sugg({"1", "1 1.5 2"})
                    .hint("Standard deviation within each group (assumed equal)."));
  }
  return v;
}

Procedure two_sample(bool welch) {
  Procedure p;
  p.id = welch ? "means.two.welch" : "means.two.t";
  p.name = welch ? "Two-Sample T-Tests Allowing Unequal Variance (Welch)"
                 : "Two-Sample T-Tests Assuming Equal Variance";
  p.category = "Means";
  p.subcategory = "Two Independent Means";
  p.order = welch ? 2 : 1;
  p.summary = welch ? "Sample size and power for comparing two means with Welch's unequal-variance t-test. Power uses "
                      "the noncentral t distribution with Welch–Satterthwaite degrees of freedom."
                    : "Sample size and power for comparing two independent means with the pooled-variance two-sample "
                      "t-test (or a z-test when σ is known), for any allocation ratio.";
  p.keywords = welch ? std::vector<std::string>{"welch", "unequal variance", "satterthwaite", "two-sample", "t-test",
                                                 "heteroscedastic", "behrens-fisher", "unequal standard deviations"}
                     : std::vector<std::string>{"two-sample", "independent", "pooled", "student", "t-test",
                                                "parallel group", "rct", "independent samples", "two groups",
                                                "difference in means", "cohen d", "effect size"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total Sample Size)"},
                 {"power", "Power", "power", "Power"},
                 {"mu1", "Effect Size (μ1)", "mu1", "Detectable μ1"}};
  p.params = two_mean_params(welch);
  p.run = [welch](const Scenario& s, Row& r, RunContext& ctx) { run_two_means(s, r, ctx, welch); };
  p.columns = [welch](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    c.push_back(col("mu1", "μ1", "num", "Mean of group 1"));
    c.push_back(col("mu2", "μ2", "num", "Mean of group 2"));
    c.push_back(col("delta", "δ = μ1 − μ2", "num", "Difference in means"));
    if (welch) {
      c.push_back(col("sd1", "σ1", "num", "SD of group 1"));
      c.push_back(col("sd2", "σ2", "num", "SD of group 2"));
    } else {
      c.push_back(col("sigma", "σ", "num", "Common SD"));
      c.push_back(col("es", "δ/σ", "num", "Standardized effect size (Cohen's d)"));
    }
    c.push_back(col("alpha", "Alpha", "alpha", "Significance level"));
    return c;
  };
  p.statement = [welch](const Scenario& s, const Row& r) {
    std::string sd = welch ? "standard deviations of " + fnum(r.num("sd1")) + " (group 1) and " + fnum(r.num("sd2")) +
                                 " (group 2)"
                           : "a common standard deviation of " + fnum(r.num("sigma"));
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " (total " +
           fint(r.num("n")) + ") " + achieves(s, r) + " to detect a difference of " + fnum(r.num("delta")) +
           " between the group means (μ1 = " + fnum(r.num("mu1")) + ", μ2 = " + fnum(r.num("mu2")) + ") with " + sd +
           ", using a " + alt_words(alt(s)) + " " +
           (s.is("test", "z") ? "z-test" : (welch ? "Welch t-test" : "pooled two-sample t-test")) +
           " at significance level α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"N1, N2", "Number of subjects in groups 1 and 2; N = N1 + N2."},
                   {"μ1, μ2", "Population means of groups 1 and 2."},
                   {"δ", "Difference μ1 − μ2 to detect."},
                   {welch ? "σ1, σ2" : "σ", welch ? "Group standard deviations." : "Common standard deviation."},
                   {"Alpha", "Significance level."},
                   {"Power", "Probability of rejecting H0: μ1 = μ2 when the true difference is δ."}};
  if (welch)
    p.notes = {"The degrees of freedom are Welch–Satterthwaite's at the true standard deviations: accurate from about 8 "
               "subjects per group, while with 3–5 per group the power can be overstated by several points."};
  p.references = {REF_CHOW, REF_JULIOUS, REF_MACHIN, REF_LENTH, REF_BK};
  if (welch)
    p.references.push_back("Welch, B. L. (1947). The generalization of 'Student's' problem when several different "
                           "population variances are involved. Biometrika 34, 28–35.");
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.curve_min = 2;
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  if (welch)
    p.example = {{"power", "0.8 0.9"}, {"mu1", "1"}, {"mu2", "0"}, {"sd1", "1"}, {"sd2", "1 1.5 2"}};
  else
    p.example = {{"power", "0.8 0.9"}, {"mu1", "0.5 1 1.5"}, {"mu2", "0"}, {"sigma", "1.2"}};
  return p;
}

Procedure two_sample_pooled() { return two_sample(false); }
Procedure two_sample_welch() { return two_sample(true); }
SSC_REGISTER(two_sample_pooled);
SSC_REGISTER(two_sample_welch);

}  // namespace
}  // namespace ssc::procs
