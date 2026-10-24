// Nonparametric tests (via asymptotic relative efficiency, Noether's formula
// and the exact sign test) and confidence-interval precision for means.
#include "procedures/binomial_tools.hpp"
#include "procedures/common.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// Asymptotic relative efficiency of the rank tests versus the t-test (or F-test)
// for location shifts in the chosen parent distribution.
double are_of(const Scenario& s) {
  const std::string& d = s.s("dist");
  if (d == "uniform") return 1.0;
  if (d == "logistic") return stats::PI * stats::PI / 9;
  if (d == "laplace") return 1.5;
  return 3 / stats::PI;  // normal
}

ParamDef dist_param() {
  return select("dist", "Parent Distribution (for the ARE)",
                {{"normal", "Normal (ARE = 3/π ≈ 0.955)", {}},
                 {"uniform", "Uniform (ARE = 1)", {}},
                 {"logistic", "Logistic (ARE = π²/9 ≈ 1.097)", {}},
                 {"laplace", "Double exponential (ARE = 1.5)", {}}},
                "normal")
      .grp(G_TEST)
      .hint("The rank test's sample size is the t-test's divided by the asymptotic relative efficiency (ARE).");
}

// =====================================================================================
// Wilcoxon signed-rank (one sample or paired)
// =====================================================================================
Procedure signed_rank() {
  Procedure p;
  p.id = "nonpar.signedrank";
  p.name = "Wilcoxon Signed-Rank Test (One Sample or Paired)";
  p.category = "Nonparametric";
  p.subcategory = "One Sample / Paired";
  p.order = 1;
  p.summary =
      "Sample size and power for the Wilcoxon signed-rank test of a location shift, obtained from the one-sample "
      "t-test with the sample size scaled by the asymptotic relative efficiency for the assumed distribution.";
  p.keywords = {"wilcoxon", "signed rank", "nonparametric", "paired", "median", "rank test", "wilcoxon signed-rank", "paired nonparametric", "one-sample nonparametric"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}};
  p.params = {alternative("Two-sided", "One-sided (shift > 0)", "One-sided (shift < 0)"),
              dist_param(),
              power_param(),
              alpha_param(),
              n_param("n", "N (Sample Size / Pairs)", "30", 2)
                  .sym("N")
                  .hint("Number of subjects, or of pairs when the test is applied to within-pair differences."),
              dropout_param(),
              num("delta", "δ (Location Shift from H0)", "0.5")
                  .grp(G_EFFECT)
                  .sym("δ")
                  .sugg({"0.5", "0.25 0.5 0.75"})
                  .hint("How far the center of the observations (or of the paired differences) truly lies from the "
                        "value under H0, in the units of the outcome."),
              num("sigma", "σ (Standard Deviation)", "1")
                  .grp(G_EFFECT)
                  .sym("σ")
                  .above(0)
                  .sugg({"1"})
                  .hint("Standard deviation of the observations (or of the paired differences); take it from "
                        "earlier studies or a pilot.")};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double are = are_of(s), d = s("delta"), sigma = s("sigma"), alpha = s("alpha");
    const Alt a = alt(s);
    auto P = [&](double n) {
      const double ne = n * are;
      return ne > 1 ? t_power(d * std::sqrt(ne) / sigma, ne - 1, alpha, a) : 0.0;
    };
    double n;
    if (s.solve == "n") {
      if (d == 0) throw CalcError("δ = 0: there is no shift to detect");
      if ((a == Alt::Upper && d < 0) || (a == Alt::Lower && d > 0))
        throw CalcError("δ has the sign opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = 2;
      o.guess = (z_factor(alpha, s("power"), a) * sigma * sigma / (d * d) + 1) / are;
      n = solve::min_n(P, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    r.set("n", n);
    put_power(s, r, P(n));
    r.set("are", are);
    r.set("delta", d);
    r.set("sigma", sigma);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("are", "ARE", "num"), col("delta", "δ", "num"), col("sigma", "σ", "num"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) +
           " to detect a location shift of " + fnum(r.num("delta")) + " (σ = " + fnum(r.num("sigma")) +
           ") with a " + alt_words(alt(s)) + " Wilcoxon signed-rank test at α = " + fnum(r.num("alpha")) +
           ", assuming an asymptotic relative efficiency of " + fnum(r.num("are"), 3) + ".";
  };
  p.definitions = {{"ARE", "Asymptotic relative efficiency of the signed-rank test relative to the t-test."},
                   {"δ", "Shift in location (mean/median) from the null value."}};
  p.references = {"Lehmann, E. L. (1975). Nonparametrics: Statistical Methods Based on Ranks. Holden-Day.",
                  "Al-Sunduqchi, M. S. (1990). Determining the Appropriate Sample Size for Inferences Based on the "
                  "Wilcoxon Statistics. PhD thesis, University of Wyoming."};
  p.notes = {"The asymptotic relative efficiency is a large-sample result. For normal data the power agrees with "
             "simulation of the exact test; for flat or heavy-tailed data at small N it can be clearly lower "
             "(uniform data, δ = 0.5, N = 34: 0.81, against 0.73 by simulation; Laplace data, N = 7: 0.83, against "
             "0.64)."};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"delta", "0.5"}, {"dist", "normal"}};
  return p;
}
SSC_REGISTER(signed_rank);

// =====================================================================================
// Mann–Whitney / Wilcoxon rank-sum
// =====================================================================================
Procedure mann_whitney() {
  Procedure p;
  p.id = "nonpar.mannwhitney";
  p.name = "Mann–Whitney U / Wilcoxon Rank-Sum Test";
  p.category = "Nonparametric";
  p.subcategory = "Two Independent Samples";
  p.order = 1;
  p.summary =
      "Sample size and power for the Mann–Whitney (Wilcoxon rank-sum) test, either from the two-sample t-test with "
      "an ARE adjustment for a location shift, or with Noether's formula from P(Y > X).";
  p.keywords = {"mann-whitney", "wilcoxon rank sum", "nonparametric", "ordinal", "rank test", "noether", "u test", "two-sample nonparametric", "ordinal outcome"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {alternative("Two-sided", "One-sided (group 1 larger)", "One-sided (group 1 smaller)"),
              select("method", "Method",
                     {{"are", "Location shift with ARE adjustment", {}}, {"noether", "Noether: P(X1 > X2)", {}}}, "are")
                  .grp(G_TEST)
                  .hint("Describe the effect as a shift between the two distributions, adjusted for their assumed "
                        "shape, or (Noether) as the chance that a group 1 value exceeds a group 2 value."),
              dist_param().show(when(is("method", {"are"}))), power_param(), alpha_param()};
  for (auto& q : two_group_params(2, "30")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  p.params.push_back(num("delta", "δ (Location Shift μ1 − μ2)", "0.5")
                         .grp(G_EFFECT)
                         .sym("δ")
                         .sugg({"0.5", "0.25 0.5 0.75"})
                         .hint("The true difference between the centers of group 1 and group 2, in the units of the "
                               "outcome.")
                         .show(when(is("method", {"are"}))));
  p.params.push_back(num("sigma", "σ (Standard Deviation)", "1")
                         .grp(G_EFFECT)
                         .sym("σ")
                         .above(0)
                         .hint("Standard deviation of the outcome within each group, assumed the same in both groups.")
                         .show(when(is("method", {"are"}))));
  p.params.push_back(num("pgt", "P(X1 > X2)", "0.64")
                         .grp(G_EFFECT)
                         .sym("P(X1>X2)")
                         .range(0, 1, true, true)
                         .sugg({"0.6", "0.64", "0.7"})
                         .hint("The chance that a random subject from group 1 has a larger value than one from group "
                               "2; 0.5 means no difference.")
                         .show(when(is("method", {"noether"}))));
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double alpha = s("alpha");
    const Alt a = alt(s);
    const bool noether = s.is("method", "noether");
    std::function<double(double, double)> P;
    double guess = NAN;
    // with no effect the power is alpha, which no sample size raises
    if (noether) {
      const double pg = s("pgt");
      if (pg == 0.5 && s.solve == "n") throw CalcError("P(X1 > X2) = 0.5 means no effect");
      if (s.solve == "n" && ((a == Alt::Upper && (pg - 0.5) < 0) || (a == Alt::Lower && (pg - 0.5) > 0)))
        throw CalcError("P(X1 > X2) lies on the side opposite to the one-sided alternative");
      P = [=](double n1, double n2) {
        const double N = n1 + n2, c = n1 / N;
        return z_power(std::sqrt(12 * N * c * (1 - c)) * (pg - 0.5), alpha, a);
      };
      if (s.solve == "n") guess = z_factor(alpha, s("power"), a) / (12 * 0.25 * (pg - 0.5) * (pg - 0.5)) / 2;
    } else {
      const double are = are_of(s), d = s("delta"), sigma = s("sigma");
      if (d == 0 && s.solve == "n") throw CalcError("δ = 0: there is no shift to detect");
      if (s.solve == "n" && ((a == Alt::Upper && d < 0) || (a == Alt::Lower && d > 0)))
        throw CalcError("δ has the sign opposite to the one-sided alternative");
      P = [=](double n1, double n2) {
        const double e1 = n1 * are, e2 = n2 * are;
        return t_power(d / (sigma * std::sqrt(1 / e1 + 1 / e2)), e1 + e2 - 2, alpha, a);
      };
      if (s.solve == "n") guess = z_factor(alpha, s("power"), a) * 2 * sigma * sigma / (d * d) / are + 1;
      r.set("are", are);
      r.set("delta", d);
      r.set("sigma", sigma);
    }
    TwoN n = s.solve == "n" ? solve_two_group_n(s, P, s("power"), 2, guess, true, ctx) : two_group_n(s);
    put_two_group(r, n);
    put_power(s, r, P(n.n1, n.n2));
    if (noether) r.set("pgt", s("pgt"));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    if (s.is("method", "noether")) {
      c.push_back(col("pgt", "P(X1>X2)", "num"));
    } else {
      for (auto x : {col("are", "ARE", "num"), col("delta", "δ", "num"), col("sigma", "σ", "num")}) c.push_back(x);
    }
    c.push_back(col("alpha", "Alpha", "alpha"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    std::string eff = s.is("method", "noether") ? "a probability P(X1 > X2) of " + fnum(r.num("pgt"))
                                                : "a location shift of " + fnum(r.num("delta")) + " (σ = " +
                                                      fnum(r.num("sigma")) + ", ARE " + fnum(r.num("are"), 3) + ")";
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
           " to detect " + eff + " with a " + alt_words(alt(s)) + " Mann–Whitney test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"P(X1>X2)", "Probability that a random observation from group 1 exceeds one from group 2."},
                   {"ARE", "Asymptotic relative efficiency relative to the two-sample t-test."}};
  p.references = {"Noether, G. E. (1987). Sample size determination for some common nonparametric tests. Journal of "
                  "the American Statistical Association 82, 645–647.",
                  "Lehmann, E. L. (1975). Nonparametrics: Statistical Methods Based on Ranks. Holden-Day."};
  p.notes = {"Both methods are large-sample results. With the ARE adjustment the power agrees with simulation for "
             "normal data and can be lower for other shapes (uniform data, δ = 0.5, 64 per group: 0.80, against "
             "0.77 by simulation); Noether's formula is conservative for large effects (P(X1 > X2) = 0.85, 15 per "
             "group: 0.91, against 0.96). The simulation method for skewed or heavy-tailed data checks a design "
             "directly."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"delta", "0.5"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(mann_whitney);

// =====================================================================================
// Kruskal–Wallis (ARE-adjusted one-way ANOVA)
// =====================================================================================
Procedure kruskal() {
  Procedure p;
  p.id = "nonpar.kruskal";
  p.name = "Kruskal–Wallis Test";
  p.category = "Nonparametric";
  p.subcategory = "Several Independent Samples";
  p.order = 1;
  p.summary = "Sample size and power for the Kruskal–Wallis test from the one-way ANOVA F-test with the per-group "
              "sample size scaled by the asymptotic relative efficiency.";
  p.keywords = {"kruskal-wallis", "nonparametric anova", "rank test", "k groups", "kruskal", "anova on ranks", "several groups"};
  p.solve_for = {{"n", "Sample Size", "n", "n (Per Group)"}, {"power", "Power", "power", "Power"}};
  p.params = {dist_param(),
              power_param(),
              alpha_param(),
              n_param("n", "n (Per-Group Sample Size)", "20", 2)
                  .sym("n")
                  .hint("Number of subjects in each group; all groups have the same size."),
              dropout_param(),
              vec("means", "Group Means (Location Parameters)", "0 0.5 1")
                  .grp(G_EFFECT)
                  .sym("Means")
                  .sugg({"0 0.5 1", "0 0 1"})
                  .hint("The expected center of each group, one value per group; their spread relative to σ sets "
                        "the effect. Separate alternative sets with ';'."),
              num("sigma", "σ (Standard Deviation)", "1")
                  .grp(G_EFFECT)
                  .sym("σ")
                  .above(0)
                  .hint("Standard deviation of the outcome within each group, assumed the same in all groups.")};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const auto& mu = s.v("means");
    const int k = static_cast<int>(mu.size());
    if (k < 2) throw CalcError("enter at least two group means");
    // from the differences from the first mean: a mean of the means themselves rounds when they are large next to
    // their spread, and leaves a spread that is not there
    double mbar = 0;
    for (double m : mu) mbar += (m - mu[0]) / k;
    double ss = 0;
    for (double m : mu) ss += (m - mu[0] - mbar) * (m - mu[0] - mbar);
    const double sigma = s("sigma"), alpha = s("alpha"), are = are_of(s);
    const double f2 = ss / k / (sigma * sigma);
    // with no effect the power is alpha, which no sample size raises
    if (!(f2 > 0) && s.solve == "n") throw CalcError("the group means are all equal");
    if (!std::isfinite(f2))
      throw CalcError("the effect size is too large to compute with (σ is tiny next to the spread of the means)");
    auto P = [&](double n) {
      const double ne = n * are, N = ne * k;
      const double df2 = N - k;
      if (!(df2 > 0)) return 0.0;
      return stats::pnf(stats::qf(alpha, k - 1, df2, false), k - 1, df2, N * f2, false);
    };
    double n;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = 2;
      o.guess = (z_factor(alpha, s("power"), Alt::Two) + 2.0 * (k - 1)) / (k * f2) / are + 1;
      n = solve::min_n(P, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    r.set("n", n);
    r.set("N", n * k);
    put_power(s, r, P(n));
    r.set("k", k);
    r.set("f", std::sqrt(f2));
    r.set("are", are);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "n", "int", "Per group"), col("N", "N", "int"), col("k", "k", "int"),
                   col("f", "f", "num", "Cohen's f of the location parameters"), col("are", "ARE", "num"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return fint(r.num("n")) + " subjects per group (" + fint(r.num("N")) + " in total, " + fint(r.num("k")) +
           " groups) " + achieves(s, r) + " to detect differences with effect size f = " + fnum(r.num("f")) +
           " using the Kruskal–Wallis test at α = " + fnum(r.num("alpha")) + " (ARE " + fnum(r.num("are"), 3) + ").";
  };
  p.references = {"Lehmann, E. L. (1975). Nonparametrics: Statistical Methods Based on Ranks. Holden-Day.",
                  "Mahoney, M. & Magel, R. (1996). Estimation of the power of the Kruskal-Wallis test. Biometrical "
                  "Journal 38, 613–630."};
  p.notes = {"The asymptotic relative efficiency is a large-sample result. For normal data the power agrees with "
             "simulation; for heavy-tailed data at small n it can be clearly lower (Laplace data, means 0 0 1, 11 "
             "per group: 0.83, against 0.72 by simulation)."};
  p.curve_param = "n";
  p.curve_label = "n (Per-Group Sample Size)";
  p.dropout_groups = {"N"};
  p.example_solve = "n";
  p.example = {{"power", "0.8"}, {"means", "0 0.5 1; 0 0 1"}};
  return p;
}
SSC_REGISTER(kruskal);

// =====================================================================================
// Sign test (exact binomial)
// =====================================================================================
double sign_power(int n, double p, double alpha, Alt a, double* actual = nullptr) {
  const double side = tail_alpha(alpha, a);
  const int cu = a == Alt::Lower ? n + 1 : upper_critical(n, 0.5, side);
  const int cl = a == Alt::Upper ? -1 : lower_critical(n, 0.5, side);
  if (actual) *actual = two_tail_prob(n, 0.5, cl, cu);
  return two_tail_prob(n, p, cl, cu);
}

Procedure sign_test() {
  Procedure p;
  p.id = "nonpar.sign";
  p.name = "Sign Test (Exact Binomial)";
  p.category = "Nonparametric";
  p.subcategory = "One Sample / Paired";
  p.order = 2;
  p.summary = "Exact power and sample size for the sign test of a median (or of paired differences), from the "
              "probability that an observation exceeds the hypothesized median.";
  p.keywords = {"sign test", "median", "binomial", "paired", "nonparametric", "exact", "paired sign test"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}};
  p.params = {alternative("Two-sided (H1: P ≠ 0.5)", "One-sided (H1: P > 0.5)", "One-sided (H1: P < 0.5)"),
              power_param(),
              alpha_param(),
              n_param("n", "N (Sample Size / Pairs)", "30", 1)
                  .sym("N")
                  .hint("Number of subjects, or of pairs when the signs of within-pair differences are counted."),
              dropout_param(),
              select("es", "Effect Size Input",
                     {{"p", "P(X > median under H0)", {}}, {"shift", "Normal shift δ/σ", {}}}, "p")
                  .grp(G_EFFECT)
                  .hint("Give the chance of an observation above the median under H0 directly, or derive it from a "
                        "standardized shift of normally distributed data."),
              num("p", "P (Probability of a Positive Sign)", "0.7")
                  .grp(G_EFFECT)
                  .sym("P")
                  .range(0, 1, true, true)
                  .sugg({"0.65", "0.7", "0.75"})
                  .hint("The chance that an observation (or a paired difference) lies above the median under H0; "
                        "0.5 means no effect.")
                  .show(when(is("es", {"p"}))),
              num("es_d", "δ/σ (Standardized Shift)", "0.5")
                  .grp(G_EFFECT)
                  .sym("δ/σ")
                  .sugg({"0.5"})
                  .hint("How far the center lies from the median under H0, in standard deviations, for normal data; "
                        "0.5 gives P ≈ 0.69.")
                  .show(when(is("es", {"shift"})))};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double pp = s.is("es", "shift") ? stats::pnorm(s("es_d")) : s("p");
    const double alpha = s("alpha");
    const Alt a = alt(s);
    int n;
    if (s.solve == "n") {
      if (pp == 0.5) throw CalcError("P = 0.5 means no effect");
      if ((a == Alt::Upper && (pp - 0.5) < 0) || (a == Alt::Lower && (pp - 0.5) > 0))
        throw CalcError("P lies on the side opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = 1;
      o.monotone = false;
      const double za = stats::qnorm(tail_alpha(alpha, a), false), zb = stats::qnorm(s("power"));
      o.guess = std::ceil(std::pow(za * 0.5 + zb * std::sqrt(pp * (1 - pp)), 2) / ((pp - 0.5) * (pp - 0.5)));
      n = static_cast<int>(solve::min_n([&](double x) { return sign_power(static_cast<int>(x), pp, alpha, a); }, s("power"), o, &ctx));
    } else {
      n = static_cast<int>(s("n"));
    }
    double act = 0;
    const double pw = sign_power(n, pp, alpha, a, &act);
    r.set("n", n);
    put_power(s, r, pw);
    r.set("p", pp);
    r.set("alpha", alpha);
    r.set("actual_alpha", act);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("p", "P", "num", "Probability of a positive sign under H1"),
                   col("alpha", "Target Alpha", "alpha"), col("actual_alpha", "Actual Alpha", "prob")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) +
           " to detect a probability of a positive sign of " + fnum(r.num("p")) + " (versus 0.5) with a " +
           alt_words(alt(s)) + " exact sign test at a target α of " + fnum(r.num("alpha")) + " (actual " +
           text::format_fixed(r.num("actual_alpha"), 4) + ").";
  };
  p.references = {"Noether, G. E. (1987). Sample size determination for some common nonparametric tests. Journal of "
                  "the American Statistical Association 82, 645–647.",
                  REF_MACHIN};
  p.notes = {"Exact power is saw-toothed in N; the reported N is the first that reaches the target."};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.curve_min = 1;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"p", "0.65 0.7 0.75"}};
  return p;
}
SSC_REGISTER(sign_test);

// =====================================================================================
// Confidence interval for one mean / for the difference of two means
// =====================================================================================
// Half-width of the interval for a given n (and the probability that the
// realized half-width is at most the target when a tolerance probability is used).
struct CiMean {
  double halfwidth;  // expected / planning half-width
  double prob;       // P(realized half-width ≤ target)
};

CiMean ci_mean(const Scenario& s, double se_factor, double df, double target) {
  using namespace stats;
  const double conf = s("conf"), sigma = s("sigma");
  const bool two = s.is("sides", "two");
  const double a = two ? (1 - conf) / 2 : 1 - conf;
  CiMean c;
  if (s.is("known", "z")) {
    c.halfwidth = qnorm(a, false) * sigma * se_factor;
    c.prob = c.halfwidth <= target ? 1.0 : 0.0;
    return c;
  }
  const double t = qt(a, df, false);
  c.halfwidth = t * sigma * se_factor;
  // realized half-width t·S·f ≤ target  ⇔  χ²_df ≤ df·(target/(t σ f))²
  c.prob = pchisq(df * std::pow(target / (t * sigma * se_factor), 2), df);
  return c;
}

void run_ci_mean(const Scenario& s, Row& r, RunContext& ctx, bool two_groups) {
  if (!s.is("sides", "two")) check_one_sided_conf(s("conf"));
  const double gamma = s.get("tolprob", 0);
  // group 2: n1 × the ratio, rounded up (at least 1: a tiny ratio left it empty, and the interval infinitely wide)
  auto n2_of = [&](double n) { return solve::round_up(s("r") * n); };
  auto factor = [&](double n) {
    if (!two_groups) return std::sqrt(1 / n);
    return std::sqrt(1 / n + 1 / n2_of(n));
  };
  auto dfof = [&](double n) { return two_groups ? n + n2_of(n) - 2 : n - 1; };
  double n;
  if (s.solve == "n") {
    const double target = s("e");
    auto ok = [&](double x) {
      const CiMean c = ci_mean(s, factor(x), dfof(x), target);
      if (gamma > 0 && !s.is("known", "z")) return c.prob >= gamma ? 1.0 : 0.0;
      return c.halfwidth <= target ? 1.0 : 0.0;
    };
    solve::NSearch o;
    o.lo = 2;
    const double z = stats::qnorm((1 - s("conf")) / (s.is("sides", "two") ? 2 : 1), false);
    o.guess = z * z * s("sigma") * s("sigma") / (target * target) * (two_groups ? 1 + 1 / s("r") : 1) + 2;
    o.goal = "the target half-width";
    n = solve::min_n(ok, 0.5, o, &ctx);
    r.set("target_e", target);
  } else {
    n = s("n1");
  }
  const CiMean c = ci_mean(s, factor(n), dfof(n), s.solve == "n" ? s("e") : NAN);
  if (two_groups) {
    r.set("n1", n);
    r.set("n2", n2_of(n));
    r.set("n", n + n2_of(n));
  } else {
    r.set("n", n);
  }
  r.set("e", c.halfwidth);
  if (s.solve == "n" && gamma > 0 && !s.is("known", "z")) r.set("prob", c.prob);
  r.set("sigma", s("sigma"));
  r.set("conf", s("conf"));
}

std::vector<ParamDef> ci_mean_params(bool two) {
  std::vector<ParamDef> v = {
      num("conf", "Confidence Level (1 − α)", "0.95")
          .grp(G_TEST)
          .sym("Confidence")
          .range(0, 1, true, true)
          .sugg({"0.95", "0.9 0.95 0.99"})
          .hint("How often intervals computed this way contain the true value; 0.95 gives a 95% interval."),
      select("sides", "Interval Type", {{"two", "Two-sided", {}}, {"one", "One-sided (upper or lower bound)", {}}},
             "two")
          .grp(G_TEST)
          .hint("A two-sided interval has a lower and an upper limit; a one-sided interval gives only an upper or "
                "only a lower bound."),
      select("known", "Standard Deviation", {{"t", "Estimated from the data (t)", {}}, {"z", "Known (z)", {}}}, "t")
          .grp(G_TEST)
          .hint("Whether the interval will use the standard deviation estimated from the data (a t interval, the "
                "usual case) or a known σ (a z interval)."),
      num(two ? "n1" : "n1", two ? "N1 (Group 1 Sample Size)" : "N (Sample Size)", "50")
          .grp(G_N)
          .sym(two ? "N1" : "N")
          .at_least(2)
          .whole()
          .sugg({"50", "20 to 100 by 20"})
          .hint(two ? "Number of subjects in group 1; group 2 gets R × N1, rounded up."
                    : "Number of subjects whose outcome is measured.")
          .show(when(not_solving({"n"}))),
  };
  if (two)
    v.push_back(num("r", "R (Allocation Ratio N2/N1)", "1")
                    .grp(G_N)
                    .sym("R")
                    .above(0)
                    .sugg({"1", "2"})
                    .hint("Size of group 2 relative to group 1 (N2 = R × N1, rounded up); 1 gives equal groups."));
  v.push_back(num("e", "E (Half-Width: Distance from Estimate to Limit)", "0.5")
                  .grp(G_EFFECT)
                  .sym("E")
                  .above(0)
                  .sugg({"0.5", "0.25 0.5 1"})
                  .hint(std::string("The precision you want: the distance from the estimated ") +
                        (two ? "difference in means" : "mean") +
                        " to the limit (half the width of a two-sided interval), in the units of the outcome.")
                  .show(when(solving({"n"}))));
  v.push_back(num("sigma", two ? "σ (Common Standard Deviation)" : "σ (Standard Deviation)", "1")
                  .grp(G_EFFECT)
                  .sym("σ")
                  .above(0)
                  .sugg({"1", "2"})
                  .hint(two ? "Standard deviation of the outcome within each group, assumed the same in both groups."
                            : "Standard deviation of the outcome; take it from earlier studies or a pilot."));
  v.push_back(num("tolprob", "Tolerance Probability (0 = use the expected width)", "0")
                  .grp(G_EFFECT)
                  .sym("γ")
                  .range(0, 1, false, true)
                  .sugg({"0", "0.8", "0.9"})
                  .hint("With 0, N gives the half-width E when the sample SD equals σ. A positive value, such as 0.9, "
                        "makes the realized half-width at most E with that probability, as the sample SD varies.")
                  .show(when({solving({"n"}), is("known", {"t"})})));
  return v;
}

/// "; the half-width is at most 0.5 with probability 90.5% (target 90%)" when the size meets a tolerance
/// probability rather than the expected width.
std::string tolerance_words(const Scenario& s, const Row& r) {
  if (!(s.solve == "n" && s.get("tolprob", 0) > 0 && std::isfinite(r.num("prob")))) return "";
  return "; the half-width is at most " + fnum(r.num("target_e")) + " with probability " + power_text(r.num("prob")) +
         " (target " + level_text(s("tolprob")) + ")";
}

Procedure ci_one_mean() {
  Procedure p;
  p.id = "means.one.ci";
  p.name = "Confidence Interval for One Mean";
  p.category = "Confidence Intervals";
  p.subcategory = "Means";
  p.order = 1;
  p.summary = "Sample size for a confidence interval of a mean with a desired half-width (precision), optionally "
              "with a tolerance probability that the realized half-width does not exceed the target.";
  p.keywords = {"precision", "margin of error", "half-width", "confidence interval", "mean", "estimate a mean"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"width", "Half-Width", "e", "Half-Width"}};
  p.params = ci_mean_params(false);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) { run_ci_mean(s, r, ctx, false); };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c = {col("conf", "Confidence", "num"), col("n", "N", "int")};
    if (s.solve == "n") c.push_back(col("target_e", "Target E", "num"));
    c.push_back(col("e", "E (Half-Width)", "num", "The half-width when the sample standard deviation equals σ"));
    if (s.solve == "n" && s.get("tolprob", 0) > 0) c.push_back(col("prob", "P(E ≤ Target)", "prob"));
    c.push_back(col("sigma", "σ", "num"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A sample of " + fint(r.num("n")) + " gives a " + (s.is("sides", "two") ? "two-sided " : "one-sided ") +
           level_text(r.num("conf")) + " confidence interval for the mean with a half-width of " + fnum(r.num("e")) +
           " when the standard deviation is " + fnum(r.num("sigma")) + tolerance_words(s, r) + ".";
  };
  p.definitions = {{"E", "Distance from the sample mean to a confidence limit (half the width of a two-sided interval)."},
                   {"P(E ≤ Target)", "Probability that the realized half-width is no larger than the target."}};
  p.references = {"Kupper, L. L. & Hafner, K. B. (1989). How appropriate are popular sample size formulas? The American "
                  "Statistician 43, 101–105.",
                  REF_MACHIN};
  p.example_solve = "n";
  p.example = {{"e", "0.25 0.5"}, {"sigma", "1 2"}};
  return p;
}
SSC_REGISTER(ci_one_mean);

Procedure ci_two_means() {
  Procedure p;
  p.id = "means.two.ci";
  p.name = "Confidence Interval for the Difference of Two Means";
  p.category = "Confidence Intervals";
  p.subcategory = "Means";
  p.order = 2;
  p.summary = "Sample size for a confidence interval of the difference between two independent means with a desired "
              "half-width, assuming a common standard deviation.";
  p.keywords = {"precision", "difference in means", "confidence interval", "two groups", "mean difference interval"};
  p.solve_for = {{"n", "Sample Size", "n1", "N1"}, {"width", "Half-Width", "e", "Half-Width"}};
  p.params = ci_mean_params(true);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) { run_ci_mean(s, r, ctx, true); };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c = {col("conf", "Confidence", "num"), col("n1", "N1", "int"), col("n2", "N2", "int"), col("n", "N", "int")};
    if (s.solve == "n") c.push_back(col("target_e", "Target E", "num"));
    c.push_back(col("e", "E (Half-Width)", "num", "The half-width when the sample standard deviation equals σ"));
    if (s.solve == "n" && s.get("tolprob", 0) > 0) c.push_back(col("prob", "P(E ≤ Target)", "prob"));
    c.push_back(col("sigma", "σ", "num"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " give a " +
           (s.is("sides", "two") ? "two-sided " : "one-sided ") + level_text(r.num("conf")) +
           " confidence interval for the difference in means with a half-width of " + fnum(r.num("e")) +
           " when the common standard deviation is " + fnum(r.num("sigma")) + tolerance_words(s, r) + ".";
  };
  p.references = {"Kupper, L. L. & Hafner, K. B. (1989). How appropriate are popular sample size formulas? The American "
                  "Statistician 43, 101–105.",
                  REF_MACHIN};
  p.example_solve = "n";
  p.example = {{"e", "0.25 0.5"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(ci_two_means);

}  // namespace
}  // namespace ssc::procs
