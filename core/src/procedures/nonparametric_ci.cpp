// Nonparametric tests (via asymptotic relative efficiency, Noether's formula
// and the exact sign test) and confidence-interval precision for means.
#include <array>

#include "procedures/binomial_tools.hpp"
#include "procedures/common.hpp"
#include "procedures/derive_misc.hpp"

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

// The ARE as a step: its value for the parent distribution.
Step are_step(const Scenario& s) {
  using namespace derive;
  const std::string& d = s.s("dist");
  const double are = are_of(s);
  if (d == "uniform") return {"Asymptotic relative efficiency", "ARE = 1 for uniform data", "ARE = " + dv(are)};
  if (d == "logistic") return {"Asymptotic relative efficiency", "ARE = π²/9 for logistic data", "π²/9 = " + dv(are)};
  if (d == "laplace")
    return {"Asymptotic relative efficiency", "ARE = 3/2 for double exponential data", "3/2 = " + dv(are)};
  return {"Asymptotic relative efficiency", "ARE = 3/π for normal data", "3/π = " + dv(are)};
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
// The signed-rank test's power: the one-sample t-test's with N × ARE subjects.
double signed_rank_power(double n, double d, double sigma, double are, double alpha, Alt a) {
  const double ne = n * are;
  return ne > 1 ? t_power(d * std::sqrt(ne) / sigma, ne - 1, alpha, a) : 0.0;
}

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
    auto P = [&](double n) { return signed_rank_power(n, d, sigma, are, alpha, a); };
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
  // The one-sample t-test's power with N × ARE subjects: solving for N, at the size found and one below, beside the
  // normal closed form; solving for power, at the size given.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const Alt a = alt(s);
    const bool solve_n = s.solve == "n";
    const double n = r.num("n"), are = r.num("are"), d = r.num("delta"), sigma = r.num("sigma"), alpha = s("alpha");
    const double ne = n * are, ncp = d * std::sqrt(ne) / sigma, pw = r.num("power");
    const std::string N = size(n, solve_n);
    Derivation out;
    out.find = solve_n ? find_size("N", s("power")) : find_power("with N = " + N);
    out.steps.push_back(are_step(s));
    if (solve_n) {  // the t-test's normal closed form divided by the ARE, beside the exact answer below
      const double za = stats::qnorm(tail_alpha(alpha, a), false), zb = stats::qnorm(s("power"));
      out.steps.push_back({"Normal approximation",
                           "N ≈ " + z_sum_symbol(a) + "² σ² / (δ² ARE)" + FOR_COMPARISON,
                           z_sum_numbers(alpha, s("power"), a) + "² × " + gv(sigma) + "² / (" + gvs(d) + "² × " +
                               dv(are) + ") = " + dv((za + zb) * (za + zb) * sigma * sigma / (d * d * are))});
    }
    out.steps.push_back({"Equivalent t-test size", "Ne = N × ARE", N + " × " + dv(are) + " = " + dv(ne)});
    out.steps.push_back(
        {"Noncentrality", "λ = δ √Ne / σ", gvs(d) + " × √" + dv(ne) + " / " + gv(sigma) + " = " + dv(ncp)});
    if (ne > 1) {
      out.steps.push_back(t_critical(alpha, ne - 1, a, "Ne − 1"));
      out.steps.push_back(t_power_step(ncp, ne - 1, alpha, a));
    } else {
      out.steps.push_back({"Power", "0 when Ne ≤ 1 (no degrees of freedom)",
                           "Ne = " + dv(ne) + " ≤ 1: power = " + dvp(pw)});
    }
    if (solve_n) {
      out.steps.push_back(
          n_search([&](double x) { return signed_rank_power(x, d, sigma, are, alpha, a); }, n, s("power"), 2));
      out.answer = "N = " + dvi(n) + ", with power " + dvp(pw, s("power"));
    } else {
      out.answer = "power = " + dvp(pw);
    }
    return out;
  };
  p.definitions = {{"ARE", "Asymptotic relative efficiency of the signed-rank test relative to the t-test."},
                   {"δ", "Shift in location (mean/median) from the null value."}};
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
// The Mann–Whitney test's power for a location shift: the two-sample t-test's with each group × ARE.
double mw_shift_power(double n1, double n2, double d, double sigma, double are, double alpha, Alt a) {
  const double e1 = n1 * are, e2 = n2 * are;
  return t_power(d / (sigma * std::sqrt(1 / e1 + 1 / e2)), e1 + e2 - 2, alpha, a);
}

// Noether's power from P(X1 > X2) = pg.
double mw_noether_power(double n1, double n2, double pg, double alpha, Alt a) {
  const double N = n1 + n2, c = n1 / N;
  return z_power(std::sqrt(12 * N * c * (1 - c)) * (pg - 0.5), alpha, a);
}

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
      P = [=](double n1, double n2) { return mw_noether_power(n1, n2, pg, alpha, a); };
      if (s.solve == "n") guess = z_factor(alpha, s("power"), a) / (12 * 0.25 * (pg - 0.5) * (pg - 0.5)) / 2;
    } else {
      const double are = are_of(s), d = s("delta"), sigma = s("sigma");
      if (d == 0 && s.solve == "n") throw CalcError("δ = 0: there is no shift to detect");
      if (s.solve == "n" && ((a == Alt::Upper && d < 0) || (a == Alt::Lower && d > 0)))
        throw CalcError("δ has the sign opposite to the one-sided alternative");
      P = [=](double n1, double n2) { return mw_shift_power(n1, n2, d, sigma, are, alpha, a); };
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
  // A shift: the two-sample t-test's power with each group × ARE; Noether: a z-test. Solving for the group sizes,
  // the power at the sizes found and one design below, beside the closed form where group 2 follows group 1.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const Alt a = alt(s);
    const bool solve_n = s.solve == "n", noether = s.is("method", "noether");
    const blocks::TwoN n{r.num("n1"), r.num("n2")};
    const double n1 = n.n1, n2 = n.n2, alpha = s("alpha"), pw = r.num("power");
    const auto [N1, N2] = group_sizes(s, n, solve_n);
    const std::string& alloc = s.s("alloc");
    Derivation out;
    out.find = solve_n ? find_size(two_group_unknown(s), s("power"))
                       : find_power("with " + group_sizes_text(s, n, false));
    // the closed form where group 2 follows group 1 (N2 = r N1), beside the exact answer below
    const bool closed = solve_n && (alloc == "equal" || alloc == "ratio");
    const double k = alloc == "ratio" ? s("r") : 1;  // N2 = k N1
    const double za = stats::qnorm(tail_alpha(alpha, a), false), zb = solve_n ? stats::qnorm(s("power")) : NAN;
    std::function<double(double, double)> power;
    if (noether) {
      const double pg = r.num("pgt"), N = n1 + n2, c = n1 / N, ncp = std::sqrt(12 * N * c * (1 - c)) * (pg - 0.5);
      if (closed) {
        const double c0 = 1 / (1 + k);
        out.steps.push_back(
            {"Normal approximation",
             "N ≈ " + z_sum_symbol(a) + "² / (12 c (1 − c) (P(X1 > X2) − 0.5)²), c = " +
                 (alloc == "ratio" ? "1/(1 + r)" : "½") + " the share of group 1, Noether" + FOR_COMPARISON,
             (alloc == "ratio" ? "c = 1 / (1 + " + gv(k) + ") = " + dv(c0) + "; " : std::string()) +
                 z_sum_numbers(alpha, s("power"), a) + "² / (12 × " + dv(c0) + " × " + dv(1 - c0) + " × (" + gv(pg) +
                 " − 0.5)²) = " + dv((za + zb) * (za + zb) / (12 * c0 * (1 - c0) * (pg - 0.5) * (pg - 0.5)))});
      }
      out.steps.push_back({"Share of group 1", "c = N1 / N", N1 + " / " + dvi(N) + " = " + dv(c)});
      out.steps.push_back({"Noncentrality", "λ = √(12 N c (1 − c)) (P(X1 > X2) − 0.5), Noether",
                           "√(12 × " + dvi(N) + " × " + dv(c) + " × " + dv(1 - c) + ") × (" + gv(pg) + " − 0.5) = " +
                               dv(ncp)});
      out.steps.push_back(z_critical(alpha, a));
      out.steps.push_back(z_power_step(ncp, alpha, a));
      power = [=](double x1, double x2) { return mw_noether_power(x1, x2, pg, alpha, a); };
    } else {
      const double are = r.num("are"), d = r.num("delta"), sigma = r.num("sigma");
      const double e1 = n1 * are, e2 = n2 * are, se = sigma * std::sqrt(1 / e1 + 1 / e2), ncp = d / se;
      out.steps.push_back(are_step(s));
      if (closed)
        out.steps.push_back(
            {"Normal approximation",
             std::string(alloc == "equal" ? "N1 = N2" : "N1") + " ≈ " + z_sum_symbol(a) + "² σ²" +
                 (alloc == "ratio" ? "(1 + 1/r)" : " × 2") + " / (δ² ARE)" + FOR_COMPARISON,
             z_sum_numbers(alpha, s("power"), a) + "² × " + gv(sigma) + "² × " +
                 (alloc == "ratio" ? "(1 + 1/" + gv(k) + ")" : std::string("2")) + " / (" + gvs(d) + "² × " +
                 dv(are) + ") = " + dv((za + zb) * (za + zb) * sigma * sigma * (1 + 1 / k) / (d * d * are))});
      out.steps.push_back({"Equivalent t-test sizes", "N1e = N1 × ARE, N2e = N2 × ARE",
                           N1 + " × " + dv(are) + " = " + dv(e1) + ", " + N2 + " × " + dv(are) + " = " + dv(e2)});
      out.steps.push_back({"Standard error", "SE = σ √(1/N1e + 1/N2e)",
                           gv(sigma) + " × √(1/" + dv(e1) + " + 1/" + dv(e2) + ") = " + dv(se)});
      out.steps.push_back({"Noncentrality", "λ = δ / SE", gvs(d) + " / " + dv(se) + " = " + dv(ncp)});
      out.steps.push_back(t_critical(alpha, e1 + e2 - 2, a, "N1e + N2e − 2"));
      out.steps.push_back(t_power_step(ncp, e1 + e2 - 2, alpha, a));
      power = [=](double x1, double x2) { return mw_shift_power(x1, x2, d, sigma, are, alpha, a); };
    }
    if (solve_n) {
      out.steps.push_back(two_group_search(s, power, n, s("power"), 2));
      out.answer =
          group_sizes_text(s, n, true) + " (N = " + dvi(n.total()) + "), with power " + dvp(pw, s("power"));
    } else {
      out.answer = "power = " + dvp(pw);
    }
    return out;
  };
  p.definitions = {{"P(X1>X2)", "Probability that a random observation from group 1 exceeds one from group 2."},
                   {"ARE", "Asymptotic relative efficiency relative to the two-sample t-test."}};
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
// The spread of the location parameters: their mean less the first one (mbar), the sum of squares about the mean
// and Cohen's f² = ss / (k σ²).  From the differences from the first mean: a mean of the means themselves rounds
// when they are large next to their spread, and leaves a spread that is not there.
struct KwEffect {
  double mbar, ss, f2;
};
KwEffect kw_effect(const std::vector<double>& mu, double sigma) {
  const int k = static_cast<int>(mu.size());
  double mbar = 0;
  for (double m : mu) mbar += (m - mu[0]) / k;
  double ss = 0;
  for (double m : mu) ss += (m - mu[0] - mbar) * (m - mu[0] - mbar);
  return {mbar, ss, ss / k / (sigma * sigma)};
}

// The Kruskal–Wallis test's power: the one-way ANOVA F-test's with n × ARE subjects in each of the k groups.
double kw_power(double n, int k, double f2, double are, double alpha) {
  const double ne = n * are, N = ne * k;
  const double df2 = N - k;
  if (!(df2 > 0)) return 0.0;
  return stats::pnf(stats::qf(alpha, k - 1, df2, false), k - 1, df2, N * f2, false);
}

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
    const double sigma = s("sigma"), alpha = s("alpha"), are = are_of(s);
    const double f2 = kw_effect(mu, sigma).f2;
    // with no effect the power is alpha, which no sample size raises
    if (!(f2 > 0) && s.solve == "n") throw CalcError("the group means are all equal");
    if (!std::isfinite(f2))
      throw CalcError("the effect size is too large to compute with (σ is tiny next to the spread of the means)");
    auto P = [&](double n) { return kw_power(n, k, f2, are, alpha); };
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const auto& mu = s.v("means");
    const int k = static_cast<int>(mu.size());
    const double n = r.num("n"), are = r.num("are"), sigma = s("sigma"), alpha = s("alpha");
    const KwEffect e = kw_effect(mu, sigma);
    const double mean = mu[0] + e.mbar, ne = n * are, N = ne * k, df2 = N - k, ncp = N * e.f2, pw = r.num("power");
    const bool solve_n = s.solve == "n";
    const std::string n_text = size(n, solve_n), K = dvi(k);
    std::string squares;
    for (std::size_t i = 0; i < mu.size(); ++i)
      squares += (i ? " + " : "") + std::string("(") + gv(mu[i]) + " − " + dvs(mean) + ")²";
    // no closed form for the noncentral F: the power at the size found and one below shows the answer
    Derivation out;
    out.find = solve_n ? find_size("n per group", s("power")) : find_power("with n = " + n_text + " per group");
    out.steps.push_back({"Effect size", "f² = Σ(μᵢ − μ̄)² / (k σ²), k the number of groups",
                         "k = " + K + ", μ̄ = " + dv(mean) + "; " + squares + " = " + dv(e.ss) + "; " + dv(e.ss) +
                             " / (" + K + " × " + gv(sigma) + "²) = " + dv(e.f2)});
    out.steps.push_back(are_step(s));
    out.steps.push_back({"Equivalent ANOVA size", "Ne = n × ARE × k",
                         n_text + " × " + dv(are) + " × " + K + " = " + dv(N)});
    out.steps.push_back({"Noncentrality", "λ = Ne f²", dv(N) + " × " + dv(e.f2) + " = " + dv(ncp)});
    if (df2 > 0) {
      const double fc = stats::qf(alpha, k - 1, df2, false);
      out.steps.push_back({"Critical value", "F = F(1 − α; k − 1, Ne − k), the F quantile",
                           "F(1 − " + gv(alpha) + "; " + dvi(k - 1) + ", " + dv(df2) + ") = " + dv(fc)});
      out.steps.push_back({"Power", "P(F′ > F), F′ ~ noncentral F(k − 1, Ne − k, λ)",
                           "P(F′ > " + dv(fc) + ") = " + dvp(pw)});
    } else {
      out.steps.push_back({"Power", "0 when Ne ≤ k (no error degrees of freedom)",
                           "Ne − k = " + dv(df2) + ": power = " + dvp(pw)});
    }
    if (solve_n) {
      out.steps.push_back(
          n_search([&](double x) { return kw_power(x, k, e.f2, are, alpha); }, n, s("power"), 2, "n"));
      out.answer = "n = " + dvi(n) + " per group (N = " + dvi(r.num("N")) + "), with power " + dvp(pw, s("power"));
    } else {
      out.answer = "power = " + dvp(pw);
    }
    return out;
  };
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const Alt a = alt(s);
    const int n = static_cast<int>(r.num("n"));
    const double pp = r.num("p"), alpha = s("alpha"), side = tail_alpha(alpha, a);
    // the critical counts of sign_power
    const int cu = a == Alt::Lower ? n + 1 : upper_critical(n, 0.5, side);
    const int cl = a == Alt::Upper ? -1 : lower_critical(n, 0.5, side);
    const auto ge = [&](int c, double p) { return stats::pbinom(c - 1, n, p, false); };  // P(X ≥ c)
    const auto le = [&](int c, double p) { return stats::pbinom(c, n, p, true); };       // P(X ≤ c)
    const bool solve_n = s.solve == "n", shift = s.is("es", "shift");
    const std::string level = a == Alt::Two ? "α/2" : "α", N = size(n, solve_n), P = shift ? dv(pp) : gv(pp);
    const std::string lvl = gv(alpha) + (a == Alt::Two ? "/2" : "");  // the level with α in place
    const std::string at = N + ", ½)";
    Derivation out;
    out.find = solve_n ? find_size("N", s("power")) : find_power("with N = " + N);
    if (shift)
      out.steps.push_back(
          {"Probability of a positive sign", "P = Φ(δ/σ)", "Φ(" + gv(s("es_d")) + ") = " + dv(pp)});
    if (solve_n) {  // the normal approximation to the binomial test, beside the exact answer below
      const double za = stats::qnorm(side, false), zb = stats::qnorm(s("power"));
      out.steps.push_back({"Normal approximation",
                           std::string("N ≈ (z(") + level_symbol(a) + ") × ½ + z(1 − β) √(P(1 − P)))² / (P − ½)²" +
                               FOR_COMPARISON,
                           "(" + dv(za) + " × 0.5 + " + dv(zb) + " × √(" + P + " × (1 − " + P + ")))² / (" + P +
                               " − 0.5)² = " +
                               dv(std::pow(za * 0.5 + zb * std::sqrt(pp * (1 - pp)), 2) / ((pp - 0.5) * (pp - 0.5)))});
    }
    std::string cf, cn;
    if (a != Alt::Lower) {
      cf = "cU = the smallest count with P(X ≥ cU) ≤ " + level;
      cn = cu > n ? "no count: P(X ≥ " + N + " | " + at + " = " + dvp(ge(n, 0.5), side) + " > " + lvl
                  : "cU = " + dvi(cu) + ": P(X ≥ " + dvi(cu) + " | " + at + " = " + dvp(ge(cu, 0.5), side) + " ≤ " +
                        lvl + " < P(X ≥ " + dvi(cu - 1) + " | " + at + " = " + dvp(ge(cu - 1, 0.5), side);
    }
    if (a != Alt::Upper) {
      const std::string sep = a == Alt::Two ? "; " : "";
      cf += sep + "cL = the largest count with P(X ≤ cL) ≤ " + level;
      cn += sep + (cl < 0 ? "no count: P(X ≤ 0 | " + at + " = " + dvp(le(0, 0.5), side) + " > " + lvl
                          : "cL = " + dvi(cl) + ": P(X ≤ " + dvi(cl) + " | " + at + " = " + dvp(le(cl, 0.5), side) +
                                " ≤ " + lvl + " < P(X ≤ " + dvi(cl + 1) + " | " + at + " = " +
                                dvp(le(cl + 1, 0.5), side));
    }
    out.steps.push_back({"Critical counts", cf + ", X ~ Bin(N, ½) the positive signs under H0", cn});
    // the rejection region's probability at p (as two_tail_prob), term by term
    const auto region = [&](double p) {
      std::string f, x, parts;
      if (a != Alt::Lower) {
        f = "P(X ≥ cU)";
        x = cu > n ? "0" : "P(X ≥ " + dvi(cu) + ")";
        parts = cu > n ? "0" : dvp(ge(cu, p));
      }
      if (a != Alt::Upper) {
        const std::string plus = a == Alt::Two ? " + " : "";
        f += plus + "P(X ≤ cL)";
        x += plus + (cl < 0 ? "0" : "P(X ≤ " + dvi(cl) + ")");
        parts += plus + (cl < 0 ? "0" : dvp(le(cl, p)));
      }
      // the terms one by one, when there are two and one of them is not 0 already
      const bool terms = a == Alt::Two && x.find("P(") != std::string::npos;
      return std::array<std::string, 3>{f, x, terms ? " = " + parts : ""};
    };
    const auto h0 = region(0.5), h1 = region(pp);
    out.steps.push_back(
        {"Actual alpha", h0[0] + ", X ~ Bin(N, ½)", h0[1] + h0[2] + " = " + dvp(r.num("actual_alpha"))});
    out.steps.push_back({"Power", h1[0] + ", X ~ Bin(N, P)",
                         h1[1] + " with X ~ Bin(" + N + ", " + P + ")" + h1[2] + " = " + dvp(r.num("power"))});
    if (solve_n) {
      out.steps.push_back(saw_search([&](double x) { return sign_power(static_cast<int>(x), pp, alpha, a); }, n,
                                     s("power"), 1));
      out.answer = "N = " + dvi(n) + ", with power " + dvp(r.num("power"), s("power")) + " (actual α " +
                   dvp(r.num("actual_alpha")) + ")";
    } else {
      out.answer = "power = " + dvp(r.num("power")) + " (actual α " + dvp(r.num("actual_alpha")) + ")";
    }
    return out;
  };
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

// Group 2: n1 × the ratio, rounded up (at least 1: a tiny ratio left it empty, and the interval infinitely wide).
double ci_n2(const Scenario& s, double n) { return solve::round_up(s("r") * n); }

// The standard error of the estimate in units of σ, and its degrees of freedom, for n (group 1's size with two).
double ci_factor(const Scenario& s, double n, bool two_groups) {
  if (!two_groups) return std::sqrt(1 / n);
  return std::sqrt(1 / n + 1 / ci_n2(s, n));
}
double ci_df(const Scenario& s, double n, bool two_groups) { return two_groups ? n + ci_n2(s, n) - 2 : n - 1; }

// Whether a size is sized by the chance that the realized half-width is at most the target (rather than by the
// expected half-width).
bool ci_tolerance(const Scenario& s) { return s.solve == "n" && s.get("tolprob", 0) > 0 && !s.is("known", "z"); }

void run_ci_mean(const Scenario& s, Row& r, RunContext& ctx, bool two_groups) {
  if (!s.is("sides", "two")) check_one_sided_conf(s("conf"));
  const double gamma = s.get("tolprob", 0);
  auto n2_of = [&](double n) { return ci_n2(s, n); };
  auto factor = [&](double n) { return ci_factor(s, n, two_groups); };
  auto dfof = [&](double n) { return ci_df(s, n, two_groups); };
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
  if (ci_tolerance(s)) r.set("prob", c.prob);
  r.set("sigma", s("sigma"));
  r.set("conf", s("conf"));
}

// How a confidence interval's half-width comes about. Solving for N, the smallest size whose half-width (or the
// chance that the realized one) meets the target, shown at the size found and one below, beside the closed form for
// a known σ; solving for the half-width, the half-width at the size given.
derive::Derivation ci_mean_steps(const Scenario& s, const Row& r, bool two_groups) {
  using namespace derive;
  const bool solve_n = s.solve == "n";
  const double n = r.num(two_groups ? "n1" : "n"), sigma = r.num("sigma"), conf = r.num("conf");
  const bool two_sided = s.is("sides", "two"), known = s.is("known", "z"), tol = ci_tolerance(s);
  const double a = two_sided ? (1 - conf) / 2 : 1 - conf, f = ci_factor(s, n, two_groups), df = ci_df(s, n, two_groups);
  const double q = known ? stats::qnorm(a, false) : stats::qt(a, df, false), e = r.num("e");
  const std::string level = two_sided ? "1 − (1 − C)/2" : "C", C = gv(conf);
  const std::string level_numbers = two_sided ? "1 − (1 − " + C + ")/2" : C;
  const std::string n_sym = two_groups ? "N1" : "N", N = size(n, solve_n);
  const std::string N2 = two_groups ? dvi(ci_n2(s, n)) : "", R = two_groups ? gv(s("r")) : "";
  const std::string root = two_groups ? "√(1/" + N + " + 1/" + N2 + ")" : "√(1/" + N + ")";
  const std::string sym = known ? "z" : "t", unknown = two_groups ? "N1, with N2 = ⌈R × N1⌉ (R = " + R + ")," : "N";
  const double target = solve_n ? s("e") : NAN, gamma = s.get("tolprob", 0);
  Derivation out;
  if (!solve_n)
    out.find =
        "the half-width E with " + (two_groups ? "N1 = " + N + " and N2 = ⌈R × N1⌉ (R = " + R + ")" : "N = " + N);
  else if (tol)
    out.find = "the smallest " + unknown + " at which the realized half-width is at most E* = " + gv(target) +
               " with probability at least γ = " + gv(gamma);
  else
    out.find = "the smallest " + unknown + " whose half-width E is at most E* = " + gv(target);
  if (solve_n && !tol) {  // the closed form with σ known, beside the answer below
    const double z = stats::qnorm(a, false), k = two_groups ? 1 + 1 / s("r") : 1;
    out.steps.push_back({"Normal approximation",
                         n_sym + " ≈ z² σ²" + (two_groups ? " (1 + 1/R)" : "") + " / E*², z = Φ⁻¹(" + level +
                             "), the size with σ known" + FOR_COMPARISON,
                         dv(z) + "² × " + gv(sigma) + "²" + (two_groups ? " × (1 + 1/" + R + ")" : "") + " / " +
                             gv(target) + "² = " + dv(z * z * sigma * sigma * k / (target * target))});
  }
  if (two_groups) out.steps.push_back({"Group 2", "N2 = ⌈R × N1⌉", "⌈" + R + " × " + N + "⌉ = " + N2});
  if (known)
    out.steps.push_back({"Critical value", "z = Φ⁻¹(" + level + "), C the confidence level",
                         "Φ⁻¹(" + level_numbers + ") = " + dv(q)});
  else
    out.steps.push_back({"Critical value",
                         "t = t(" + level + "; ν), C the confidence level, the t quantile with ν = " +
                             (two_groups ? "N1 + N2 − 2" : "N − 1") + " degrees of freedom",
                         "t(" + level_numbers + "; " + dvi(df) + ") = " + dv(q)});
  out.steps.push_back({"Standard error", std::string("SE = σ ") + (two_groups ? "√(1/N1 + 1/N2)" : "√(1/N)"),
                       gv(sigma) + " × " + root + " = " + dv(sigma * f)});
  out.steps.push_back({"Half-width", "E = " + sym + " SE", dv(q) + " × " + dv(sigma * f) + " = " + dv(e)});
  if (!solve_n) {
    out.answer = "E = " + dv(e);
    return out;
  }
  if (tol)
    out.steps.push_back({"Chance of the target half-width",
                         "P(E′ ≤ E*) = P(χ²(ν) ≤ ν (E*/E)²), with E′ = " + sym + " S " +
                             (two_groups ? "√(1/N1 + 1/N2)" : "√(1/N)") +
                             " the half-width with the sample SD S and E* the target",
                         "P(χ²(" + dvi(df) + ") ≤ " + dvi(df) + " × (" + gv(target) + " / " + dv(e) + ")²) = " +
                             dvp(r.num("prob"), gamma)});
  const std::string rule = tol ? n_sym + " = the smallest whole number at which P(E′ ≤ E*) reaches γ"
                               : n_sym + " = the smallest whole number whose half-width E is at most the target E*";
  const std::string with = two_groups ? ", with N2 = ⌈R × N1⌉" : "";
  const CiMean at = ci_mean(s, f, df, target);
  if (n - 1 < 2) {
    out.steps.push_back({"Smallest size", rule + with,
                         n_sym + " = " + N + " is the smallest size allowed, and its " +
                             (tol ? "P(E′ ≤ " + gv(target) + ") is already " + dvp(at.prob, gamma) + " ≥ " + gv(gamma)
                                  : "half-width is already " + dvp(at.halfwidth, target) + " ≤ " + gv(target))});
  } else {
    const CiMean below = ci_mean(s, ci_factor(s, n - 1, two_groups), ci_df(s, n - 1, two_groups), target);
    if (tol)
      out.steps.push_back({"Smallest size", rule + with,
                           "P(E′ ≤ " + gv(target) + ") at " + n_sym + " = " + dvi(n - 1) + ": " +
                               dvp(below.prob, gamma) + " < " + gv(gamma) + " ≤ at " + n_sym + " = " + N + ": " +
                               dvp(at.prob, gamma)});
    else
      out.steps.push_back({"Smallest size", rule + with,
                           "half-width at " + n_sym + " = " + dvi(n - 1) + ": " + dvp(below.halfwidth, target) +
                               " > " + gv(target) + " ≥ half-width at " + n_sym + " = " + N + ": " +
                               dvp(at.halfwidth, target)});
  }
  out.answer = (two_groups ? "N1 = " + N + ", N2 = " + N2 + " (N = " + dvi(r.num("n")) + ")" : "N = " + N) +
               ", with half-width " + dv(e) +
               (tol ? " and P(E′ ≤ " + gv(target) + ") = " + dvp(r.num("prob"), gamma) : std::string());
  return out;
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
  p.derivation = [](const Scenario& s, const Row& r) { return ci_mean_steps(s, r, false); };
  p.definitions = {{"E", "Distance from the sample mean to a confidence limit (half the width of a two-sided interval)."},
                   {"P(E ≤ Target)", "Probability that the realized half-width is no larger than the target."}};
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
  p.derivation = [](const Scenario& s, const Row& r) { return ci_mean_steps(s, r, true); };
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
  p.example_solve = "n";
  p.example = {{"e", "0.25 0.5"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(ci_two_means);

}  // namespace
}  // namespace ssc::procs
