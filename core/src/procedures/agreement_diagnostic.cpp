// Agreement (intraclass correlation, Cohen's kappa) and diagnostic accuracy
// (area under the ROC curve, sensitivity and specificity of one test, and two
// tests compared on the same subjects).
#include <unordered_map>

#include "procedures/apart.hpp"
#include "procedures/binomial_tools.hpp"
#include "procedures/common.hpp"
#include "procedures/mcnemar.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// =====================================================================================
// Intraclass correlation (one-way random effects, k ratings per subject)
// =====================================================================================
// With n subjects rated k times, F = MSB/MSW equals c(ρ)·F* where
// c(ρ) = (1 + (k − 1)ρ)/(1 − ρ) and F* ~ F(n − 1, n(k − 1)); the test of
// H0: ρ = ρ0 compares F with c(ρ0) times an F quantile, so power is exact.
// 1 + (k − 1)ρ is taken with one rounding (fma): next to ρ = −1/(k − 1) it is the difference of two numbers near 1,
// and the rounding of (k − 1)ρ alone, up to 1e-16, left it a thousandth off at ρ = −1/(k − 1) + 1e-13.
double icc_c(double rho, double k) { return std::fma(k - 1, rho, 1) / (1 - rho); }

// The power from the ratio c(ρ0)/c(ρ1) of the two multipliers.
double icc_power_ratio(double n, double k, double ratio, double alpha, Alt a) {
  using namespace stats;
  const double df1 = n - 1, df2 = n * (k - 1);
  if (!(df1 > 0 && df2 > 0)) return NAN;
  const double side = tail_alpha(alpha, a);
  double pw = 0;
  if (a != Alt::Lower) pw += pf(qf(side, df1, df2, false) * ratio, df1, df2, false);
  if (a != Alt::Upper) pw += pf(qf(side, df1, df2, true) * ratio, df1, df2, true);
  return stats::clamp01(pw);
}

double icc_power(double n, double k, double rho0, double rho1, double alpha, Alt a) {
  return icc_power_ratio(n, k, icc_c(rho0, k) / icc_c(rho1, k), alpha, a);
}

// The ICC whose multiplier is c = c0·e^u and whose distance from ρ0 is d (accurate as
// k·c0·(e^u − 1)/((c + k − 1)(c0 + k − 1))), kept to a few doubles near both limits, where the power changes much
// from one double to the next: (c − 1)/(c + k − 1), written as 1 − k/(c + k − 1) near 1 and as
// −1/(k − 1) + kc/((k − 1)(c + k − 1)) near −1/(k − 1), or ρ0 + d next to ρ0.
double icc_of(double u, double rho0, double c0, double d, double k) {
  const double c = c0 * std::exp(u);
  if (c > k + 1) return 1 - k / (c + k - 1);  // above 0.5
  if (c < (k - 1) / (2 * k - 1)) return -1 / (k - 1) + k * c / ((k - 1) * (c + k - 1));  // below −1/(2(k − 1))
  if (std::fabs(d) < 0.5 * std::fabs(rho0 + d)) return rho0 + d;
  return (c - 1) / (c + k - 1);
}

Procedure icc() {
  Procedure p;
  p.id = "agree.icc";
  p.name = "Tests for an Intraclass Correlation (Reliability)";
  p.category = "Agreement";
  p.subcategory = "Intraclass Correlation";
  p.order = 1;
  p.summary =
      "Number of subjects for testing H0: ρ = ρ0 for the intraclass correlation of a one-way random-effects model in "
      "which every subject is measured by k raters (or k times). Power is exact, from the F distribution (Walter, "
      "Eliasziw & Donner, 1998).";
  p.keywords = {"intraclass correlation", "icc", "reliability", "raters", "agreement", "repeatability", "test-retest", "inter-rater reliability"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Subjects)"},
                 {"power", "Power", "power", "Power"},
                 {"rho1", "Effect Size (ρ1)", "rho1", "Detectable ρ1"}};
  p.params = {
      alternative("Two-sided (H1: ρ ≠ ρ0)", "One-sided (H1: ρ > ρ0)", "One-sided (H1: ρ < ρ0)", "upper"),
      power_param(),
      alpha_param(),
      n_param("n", "N (Number of Subjects)", "30", 2)
          .sym("N")
          .hint("Number of subjects, each rated by all k raters (or measured k times)."),
      dropout_param(),
      num("k", "k (Raters or Measurements per Subject)", "2")
          .grp(G_N)
          .sym("k")
          .at_least(2)
          .whole()
          .sugg({"2", "3", "2 3 4"})
          .hint("How many raters rate each subject, or how many times each subject is measured."),
      num("rho0", "ρ0 (ICC under H0)", "0.4")
          .grp(G_EFFECT)
          .sym("ρ0")
          .range(-1, 1, true, true)
          .sugg({"0", "0.4", "0.6"})
          .hint("The intraclass correlation (reliability) under the null hypothesis, typically a level that the "
                "true reliability should be shown to exceed."),
      num("rho1", "ρ1 (ICC under H1)", "0.7")
          .grp(G_EFFECT)
          .sym("ρ1")
          .range(-1, 1, true, true)
          .sugg({"0.7", "0.6 0.7 0.8"})
          .hint("The intraclass correlation you expect to be true; values near 1 mean that the ratings of a "
                "subject agree closely.")
          .show(when(not_solving({"rho1"}))),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double k = s("k"), rho0 = s("rho0"), alpha = s("alpha");
    const Alt a = alt(s);
    if (!(rho0 > -1 / (k - 1))) throw CalcError("ρ0 must exceed −1/(k − 1)");
    double n, rho1;
    if (s.solve == "rho1") {
      n = s("n");
      const double target = s("power");
      const bool up = a != Alt::Lower;
      auto P = [&](double x) { return icc_power(n, k, rho0, x, alpha, a); };
      // the search reaches the last ICCs below 1 and above −1/(k − 1) (where 1 + (k − 1)ρ is still positive, and
      // which a given ρ1 may be), so that a ρ0 within 1e-9 of either keeps room for ρ1
      const double hi_lim = std::nextafter(1.0, 0.0);
      double lo_lim = std::nextafter(-1 / (k - 1), 0.0);
      while (!(std::fma(k - 1, lo_lim, 1) > 0)) lo_lim = std::nextafter(lo_lim, 0.0);
      if (up ? rho0 >= hi_lim : rho0 <= lo_lim)
        throw CalcError("ρ1 cannot be told apart from ρ0 = " + text::format_roundtrip(rho0) +
                        ": no number lies between it and " + (up ? "1" : "−1/(k − 1)"));
      // ρ1 is found as u = log(c(ρ1)/c(ρ0)), the log of the factor by which it multiplies F: the power is smooth in
      // u, which runs over all numbers as ρ1 nears either limit, and a tolerance relative to u is relative to the
      // effect near ρ0 and to the distance from the limit near 1 or −1/(k − 1) (where one relative to ρ1 − ρ0 left
      // answers far from the target power: 0.49 for 0.5)
      const double c0 = icc_c(rho0, k);
      auto Pu = [&](double u) { return icc_power_ratio(n, k, std::exp(-u), alpha, a); };
      auto u_of = [&](double x) { return std::log(icc_c(x, k) / c0); };
      const double u = up ? solve::solve_increasing(Pu, target, 0, u_of(rho0 + 0.5 * (1 - rho0)), u_of(hi_lim), 1e-12,
                                                    "ρ1", "ρ0")
                          : solve::solve_decreasing(Pu, target, 0, u_of(rho0 - 0.5 * (rho0 - lo_lim)), u_of(lo_lim),
                                                    1e-12, "ρ1", "ρ0");
      // a double next to the root, which must keep the difference from ρ0 to four digits; near the limits the
      // doubles lie too far apart to hold the root, and the answer is the first of them that reaches the target
      const double c = c0 * std::exp(u), d = k * c0 * std::expm1(u) / ((c + k - 1) * (c0 + k - 1));
      rho1 = std::fmin(hi_lim, std::fmax(lo_lim, icc_of(u, rho0, c0, d, k)));
      rho1 = solve::reach(P, rho1, up ? hi_lim : lo_lim, up, target, rho0, d, "ρ1", "ρ0",
                          "the target power cannot be reached for any admissible value of ρ1");
    } else {
      rho1 = s("rho1");
      if (!(rho1 > -1 / (k - 1))) throw CalcError("ρ1 must exceed −1/(k − 1)");
      // with no effect the power is alpha, which no sample size raises
      if (rho1 == rho0 && s.solve == "n") throw CalcError("ρ1 equals ρ0");
      if ((a == Alt::Upper && rho1 < rho0) || (a == Alt::Lower && rho1 > rho0))
        throw CalcError("ρ1 lies on the side opposite to the one-sided alternative");
      if (s.solve == "n") {
        solve::NSearch o;
        o.lo = 2;
        // Fisher-type normal approximation for a starting value
        const double z = stats::qnorm(tail_alpha(alpha, a), false) + stats::qnorm(s("power"));
        const double d = std::log(icc_c(rho1, k) / icc_c(rho0, k));
        o.guess = 2 * k * z * z / ((k - 1) * d * d) + 2;
        n = solve::min_n([&](double x) { return icc_power(x, k, rho0, rho1, alpha, a); }, s("power"), o, &ctx);
      } else {
        n = s("n");
      }
    }
    r.set("n", n);
    put_power(s, r, icc_power(n, k, rho0, rho1, alpha, a));
    r.set("k", k);
    r.set("total", n * k);
    r.set("rho0", rho0);
    r.set("rho1", rho1);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int", "Number of subjects"), col("k", "k", "int", "Ratings per subject"),
                   col("total", "N × k", "int", "Total number of ratings"), col("rho0", "ρ0", "num"),
                   col("rho1", "ρ1", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    // apart from each other and from the limits 1 and −1/(k − 1)
    const auto [r1, r0] = fnum_apart(r.num("rho1"), r.num("rho0"), {1.0, -1 / (r.num("k") - 1)});
    return "A sample of " + fint(r.num("n")) + " subjects, each rated " + fint(r.num("k")) + " times, " +
           achieves_singular(s, r) + " to detect an intraclass correlation of " + r1 +
           " when the value under the null hypothesis is " + r0 + ", using a " +
           alt_words(alt(s)) + " F-test at significance level α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"N", "Number of subjects."},
                   {"k", "Number of raters (or repeated measurements) per subject."},
                   {"ρ0, ρ1", "Intraclass correlation under the null and the alternative hypotheses."}};
  p.references = {"Walter, S. D., Eliasziw, M. & Donner, A. (1998). Sample size and optimal designs for reliability "
                  "studies. Statistics in Medicine 17, 101–110.",
                  "Shrout, P. E. & Fleiss, J. L. (1979). Intraclass correlations: uses in assessing rater reliability. "
                  "Psychological Bulletin 86, 420–428."};
  p.curve_param = "n";
  p.curve_label = "N (Subjects)";
  p.curve_min = 2;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"k", "2 3"}, {"rho0", "0.4"}, {"rho1", "0.7"}};
  return p;
}
SSC_REGISTER(icc);

// =====================================================================================
// Cohen's kappa: two raters, binary ratings
// =====================================================================================
// Asymptotic variance of the estimated kappa (Fleiss, Cohen & Everitt, 1969) for
// a 2×2 table with equal marginal probabilities π for both raters, as used by
// Cantor (1996): n·Var(κ̂) = Q(κ).
double kappa_q(double kappa, double pi) {
  const double a = pi, b = 1 - pi;
  const double pe = a * a + b * b;
  const double p11 = a * a + kappa * a * b, p22 = b * b + kappa * a * b, p12 = (1 - kappa) * a * b;
  const double t1 = p11 * std::pow(1 - 2 * a * (1 - kappa), 2) + p22 * std::pow(1 - 2 * b * (1 - kappa), 2);
  const double t2 = (1 - kappa) * (1 - kappa) * 2 * p12;  // p.i + pj. = 1 for equal marginals
  const double t3 = std::pow(kappa - pe * (1 - kappa), 2);
  return (t1 + t2 - t3) / ((1 - pe) * (1 - pe));
}

// Smallest and largest kappa for which all four cells are probabilities.
double kappa_min(double pi) { return -std::min(pi, 1 - pi) / std::max(pi, 1 - pi); }

double kappa_power(double n, double k0, double k1, double pi, double alpha, Alt a) {
  using namespace stats;
  const double s0 = std::sqrt(kappa_q(k0, pi) / n), s1 = std::sqrt(kappa_q(k1, pi) / n);
  const double side = tail_alpha(alpha, a);
  const double z = qnorm(side, false), d = k1 - k0;
  double pw = 0;
  if (a != Alt::Lower) pw += pnorm((d - z * s0) / s1);
  if (a != Alt::Upper) pw += pnorm((-d - z * s0) / s1);
  return stats::clamp01(pw);
}

Procedure kappa() {
  Procedure p;
  p.id = "agree.kappa";
  p.name = "Tests for Cohen's Kappa (Two Raters, Binary Ratings)";
  p.category = "Agreement";
  p.subcategory = "Kappa";
  p.order = 1;
  p.summary =
      "Number of subjects for testing H0: κ = κ0 for the agreement of two raters who classify each subject as "
      "positive or negative, given the proportion of positive ratings; uses the large-sample variance of the "
      "estimated kappa under the null and the alternative (Fleiss, Cohen & Everitt, 1969; Cantor, 1996).";
  p.keywords = {"kappa", "cohen", "agreement", "raters", "inter-rater", "reliability", "binary", "cohen's kappa", "inter-observer", "inter-rater agreement"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Subjects)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: κ ≠ κ0)", "One-sided (H1: κ > κ0)", "One-sided (H1: κ < κ0)", "upper"),
      power_param(),
      alpha_param(),
      n_param("n", "N (Number of Subjects)", "100", 2)
          .sym("N")
          .hint("Number of subjects, each classified as positive or negative by both raters."),
      dropout_param(),
      num("kappa0", "κ0 (Kappa under H0)", "0.4")
          .grp(G_EFFECT)
          .sym("κ0")
          .range(-1, 1, false, true)
          .sugg({"0", "0.4", "0.6"})
          .hint("The kappa under the null hypothesis, typically a level of agreement that the true kappa should be "
                "shown to exceed (0 is agreement by chance alone)."),
      num("kappa1", "κ1 (Kappa under H1)", "0.7")
          .grp(G_EFFECT)
          .sym("κ1")
          .range(-1, 1, false, true)
          .sugg({"0.7", "0.6 0.7 0.8"})
          .hint("The kappa you expect to be true; 1 is perfect agreement and 0 is agreement no better than "
                "chance."),
      num("prev", "π (Proportion of Positive Ratings)", "0.3")
          .grp(G_EFFECT)
          .sym("π")
          .range(0, 1, true, true)
          .sugg({"0.3", "0.5", "0.1 0.3 0.5"})
          .hint("Marginal probability that a rater classifies a subject as positive (the same for both raters)."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double k0 = s("kappa0"), k1 = s("kappa1"), pi = s("prev"), alpha = s("alpha");
    const Alt a = alt(s);
    const double kmin = kappa_min(pi);
    if (!(k0 >= kmin && k1 >= kmin))
      throw CalcError("with π = " + fnum(pi) + ", kappa cannot be below " + fnum(kmin, 4));
    // with no effect the power is alpha, which no sample size raises
    if (k1 == k0 && s.solve == "n") throw CalcError("κ1 equals κ0");
    if ((a == Alt::Upper && k1 < k0) || (a == Alt::Lower && k1 > k0))
      throw CalcError("κ1 lies on the side opposite to the one-sided alternative");
    double n;
    if (s.solve == "n") {
      const double z0 = stats::qnorm(tail_alpha(alpha, a), false), z1 = stats::qnorm(s("power"));
      const double guess = std::pow(z0 * std::sqrt(kappa_q(k0, pi)) + z1 * std::sqrt(kappa_q(k1, pi)), 2) /
                           ((k1 - k0) * (k1 - k0));
      solve::NSearch o;
      o.lo = 2;
      o.guess = guess;
      n = solve::min_n([&](double x) { return kappa_power(x, k0, k1, pi, alpha, a); }, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    r.set("n", n);
    put_power(s, r, kappa_power(n, k0, k1, pi, alpha, a));
    r.set("kappa0", k0);
    r.set("kappa1", k1);
    r.set("prev", pi);
    r.set("po1", 1 - 2 * (1 - k1) * pi * (1 - pi));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int", "Number of subjects"), col("kappa0", "κ0", "num"), col("kappa1", "κ1", "num"),
                   col("prev", "π", "num", "Proportion of positive ratings"),
                   col("po1", "Agreement under H1", "num", "Proportion of subjects both raters classify alike"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const auto [k1, k0] = fnum_apart(r.num("kappa1"), r.num("kappa0"));
    return "A sample of " + fint(r.num("n")) + " subjects rated by two raters " + achieves_singular(s, r) +
           " to detect a kappa of " + k1 + " when the value under the null hypothesis is " + k0 +
           ", assuming a proportion of positive ratings of " + fnum_fine(r.num("prev")) + ", using a " +
           alt_words(alt(s)) + " z-test at significance level α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"κ0, κ1", "Cohen's kappa under the null and the alternative hypotheses."},
                   {"π", "Probability that a rater classifies a subject as positive (equal for both raters)."},
                   {"Agreement under H1", "Probability that the two raters agree when κ = κ1."}};
  p.references = {"Fleiss, J. L., Cohen, J. & Everitt, B. S. (1969). Large sample standard errors of kappa and "
                  "weighted kappa. Psychological Bulletin 72, 323–327.",
                  "Cantor, A. B. (1996). Sample-size calculations for Cohen's kappa. Psychological Methods 1, 150–153."};
  p.curve_param = "n";
  p.curve_label = "N (Subjects)";
  p.curve_min = 2;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"kappa0", "0.4"}, {"kappa1", "0.6 0.7"}, {"prev", "0.3"}};
  return p;
}
SSC_REGISTER(kappa);

// =====================================================================================
// Area under the ROC curve
// =====================================================================================
// Hanley & McNeil (1982): Var(AUC) = [A(1 − A) + (n1 − 1)(Q1 − A²) + (n0 − 1)(Q2 − A²)]/(n1 n0)
// with Q1 = A/(2 − A) and Q2 = 2A²/(1 + A); n1 cases (diseased), n0 controls.
double auc_var(double A, double n1, double n0) {
  const double q1 = A / (2 - A), q2 = 2 * A * A / (1 + A);
  return (A * (1 - A) + (n1 - 1) * (q1 - A * A) + (n0 - 1) * (q2 - A * A)) / (n1 * n0);
}

double auc_power(double n1, double n0, double a0, double a1, double alpha, Alt a) {
  using namespace stats;
  const double s0 = std::sqrt(auc_var(a0, n1, n0)), s1 = std::sqrt(auc_var(a1, n1, n0));
  const double side = tail_alpha(alpha, a);
  const double z = qnorm(side, false), d = a1 - a0;
  double pw = 0;
  if (a != Alt::Lower) pw += pnorm((d - z * s0) / s1);
  if (a != Alt::Upper) pw += pnorm((-d - z * s0) / s1);
  return stats::clamp01(pw);
}

Procedure auc() {
  Procedure p;
  p.id = "diag.auc";
  p.name = "Tests for the Area Under a ROC Curve";
  p.category = "Diagnostic Accuracy";
  p.subcategory = "ROC Curves";
  p.order = 1;
  p.summary =
      "Numbers of cases (with the condition) and controls for testing whether the area under the ROC curve of a "
      "diagnostic test differs from a reference value (0.5 for an uninformative test), using the variance of the "
      "estimated AUC of Hanley & McNeil (1982).";
  p.keywords = {"roc", "auc", "area under the curve", "diagnostic", "c-statistic", "discrimination", "biomarker", "receiver operating characteristic", "hanley-mcneil"};
  p.solve_for = {{"n", "Sample Size", "n1", "Cases"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: AUC ≠ AUC0)", "One-sided (H1: AUC > AUC0)", "One-sided (H1: AUC < AUC0)", "upper"),
      power_param(),
      alpha_param(),
      n_param("n1", "Number of Cases (With the Condition)", "50", 2)
          .sym("Cases")
          .hint("Number of subjects with the condition (cases); the controls number this times the controls per "
                "case, rounded up."),
      num("ratio", "Controls per Case", "1")
          .grp(G_N)
          .sym("Controls/Case")
          .above(0)
          .sugg({"1", "2", "1 2 3"})
          .hint("Number of subjects without the condition (controls) per case; 1 gives equal numbers."),
      num("auc0", "AUC0 (Area under H0)", "0.5")
          .grp(G_EFFECT)
          .sym("AUC0")
          .range(0, 1, true, true)
          .sugg({"0.5", "0.7"})
          .hint("The area under the ROC curve under the null hypothesis; 0.5 is a test no better than chance."),
      num("auc1", "AUC1 (Area under H1)", "0.7")
          .grp(G_EFFECT)
          .sym("AUC1")
          .range(0, 1, true, true)
          .sugg({"0.7", "0.65 0.7 0.75 0.8"})
          .hint("The area under the ROC curve you expect: the chance that a random case scores higher on the test "
                "than a random control."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double a0 = s("auc0"), a1 = s("auc1"), ratio = s("ratio"), alpha = s("alpha");
    const Alt a = alt(s);
    // with no effect the power is alpha, which no sample size raises
    if (a1 == a0 && s.solve == "n") throw CalcError("AUC1 equals AUC0");
    if ((a == Alt::Upper && a1 < a0) || (a == Alt::Lower && a1 > a0))
      throw CalcError("AUC1 lies on the side opposite to the one-sided alternative");
    auto controls = [&](double n1) { return solve::round_up(n1 * ratio); };
    double n1;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = solve::first_whole([&](double x) { return controls(x) >= 2; }, 2, o.hi,
                                "the ratio of controls to cases is too small: there would be fewer than 2 controls");
      n1 = solve::min_n([&](double x) { return auc_power(x, controls(x), a0, a1, alpha, a); }, s("power"), o, &ctx);
    } else {
      n1 = s("n1");
    }
    const double n0 = controls(n1);
    if (n0 < 2)
      throw CalcError("there would be " + fint(n0) + (n0 == 1 ? " control" : " controls") +
                      " (cases × the ratio, rounded up); at least 2 are needed");
    r.set("n1", n1);
    r.set("n0", n0);
    r.set("n", n1 + n0);
    put_power(s, r, auc_power(n1, n0, a0, a1, alpha, a));
    r.set("auc0", a0);
    r.set("auc1", a1);
    r.set("se1", std::sqrt(auc_var(a1, n1, n0)));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n1", "Cases", "int", "Subjects with the condition"),
                   col("n0", "Controls", "int", "Subjects without the condition"), col("n", "N", "int", "Total"),
                   col("auc0", "AUC0", "num"), col("auc1", "AUC1", "num"),
                   col("se1", "SE(AUC) under H1", "num", "Standard error of the estimated AUC when AUC = AUC1"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const auto [a1, a0] = fnum_apart(r.num("auc1"), r.num("auc0"));
    return "Samples of " + fint(r.num("n1")) + " cases and " + fint(r.num("n0")) + " controls " +
           achieves(s, r) + " to detect an area under the ROC curve of " + a1 +
           " when the value under the null hypothesis is " + a0 + ", using a " + alt_words(alt(s)) +
           " z-test at significance level α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"Cases, Controls", "Subjects with and without the condition (by the reference standard)."},
                   {"AUC", "Area under the ROC curve: the probability that a random case scores higher than a random "
                           "control."}};
  p.references = {"Hanley, J. A. & McNeil, B. J. (1982). The meaning and use of the area under a receiver operating "
                  "characteristic (ROC) curve. Radiology 143, 29–36.",
                  "Obuchowski, N. A. (1998). Sample size calculations in studies of test accuracy. Statistical Methods "
                  "in Medical Research 7, 371–392."};
  p.notes = {"The variance uses Hanley and McNeil's exponential approximation; for binormal data it is close but "
             "slightly conservative."};
  p.curve_param = "n1";
  p.curve_label = "Number of Cases";
  p.curve_min = 2;
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"auc0", "0.5"}, {"auc1", "0.65 0.7 0.75"}, {"ratio", "1"}};
  return p;
}
SSC_REGISTER(auc);

// =====================================================================================
// Sensitivity and specificity (exact binomial tests)
// =====================================================================================
// Smallest c with P(X ≥ c | n, p0) ≤ α: the exact one-sided test of H0: p ≤ p0.
double exact_upper_power(int n, double p0, double p1, double alpha) {
  if (n < 1) return 0;
  const int c = upper_critical(n, p0, alpha);
  return c > n ? 0.0 : stats::pbinom(c - 1, n, p1, false);
}

Procedure sens_spec() {
  Procedure p;
  p.id = "diag.sens";
  p.name = "Tests for Sensitivity and Specificity (One Diagnostic Test)";
  p.category = "Diagnostic Accuracy";
  p.subcategory = "Sensitivity and Specificity";
  p.order = 1;
  p.summary =
      "Total number of subjects for showing that a diagnostic test's sensitivity exceeds Se0 and/or its specificity "
      "exceeds Sp0, with exact one-sided binomial tests among the subjects with and without the condition, given the "
      "prevalence of the condition in the sampled population.";
  p.keywords = {"sensitivity", "specificity", "diagnostic accuracy", "screening", "prevalence", "binomial", "true positive rate", "true negative rate"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      select("goal", "Hypotheses",
             {{"both", "Sensitivity and specificity (both tests significant)", {}},
              {"se", "Sensitivity only", {}},
              {"sp", "Specificity only", {}}},
             "both")
          .grp(G_TEST)
          .hint("What the study must show: both a high enough sensitivity and specificity (both tests must "
                "succeed), or only one of them."),
      power_param(),
      alpha_param("0.05", "Alpha (One-Sided, Each Test)").range(0, 0.5, true, true),
      n_param("n", "N (Total Number of Subjects)", "200", 2)
          .sym("N")
          .hint("Total number of subjects tested, with and without the condition; the prevalence splits them."),
      num("prev", "Prevalence of the Condition", "0.2")
          .grp(G_N)
          .sym("Prevalence")
          .range(0, 1, true, true)
          .sugg({"0.1", "0.2", "0.3"})
          .hint("The share of the sampled subjects who have the condition, between 0 and 1; it sets how many "
                "subjects with and without it to expect."),
      num("se0", "Se0 (Sensitivity under H0)", "0.7")
          .grp(G_EFFECT)
          .sym("Se0")
          .range(0, 1, true, true)
          .hint("The sensitivity (share of subjects with the condition who test positive) that the test must be "
                "shown to exceed.")
          .show(when(is_not("goal", {"sp"}))),
      num("se1", "Se1 (Sensitivity under H1)", "0.85")
          .grp(G_EFFECT)
          .sym("Se1")
          .range(0, 1, true, true)
          .hint("The sensitivity you expect the test to have; it must exceed Se0.")
          .show(when(is_not("goal", {"sp"}))),
      num("sp0", "Sp0 (Specificity under H0)", "0.8")
          .grp(G_EFFECT)
          .sym("Sp0")
          .range(0, 1, true, true)
          .hint("The specificity (share of subjects without the condition who test negative) that the test must "
                "be shown to exceed.")
          .show(when(is_not("goal", {"se"}))),
      num("sp1", "Sp1 (Specificity under H1)", "0.9")
          .grp(G_EFFECT)
          .sym("Sp1")
          .range(0, 1, true, true)
          .hint("The specificity you expect the test to have; it must exceed Sp0.")
          .show(when(is_not("goal", {"se"}))),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double prev = s("prev"), alpha = s("alpha");
    const bool use_se = !s.is("goal", "sp"), use_sp = !s.is("goal", "se");
    const double se0 = use_se ? s("se0") : NAN, se1 = use_se ? s("se1") : NAN;
    const double sp0 = use_sp ? s("sp0") : NAN, sp1 = use_sp ? s("sp1") : NAN;
    // the tests are one-sided; at Se1 = Se0 (or Sp1 = Sp0) the power is the exact test's size, which no sample
    // size raises above alpha
    const bool n_mode = s.solve == "n";
    if (use_se && (n_mode ? !(se1 > se0) : se1 < se0))
      throw CalcError(n_mode ? "Se1 must exceed Se0" : "Se1 cannot be smaller than Se0");
    if (use_sp && (n_mode ? !(sp1 > sp0) : sp1 < sp0))
      throw CalcError(n_mode ? "Sp1 must exceed Sp0" : "Sp1 cannot be smaller than Sp0");
    // subjects with / without the condition: whole numbers not exceeding the expected counts
    auto nd = [&](double n) { return std::floor(n * prev + 1e-9); };
    auto nn = [&](double n) { return std::floor(n * (1 - prev) + 1e-9); };
    auto power_parts = [&](double n, double& pse, double& psp) {
      pse = use_se ? exact_upper_power(static_cast<int>(nd(n)), se0, se1, alpha) : 1.0;
      psp = use_sp ? exact_upper_power(static_cast<int>(nn(n)), sp0, sp1, alpha) : 1.0;
      return pse * psp;
    };
    double n;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = 2;
      o.monotone = false;  // exact binomial power is saw-toothed
      // a tooth is about 1/min(p0, 1 − p0) subjects with (or without) the condition; the first one that reaches the
      // target can end two teeth before the next begins (power 0.83 at N = 161, then 0.80 only from N = 242)
      if (use_se) o.period = 1 / (prev * std::fmin(se0, 1 - se0));
      if (use_sp) o.period = std::fmax(o.period, 1 / ((1 - prev) * std::fmin(sp0, 1 - sp0)));
      const double z = stats::qnorm(alpha, false) + stats::qnorm(s("power"));
      double guess = 2;
      if (use_se) guess = std::max(guess, z * z * se1 * (1 - se1) / ((se1 - se0) * (se1 - se0)) / prev);
      if (use_sp) guess = std::max(guess, z * z * sp1 * (1 - sp1) / ((sp1 - sp0) * (sp1 - sp0)) / (1 - prev));
      o.guess = guess;
      double a, b;
      n = solve::min_n([&](double x) { return power_parts(x, a, b); }, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    double pse, psp;
    const double pw = power_parts(n, pse, psp);
    r.set("n", n);
    r.set("nd", nd(n));
    r.set("nn", nn(n));
    put_power(s, r, pw);
    if (use_se) {
      r.set("pse", pse);
      r.set("se0", se0);
      r.set("se1", se1);
    }
    if (use_sp) {
      r.set("psp", psp);
      r.set("sp0", sp0);
      r.set("sp1", sp1);
    }
    r.set("prev", prev);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    c.push_back(col("n", "N", "int", "Total number of subjects"));
    c.push_back(col("nd", "With Condition", "int", "Subjects with the condition (N × prevalence, rounded down)"));
    c.push_back(col("nn", "Without", "int", "Subjects without the condition (N × (1 − prevalence), rounded down)"));
    if (!s.is("goal", "sp")) {
      c.push_back(col("pse", "Power (Se)", "prob", "Power of the sensitivity test"));
      c.push_back(col("se0", "Se0", "num"));
      c.push_back(col("se1", "Se1", "num"));
    }
    if (!s.is("goal", "se")) {
      c.push_back(col("psp", "Power (Sp)", "prob", "Power of the specificity test"));
      c.push_back(col("sp0", "Sp0", "num"));
      c.push_back(col("sp1", "Sp1", "num"));
    }
    c.push_back(col("prev", "Prevalence", "num"));
    c.push_back(col("alpha", "Alpha", "alpha"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    std::string what;
    if (!s.is("goal", "sp")) {
      const auto [e1, e0] = fnum_apart(r.num("se1"), r.num("se0"));
      what = "a sensitivity of " + e1 + " (null value " + e0 + ")";
    }
    if (s.is("goal", "both")) what += " and ";
    if (!s.is("goal", "se")) {
      const auto [p1, p0] = fnum_apart(r.num("sp1"), r.num("sp0"));
      what += "a specificity of " + p1 + " (null value " + p0 + ")";
    }
    // both counts are rounded down, so they can add up to one less than N
    return "A total of " + fint(r.num("n")) + " subjects (" + fint(r.num("nd")) + " with and " + fint(r.num("nn")) +
           " without the condition: the expected numbers at a prevalence of " + fnum_fine(r.num("prev")) +
           ", rounded down) " + achieves_singular(s, r) +
           " to show " + what + " with exact one-sided binomial tests at α = " + fnum(r.num("alpha")) +
           (s.is("goal", "both") ? " each (power: both tests significant)." : ".");
  };
  p.definitions = {{"Se, Sp", "Sensitivity (true positive rate) and specificity (true negative rate)."},
                   {"Power", "Probability that the test(s) reject the null hypotheses when Se = Se1 and Sp = Sp1; for "
                             "both, the product of the two powers (independent subjects)."}};
  p.references = {"Flahault, A., Cadilhac, M. & Thomas, G. (2005). Sample size calculation should be performed for "
                  "design accuracy in diagnostic test studies. Journal of Clinical Epidemiology 58, 859–862.",
                  "Pepe, M. S. (2003). The Statistical Evaluation of Medical Tests for Classification and Prediction. "
                  "Oxford University Press."};
  p.notes = {"Exact binomial power is not monotone in N; the smallest N reaching the target is reported."};
  p.curve_param = "n";
  p.curve_label = "N (Total)";
  p.curve_min = 2;
  p.curve_integer = true;
  p.example_solve = "n";
  p.example = {{"goal", "both"}, {"power", "0.8"}, {"prev", "0.2 0.3"}, {"se0", "0.7"}, {"se1", "0.85"},
               {"sp0", "0.8"}, {"sp1", "0.9"}};
  return p;
}
SSC_REGISTER(sens_spec);

// =====================================================================================
// Two diagnostic tests on the same subjects (paired design)
// =====================================================================================
// Among the subjects with the condition each one's two results form a pair, and the sensitivities differ when the
// pairs positive on test 1 only (P10) and on test 2 only (P01) differ in frequency: McNemar's test.  Among those
// without the condition the same holds for the specificities, with "correct" in place of "positive".
struct PairedCells {
  double p10 = 0, p01 = 0;  // test 1 correct only, test 2 correct only
};

// The discordant cells from the two tests' accuracies a and b (sensitivities or specificities) and their dependence:
// conditionally independent, a correlation between the results, or the share of discordant pairs.
PairedCells paired_cells(const Scenario& s, double a, double b, const char* rho_key, const char* disc_key,
                         const char* what) {
  double both;  // the share correct on both tests
  if (s.is("dep", "corr"))
    both = a * b + s(rho_key) * std::sqrt(a * (1 - a) * b * (1 - b));
  else if (s.is("dep", "disc"))
    both = (a + b - s(disc_key)) / 2;
  else
    both = a * b;
  const PairedCells c{a - both, b - both};
  // all four cells must be probabilities: both correct, test 1 only, test 2 only, neither
  if (!(c.p10 >= -1e-12 && c.p01 >= -1e-12 && both >= -1e-12 && 1 - both - c.p10 - c.p01 >= -1e-12))
    throw CalcError(std::string("the ") + (s.is("dep", "disc") ? "share of discordant results" : "correlation") +
                    " is impossible for these " + what + ": the two tests' results would not form valid shares");
  return {std::max(0.0, c.p10), std::max(0.0, c.p01)};
}

Procedure paired_diagnostic() {
  Procedure p;
  p.id = "diag.paired";
  p.name = "Comparing Two Diagnostic Tests on the Same Subjects (Paired Sensitivity and Specificity)";
  p.category = "Diagnostic Accuracy";
  p.subcategory = "Sensitivity and Specificity";
  p.order = 2;
  p.summary =
      "Total number of subjects for showing that two diagnostic tests, both applied to every subject, differ in "
      "sensitivity and/or specificity, with McNemar's test among the subjects with and without the condition "
      "(exact conditional or Connor's normal approximation). The dependence between the tests is given as "
      "conditional independence, a correlation or the share of discordant results.";
  p.keywords = {"compare two diagnostic tests", "paired diagnostic", "comparative accuracy", "sensitivity",
                "specificity", "mcnemar", "screening", "paired design", "same patients", "new test versus standard"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      select("goal", "Hypotheses",
             {{"se", "Sensitivities differ", {}},
              {"sp", "Specificities differ", {}},
              {"both", "Both differ (both tests significant)", {}}},
             "se")
          .grp(G_TEST)
          .hint("What the study must show: a difference in sensitivity, in specificity, or in both (then both "
                "McNemar tests must succeed)."),
      alternative("Two-sided (H1: the tests differ)", "One-sided (H1: test 1 is more accurate)",
                  "One-sided (H1: test 1 is less accurate)"),
      select("method", "Test / Power Calculation",
             {{"exact", "Exact conditional test, exact power", {}}, {"normal", "Normal approximation (Connor)", {}}},
             "exact")
          .grp(G_TEST)
          .hint("The exact conditional (binomial) McNemar test with its exact power, or Connor's normal "
                "approximation."),
      power_param(),
      alpha_param("0.05", "Alpha (Each Test)"),
      n_param("n", "N (Total Number of Subjects)", "500", 2)
          .sym("N")
          .hint("Total number of subjects, each given both tests; the prevalence splits them into those with and "
                "without the condition."),
      num("prev", "Prevalence of the Condition", "0.2")
          .grp(G_N)
          .sym("Prevalence")
          .range(0, 1, true, true)
          .sugg({"0.1", "0.2", "0.3"})
          .hint("The share of the sampled subjects who have the condition, between 0 and 1."),
      num("se1", "Se1 (Sensitivity of Test 1)", "0.9")
          .grp(G_EFFECT)
          .sym("Se1")
          .range(0, 1, true, true)
          .sugg({"0.9", "0.85 0.9"})
          .hint("The share of subjects with the condition whom test 1 (e.g. the new test) finds positive.")
          .show(when(is_not("goal", {"sp"}))),
      num("se2", "Se2 (Sensitivity of Test 2)", "0.8")
          .grp(G_EFFECT)
          .sym("Se2")
          .range(0, 1, true, true)
          .sugg({"0.8", "0.75 0.8"})
          .hint("The share of subjects with the condition whom test 2 (e.g. the standard test) finds positive.")
          .show(when(is_not("goal", {"sp"}))),
      num("sp1", "Sp1 (Specificity of Test 1)", "0.85")
          .grp(G_EFFECT)
          .sym("Sp1")
          .range(0, 1, true, true)
          .sugg({"0.85", "0.8 0.85"})
          .hint("The share of subjects without the condition whom test 1 finds negative.")
          .show(when(is_not("goal", {"se"}))),
      num("sp2", "Sp2 (Specificity of Test 2)", "0.9")
          .grp(G_EFFECT)
          .sym("Sp2")
          .range(0, 1, true, true)
          .sugg({"0.9", "0.9 0.95"})
          .hint("The share of subjects without the condition whom test 2 finds negative.")
          .show(when(is_not("goal", {"se"}))),
      select("dep", "Dependence between the Tests",
             {{"indep", "Conditionally independent", {}},
              {"corr", "Correlation between the results", {}},
              {"disc", "Share of discordant results", {}}},
             "indep")
          .grp(G_EFFECT)
          .hint("How the two tests' results go together in the same subject. Tests that look at the same thing tend "
                "to agree (positive correlation), which leaves fewer discordant results and needs fewer subjects."),
      num("rse", "ρ (Correlation among Subjects with the Condition)", "0.3")
          .grp(G_EFFECT)
          .sym("ρ (Se)")
          .range(-1, 1, true, true)
          .sugg({"0.3", "0 0.3 0.5"})
          .hint("Correlation between the two tests' results (positive or not) among subjects with the condition.")
          .show(when({is("dep", {"corr"}), is_not("goal", {"sp"})})),
      num("rsp", "ρ (Correlation among Subjects without the Condition)", "0.3")
          .grp(G_EFFECT)
          .sym("ρ (Sp)")
          .range(-1, 1, true, true)
          .sugg({"0.3", "0 0.3 0.5"})
          .hint("Correlation between the two tests' results among subjects without the condition.")
          .show(when({is("dep", {"corr"}), is_not("goal", {"se"})})),
      num("dse", "ψ (Discordant Share among Subjects with the Condition)", "0.15")
          .grp(G_EFFECT)
          .sym("ψ (Se)")
          .range(0, 1, true, true)
          .sugg({"0.15", "0.1 0.15 0.2"})
          .hint("The share of subjects with the condition whose two results differ (positive on one test only); "
                "at least |Se1 − Se2|.")
          .show(when({is("dep", {"disc"}), is_not("goal", {"sp"})})),
      num("dsp", "ψ (Discordant Share among Subjects without the Condition)", "0.15")
          .grp(G_EFFECT)
          .sym("ψ (Sp)")
          .range(0, 1, true, true)
          .sugg({"0.15", "0.1 0.15 0.2"})
          .hint("The share of subjects without the condition whose two results differ; at least |Sp1 − Sp2|.")
          .show(when({is("dep", {"disc"}), is_not("goal", {"se"})})),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double prev = s("prev"), alpha = s("alpha");
    const bool use_se = !s.is("goal", "sp"), use_sp = !s.is("goal", "se");
    const Alt a = alt(s);
    PairedCells se, sp;
    auto check = [&](double x1, double x2, const char* name) {
      // with no difference the power is the test's size, which no sample size raises above alpha
      if (x1 == x2 && s.solve == "n")
        throw CalcError(std::string("the two ") + name + " are equal; there is no difference to detect");
      if ((a == Alt::Upper && x1 < x2) || (a == Alt::Lower && x1 > x2))
        throw CalcError(std::string("the difference in ") + name + " is on the side opposite to the one-sided "
                        "alternative");
    };
    if (use_se) {
      check(s("se1"), s("se2"), "sensitivities");
      se = paired_cells(s, s("se1"), s("se2"), "rse", "dse", "sensitivities");
    }
    if (use_sp) {
      check(s("sp1"), s("sp2"), "specificities");
      sp = paired_cells(s, s("sp1"), s("sp2"), "rsp", "dsp", "specificities");
    }
    // subjects with / without the condition: whole numbers not exceeding the expected counts
    auto nd = [&](double n) { return std::floor(n * prev + 1e-9); };
    auto nn = [&](double n) { return std::floor(n * (1 - prev) + 1e-9); };
    McNemarCache cse, csp;  // the exact test's rejection probabilities, kept between sample sizes
    // and each group's power by its number of pairs: at a prevalence of 0.001 a thousand totals in a row have the
    // same number of subjects with the condition, and the search's scan of totals would compute it again for each
    std::unordered_map<int, double> mse, msp;
    auto part = [&](std::unordered_map<int, double>& memo, double pairs, const PairedCells& c, McNemarCache& cache) {
      if (pairs < 1) return 0.0;
      const int k = static_cast<int>(pairs);
      const auto it = memo.find(k);
      return it != memo.end() ? it->second : memo.emplace(k, mcnemar_power(s, k, c.p10, c.p01, cache, ctx)).first->second;
    };
    auto power_parts = [&](double n, double& pse, double& psp) {
      pse = !use_se ? 1.0 : part(mse, nd(n), se, cse);
      psp = !use_sp ? 1.0 : part(msp, nn(n), sp, csp);
      return pse * psp;
    };
    // the exact test is computed for up to a million pairs in each group
    const bool normal = s.is("method", "normal");
    // (and no more than the 1e9 subjects any count may have: a prevalence of 1e-303 would give infinity)
    const double n_max =
        normal ? 1e8 : std::fmin(1e9, std::floor(1e6 / std::max(use_se ? prev : 0.0, use_sp ? 1 - prev : 0.0)));
    double n;
    if (s.solve == "n") {
      const double za = stats::qnorm(tail_alpha(alpha, a), false), zb = stats::qnorm(s("power"));
      auto connor = [&](const PairedCells& c) {  // pairs needed by Connor's formula
        const double psi = c.p10 + c.p01, d = c.p10 - c.p01;
        return std::pow(za * std::sqrt(psi) + zb * std::sqrt(psi - d * d), 2) / (d * d);
      };
      solve::NSearch o;
      o.lo = 2;
      o.hi = n_max;
      o.monotone = normal;
      double guess = 2;
      if (use_se) guess = std::max(guess, connor(se) / prev);
      if (use_sp) guess = std::max(guess, connor(sp) / (1 - prev));
      o.guess = std::ceil(guess);
      double x, y;
      n = solve::min_n([&](double m) { return power_parts(m, x, y); }, s("power"), o, &ctx);
    } else {
      n = s("n");
      if (n > n_max)
        throw CalcError("the exact test is computed for up to 1,000,000 subjects with (or without) the condition; use "
                        "the normal approximation for more");
    }
    double pse, psp;
    const double pw = power_parts(n, pse, psp);
    r.set("n", n);
    r.set("nd", nd(n));
    r.set("nn", nn(n));
    put_power(s, r, pw);
    if (use_se) {
      r.set("pse", pse);
      r.set("se1", s("se1"));
      r.set("se2", s("se2"));
      r.set("psise", se.p10 + se.p01);
    }
    if (use_sp) {
      r.set("psp", psp);
      r.set("sp1", s("sp1"));
      r.set("sp2", s("sp2"));
      r.set("psisp", sp.p10 + sp.p01);
    }
    r.set("prev", prev);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    c.push_back(col("n", "N", "int", "Total number of subjects"));
    c.push_back(col("nd", "With Condition", "int", "Subjects with the condition (N × prevalence, rounded down)"));
    c.push_back(col("nn", "Without", "int", "Subjects without the condition (N × (1 − prevalence), rounded down)"));
    if (!s.is("goal", "sp")) {
      c.push_back(col("pse", "Power (Se)", "prob", "Power of the McNemar test of the sensitivities"));
      c.push_back(col("se1", "Se1", "num"));
      c.push_back(col("se2", "Se2", "num"));
      c.push_back(col("psise", "ψ (Se)", "num", "Share of subjects with the condition with discordant results"));
    }
    if (!s.is("goal", "se")) {
      c.push_back(col("psp", "Power (Sp)", "prob", "Power of the McNemar test of the specificities"));
      c.push_back(col("sp1", "Sp1", "num"));
      c.push_back(col("sp2", "Sp2", "num"));
      c.push_back(col("psisp", "ψ (Sp)", "num", "Share of subjects without the condition with discordant results"));
    }
    c.push_back(col("prev", "Prevalence", "num"));
    c.push_back(col("alpha", "Alpha", "alpha"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    // the two tests' values apart from each other (0.90001 and 0.9, not 0.9 and 0.9), a prevalence never as 0 or 1
    std::string what;
    if (!s.is("goal", "sp")) {
      const auto [se1, se2] = fnum_apart(r.num("se1"), r.num("se2"));
      what = "the sensitivities " + se1 + " and " + se2;
    }
    if (s.is("goal", "both")) what += " and ";
    if (!s.is("goal", "se")) {
      const auto [sp1, sp2] = fnum_apart(r.num("sp1"), r.num("sp2"));
      what += "the specificities " + sp1 + " and " + sp2;
    }
    return "A total of " + fint(r.num("n")) + " subjects (" + fint(r.num("nd")) + " with and " + fint(r.num("nn")) +
           " without the condition: the expected numbers at a prevalence of " + fnum_fine(r.num("prev")) +
           ", rounded down), each given both tests, " + achieves_singular(s, r) + " to detect the difference between " +
           what + " of the two tests using " + alt_words(alt(s)) + " McNemar " +
           (s.is("goal", "both") ? "tests" : "test") +
           (s.is("method", "exact") ? " (exact conditional)" : " (normal approximation)") + " at α = " +
           fnum(r.num("alpha")) + (s.is("goal", "both") ? " each (power: both tests significant)." : ".");
  };
  p.definitions = {{"Se1, Se2", "Sensitivities of tests 1 and 2: shares of subjects with the condition testing positive."},
                   {"Sp1, Sp2", "Specificities: shares of subjects without the condition testing negative."},
                   {"ψ", "Share of subjects whose two test results differ (only they inform McNemar's test)."},
                   {"Power", "For both, the product of the two tests' powers (different subjects)."}};
  p.notes = {"With exact power, power is not monotone in N; the smallest N reaching the target is reported.",
             "The dependence between the tests matters: conditional independence is a conservative guess when both "
             "tests look at the same feature, and it understates the sample size when their results tend to "
             "disagree."};
  p.references = {"Connor, R. J. (1987). Sample size for testing differences in proportions for the paired-sample "
                  "design. Biometrics 43, 207–211.",
                  "Alonzo, T. A., Pepe, M. S. & Moskowitz, C. S. (2002). Sample size calculations for comparative "
                  "studies of medical tests for detecting presence of disease. Statistics in Medicine 21, 835–852.",
                  "Pepe, M. S. (2003). The Statistical Evaluation of Medical Tests for Classification and Prediction. "
                  "Oxford University Press."};
  p.curve_param = "n";
  p.curve_label = "N (Total)";
  p.curve_min = 2;
  p.example_solve = "n";
  p.example = {{"goal", "se"}, {"power", "0.8"}, {"prev", "0.2 0.3"}, {"se1", "0.9"}, {"se2", "0.8"}};
  return p;
}
SSC_REGISTER(paired_diagnostic);

}  // namespace
}  // namespace ssc::procs
