// Non-inferiority / superiority-by-a-margin and equivalence (TOST) tests for
// one mean (or paired differences) and for two independent means.
#include "procedures/common.hpp"
#include "procedures/integrals.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// ------------------------------------------------------------------ margins
// H0: δ ≤ δ0 vs H1: δ > δ0 (upper) or H0: δ ≥ δ0 vs H1: δ < δ0 (lower).
struct MarginTest {
  double delta0;
  bool upper;
};

MarginTest margin_test(const Scenario& s) {
  const double m = std::fabs(s("margin"));
  const bool higher = s.is("direction", "higher");
  const bool ni = s.is("testtype", "ni");
  if (ni) return {higher ? -m : m, higher};
  return {higher ? m : -m, higher};
}

double one_sided_power(double ncp, double df, double alpha, bool upper, bool z) {
  using namespace stats;
  if (z) {
    const double zc = qnorm(alpha, false);
    return upper ? pnorm(ncp - zc) : pnorm(-ncp - zc);
  }
  const double t = qt(alpha, df, false);
  return upper ? pnt(t, df, ncp, false) : pnt(-t, df, ncp, true);
}

std::vector<ParamDef> margin_params(const std::string& diff_label, const std::string& diff_sym) {
  return {
      select("testtype", "Test Type",
             {{"ni", "Non-inferiority", {}}, {"sup", "Superiority by a margin", {}}}, "ni")
          .grp(G_TEST)
          .hint("Non-inferiority: show the new treatment is not worse by more than the margin. Superiority by a "
                "margin: show it is better by more than the margin."),
      select("direction", "Higher Values Are",
             {{"higher", "Better (e.g. response, efficacy)", {}}, {"lower", "Worse (e.g. pain, adverse events)", {}}},
             "higher")
          .grp(G_TEST)
          .hint("Whether a larger value of the outcome is a better result (e.g. a response score) or a worse one "
                "(e.g. a pain score); it sets the side of the margin."),
      select("test", "Test Statistic", {{"t", "T-test (σ estimated)", {}}, {"z", "Z-test (σ known)", {}}}, "t")
          .grp(G_TEST)
          .hint("The t-test uses the sample standard deviation; the z-test assumes σ is known."),
      power_param(),
      alpha_param("0.025", "Alpha (One-Sided Significance Level)").range(0, 0.5, true, true),
      num("margin", "Margin (|NIM| or |SM|)", "0.5")
          .grp(G_EFFECT)
          .sym("Margin")
          .above(0)
          .sugg({"0.5", "0.25 0.5", "1"})
          .hint("Magnitude of the non-inferiority or superiority margin, on the scale of the difference."),
      num("delta", diff_label, "0")
          .grp(G_EFFECT)
          .sym(diff_sym)
          .sugg({"0", "0.1", "1", "-1", "-0.1 0 0.1"})
          .hint("True difference (treatment − reference) at which power is computed."),
  };
}

std::string margin_hypothesis(const Scenario& s, const MarginTest& mt, const std::string& d) {
  return "H0: " + d + (mt.upper ? " ≤ " : " ≥ ") + fnum(mt.delta0) + " versus H1: " + d + (mt.upper ? " > " : " < ") +
         fnum(mt.delta0);
}

// ------------------------------------------------------------------ one mean NI
Procedure one_ni() {
  Procedure p;
  p.id = "means.one.ni";
  p.name = "Non-Inferiority Tests for One Mean (or Paired Means)";
  p.category = "Means";
  p.subcategory = "One Mean";
  p.order = 3;
  p.summary =
      "One-sided t-test (or z-test) of a single mean, or of the mean of paired differences, against a "
      "non-inferiority or superiority margin. Power is exact from the noncentral t distribution.";
  p.keywords = {"non-inferiority", "noninferiority", "superiority by a margin", "one-sided", "paired", "margin", "single-arm", "paired differences"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}};
  p.params = margin_params("δ (True Difference from Reference)", "δ");
  p.params.insert(p.params.begin() + 5,
                  n_param("n", "N (Sample Size)", "30", 2)
                      .sym("N")
                      .hint("Number of subjects, or of pairs when the outcome is a within-pair difference."));
  p.params.insert(p.params.begin() + 6, dropout_param());
  p.params.push_back(num("sigma", "σ (SD of Responses or of Differences)", "1")
                         .grp(G_EFFECT)
                         .sym("σ")
                         .above(0)
                         .sugg({"1", "0.8 1 1.2"})
                         .hint("Standard deviation of the responses, or of the within-pair differences in a paired "
                               "design; take it from earlier studies or a pilot."));
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const MarginTest mt = margin_test(s);
    const double sigma = s("sigma"), alpha = s("alpha"), delta = s("delta");
    const bool z = s.is("test", "z");
    auto P = [&](double n) { return one_sided_power((delta - mt.delta0) * std::sqrt(n) / sigma, n - 1, alpha, mt.upper, z); };
    double n;
    if (s.solve == "n") {
      const double eff = mt.upper ? delta - mt.delta0 : mt.delta0 - delta;
      if (eff <= 0) throw CalcError("the true difference does not lie inside the alternative hypothesis");
      solve::NSearch o;
      o.lo = 2;
      o.guess = z_factor(alpha, s("power"), Alt::Upper) * sigma * sigma / (eff * eff) + 1;
      n = solve::min_n(P, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    r.set("n", n);
    put_power(s, r, P(n));
    r.set("margin", std::fabs(s("margin")));
    r.set("delta0", mt.delta0);
    r.set("delta", delta);
    r.set("sigma", sigma);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("margin", "Margin", "num"), col("delta0", "δ0 (H0 bound)", "num"),
                   col("delta", "δ", "num", "True difference"), col("sigma", "σ", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const MarginTest mt = margin_test(s);
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) + " for the one-sided test of " +
           margin_hypothesis(s, mt, "δ") + " (" + (s.is("testtype", "ni") ? "non-inferiority" : "superiority") +
           " margin " + fnum(r.num("margin")) + ") when the true difference is δ = " + fnum(r.num("delta")) +
           " and the standard deviation is " + fnum(r.num("sigma")) + ", at significance level " +
           fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"Margin", "Non-inferiority (NIM) or superiority (SM) margin, entered as a positive number."},
                   {"δ0", "Boundary of the null hypothesis implied by the margin and the direction."},
                   {"δ", "True mean (or mean difference) relative to the reference value."},
                   {"Alpha", "One-sided significance level."}};
  p.references = {REF_CHOW, REF_JULIOUS,
                  "Blackwelder, W. C. (1982). Proving the null hypothesis in clinical trials. Controlled Clinical "
                  "Trials 3, 345–353."};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"margin", "0.5"}, {"delta", "0 0.1"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(one_ni);

// ------------------------------------------------------------------ two means NI
Procedure two_ni() {
  Procedure p;
  p.id = "means.two.ni";
  p.name = "Non-Inferiority Tests for Two Means";
  p.category = "Means";
  p.subcategory = "Two Independent Means";
  p.order = 5;
  p.summary =
      "One-sided two-sample t-test (or z-test) of the difference between a treatment and a reference mean against a "
      "non-inferiority or superiority margin, for any allocation ratio.";
  p.keywords = {"non-inferiority", "noninferiority", "superiority by a margin", "two-sample", "margin", "active control", "active comparator"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = margin_params("δ (True Difference μ1 − μ2)", "δ");
  {
    auto alloc = two_group_params(2);
    p.params.insert(p.params.begin() + 5, alloc.begin(), alloc.end());
    p.params.insert(p.params.begin() + 5 + static_cast<long>(alloc.size()), dropout_param());
  }
  p.params.push_back(num("sigma", "σ (Common Standard Deviation)", "1")
                         .grp(G_EFFECT)
                         .sym("σ")
                         .above(0)
                         .sugg({"1", "0.8 1 1.2"})
                         .hint("Standard deviation of the outcome within each group, assumed the same in both "
                               "groups."));
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const MarginTest mt = margin_test(s);
    const double sigma = s("sigma"), alpha = s("alpha"), delta = s("delta");
    const bool z = s.is("test", "z");
    auto P = [&](double n1, double n2) {
      const double se = sigma * std::sqrt(1 / n1 + 1 / n2);
      return one_sided_power((delta - mt.delta0) / se, n1 + n2 - 2, alpha, mt.upper, z);
    };
    TwoN n;
    if (s.solve == "n") {
      const double eff = mt.upper ? delta - mt.delta0 : mt.delta0 - delta;
      if (eff <= 0) throw CalcError("the true difference does not lie inside the alternative hypothesis");
      const double guess = z_factor(alpha, s("power"), Alt::Upper) * 2 * sigma * sigma / (eff * eff) + 1;
      n = solve_two_group_n(s, P, s("power"), 2, guess, true, ctx);
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    put_power(s, r, P(n.n1, n.n2));
    r.set("margin", std::fabs(s("margin")));
    r.set("delta0", mt.delta0);
    r.set("delta", delta);
    r.set("sigma", sigma);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("margin", "Margin", "num"), col("delta0", "δ0 (H0 bound)", "num"),
                   col("delta", "δ", "num", "True difference μ1 − μ2"), col("sigma", "σ", "num"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const MarginTest mt = margin_test(s);
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
           " for the one-sided test of " + margin_hypothesis(s, mt, "μ1 − μ2") + " (" +
           (s.is("testtype", "ni") ? "non-inferiority" : "superiority") + " margin " + fnum(r.num("margin")) +
           ") when the true difference is " + fnum(r.num("delta")) + " and the common standard deviation is " +
           fnum(r.num("sigma")) + ", at significance level " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"N1, N2", "Sample sizes of the treatment (1) and reference (2) groups."},
                   {"Margin", "Non-inferiority or superiority margin (positive)."},
                   {"δ0", "Null boundary for μ1 − μ2 implied by the margin and direction."},
                   {"δ", "True difference μ1 − μ2."}};
  p.references = {REF_CHOW, REF_JULIOUS,
                  "Blackwelder, W. C. (1982). Proving the null hypothesis in clinical trials. Controlled Clinical "
                  "Trials 3, 345–353."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"margin", "0.5"}, {"delta", "0"}, {"sigma", "1 1.5"}};
  return p;
}
SSC_REGISTER(two_ni);

// ------------------------------------------------------------------ equivalence
std::vector<ParamDef> limit_params(const std::string& d) {
  return {
      select("limits", "Equivalence Limits",
             {{"sym", "Symmetric: −E and +E", {}}, {"asym", "Lower and upper limits", {}}}, "sym")
          .grp(G_EFFECT)
          .hint("Use limits placed symmetrically around zero (−E and +E), or enter a separate lower and upper "
                "limit."),
      num("e", "E (Equivalence Margin)", "0.5")
          .grp(G_EFFECT)
          .sym("E")
          .above(0)
          .sugg({"0.5", "0.25 0.5", "1"})
          .hint("Equivalence is concluded when the difference lies within (−E, E).")
          .show(when(is("limits", {"sym"}))),
      num("el", "EL (Lower Equivalence Limit)", "-0.5")
          .grp(G_EFFECT)
          .sym("EL")
          .sugg({"-0.5", "-1"})
          .hint("Lower end of the range of differences that count as equivalent (usually negative); it must be "
                "below EU.")
          .show(when(is("limits", {"asym"}))),
      num("eu", "EU (Upper Equivalence Limit)", "0.5")
          .grp(G_EFFECT)
          .sym("EU")
          .sugg({"0.5", "1"})
          .hint("Upper end of the range of differences that count as equivalent (usually positive); it must be "
                "above EL.")
          .show(when(is("limits", {"asym"}))),
      num("delta", "δ (True Difference " + d + ")", "0")
          .grp(G_EFFECT)
          .sym("δ")
          .sugg({"0", "0.1", "0 0.1 0.2"})
          .hint("True difference at which power is computed; must lie between the limits for useful power."),
  };
}

std::pair<double, double> limits(const Scenario& s) {
  if (s.is("limits", "asym")) {
    const double el = s("el"), eu = s("eu");
    if (!(el < eu)) throw CalcError("the lower equivalence limit must be below the upper limit");
    return {el, eu};
  }
  return {-s("e"), s("e")};
}

Procedure one_equiv() {
  Procedure p;
  p.id = "means.one.equiv";
  p.name = "Equivalence Tests for One Mean (or Paired Means)";
  p.category = "Means";
  p.subcategory = "One Mean";
  p.order = 4;
  p.summary =
      "Two one-sided tests (TOST) showing that a mean, or the mean of paired differences, lies within equivalence "
      "limits. Power is computed exactly by integrating over the distribution of the sample standard deviation.";
  p.keywords = {"equivalence", "tost", "two one-sided tests", "paired", "one-sample", "paired equivalence"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}};
  p.params = {select("test", "Test Statistic", {{"t", "T-tests (σ estimated)", {}}, {"z", "Z-tests (σ known)", {}}}, "t")
                  .grp(G_TEST)
                  .hint("The t-tests use the sample standard deviation; the z-tests assume σ is known."),
              power_param(),
              alpha_param("0.05", "Alpha (Each One-Sided Test)").range(0, 0.5, true, true),
              n_param("n", "N (Sample Size)", "30", 2)
                  .sym("N")
                  .hint("Number of subjects, or of pairs when the outcome is a within-pair difference."),
              dropout_param()};
  for (auto& q : limit_params("from Reference")) p.params.push_back(q);
  p.params.push_back(num("sigma", "σ (Standard Deviation)", "1")
                         .grp(G_EFFECT)
                         .sym("σ")
                         .above(0)
                         .sugg({"1", "0.8 1"})
                         .hint("Standard deviation of the responses, or of the within-pair differences in a paired "
                               "design; take it from earlier studies or a pilot."));
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const auto lim = limits(s);
    const double sigma = s("sigma"), alpha = s("alpha"), delta = s("delta");
    const bool z = s.is("test", "z");
    if (!(delta > lim.first && delta < lim.second))
      throw CalcError("the true difference must lie strictly between the equivalence limits");
    auto P = [&](double n) { return tost_power(lim.first, lim.second, delta, sigma / std::sqrt(n), z ? INFINITY : n - 1, alpha); };
    double n;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = 2;
      const double gap = std::min(delta - lim.first, lim.second - delta);
      const double beta = 1 - s("power");
      const double zf = stats::qnorm(alpha, false) + stats::qnorm(std::fabs(delta - (lim.first + lim.second) / 2) < 1e-12 ? 1 - beta / 2 : 1 - beta);
      o.guess = zf * zf * sigma * sigma / (gap * gap) + 1;
      n = solve::min_n(P, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    r.set("n", n);
    put_power(s, r, P(n));
    r.set("el", lim.first);
    r.set("eu", lim.second);
    r.set("delta", delta);
    r.set("sigma", sigma);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("el", "EL", "num", "Lower limit"), col("eu", "EU", "num", "Upper limit"),
                   col("delta", "δ", "num", "True difference"), col("sigma", "σ", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) +
           " to conclude equivalence with two one-sided " + (s.is("test", "z") ? "z-tests" : "t-tests") +
           " at level " + fnum(r.num("alpha")) + " each, when the equivalence limits are " + fnum(r.num("el")) +
           " and " + fnum(r.num("eu")) + ", the true difference is " + fnum(r.num("delta")) +
           " and the standard deviation is " + fnum(r.num("sigma")) + ".";
  };
  p.definitions = {{"EL, EU", "Equivalence limits for the difference."},
                   {"δ", "True difference (mean minus reference value, or mean paired difference)."},
                   {"Alpha", "Level of each one-sided test (the TOST procedure has size ≤ alpha)."}};
  p.references = {REF_CHOW, REF_JULIOUS,
                  "Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and the power approach "
                  "for assessing the equivalence of average bioavailability. Journal of Pharmacokinetics and "
                  "Biopharmaceutics 15, 657–680.",
                  "Owen, D. B. (1965). A special case of a bivariate non-central t-distribution. Biometrika 52, 437–446."};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"e", "0.5"}, {"delta", "0 0.1"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(one_equiv);

Procedure two_equiv() {
  Procedure p;
  p.id = "means.two.equiv";
  p.name = "Equivalence Tests for Two Means (TOST)";
  p.category = "Means";
  p.subcategory = "Two Independent Means";
  p.order = 6;
  p.summary =
      "Two one-sided two-sample tests (TOST) for showing that the difference between two means lies within "
      "equivalence limits. Exact power integrates over the distribution of the pooled standard deviation.";
  p.keywords = {"equivalence", "tost", "two one-sided tests", "parallel", "two-sample", "parallel-group equivalence"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {select("test", "Test Statistic", {{"t", "T-tests (σ estimated)", {}}, {"z", "Z-tests (σ known)", {}}}, "t")
                  .grp(G_TEST)
                  .hint("The t-tests use the sample standard deviation; the z-tests assume σ is known."),
              power_param(), alpha_param("0.05", "Alpha (Each One-Sided Test)").range(0, 0.5, true, true)};
  for (auto& q : two_group_params(2)) p.params.push_back(q);
  p.params.push_back(dropout_param());
  for (auto& q : limit_params("μ1 − μ2")) p.params.push_back(q);
  p.params.push_back(num("sigma", "σ (Common Standard Deviation)", "1")
                         .grp(G_EFFECT)
                         .sym("σ")
                         .above(0)
                         .sugg({"1", "0.8 1"})
                         .hint("Standard deviation of the outcome within each group, assumed the same in both "
                               "groups."));
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const auto lim = limits(s);
    const double sigma = s("sigma"), alpha = s("alpha"), delta = s("delta");
    const bool z = s.is("test", "z");
    if (!(delta > lim.first && delta < lim.second))
      throw CalcError("the true difference must lie strictly between the equivalence limits");
    auto P = [&](double n1, double n2) {
      return tost_power(lim.first, lim.second, delta, sigma * std::sqrt(1 / n1 + 1 / n2), z ? INFINITY : n1 + n2 - 2,
                        alpha);
    };
    TwoN n;
    if (s.solve == "n") {
      const double gap = std::min(delta - lim.first, lim.second - delta);
      const double beta = 1 - s("power");
      const double zf = stats::qnorm(alpha, false) +
                        stats::qnorm(std::fabs(delta - (lim.first + lim.second) / 2) < 1e-12 ? 1 - beta / 2 : 1 - beta);
      n = solve_two_group_n(s, P, s("power"), 2, 2 * zf * zf * sigma * sigma / (gap * gap) + 1, true, ctx);
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    put_power(s, r, P(n.n1, n.n2));
    r.set("el", lim.first);
    r.set("eu", lim.second);
    r.set("delta", delta);
    r.set("sigma", sigma);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("el", "EL", "num", "Lower limit"), col("eu", "EU", "num", "Upper limit"),
                   col("delta", "δ", "num", "True difference μ1 − μ2"), col("sigma", "σ", "num"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
           " to conclude equivalence with two one-sided " + (s.is("test", "z") ? "z-tests" : "t-tests") +
           " at level " + fnum(r.num("alpha")) + " each, when the equivalence limits for μ1 − μ2 are " +
           fnum(r.num("el")) + " and " + fnum(r.num("eu")) + ", the true difference is " + fnum(r.num("delta")) +
           " and the common standard deviation is " + fnum(r.num("sigma")) + ".";
  };
  p.definitions = {{"EL, EU", "Equivalence limits for μ1 − μ2."},
                   {"δ", "True difference μ1 − μ2."},
                   {"Alpha", "Level of each one-sided test."}};
  p.references = {REF_CHOW, REF_JULIOUS,
                  "Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and the power approach "
                  "for assessing the equivalence of average bioavailability. Journal of Pharmacokinetics and "
                  "Biopharmaceutics 15, 657–680."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"e", "0.5"}, {"delta", "0 0.1"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(two_equiv);

}  // namespace
}  // namespace ssc::procs
