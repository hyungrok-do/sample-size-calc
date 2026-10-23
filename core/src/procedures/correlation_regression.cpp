// Correlation (one and two correlations, CI precision) and regression
// (multiple regression R², simple linear regression slope, logistic and Cox
// regression with one covariate of interest).
#include "procedures/apart.hpp"
#include "procedures/common.hpp"
#include "procedures/integrals.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// The exact distribution of the sample correlation is stats::pcorr (stats::pcorr_z on Fisher's scale).

// Root of h(z) = 0 for h increasing in Fisher's z = atanh(r) on [zmin, zmax], to within tol + rtol·|z|.  The
// bracket starts narrow around zg, where the normal-theory guess is accurate to O(1/n), and widens only if needed:
// few exact CDF evaluations, and none near r = ±1 where they are expensive.
template <class H>
double corr_root_z(H&& h, double zg, double zmin, double zmax, double tol = 1e-12, double rtol = 0) {
  zg = std::min(zmax, std::max(zmin, zg));
  double w = 0.05;
  double lo = std::max(zmin, zg - w), hi = std::min(zmax, zg + w);
  double hlo = h(lo), hhi = h(hi);
  while (hlo > 0 && lo > zmin) {
    hi = lo;
    hhi = hlo;
    w *= 3;
    lo = std::max(zmin, zg - w);
    hlo = h(lo);
  }
  while (hhi < 0 && hi < zmax) {
    lo = hi;
    hlo = hhi;
    w *= 3;
    hi = std::min(zmax, zg + w);
    hhi = h(hi);
  }
  if ((hlo > 0) == (hhi > 0)) return hlo > 0 ? lo : hi;
  return numerics::brent(h, lo, hi, hlo, hhi, tol, 300, rtol);
}

// The critical values are found and kept on Fisher's scale: the doubles next to a critical value near ±1 are too
// coarse to hold it (those next to 1 − 1e-14 are a hundredth of 1e-14 apart).  |z| up to 60 covers an α of 1e-30
// with 4 subjects at the last ρ0 below 1.
constexpr double Z_LIMIT = 60;

// Critical value z_c = atanh(r_c) with P(r ≥ r_c | rho0) = a (upper) under the exact null distribution.
double zcrit_upper(double n, double rho0, double a) {
  if (rho0 == 0) {
    // r = t/√(n − 2 + t²), whose atanh is asinh(t/√(n − 2)): r itself rounds to 1 for a tiny alpha with few
    // subjects (t of 1e8 or more), where atanh(r) was infinite and the side was dropped
    return std::asinh(stats::qt(a, n - 2, false) / std::sqrt(n - 2));
  }
  const double zeta = std::atanh(rho0);
  return corr_root_z([&](double z) { return a - stats::pcorr_z(z, n, zeta, false); },
                     zeta + stats::qnorm(a, false) / std::sqrt(n - 3), -Z_LIMIT, Z_LIMIT);
}
double zcrit_lower(double n, double rho0, double a) {
  if (rho0 == 0) return -zcrit_upper(n, 0, a);
  const double zeta = std::atanh(rho0);
  return corr_root_z([&](double z) { return stats::pcorr_z(z, n, zeta, true) - a; },
                     zeta - stats::qnorm(a, false) / std::sqrt(n - 3), -Z_LIMIT, Z_LIMIT);
}

// Exact critical values of the sample correlation on Fisher's scale (NaN for an unused side).
struct CorrCrit {
  double up = NAN, lo = NAN;
};

CorrCrit corr_crit(const Scenario& s, double n) {
  CorrCrit c;
  if (s.is("method", "fisher")) return c;
  const double rho0 = s("rho0"), alpha = s("alpha");
  const Alt a = alt(s);
  const double side = tail_alpha(alpha, a);
  if (a != Alt::Lower) c.up = zcrit_upper(n, rho0, side);
  if (a != Alt::Upper) c.lo = zcrit_lower(n, rho0, side);
  return c;
}

// The power at ζ1 = atanh(ρ1), at the distance t = ζ1 − atanh(ρ0) on Fisher's scale (the one the approximation
// takes, the other the exact distribution).
double corr_power_z(const Scenario& s, double n, double zeta1, double t, const CorrCrit& c) {
  const double alpha = s("alpha");
  const Alt a = alt(s);
  if (s.is("method", "fisher")) return z_power(t * std::sqrt(n - 3), alpha, a);
  double pw = 0;
  if (std::isfinite(c.up)) pw += stats::pcorr_z(c.up, n, zeta1, false);
  if (std::isfinite(c.lo)) pw += stats::pcorr_z(c.lo, n, zeta1, true);
  return stats::clamp01(pw);
}

double corr_power(const Scenario& s, double n, double rho1, const CorrCrit& c) {
  const double zeta1 = std::atanh(rho1);
  return corr_power_z(s, n, zeta1, zeta1 - std::atanh(s("rho0")), c);
}

double corr_power(const Scenario& s, double n, double rho1) { return corr_power(s, n, rho1, corr_crit(s, n)); }

Procedure one_corr() {
  Procedure p;
  p.id = "corr.one";
  p.name = "Tests for One Correlation (Pearson)";
  p.category = "Correlation";
  p.subcategory = "One Correlation";
  p.order = 1;
  p.summary =
      "Sample size and power for testing H0: ρ = ρ0 for a Pearson correlation (bivariate normal data). The exact "
      "method uses the exact distribution of the sample correlation for any ρ0; Fisher's z approximation is also "
      "available.";
  p.keywords = {"pearson", "correlation", "association", "fisher z", "bivariate normal", "correlation coefficient", "pearson r"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}, {"rho1", "Effect Size (ρ1)", "rho1", "Detectable ρ1"}};
  p.params = {
      alternative("Two-sided (H1: ρ ≠ ρ0)", "One-sided (H1: ρ > ρ0)", "One-sided (H1: ρ < ρ0)"),
      select("method", "Calculation Method",
             {{"exact", "Exact distribution of r", {}}, {"fisher", "Fisher's z transformation", {}}}, "exact")
          .grp(G_TEST)
          .hint("The exact method uses the exact distribution of the sample correlation; Fisher's z "
                "transformation is a quick normal approximation, a little less accurate in small samples."),
      power_param(),
      alpha_param(),
      n_param("n", "N (Sample Size)", "50", 4).sym("N").hint("Number of subjects on whom both variables are measured."),
      dropout_param(),
      num("rho0", "ρ0 (Correlation under H0)", "0")
          .grp(G_EFFECT)
          .sym("ρ0")
          .range(-1, 1, true, true)
          .sugg({"0", "0.3"})
          .hint("The correlation assumed by the null hypothesis, between −1 and 1; usually 0 (no linear "
                "association)."),
      num("rho1", "ρ1 (Correlation under H1)", "0.3")
          .grp(G_EFFECT)
          .sym("ρ1")
          .range(-1, 1, true, true)
          .sugg({"0.3", "0.1 0.3 0.5"})
          .hint("The true correlation you expect, between −1 and 1; the study is sized to tell it apart from ρ0.")
          .show(when(not_solving({"rho1"}))),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double rho0 = s("rho0");
    const Alt a = alt(s);
    double n, rho1;
    if (s.solve == "n") {
      rho1 = s("rho1");
      if (rho1 == rho0) throw CalcError("ρ1 equals ρ0");
      if ((a == Alt::Upper && rho1 < rho0) || (a == Alt::Lower && rho1 > rho0))
        throw CalcError("ρ1 lies on the side opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = 4;
      const double dz = std::atanh(rho1) - std::atanh(rho0);
      o.guess = z_factor(s("alpha"), s("power"), a) / (dz * dz) + 3;
      n = solve::min_n([&](double x) { return corr_power(s, x, rho1); }, s("power"), o, &ctx);
    } else {
      n = s("n");
      if (s.solve == "power") {
        rho1 = s("rho1");
      } else {
        const bool up = a != Alt::Lower;
        const double target = s("power"), alpha = s("alpha");
        if (n < 4) throw CalcError("N must be at least 4");
        // Fisher-z guess, then the exact power on the side of the alternative
        const double z = stats::qnorm(tail_alpha(alpha, a), false) + stats::qnorm(target);
        const double guess = std::tanh(std::atanh(rho0) + (up ? z : -z) / std::sqrt(n - 3));
        const CorrCrit crit = corr_crit(s, n);
        if (corr_power(s, n, rho0, crit) >= target)  // at ρ1 = ρ0 the power is the significance level
          throw CalcError("the target power must exceed the power with no effect (the significance level)");
        // up to the last correlation below 1 (or above −1), so that a ρ0 near it keeps room for ρ1
        const double last = up ? std::nextafter(1.0, 0.0) : std::nextafter(-1.0, 0.0);
        auto power = [&](double x) { return corr_power(s, n, x, crit); };
        if (up ? rho0 >= last : rho0 <= last)
          throw CalcError("ρ1 cannot be told apart from ρ0 = " + text::format_roundtrip(rho0) +
                          ": no number lies between it and " + (up ? "1" : "−1"));
        const std::string unreachable = "the target power cannot be reached with this N";
        // the search is for the distance t = z1 − z0 on Fisher's scale, to a tolerance relative to it (one of 1e-12
        // on z left the power 5e-9 off the target with 1e9 subjects), with the power taken at t itself: at the
        // double ρ1 = tanh(z0 + t) it changes in steps near ±1, and the root at a step's edge stood half a double
        // off the difference it needs
        const double z0 = std::atanh(rho0), tg = std::atanh(guess) - z0, t_last = std::atanh(last) - z0;
        auto power_t = [&](double x) { return corr_power_z(s, n, z0 + x, x, crit); };
        const double t = up ? corr_root_z([&](double x) { return power_t(x) - target; }, tg, 0, t_last, 0, 1e-12)
                            : corr_root_z([&](double x) { return target - power_t(x); }, tg, t_last, 0, 0, 1e-12);
        // a search that ends at the last correlation found no root: the target is out of reach, whatever the
        // difference from ρ0 would have been
        if (t == t_last && power(last) < target - solve::power_band(target)) throw CalcError(unreachable);
        // a double next to the root, which must keep the difference tanh(z1) − ρ0 = sinh(z1 − z0)/(cosh z1 cosh z0)
        // to four digits; near ±1 the power changes by more than 1e-9 from one double to the next, and the root can
        // lie between two of them: the answer is the first that reaches the target
        const double z1 = z0 + t, d = std::sinh(t) / (std::cosh(z1) * std::cosh(z0));
        rho1 = up ? std::fmin(last, std::tanh(z1)) : std::fmax(last, std::tanh(z1));
        rho1 = solve::reach(power, rho1, last, up, target, rho0, d, "ρ1", "ρ0", unreachable);
      }
    }
    r.set("n", n);
    put_power(s, r, corr_power(s, n, rho1));
    r.set("rho0", rho0);
    r.set("rho1", rho1);
    r.set("alpha", s("alpha"));
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("rho0", "ρ0", "num"), col("rho1", "ρ1", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const auto [r1, r0] = fnum_apart(r.num("rho1"), r.num("rho0"));
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) +
           " to detect a correlation of ρ1 = " + r1 + " when the null correlation is ρ0 = " + r0 + ", using a " +
           alt_words(alt(s)) + " test" +
           (s.is("method", "fisher") ? " based on Fisher's z transformation" : " based on the exact distribution of r") +
           " at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"ρ0", "Population correlation under H0."}, {"ρ1", "Population correlation under H1."}};
  p.references = {"Guenther, W. C. (1977). Desk calculation of probabilities for the distribution of the sample "
                  "correlation coefficient. The American Statistician 31, 45–48.",
                  "Fisher, R. A. (1921). On the 'probable error' of a coefficient of correlation deduced from a small "
                  "sample. Metron 1, 3–32.",
                  REF_COHEN, REF_MACHIN};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.curve_min = 4;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"rho0", "0"}, {"rho1", "0.2 0.3 0.4"}};
  return p;
}
SSC_REGISTER(one_corr);

Procedure two_corr() {
  Procedure p;
  p.id = "corr.two";
  p.name = "Tests for Two Independent Correlations";
  p.category = "Correlation";
  p.subcategory = "Two Correlations";
  p.order = 1;
  p.summary = "Compares Pearson correlations from two independent samples using Fisher's z transformation.";
  p.keywords = {"two correlations", "compare correlations", "fisher z", "independent samples", "difference between correlations"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {alternative("Two-sided (H1: ρ1 ≠ ρ2)", "One-sided (H1: ρ1 > ρ2)", "One-sided (H1: ρ1 < ρ2)"),
              power_param(), alpha_param()};
  for (auto& q : two_group_params(4, "100")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  p.params.push_back(num("rho1", "ρ1 (Correlation in Group 1)", "0.5")
                         .grp(G_EFFECT)
                         .sym("ρ1")
                         .range(-1, 1, true, true)
                         .sugg({"0.5"})
                         .hint("The correlation between the two variables that you expect in group 1, between −1 "
                               "and 1."));
  p.params.push_back(num("rho2", "ρ2 (Correlation in Group 2)", "0.3")
                         .grp(G_EFFECT)
                         .sym("ρ2")
                         .range(-1, 1, true, true)
                         .sugg({"0.3", "0.1 0.2 0.3"})
                         .hint("The correlation between the same two variables that you expect in group 2; the test "
                               "compares it with ρ1."));
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double r1 = s("rho1"), r2 = s("rho2"), alpha = s("alpha");
    const Alt a = alt(s);
    const double dz = std::atanh(r1) - std::atanh(r2);
    auto P = [&](double n1, double n2) { return z_power(dz / std::sqrt(1 / (n1 - 3) + 1 / (n2 - 3)), alpha, a); };
    TwoN n;
    if (s.solve == "n") {
      if (dz == 0) throw CalcError("ρ1 equals ρ2");
      if ((a == Alt::Upper && dz < 0) || (a == Alt::Lower && dz > 0))
        throw CalcError("ρ1 − ρ2 has the sign opposite to the one-sided alternative");
      n = solve_two_group_n(s, P, s("power"), 4, 2 * z_factor(alpha, s("power"), a) / (dz * dz) + 3, true, ctx);
    } else {
      n = two_group_n(s);
      if (n.n1 < 4 || n.n2 < 4) throw CalcError("each group needs at least 4 subjects");
    }
    put_two_group(r, n);
    put_power(s, r, P(n.n1, n.n2));
    r.set("rho1", r1);
    r.set("rho2", r2);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("rho1", "ρ1", "num"), col("rho2", "ρ2", "num"), col("alpha", "Alpha", "alpha")}) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const auto [r1, r2] = fnum_apart(r.num("rho1"), r.num("rho2"));
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
           " to detect a difference between correlations of " + r1 + " and " + r2 +
           " using a " + alt_words(alt(s)) + " Fisher z-test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"ρ1, ρ2", "Population correlations in the two groups."}};
  p.references = {REF_COHEN, "Fisher, R. A. (1921). On the 'probable error' of a coefficient of correlation deduced "
                             "from a small sample. Metron 1, 3–32."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.curve_min = 4;
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8"}, {"rho1", "0.5"}, {"rho2", "0.1 0.2 0.3"}};
  return p;
}
SSC_REGISTER(two_corr);

Procedure corr_ci() {
  Procedure p;
  p.id = "corr.ci";
  p.name = "Confidence Interval for a Correlation";
  p.category = "Confidence Intervals";
  p.subcategory = "Correlation";
  p.order = 1;
  p.summary = "Sample size giving a Fisher-z confidence interval for a Pearson correlation of a desired width.";
  p.keywords = {"correlation", "confidence interval", "precision", "fisher z", "fisher z interval", "precision of a correlation"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"width", "Interval Width", "width", "Width"}};
  p.params = {num("conf", "Confidence Level", "0.95")
                  .grp(G_TEST)
                  .sym("Confidence")
                  .range(0, 1, true, true)
                  .sugg({"0.95", "0.9 0.95 0.99"})
                  .hint("How often intervals computed this way contain the true correlation; 0.95 gives a 95% "
                        "interval."),
              n_param("n", "N (Sample Size)", "100", 4)
                  .sym("N")
                  .hint("Number of subjects on whom both variables are measured.")
                  .show(when(not_solving({"n"}))),
              num("width", "Width (Upper − Lower Limit)", "0.2")
                  .grp(G_EFFECT)
                  .sym("Width")
                  .range(0, 2, true, true)
                  .sugg({"0.2", "0.1 0.2 0.3"})
                  .hint("The precision you want: the distance from the lower to the upper limit of the two-sided "
                        "interval for the correlation.")
                  .show(when(solving({"n"}))),
              num("rho", "r (Anticipated Sample Correlation)", "0.5")
                  .grp(G_EFFECT)
                  .sym("r")
                  .range(-1, 1, true, true)
                  .sugg({"0.3", "0.5", "0.7"})
                  .hint("The correlation you expect to observe. The interval is widest near 0 and narrows toward −1 "
                        "or 1, so a value closer to 0 is the cautious choice.")};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double rho = s("rho"), z = stats::qnorm((1 - s("conf")) / 2, false);
    auto width = [&](double n) {
      const double h = z / std::sqrt(n - 3);
      return std::tanh(std::atanh(rho) + h) - std::tanh(std::atanh(rho) - h);
    };
    double n;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = 4;
      const double target = s("width");
      o.goal = "the target width";
      n = solve::min_n([&](double x) { return width(x) <= target ? 1.0 : 0.0; }, 0.5, o, &ctx);
      r.set("target_width", target);
    } else {
      n = s("n");
    }
    const double h = z / std::sqrt(n - 3);
    r.set("n", n);
    r.set("width", width(n));
    r.set("lower", std::tanh(std::atanh(rho) - h));
    r.set("upper", std::tanh(std::atanh(rho) + h));
    r.set("rho", rho);
    r.set("conf", s("conf"));
  };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c = {col("conf", "Confidence", "num"), col("n", "N", "int")};
    if (s.solve == "n") c.push_back(col("target_width", "Target Width", "num"));
    for (auto x : {col("width", s.solve == "n" ? "Actual Width" : "Width", "num"), col("rho", "r", "num"),
                   col("lower", "Lower Limit", "num"), col("upper", "Upper Limit", "num")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const auto [lower, upper] = fnum_apart(r.num("lower"), r.num("upper"), {-1.0, 1.0});
    return "A sample of " + fint(r.num("n")) + " gives a " + level_text(r.num("conf")) +
           " Fisher-z confidence interval of width " + fnum(r.num("width")) + " (" + lower + " to " + upper +
           ") when the sample correlation is " + fnum_fine(r.num("rho")) + ".";
  };
  p.references = {REF_MACHIN};
  p.example_solve = "n";
  p.example = {{"width", "0.1 0.2"}, {"rho", "0.3 0.5 0.7"}};
  return p;
}
SSC_REGISTER(corr_ci);

// =====================================================================================
// Multiple regression (test of R² or of an R² increase)
// =====================================================================================
double mreg_power(double n, double f2, double b, double a, double alpha) {
  const double df2 = n - a - b - 1;
  if (!(df2 > 0)) return NAN;
  const double fc = stats::qf(alpha, b, df2, false);
  return stats::pnf(fc, b, df2, n * f2, false);
}

Procedure multiple_regression() {
  Procedure p;
  p.id = "reg.multiple";
  p.name = "Multiple Regression (R² Increase F-Test)";
  p.category = "Regression";
  p.subcategory = "Linear Regression";
  p.order = 1;
  p.summary =
      "Sample size and power for the F-test that a set of B predictors adds to a regression already containing A "
      "predictors (A = 0 tests the overall R²), using Cohen's f² = ΔR²/(1 − R²full) and the noncentral F "
      "distribution with λ = N f².";
  p.keywords = {"multiple regression", "r-squared", "f-test", "cohen f2", "linear model", "predictors", "ancova", "r2 change", "hierarchical regression", "partial f-test", "covariates"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}, {"r2t", "Effect Size (ΔR²)", "r2t", "Detectable ΔR²"}};
  p.params = {
      power_param(),
      alpha_param(),
      n_param("n", "N (Sample Size)", "100", 3)
          .sym("N")
          .hint("Number of subjects (observations) in the regression; it must exceed A + B + 1."),
      dropout_param(),
      num("b", "B (Number of Tested Predictors)", "3")
          .grp(G_EFFECT)
          .sym("B")
          .at_least(1)
          .whole()
          .sugg({"1", "3", "5"})
          .hint("Number of predictors whose joint contribution is tested."),
      num("r2t", "ΔR² (R² Due to Tested Predictors)", "0.1")
          .grp(G_EFFECT)
          .sym("ΔR²")
          .range(0, 1, true, true)
          .sugg({"0.02", "0.13", "0.26", "0.05 0.1 0.15"})
          .hint("The share of the outcome's variance that the tested predictors add to the control predictors. "
                "Without control predictors, Cohen calls 0.02 small, 0.13 medium and 0.26 large.")
          .show(when(not_solving({"r2t"}))),
      num("a", "A (Number of Control Predictors)", "0")
          .grp(G_EFFECT)
          .sym("A")
          .at_least(0)
          .whole()
          .sugg({"0", "2"})
          .hint("Number of predictors adjusted for but not tested; 0 tests the overall R² of the B predictors."),
      num("r2c", "R² of Control Predictors", "0")
          .grp(G_EFFECT)
          .sym("R²c")
          .range(0, 1, false, true)
          .sugg({"0", "0.2"})
          .hint("The share of the outcome's variance already explained by the control predictors, between 0 and "
                "1 (0 when A = 0)."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double b = s("b"), a = s("a"), r2c = s("r2c"), alpha = s("alpha");
    double n, r2t;
    auto f2_of = [&](double d) { return d / (1 - r2c - d); };
    if (s.solve == "r2t") {
      n = s("n");
      if (n < a + b + 2) throw CalcError("N must exceed A + B + 1");
      // solve on the unbounded f² scale, then ΔR² = f²(1 − R²c)/(1 + f²), up to the f² of the last ΔR² that may
      // be given (R²c + ΔR² < 1)
      double last = std::nextafter(1 - r2c, 0.0);
      while (!(r2c + last < 1 && 1 - r2c - last > 0)) last = std::nextafter(last, 0.0);
      const std::string unreachable = "the target power cannot be reached for any admissible value of ΔR²";
      const double f2 = solve::solve_increasing([&](double x) { return mreg_power(n, x, b, a, alpha); }, s("power"),
                                                1e-12, 0.15, f2_of(last), 1e-12, "ΔR²");
      // next to 1 − R²c the doubles lie too far apart in f² = ΔR²/(1 − R²c − ΔR²) to hold the root: the answer is
      // the first that reaches the target, as for a correlation near 1
      r2t = f2 * (1 - r2c) / (1 + f2);
      r2t = solve::reach([&](double x) { return mreg_power(n, f2_of(x), b, a, alpha); }, std::fmin(r2t, last), last,
                         true, s("power"), 0, r2t, "ΔR²", "0", unreachable);
    } else {
      r2t = s("r2t");
      if (!(r2c + r2t < 1)) throw CalcError("R² of the control and tested predictors must sum to less than 1");
      if (s.solve == "n") {
        solve::NSearch o;
        o.lo = a + b + 2;
        const double f2 = f2_of(r2t);
        o.guess = (z_factor(alpha, s("power"), Alt::Two) + 2 * (b - 1)) / f2 + a + b + 1;
        n = solve::min_n([&](double x) { return mreg_power(x, f2, b, a, alpha); }, s("power"), o, &ctx);
      } else {
        n = s("n");
        if (n < a + b + 2) throw CalcError("N must exceed A + B + 1");
      }
    }
    const double f2 = f2_of(r2t);
    r.set("n", n);
    put_power(s, r, mreg_power(n, f2, b, a, alpha));
    r.set("b", b);
    r.set("a", a);
    r.set("r2t", r2t);
    r.set("r2c", r2c);
    r.set("f2", f2);
    r.set("lambda", n * f2);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("b", "B", "int", "Tested predictors"), col("a", "A", "int", "Control predictors"),
                   col("r2t", "ΔR²", "num"), col("r2c", "R²c", "num"), col("f2", "f²", "num", "Cohen's effect size"),
                   col("lambda", "λ", "num", "Noncentrality N f²"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) + " to detect an R² increase of " +
           fnum(r.num("r2t")) + " (f² = " + fnum(r.num("f2")) + ") attributable to " + fint(r.num("b")) +
           " tested predictor(s), adjusting for " + fint(r.num("a")) + " control predictor(s) with R² = " +
           fnum(r.num("r2c")) + ", using an F-test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"B", "Number of predictors being tested."},
                   {"A", "Number of predictors already in the model."},
                   {"ΔR²", "Increase in R² due to the B tested predictors."},
                   {"f²", "Cohen's effect size ΔR²/(1 − R²c − ΔR²)."}};
  p.references = {REF_COHEN, "Gatsonis, C. & Sampson, A. R. (1989). Multiple correlation: exact power and sample size "
                             "calculations. Psychological Bulletin 106, 516–524.", REF_BK};
  p.notes = {"The power treats the predictor values as fixed (noncentrality λ = N f², Cohen's model). When the "
             "predictors are random, as in an observational study, the exact power (Gatsonis and Sampson) is lower: "
             "for B = 3 and ΔR² = 0.1, N = 103 gives 0.79 rather than 0.80 (105 reach 0.80), and the gap grows with "
             "more predictors and smaller N."};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.curve_min = 3;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"b", "1 3 5"}, {"r2t", "0.1"}};
  return p;
}
SSC_REGISTER(multiple_regression);

// =====================================================================================
// Simple linear regression slope
// =====================================================================================
Procedure linreg_slope() {
  Procedure p;
  p.id = "reg.slope";
  p.name = "Simple Linear Regression (Test of the Slope)";
  p.category = "Regression";
  p.subcategory = "Linear Regression";
  p.order = 2;
  p.summary = "Sample size and power for the t-test of a regression slope, H0: β1 = B0, given the SD of the predictor "
              "and the residual (or response) SD.";
  p.keywords = {"linear regression", "slope", "dose response", "trend", "t-test", "simple regression", "regression coefficient", "beta coefficient"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}, {"b1", "Effect Size (β1)", "b1", "Detectable β1"}};
  p.params = {
      alternative("Two-sided (H1: β1 ≠ B0)", "One-sided (H1: β1 > B0)", "One-sided (H1: β1 < B0)"),
      power_param(),
      alpha_param(),
      n_param("n", "N (Sample Size)", "50", 3)
          .sym("N")
          .hint("Number of subjects, each with a value of the predictor X and of the outcome Y."),
      dropout_param(),
      num("b0", "B0 (Slope under H0)", "0")
          .grp(G_EFFECT)
          .sym("B0")
          .sugg({"0"})
          .hint("The slope assumed by the null hypothesis; usually 0 (the mean of Y does not change with X)."),
      num("b1", "β1 (Slope under H1)", "0.5")
          .grp(G_EFFECT)
          .sym("β1")
          .sugg({"0.5", "0.25 0.5 1"})
          .hint("The true slope you expect: the change in the mean of Y for a one-unit increase of X.")
          .show(when(not_solving({"b1"}))),
      num("sx", "σx (SD of the Predictor X)", "1")
          .grp(G_EFFECT)
          .sym("σx")
          .above(0)
          .sugg({"1", "2"})
          .hint("Standard deviation of the predictor X among the subjects, in its own units (Σ(x − x̄)² = (N − 1)σx²); "
                "a wider spread of X gives more power."),
      select("sdtype", "Variability Input",
             {{"resid", "σ: residual SD (SD of Y given X)", {}}, {"sy", "σy: total SD of Y", {}}}, "resid")
          .grp(G_EFFECT)
          .hint("Whether the SD value is the scatter of Y around the regression line (residual SD) or the total "
                "SD of Y, from which the residual SD is derived."),
      num("sigma", "SD Value", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .sugg({"1", "2"})
          .hint("The standard deviation of the kind chosen above, in the units of Y."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double b0 = s("b0"), sx = s("sx"), alpha = s("alpha");
    const Alt a = alt(s);
    auto resid_sd = [&](double b1) {
      if (s.is("sdtype", "sy")) {
        const double v = s("sigma") * s("sigma") - b1 * b1 * sx * sx;
        if (!(v > 0)) throw CalcError("β1·σx must be smaller than σy");
        return std::sqrt(v);
      }
      return s("sigma");
    };
    // from the difference d = β1 − B0, which the search for β1 keeps apart from a large B0 (which would round it away)
    auto Pd = [&](double n, double d, double b1) {
      const double ncp = d * sx * std::sqrt(n - 1) / resid_sd(b1);
      return t_power(ncp, n - 2, alpha, a);
    };
    auto P = [&](double n, double b1) { return Pd(n, b1 - b0, b1); };
    double n, b1;
    if (s.solve == "n") {
      b1 = s("b1");
      if (b1 == b0) throw CalcError("β1 equals B0");
      if ((a == Alt::Upper && (b1 - b0) < 0) || (a == Alt::Lower && (b1 - b0) > 0))
        throw CalcError("β1 − B0 has the sign opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = 3;
      const double sd = resid_sd(b1);
      o.guess = z_factor(alpha, s("power"), a) * sd * sd / ((b1 - b0) * (b1 - b0) * sx * sx) + 2;
      n = solve::min_n([&](double x) { return P(x, b1); }, s("power"), o, &ctx);
    } else {
      n = s("n");
      if (s.solve == "power") {
        b1 = s("b1");
      } else {
        const double sign = a == Alt::Lower ? -1 : 1;
        const double scale = s("sigma") / sx;
        if (s.is("sdtype", "sy") && !(s("sigma") * s("sigma") - b0 * b0 * sx * sx > 0))
          throw CalcError("B0·σx must be smaller than σy (the null slope alone would explain all of Y's variance)");
        const double d = solve::solve_increasing([&](double dd) {
          try {
            return Pd(n, sign * dd, b0 + sign * dd);
          } catch (const CalcError&) {
            return 1.0;
          }
        }, s("power"), 0, scale, std::fmax(1e9 * scale, solve::BOUND_SPAN));
        b1 = solve::offset(b0, sign * d, [&](double dd) { return Pd(n, dd, b0 + dd); }, s("power"), "β1", "B0", true);
      }
    }
    r.set("n", n);
    put_power(s, r, P(n, b1));
    r.set("b0", b0);
    r.set("b1", b1);
    r.set("sx", sx);
    r.set("sigma", resid_sd(b1));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("b0", "B0", "num"), col("b1", "β1", "num"), col("sx", "σx", "num"),
                   col("sigma", "σ (Residual)", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const auto [b1, b0] = fnum_apart(r.num("b1"), r.num("b0"));
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) + " to detect a slope of " + b1 +
           " (null value " + b0 + ") when the predictor's SD is " +
           fnum(r.num("sx")) + " and the residual SD is " + fnum(r.num("sigma")) + ", using a " + alt_words(alt(s)) +
           " t-test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"β1", "Slope under H1."}, {"σx", "Standard deviation of the predictor X."},
                   {"σ", "Residual standard deviation of Y given X."}};
  p.references = {"Dupont, W. D. & Plummer, W. D. (1998). Power and sample size calculations for studies involving "
                  "linear regression. Controlled Clinical Trials 19, 589–601.",
                  REF_MACHIN};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.curve_min = 3;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"b1", "0.25 0.5"}, {"sx", "1"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(linreg_slope);

// =====================================================================================
// Logistic regression with one covariate of interest (Hsieh, Bloch & Larsen 1998)
// =====================================================================================
Procedure logistic() {
  Procedure p;
  p.id = "reg.logistic";
  p.name = "Logistic Regression (One Covariate of Interest)";
  p.category = "Regression";
  p.subcategory = "Logistic Regression";
  p.order = 1;
  p.summary =
      "Sample size and power for the Wald test of one covariate in a logistic regression — binary or normally "
      "distributed — with a variance-inflation adjustment 1/(1 − R²) for correlation with other covariates "
      "(Hsieh, Bloch & Larsen, 1998).";
  p.keywords = {"logistic regression", "odds ratio", "binary outcome", "covariate", "hsieh", "epidemiology", "binary regression", "logit"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: OR ≠ 1)", "One-sided (H1: OR > 1)", "One-sided (H1: OR < 1)"),
      power_param(),
      alpha_param(),
      n_param("n", "N (Sample Size)", "200", 2).sym("N").hint("Total number of subjects in the logistic regression."),
      dropout_param(),
      select("xtype", "Covariate X", {{"binary", "Binary (0/1)", {}}, {"normal", "Continuous, normally distributed", {}}}, "binary")
          .grp(G_EFFECT)
          .hint("Whether the covariate of interest is a yes-or-no (0/1) variable, such as an exposure, or a "
                "continuous variable with a roughly normal distribution."),
      num("p0", "P0 (Event Probability at X = 0 / at the Mean of X)", "0.2")
          .grp(G_EFFECT)
          .sym("P0")
          .range(0, 1, true, true)
          .sugg({"0.1", "0.2", "0.3"})
          .hint("Probability of the event for a subject with X = 0 (binary X) or with X at its mean (continuous "
                "X), between 0 and 1."),
      num("orr", "OR (Odds Ratio: X = 1 vs 0, or per SD of X)", "2")
          .grp(G_EFFECT)
          .sym("OR")
          .above(0)
          .sugg({"1.5", "2", "1.5 2 2.5"})
          .hint("The odds ratio to detect: the odds of the event at X = 1 over those at X = 0 for a binary X, or "
                "the change in odds per one-SD increase of a continuous X."),
      num("px", "Proportion with X = 1", "0.5")
          .grp(G_EFFECT)
          .sym("P(X=1)")
          .range(0, 1, true, true)
          .sugg({"0.5", "0.3"})
          .hint("The share of subjects with X = 1 (e.g. exposed), between 0 and 1.")
          .show(when(is("xtype", {"binary"}))),
      num("r2", "R² of X with Other Covariates", "0")
          .grp(G_EFFECT)
          .sym("R²")
          .range(0, 1, false, true)
          .sugg({"0", "0.1", "0.3"})
          .hint("How much of the variation of X the other covariates in the model explain (from a regression of "
                "X on them); 0 when there are none."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    using namespace stats;
    const double p0 = s("p0"), orr = s("orr"), r2 = s("r2"), alpha = s("alpha");
    const Alt a = alt(s);
    if (orr == 1) throw CalcError("an odds ratio of 1 has no effect to detect");
    if ((a == Alt::Upper && orr < 1) || (a == Alt::Lower && orr > 1))
      throw CalcError("the odds ratio is on the side opposite to the one-sided alternative");
    const double za = qnorm(tail_alpha(alpha, a), false);
    std::function<double(double)> P;
    double p1 = NAN;
    if (s.is("xtype", "binary")) {
      const double B = s("px");
      const double o1 = orr * p0 / (1 - p0);
      p1 = o1 / (1 + o1);
      const double pb = (1 - B) * p0 + B * p1;
      P = [=](double n) {
        const double num = std::sqrt(n * (1 - r2) * (1 - B)) * std::fabs(p1 - p0) - za * std::sqrt(pb * (1 - pb) / B);
        return pnorm(num / std::sqrt(p0 * (1 - p0) + p1 * (1 - p1) * (1 - B) / B));
      };
    } else {
      const double b = std::log(orr);
      P = [=](double n) { return pnorm(std::sqrt(n * (1 - r2) * p0 * (1 - p0)) * std::fabs(b) - za); };
    }
    double n;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = 2;
      o.guess = 100;
      n = solve::min_n(P, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    r.set("n", n);
    put_power(s, r, P(n));
    r.set("p0", p0);
    if (std::isfinite(p1)) r.set("p1", p1);
    r.set("orr", orr);
    if (s.is("xtype", "binary")) r.set("px", s("px"));
    r.set("r2", r2);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    c.push_back(col("n", "N", "int"));
    c.push_back(col("p0", "P0", "num"));
    if (s.is("xtype", "binary")) {
      c.push_back(col("p1", "P1 (X = 1)", "num"));
      c.push_back(col("px", "P(X=1)", "num"));
    }
    for (auto x : {col("orr", "OR", "num"), col("r2", "R²", "num"), col("alpha", "Alpha", "alpha")}) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const bool bin = s.is("xtype", "binary");
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) + " to detect an odds ratio of " +
           fnum(r.num("orr")) + (bin ? " for a binary covariate (" + fpct(r.num("px"), 0) + " exposed)" : " per standard deviation of a normal covariate") +
           " when the event probability is " + fnum(r.num("p0")) + (bin ? " among the unexposed" : " at the covariate mean") +
           " and R² with other covariates is " + fnum(r.num("r2")) + ", using a " + alt_words(alt(s)) +
           " Wald test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"P0", "Event probability at X = 0 (binary X) or at the mean of X (normal X)."},
                   {"OR", "Odds ratio for X = 1 vs 0, or per one-SD increase of X."},
                   {"R²", "Squared multiple correlation of X with the other covariates (VIF = 1/(1 − R²))."}};
  p.references = {"Hsieh, F. Y., Bloch, D. A. & Larsen, M. D. (1998). A simple method of sample size calculation for "
                  "linear and logistic regression. Statistics in Medicine 17, 1623–1634.",
                  "Hsieh, F. Y. (1989). Sample size tables for logistic regression. Statistics in Medicine 8, 795–802."};
  p.notes = {"Hsieh's formula is accurate for a binary covariate. For a normal covariate it overstates the power "
             "when the odds ratio per standard deviation is large: at OR 2.5 (P0 = 0.3, N = 45) the Wald test's "
             "power is about 0.70 by simulation rather than 0.80, and at OR 1.5 about 0.77."};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"p0", "0.2"}, {"orr", "1.5 2 2.5"}, {"px", "0.5"}};
  return p;
}
SSC_REGISTER(logistic);

// =====================================================================================
// Cox regression (Hsieh & Lavori 2000)
// =====================================================================================
Procedure cox() {
  Procedure p;
  p.id = "reg.cox";
  p.name = "Cox Regression (One Covariate of Interest)";
  p.category = "Regression";
  p.subcategory = "Survival Regression";
  p.order = 1;
  p.summary =
      "Number of events and subjects for the test of one covariate in a Cox proportional-hazards model, allowing "
      "correlation with other covariates (Hsieh & Lavori, 2000): D = (z_α + z_β)² / ((1 − R²) σ² (ln HR)²).";
  p.keywords = {"cox", "proportional hazards", "hazard ratio", "survival", "events", "time to event", "multivariable survival", "covariate adjustment"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: HR ≠ 1)", "One-sided (H1: HR > 1)", "One-sided (H1: HR < 1)"),
      power_param(),
      alpha_param(),
      n_param("n", "N (Sample Size)", "200", 1)
          .sym("N")
          .hint("Total number of subjects; the power depends on how many of them have the event."),
      dropout_param(),
      num("hr", "HR (Hazard Ratio per Unit of X)", "1.5")
          .grp(G_EFFECT)
          .sym("HR")
          .above(0)
          .sugg({"0.7", "1.5", "1.3 1.5 2"})
          .hint("The hazard ratio to detect for a one-unit increase of X; for a binary X, the hazard at X = 1 "
                "relative to X = 0."),
      select("xtype", "Covariate X", {{"binary", "Binary (0/1)", {}}, {"cont", "Continuous (give its SD)", {}}},
             "binary")
          .grp(G_EFFECT)
          .hint("Whether the covariate of interest is a yes-or-no (0/1) variable or a continuous one, whose SD "
                "you then give."),
      num("px", "Proportion with X = 1", "0.5")
          .grp(G_EFFECT)
          .sym("P(X=1)")
          .range(0, 1, true, true)
          .sugg({"0.5"})
          .hint("The share of subjects with X = 1, between 0 and 1.")
          .show(when(is("xtype", {"binary"}))),
      num("sx", "σx (SD of X)", "1")
          .grp(G_EFFECT)
          .sym("σx")
          .above(0)
          .sugg({"1"})
          .hint("Standard deviation of X among the subjects, in the units to which the hazard ratio refers.")
          .show(when(is("xtype", {"cont"}))),
      num("pevent", "Overall Probability of an Event", "0.5")
          .grp(G_EFFECT)
          .sym("P(Event)")
          .range(0, 1, true, false)
          .sugg({"0.3", "0.5", "0.7"})
          .hint("The share of subjects expected to have the event (e.g. death or relapse) during the study, "
                "between 0 and 1."),
      num("r2", "R² of X with Other Covariates", "0")
          .grp(G_EFFECT)
          .sym("R²")
          .range(0, 1, false, true)
          .sugg({"0", "0.2"})
          .hint("How much of the variation of X the other covariates in the model explain (from a regression of "
                "X on them); 0 when there are none."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double hr = s("hr"), pe = s("pevent"), r2 = s("r2"), alpha = s("alpha");
    const Alt a = alt(s);
    if (hr == 1) throw CalcError("a hazard ratio of 1 has no effect to detect");
    const double var = s.is("xtype", "binary") ? s("px") * (1 - s("px")) : s("sx") * s("sx");
    const double drift_per_event = std::sqrt((1 - r2) * var) * std::log(hr);
    auto P = [&](double n) { return z_power(std::sqrt(n * pe) * drift_per_event, alpha, a); };
    double n;
    if (s.solve == "n") {
      if ((a == Alt::Upper && hr < 1) || (a == Alt::Lower && hr > 1))
        throw CalcError("the hazard ratio is on the side opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = 1;
      o.guess = z_factor(alpha, s("power"), a) / (drift_per_event * drift_per_event * pe);
      n = solve::min_n(P, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    r.set("n", n);
    put_power(s, r, P(n));
    r.set("events", std::ceil(n * pe - 1e-9));
    r.set("hr", hr);
    r.set("var", var);
    r.set("pevent", pe);
    r.set("r2", r2);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("events", "Events", "int", "Expected number of events N × P(Event)"),
                   col("hr", "HR", "num"), col("var", "Var(X)", "num"), col("pevent", "P(Event)", "num"),
                   col("r2", "R²", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A sample of " + fint(r.num("n")) + " subjects with about " + fint(std::ceil(r.num("events") - 1e-9)) + " events " +
           achieves_singular(s, r) + " to detect a hazard ratio of " + fnum(r.num("hr")) +
           " for the covariate (variance " + fnum(r.num("var")) + ", R² with other covariates " + fnum(r.num("r2")) +
           ") using a " + alt_words(alt(s)) + " test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"Events", "Expected number of events: N × P(Event)."},
                   {"Var(X)", "Variance of the covariate (p(1 − p) for a binary covariate)."}};
  p.references = {"Hsieh, F. Y. & Lavori, P. W. (2000). Sample-size calculations for the Cox proportional hazards "
                  "regression model with nonbinary covariates. Controlled Clinical Trials 21, 552–560.",
                  "Schoenfeld, D. A. (1983). Sample-size formula for the proportional-hazards regression model. "
                  "Biometrics 39, 499–503."};
  p.notes = {"Hsieh and Lavori's formula is a large-sample approximation. For a binary covariate far from half "
             "exposed and a hazard ratio far from 1, it overstates the power when the larger group has the higher "
             "hazard (20% exposed, HR 0.6: 0.90, against 0.83 by simulation) and understates it in the opposite "
             "case (80% exposed: 0.94 by simulation). With fewer than about 30 events it is rough."};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.curve_min = 1;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"hr", "1.3 1.5 2"}, {"pevent", "0.5"}};
  return p;
}
SSC_REGISTER(cox);

}  // namespace
}  // namespace ssc::procs
