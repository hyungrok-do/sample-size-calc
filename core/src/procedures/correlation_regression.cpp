// Correlation (one and two correlations, CI precision) and regression
// (multiple regression R², simple linear regression slope, logistic and Cox
// regression with one covariate of interest).
#include <cstdint>
#include <cstring>

#include "procedures/apart.hpp"
#include "procedures/common.hpp"
#include "procedures/derive_models.hpp"
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

// The exact power's tails at ζ1 = atanh(ρ1): the chances of a sample correlation beyond each critical value (0 for
// an unused side), and the power they add up to.
std::pair<double, double> corr_tails_z(double n, double zeta1, const CorrCrit& c) {
  return {std::isfinite(c.up) ? stats::pcorr_z(c.up, n, zeta1, false) : 0.0,
          std::isfinite(c.lo) ? stats::pcorr_z(c.lo, n, zeta1, true) : 0.0};
}
double corr_tails_power(const std::pair<double, double>& tails) { return stats::clamp01(tails.first + tails.second); }

// The power at ζ1 = atanh(ρ1), at the distance t = ζ1 − atanh(ρ0) on Fisher's scale (the one the approximation
// takes, the other the exact distribution).
double corr_power_z(const Scenario& s, double n, double zeta1, double t, const CorrCrit& c) {
  const double alpha = s("alpha");
  const Alt a = alt(s);
  if (s.is("method", "fisher")) return z_power(t * std::sqrt(n - 3), alpha, a);
  return corr_tails_power(corr_tails_z(n, zeta1, c));
}

double corr_power(const Scenario& s, double n, double rho1, const CorrCrit& c) {
  const double zeta1 = std::atanh(rho1);
  return corr_power_z(s, n, zeta1, zeta1 - std::atanh(s("rho0")), c);
}

double corr_power(const Scenario& s, double n, double rho1) { return corr_power(s, n, rho1, corr_crit(s, n)); }

// The critical correlations of the exact test of one correlation with N subjects (written `N`): for ρ0 = 0 from the
// t-test, as zcrit_upper finds them, otherwise roots of the exact null distribution, which stats::pcorr_z evaluates
// on Fisher's scale (ζ0 = atanh ρ0, which the steps before have shown).
void exact_corr_critical_steps(const Scenario& s, double n, const std::string& N, double rho0, const CorrCrit& c,
                               derive::Steps& out) {
  using namespace derive;
  const double alpha = s("alpha"), zeta0 = std::atanh(rho0);
  const Alt a = alt(s);
  const bool hi = std::isfinite(c.up), lo = std::isfinite(c.lo);
  const std::string level = a == Alt::Two ? "α/2" : "α";
  const std::string join = hi && lo ? ", " : "";
  std::string f, t;
  if (rho0 == 0) {
    const double tc = stats::qt(tail_alpha(alpha, a), n - 2, false);
    out.push_back(t_critical(alpha, n - 2, a, "N − 2"));
    const std::string root = dv(tc) + " / √(" + N + " − 2 + " + dv(tc) + "²)";
    if (hi) {
      f = "r_U = t / √(N − 2 + t²)";
      t = "r_U = " + root + " = " + dv(std::tanh(c.up));
    }
    if (lo) {
      f += join + (hi ? "r_L = −r_U" : "r_L = −t / √(N − 2 + t²)");
      t += join + (hi ? "r_L = " : "r_L = −" + root + " = ") + dv(std::tanh(c.lo));
    }
    out.push_back({"Critical correlations", f + ", the t-test of ρ = 0 written for r", t});
    return;
  }
  if (hi) {
    f = "P(r ≥ r_U | ρ0) = " + level;
    t = "P(r ≥ " + dv(std::tanh(c.up)) + " | ρ0 = " + gv(rho0) + ") = " + dvp(stats::pcorr_z(c.up, n, zeta0, false));
  }
  if (lo) {
    f += join + "P(r ≤ r_L | ρ0) = " + level;
    t += join + "P(r ≤ " + dv(std::tanh(c.lo)) + " | ρ0 = " + gv(rho0) + ") = " +
         dvp(stats::pcorr_z(c.lo, n, zeta0, true));
  }
  out.push_back({"Critical correlations",
                 f + ", r the sample correlation of N subjects, from its exact distribution (on Fisher's scale ζ)", t});
}

// The exact power at ρ1 (written `rho1`) beyond the critical correlations, whose tails corr_tails_z gave: in symbols,
// and its terms up to the total.
std::pair<std::string, std::string> exact_corr_power_terms(const CorrCrit& c, const std::string& rho1,
                                                           const std::pair<double, double>& tails) {
  using namespace derive;
  const bool hi = std::isfinite(c.up), lo = std::isfinite(c.lo);
  const std::string plus = hi && lo ? " + " : "", at = " | ρ1 = " + rho1 + ")";
  std::string f, t;
  if (hi) {
    f = "P(r ≥ r_U | ρ1)";
    t = "P(r ≥ " + dv(std::tanh(c.up)) + at;
  }
  if (lo) {
    f += plus + "P(r ≤ r_L | ρ1)";
    t += plus + "P(r ≤ " + dv(std::tanh(c.lo)) + at;
  }
  if (hi && lo) t += " = " + dvp(tails.first) + " + " + dvp(tails.second);
  return {f + ", from the exact distribution of r", t};
}

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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const double n = r.num("n"), rho0 = r.num("rho0"), rho1 = r.num("rho1"), alpha = s("alpha");
    const Alt a = alt(s);
    const bool fisher = s.is("method", "fisher"), solved = s.solve == "n";
    const double zeta0 = std::atanh(rho0), zeta1 = std::atanh(rho1), ncp = (zeta1 - zeta0) * std::sqrt(n - 3);
    const std::string N = size(n, solved);
    // the exact power as corr_power computes it, with its critical values and tails kept for the steps
    const CorrCrit c = fisher ? CorrCrit{} : corr_crit(s, n);
    const auto tails = fisher ? std::pair<double, double>{} : corr_tails_z(n, zeta1, c);
    const double power = fisher ? corr_power(s, n, rho1, c) : corr_tails_power(tails);
    const std::string zeta = "ζ = atanh ρ = ½ ln((1 + ρ)/(1 − ρ))";
    const std::string zeta0_text = "ζ0 = atanh(" + gv(rho0) + ") = " + dv(zeta0);
    Derivation out;
    if (s.solve == "rho1") {
      const double target = s("power");
      const std::string rho1_text = derived(apart(rho1, rho0).first);
      out.find = find_effect("ρ1", target, "with N = " + gvi(n));
      out.steps.push_back({"Fisher z transform", zeta, zeta0_text});
      if (fisher) {  // the target power gives λ, λ the distance on Fisher's scale, and that ρ1
        out.steps.push_back(z_critical(alpha, a));
        out.steps.push_back(z_ncp_needed(ncp, alpha, a, target));
        out.steps.push_back({"Fisher z detected", "ζ1 = ζ0 + λ / √(N − 3)",
                             dv(zeta0) + " + " + dvs(ncp) + " / √(" + N + " − 3) = " + dv(zeta1)});
        out.steps.push_back({"Correlation detected", "ρ1 = tanh ζ1", "tanh(" + dv(zeta1) + ") = " + rho1_text});
      } else {  // the critical values follow from ρ0 alone; ρ1 is where the exact power beyond them reaches the target
        const double za = z_alpha(alpha, a), zb = stats::qnorm(target), sign = a == Alt::Lower ? -1 : 1;
        out.steps.push_back({"Fisher z approximation",
                             std::string("ρ1 ≈ tanh(ζ0 ") + (sign < 0 ? "−" : "+") + " (z(" + level_symbol(a) +
                                 ") + z(1 − β)) / √(N − 3)), for comparison; the steps below give the answer",
                             "tanh(" + dv(zeta0) + (sign < 0 ? " − " : " + ") + "(" + dv(za) + " + " + dv(zb) +
                                 ") / √(" + N + " − 3)) = " +
                                 dv(std::tanh(zeta0 + sign * (za + zb) / std::sqrt(n - 3)))});
        exact_corr_critical_steps(s, n, N, rho0, c, out.steps);
        const auto [f, t] = exact_corr_power_terms(c, rho1_text, tails);
        out.steps.push_back(effect_found("ρ1", rho1_text, f, t, power, target));
      }
      out.answer = "ρ1 = " + rho1_text;
      return out;
    }
    if (solved)
      out.find = find_size("N", s("power"));
    else
      out.find = find_power("with N = " + gvi(n));
    const std::string zeta1_text = "ζ1 = atanh(" + gv(rho1) + ") = " + dv(zeta1);
    // the exact test with ρ0 = 0 is a t-test, which needs no Fisher transform unless for the closed form
    if (fisher || solved || rho0 != 0)
      out.steps.push_back({"Fisher z transform", zeta, zeta0_text + ", " + zeta1_text});
    if (solved) {
      const double za = z_alpha(alpha, a), zb = stats::qnorm(s("power"));
      out.steps.push_back({"Fisher z approximation",
                           std::string("N ≈ ((z(") + level_symbol(a) +
                               ") + z(1 − β)) / (ζ1 − ζ0))² + 3, for comparison; the steps below give the answer",
                           "((" + dv(za) + " + " + dv(zb) + ") / (" + dv(zeta1) + " − " + dvs(zeta0) + "))² + 3 = " +
                               dv(std::pow((za + zb) / (zeta1 - zeta0), 2) + 3)});
    }
    if (fisher) {
      out.steps.push_back({"Noncentrality", "λ = (ζ1 − ζ0) √(N − 3)",
                           "(" + dv(zeta1) + " − " + dvs(zeta0) + ") × √(" + N + " − 3) = " + dv(ncp)});
      out.steps.push_back(z_critical(alpha, a));
      out.steps.push_back(z_power_step(ncp, alpha, a));
    } else {
      exact_corr_critical_steps(s, n, N, rho0, c, out.steps);
      const auto [f, t] = exact_corr_power_terms(c, gv(rho1), tails);
      out.steps.push_back({"Power", f, t + " = " + dvp(power)});
    }
    if (solved) {
      out.steps.push_back(
          n_search([&](double x) { return x == n ? power : corr_power(s, x, rho1); }, n, s("power"), 4));
      out.answer = size_answer(n, power);
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
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

// The power of Fisher's z-test of two correlations whose transforms differ by dz, with group sizes n1 and n2.
double two_corr_power(double dz, double n1, double n2, double alpha, Alt a) {
  return z_power(dz / std::sqrt(1 / (n1 - 3) + 1 / (n2 - 3)), alpha, a);
}

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
    auto P = [&](double n1, double n2) { return two_corr_power(dz, n1, n2, alpha, a); };
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const TwoN n{r.num("n1"), r.num("n2")};
    const double rho1 = r.num("rho1"), rho2 = r.num("rho2"), alpha = s("alpha");
    const Alt a = alt(s);
    const bool solved = s.solve == "n";
    const double zeta1 = std::atanh(rho1), zeta2 = std::atanh(rho2), dz = zeta1 - zeta2;
    const double se = std::sqrt(1 / (n.n1 - 3) + 1 / (n.n2 - 3)), ncp = dz / se;
    const double power = two_corr_power(dz, n.n1, n.n2, alpha, a);
    Derivation out;
    out.steps.push_back({"Fisher z transform", "ζ = atanh ρ = ½ ln((1 + ρ)/(1 − ρ))",
                         "ζ1 = atanh(" + gv(rho1) + ") = " + dv(zeta1) + ", ζ2 = atanh(" + gv(rho2) + ") = " +
                             dv(zeta2)});
    if (solved) {
      out.find = find_size(two_group_unknown(s), s("power"));
      // the closed form where group 2 follows group 1: 1/(N1 − 3) + 1/(N2 − 3) ≈ (1 + 1/r)/(N1 − 3)
      const std::string& alloc = s.s("alloc");
      if (alloc == "equal" || alloc == "ratio") {
        const double k = alloc == "ratio" ? 1 + 1 / s("r") : 2;
        const double za = z_alpha(alpha, a), zb = stats::qnorm(s("power"));
        out.steps.push_back(
            {"Fisher z approximation",
             std::string(alloc == "equal" ? "N1 = N2 ≈ 2" : "N1 ≈ (1 + 1/r)") + " ((z(" + level_symbol(a) +
                 ") + z(1 − β)) / (ζ1 − ζ2))² + 3, for comparison; the steps below give the answer",
             std::string(alloc == "equal" ? "2" : "(1 + 1/" + gv(s("r")) + ")") + " × ((" + dv(za) + " + " + dv(zb) +
                 ") / (" + dv(zeta1) + " − " + dvs(zeta2) + "))² + 3 = " + dv(k * std::pow((za + zb) / dz, 2) + 3)});
      }
    } else {
      out.find = find_power("with " + two_sizes(n, false));
    }
    out.steps.push_back({"Standard error", "SE = √(1/(N1 − 3) + 1/(N2 − 3))",
                         "√(1/(" + size(n.n1, solved) + " − 3) + 1/(" + size(n.n2, solved) + " − 3)) = " + dv(se)});
    out.steps.push_back({"Noncentrality", "λ = (ζ1 − ζ2) / SE",
                         "(" + dv(zeta1) + " − " + dvs(zeta2) + ") / " + dv(se) + " = " + dv(ncp)});
    out.steps.push_back(z_critical(alpha, a));
    out.steps.push_back(z_power_step(ncp, alpha, a));
    if (solved) {
      out.steps.push_back(two_group_search(
          s, [&](double x1, double x2) { return two_corr_power(dz, x1, x2, alpha, a); }, n, s("power"), 4));
      out.answer = sizes_answer(n, "power " + dvp(power));
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
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

// The width of the Fisher-z interval around the sample correlation rho with n subjects, z the normal quantile of
// the confidence level.
double corr_ci_width(double rho, double z, double n) {
  const double h = z / std::sqrt(n - 3);
  return std::tanh(std::atanh(rho) + h) - std::tanh(std::atanh(rho) - h);
}

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
    auto width = [&](double n) { return corr_ci_width(rho, z, n); };
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
    // the bounds told apart from each other and from the estimate (0.9999999938 next to 0.99999999, not as it)
    const auto t = fnum_apart_list({r.num("lower"), r.num("rho"), r.num("upper")}, {-1.0, 1.0});
    return "A sample of " + fint(r.num("n")) + " gives a " + level_text(r.num("conf")) +
           " Fisher-z confidence interval of width " + fnum(r.num("width")) + " (" + t[0] + " to " + t[2] +
           ") when the sample correlation is " + t[1] + ".";
  };
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const double n = r.num("n"), rho = r.num("rho"), conf = r.num("conf"), width = r.num("width");
    const double z = stats::qnorm((1 - conf) / 2, false), zeta = std::atanh(rho), h = z / std::sqrt(n - 3);
    const bool solved = s.solve == "n";
    const std::string N = size(n, solved), level = "(level C = " + gv(conf) + ")";
    // the limits told apart from each other and from the estimate, as the statement shows them
    const auto t = fnum_apart_list({r.num("lower"), rho, r.num("upper")}, {-1.0, 1.0});
    const auto printed = [](const std::string& x) { return !x.empty() && x[0] == '-' ? "−" + x.substr(1) : x; };
    const std::string lower = derived(printed(t[0])), upper = derived(printed(t[2]));
    Derivation out;
    out.steps.push_back(
        {"Fisher z transform", "ζ = atanh r = ½ ln((1 + r)/(1 − r))", "atanh(" + gv(rho) + ") = " + dv(zeta)});
    out.steps.push_back({"Critical value", "z = Φ⁻¹(1 − (1 − C)/2), C the confidence level",
                         "Φ⁻¹(1 − (1 − " + gv(conf) + ")/2) = " + dv(z)});
    if (solved) {
      // Bonett and Wright's closed form from the interval's first-order width 2z(1 − r²)/√(N − 3)
      const double target = r.num("target_width");
      out.find = "the smallest N whose confidence interval " + level + " is at most W = " + gv(target) + " wide";
      out.steps.push_back({"Approximation",
                           "N ≈ 4 z² (1 − r²)² / W² + 3, for comparison; the steps below give the answer",
                           "4 × " + dv(z) + "² × (1 − " + gvs(rho) + "²)² / " + gv(target) + "² + 3 = " +
                               dv(4 * z * z * std::pow(1 - rho * rho, 2) / (target * target) + 3)});
    } else {
      out.find = "the width of the confidence interval " + level + " with N = " + gvi(n);
    }
    out.steps.push_back(
        {"Half-width on Fisher's scale", "h = z / √(N − 3)", dv(z) + " / √(" + N + " − 3) = " + dv(h)});
    out.steps.push_back({"Limits", "L = tanh(ζ − h), U = tanh(ζ + h)",
                         "L = tanh(" + dv(zeta) + " − " + dv(h) + ") = " + lower + ", U = tanh(" + dv(zeta) + " + " +
                             dv(h) + ") = " + upper});
    out.steps.push_back({"Width", "W = U − L", upper + " − " + (r.num("lower") < 0 ? "(" + lower + ")" : lower) +
                                                   " = " + dv(width)});
    if (solved) {
      // as n_search, for the width: the smallest size whose interval is no wider than the target
      const double target = r.num("target_width");
      if (n - 1 < 4)
        out.steps.push_back({"Smallest size", "N = the smallest size whose width is at most W",
                             "N = " + dvi(n) + " is the smallest size allowed, and its width is already " +
                                 dvp(width, target) + " ≤ " + gv(target)});
      else
        out.steps.push_back({"Smallest size", "N = the smallest whole number whose width is at most W",
                             "width at N = " + dvi(n - 1) + ": " + dvp(corr_ci_width(rho, z, n - 1), target) + " > " +
                                 gv(target) + " ≥ width at N = " + dvi(n) + ": " + dvp(width, target)});
      out.answer = "N = " + dvi(n) + ", with width " + dvp(width, target) + " (from " + lower + " to " + upper + ")";
    } else {
      out.answer = "W = " + dv(width) + " (from " + lower + " to " + upper + ")";
    }
    return out;
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

// Cohen's f² = ΔR²/(1 − R²c − ΔR²), with the denominator to one rounding: 1 − R²c is not exact in doubles (R²c = 0.3),
// and next to the last ΔR² its rounding is a large part of the denominator (f² 17% off at the last double)
double mreg_f2(double d, double r2c) {
  const double rest = 1 - r2c, rest_lo = -r2c - (rest - 1);  // 1 − R²c = rest + rest_lo exactly (Dekker)
  return d / ((rest - d) + rest_lo);
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
          .range(0, 1, false, true)
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
    // 1 − R²c = rest + rest_lo exactly (Dekker), as mreg_f2 takes it
    const double rest = 1 - r2c, rest_lo = -r2c - (rest - 1);
    auto f2_of = [&](double d) { return mreg_f2(d, r2c); };
    if (s.solve == "r2t") {
      n = s("n");
      if (n < a + b + 2) throw CalcError("N must exceed A + B + 1");
      // solve on the unbounded f² scale, then ΔR² = f²(1 − R²c)/(1 + f²), up to the f² of the last ΔR² that may
      // be given (R²c + ΔR² < 1).  It is found by halving over the doubles below 1 − R²c, as one double at a time
      // would take some 1e13 steps next to R²c = 1 − 1e-14, where R²c + ΔR² rounds as 1 does
      auto admissible = [&](double d) { return r2c + d < 1 && (rest - d) + rest_lo > 0; };
      double last = rest;
      if (!admissible(last)) {
        std::uint64_t lo = 0, hi;  // as bit patterns, which order positive doubles as their values
        std::memcpy(&hi, &rest, sizeof hi);
        while (hi - lo > 1) {
          const std::uint64_t mid = lo + (hi - lo) / 2;
          double v;
          std::memcpy(&v, &mid, sizeof v);
          (admissible(v) ? lo : hi) = mid;
        }
        std::memcpy(&last, &lo, sizeof last);
      }
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
        // with no effect the power is alpha, which no sample size raises
        if (r2t == 0) throw CalcError("ΔR² = 0: there is no effect to detect");
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
    // ΔR² and R²c with the digits it takes to tell each from 1 less the other, their limit (an increase of
    // 0.699999999999886 with R²c = 0.3, not 0.7)
    const double r2t = r.num("r2t"), r2c = r.num("r2c");
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) + " to detect an R² increase of " +
           fnum_apart(r2t, 1 - r2c).first + " (f² = " + fnum(r.num("f2")) + ") attributable to " + fint(r.num("b")) +
           " tested predictor(s), adjusting for " + fint(r.num("a")) + " control predictor(s) with R² = " +
           fnum_apart(r2c, 1 - r2t).first + ", using an F-test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const double n = r.num("n"), b = r.num("b"), a = r.num("a"), r2t = r.num("r2t"), r2c = r.num("r2c");
    const double alpha = s("alpha"), f2 = mreg_f2(r2t, r2c), df2 = n - a - b - 1, lambda = n * f2;
    const double fc = stats::qf(alpha, b, df2, false), power = mreg_power(n, f2, b, a, alpha);
    const bool solved = s.solve == "n";
    const std::string N = size(n, solved);
    const Step critical{"Critical value",
                        "F_c = F(1 − α; B, ν), the F quantile with B and ν = N − A − B − 1 degrees of freedom",
                        "ν = " + N + " − " + gvi(a) + " − " + gvi(b) + " − 1 = " + dvi(df2) + ", F(1 − " + gv(alpha) +
                            "; " + gvi(b) + ", " + dvi(df2) + ") = " + dv(fc)};
    Derivation out;
    if (s.solve == "r2t") {  // the target power gives λ, λ the effect size f², and f² the R² increase
      const double target = s("power");
      out.find = find_effect("ΔR²", target, "with N = " + gvi(n));
      out.steps.push_back(critical);
      out.steps.push_back({"Noncentrality needed", "λ solves P(F′ > F_c) = 1 − β, F′ ~ noncentral F(B, ν, λ)",
                           "P(F′ > " + dv(fc) + ") = " + gv(target) + " when λ = " + dv(lambda)});
      out.steps.push_back(
          {"Effect size detected (Cohen's f²)", "f² = λ / N", dv(lambda) + " / " + N + " = " + dsig(f2)});
      // ΔR² with the digits that tell it from 1 − R²c, its limit, as the statement shows it
      const std::string delta = derived(fnum_apart(r2t, 1 - r2c).first);
      out.steps.push_back({"R² increase detected", "ΔR² = f² (1 − R²c) / (1 + f²), from f² = ΔR² / (1 − R²c − ΔR²)",
                           dsig(f2) + " × (1 − " + gsig(r2c) + ") / (1 + " + dsig(f2) + ") = " + delta});
      out.answer = "ΔR² = " + delta + " (f² = " + dsig(f2) + ")";
      return out;
    }
    out.find = solved ? find_size("N", s("power")) : find_power("with N = " + gvi(n));
    out.steps.push_back({"Effect size (Cohen's f²)", "f² = ΔR² / (1 − R²c − ΔR²)",
                         gsig(r2t) + " / (1 − " + gsig(r2c) + " − " + gsig(r2t) + ") = " + dsig(f2)});
    out.steps.push_back({"Noncentrality", "λ = N f²", N + " × " + dsig(f2) + " = " + dv(lambda)});
    out.steps.push_back(critical);
    out.steps.push_back(
        {"Power", "P(F′ > F_c), F′ ~ noncentral F(B, ν, λ)", "P(F′ > " + dv(fc) + ") = " + dvp(power)});
    if (solved) {
      out.steps.push_back(n_search([&](double x) { return x == n ? power : mreg_power(x, f2, b, a, alpha); }, n,
                                   s("power"), a + b + 2));
      out.answer = size_answer(n, power);
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
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
// The residual SD of Y given X at the slope b1: given, or derived from the total SD of Y.
double slope_resid_sd(const Scenario& s, double b1) {
  if (s.is("sdtype", "sy")) {
    const double sx = s("sx");
    const double v = s("sigma") * s("sigma") - b1 * b1 * sx * sx;
    if (!(v > 0)) throw CalcError("β1·σx must be smaller than σy");
    return std::sqrt(v);
  }
  return s("sigma");
}

// The power with n subjects at the slope b1, from the difference d = β1 − B0, which the search for β1 keeps apart
// from a large B0 (which would round it away).
double slope_power(const Scenario& s, double n, double d, double b1) {
  const double ncp = d * s("sx") * std::sqrt(n - 1) / slope_resid_sd(s, b1);
  return t_power(ncp, n - 2, s("alpha"), alt(s));
}

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
    auto resid_sd = [&](double b1) { return slope_resid_sd(s, b1); };
    auto Pd = [&](double n, double d, double b1) { return slope_power(s, n, d, b1); };
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const double n = r.num("n"), b0 = r.num("b0"), b1 = r.num("b1"), sx = r.num("sx"), alpha = s("alpha");
    const Alt a = alt(s);
    const bool total_sd = s.is("sdtype", "sy"), solved = s.solve == "n";
    const double d = b1 - b0, sd = slope_resid_sd(s, b1), ncp = d * sx * std::sqrt(n - 1) / sd, df = n - 2;
    const double power = slope_power(s, n, d, b1);
    const std::string N = size(n, solved), SX = gsig(sx), SD = total_sd ? dsig(sd) : gsig(sd);
    // the residual SD at β1 (written `b1_text`), from the total SD of Y
    const auto resid = [&](const std::string& b1_text) {
      return "σ = √(" + gsig(s("sigma")) + "² − " + b1_text + "² × " + SX + "²) = " + SD;
    };
    Derivation out;
    if (s.solve == "b1") {  // the target power gives λ, λ the difference, the difference β1
      const double target = s("power");
      const auto texts = apart(b1, b0);
      const std::string b1_text = derived(texts.first), b0_text = given(texts.second);
      out.find = find_effect("β1", target, "with N = " + gvi(n));
      out.steps.push_back(t_critical(alpha, df, a, "N − 2"));
      out.steps.push_back(t_ncp_needed(ncp, df, alpha, a, target));
      const std::string delta = dvs(ncp) + " × " + SD + " / (" + SX + " × √(" + N + " − 1)) = " + dv(d);
      const std::string slope = b0_text + " + " + dvs(d) + " = " + b1_text;
      if (!total_sd) {
        out.steps.push_back({"Difference detected", "δ = λ σ / (σx √(N − 1))", delta});
        out.steps.push_back({"Slope detected", "β1 = B0 + δ", slope});
      } else {  // σ depends on β1, so β1 solves the relation with σ at β1
        out.steps.push_back(
            {"Slope detected",
             "β1 = B0 + δ, δ = λ σ / (σx √(N − 1)), σ = √(σy² − β1² σx²) the residual SD at β1: σ changes with β1, "
             "so β1 is found by a root search",
             "at β1 = " + b1_text + ": " + resid(b1_text) + ", δ = " + delta + ", β1 = " + slope});
        out.steps.push_back(effect_found("β1", b1_text, t_power_step(ncp, df, alpha, a).formula, "", power, target));
      }
      out.answer = "β1 = " + b1_text + " (δ = " + dv(d) + ")";
      return out;
    }
    out.find = solved ? find_size("N", s("power")) : find_power("with N = " + gvi(n));
    if (total_sd) out.steps.push_back({"Residual SD", "σ = √(σy² − β1² σx²)", resid(gvs(b1))});
    out.steps.push_back({"Difference to detect", "δ = β1 − B0", gv(b1) + " − " + gvs(b0) + " = " + dv(d)});
    if (solved) {
      // the closed form a normal statistic would give, beside the exact answer below
      const double za = z_alpha(alpha, a), zb = stats::qnorm(s("power"));
      out.steps.push_back(
          {"Normal approximation",
           "N ≈ " + z_sum_formula(a) + " σ² / (δ σx)² + 1, for comparison; the steps below give the answer",
           z_sum_numbers(alpha, a, s("power")) + " × " + SD + "² / (" + dvs(d) + " × " + SX + ")² + 1 = " +
               dv((za + zb) * (za + zb) * sd * sd / (d * d * sx * sx) + 1)});
    }
    out.steps.push_back({"Noncentrality", "λ = δ σx √(N − 1) / σ",
                         dvs(d) + " × " + SX + " × √(" + N + " − 1) / " + SD + " = " + dv(ncp)});
    out.steps.push_back(t_critical(alpha, df, a, "N − 2"));
    out.steps.push_back(t_power_step(ncp, df, alpha, a));
    if (solved) {
      out.steps.push_back(n_search([&](double x) { return slope_power(s, x, d, b1); }, n, s("power"), 3));
      out.answer = size_answer(n, power);
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
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
// The power of the Wald test of the covariate: for a binary X the normal approximation with the SDs under H0 and H1
// (sd0, sd1), for a normal X that of a statistic N(shift, 1).
struct Logistic {
  bool binary;
  double p0, orr, r2, za;
  double other;  // a two-sided test also rejects on the side away from the effect: without that tail the power of
                 // no effect would be alpha/2
  double B = NAN, o1 = NAN, p1 = NAN, pb = NAN, sd0 = NAN, crit = NAN, sd1 = NAN;  // binary X
  double b = NAN;                                                                   // normal X: ln OR per SD
  explicit Logistic(const Scenario& s)
      : binary(s.is("xtype", "binary")), p0(s("p0")), orr(s("orr")), r2(s("r2")),
        za(stats::qnorm(tail_alpha(s("alpha"), alt(s)), false)), other(alt(s) == Alt::Two ? 1 : 0) {
    if (binary) {
      B = s("px");
      o1 = orr * p0 / (1 - p0);
      p1 = o1 / (1 + o1);
      pb = (1 - B) * p0 + B * p1;
      sd0 = std::sqrt(pb * (1 - pb) / B);
      crit = za * sd0;
      sd1 = std::sqrt(p0 * (1 - p0) + p1 * (1 - p1) * (1 - B) / B);
    } else {
      b = std::log(orr);
    }
  }
  double shift(double n) const {
    return binary ? std::sqrt(n * (1 - r2) * (1 - B)) * std::fabs(p1 - p0)
                  : std::sqrt(n * (1 - r2) * p0 * (1 - p0)) * std::fabs(b);
  }
  double power(double n) const {
    const double d = shift(n);
    if (binary) return stats::pnorm((d - crit) / sd1) + other * stats::pnorm((-d - crit) / sd1);
    return stats::pnorm(d - za) + other * stats::pnorm(-d - za);
  }
};

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
    const double p0 = s("p0"), orr = s("orr"), r2 = s("r2"), alpha = s("alpha");
    const Alt a = alt(s);
    // with no effect the power is alpha, which no sample size raises
    if (orr == 1 && s.solve == "n") throw CalcError("an odds ratio of 1 has no effect to detect");
    if ((a == Alt::Upper && orr < 1) || (a == Alt::Lower && orr > 1))
      throw CalcError("the odds ratio is on the side opposite to the one-sided alternative");
    const Logistic m(s);
    auto P = [&](double n) { return m.power(n); };
    const double p1 = m.p1;
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const Logistic m(s);
    const double n = r.num("n"), alpha = s("alpha"), d = m.shift(n), power = m.power(n);
    const Alt a = alt(s);
    const bool solved = s.solve == "n";
    const std::string N = size(n, solved), P0 = gsig(m.p0), R2 = gsig(m.r2);
    Derivation out;
    out.find = solved ? find_size("N", s("power")) : find_power("with N = " + gvi(n));
    // Hsieh's closed form: the N at which the shift reaches z σ0 + z(1 − β) σ1 (σ0 = σ1 = 1 for a normal X)
    const double za = m.za, zb = solved ? stats::qnorm(s("power")) : NAN;
    if (m.binary) {
      const std::string B = gsig(m.B);
      out.steps.push_back({"Event probability at X = 1", "P1 = O1 / (1 + O1), O1 = OR × P0 / (1 − P0) its odds",
                           "O1 = " + gv(m.orr) + " × " + P0 + " / (1 − " + P0 + ") = " + dsig(m.o1) + ", P1 = " +
                               dsig(m.o1) + " / (1 + " + dsig(m.o1) + ") = " + dsig(m.p1)});
      out.steps.push_back({"Overall event probability", "P̄ = (1 − B) P0 + B P1, B = P(X = 1)",
                           "(1 − " + B + ") × " + P0 + " + " + B + " × " + dsig(m.p1) + " = " + dsig(m.pb)});
      out.steps.push_back(z_critical(alpha, a));
      out.steps.push_back(
          {"Standard deviations",
           "σ0 = √(P̄ (1 − P̄) / B) under H0, σ1 = √(P0 (1 − P0) + P1 (1 − P1)(1 − B) / B) under H1",
           "σ0 = √(" + dsig(m.pb) + " × (1 − " + dsig(m.pb) + ") / " + B + ") = " + dv(m.sd0) + ", σ1 = √(" + P0 +
               " × (1 − " + P0 + ") + " + dsig(m.p1) + " × (1 − " + dsig(m.p1) + ") × (1 − " + B + ") / " + B +
               ") = " + dv(m.sd1)});
      if (solved)
        out.steps.push_back(
            {"Hsieh's formula", std::string("N ≈ (z(") + level_symbol(a) +
                                    ") σ0 + z(1 − β) σ1)² / ((1 − R²)(1 − B)(P1 − P0)²), for comparison; the steps "
                                    "below give the answer",
             "(" + dv(za) + " × " + dv(m.sd0) + " + " + dv(zb) + " × " + dv(m.sd1) + ")² / ((1 − " + R2 + ") × (1 − " +
                 B + ") × (" + dsig(m.p1) + " − " + P0 + ")²) = " +
                 dv(std::pow(za * m.sd0 + zb * m.sd1, 2) / ((1 - m.r2) * (1 - m.B) * std::pow(m.p1 - m.p0, 2)))});
      out.steps.push_back({"Shift", "δ = √(N (1 − R²)(1 − B)) |P1 − P0|",
                           "√(" + N + " × (1 − " + R2 + ") × (1 − " + B + ")) × |" + dsig(m.p1) + " − " + P0 + "| = " +
                               dv(d)});
      auto term = [&](double x) {
        return "Φ((" + dv(x) + " − " + dv(m.za) + " × " + dv(m.sd0) + ") / " + dv(m.sd1) + ")";
      };
      if (a == Alt::Two)
        out.steps.push_back({"Power", "Φ((δ − z σ0) / σ1) + Φ((−δ − z σ0) / σ1)",
                             term(d) + " + " + term(-d) + " = " + dvp(stats::pnorm((d - m.crit) / m.sd1)) + " + " +
                                 dvp(stats::pnorm((-d - m.crit) / m.sd1)) + " = " + dvp(power)});
      else
        out.steps.push_back(
            {"Power", "Φ((δ − z σ0) / σ1), on the side of the effect", term(d) + " = " + dvp(power)});
    } else {
      out.steps.push_back({"Log odds ratio per SD of X", "β = ln OR", "ln " + gv(m.orr) + " = " + dv(m.b)});
      out.steps.push_back(z_critical(alpha, a));
      if (solved)
        out.steps.push_back(
            {"Hsieh's formula",
             "N ≈ " + z_sum_formula(a) +
                 " / ((1 − R²) P0 (1 − P0) β²), for comparison; the steps below give the answer",
             z_sum_numbers(alpha, a, s("power")) + " / ((1 − " + R2 + ") × " + P0 + " × (1 − " + P0 + ") × " +
                 dvs(m.b) + "²) = " + dv((za + zb) * (za + zb) / ((1 - m.r2) * m.p0 * (1 - m.p0) * m.b * m.b))});
      out.steps.push_back({"Noncentrality", "λ = √(N (1 − R²) P0 (1 − P0)) |β|",
                           "√(" + N + " × (1 − " + R2 + ") × " + P0 + " × (1 − " + P0 + ")) × |" + dv(m.b) + "| = " +
                               dv(d)});
      // λ is taken on the side of the effect, so that a one-sided test's power is Φ(λ − z) whichever its side
      out.steps.push_back(z_power_step(d, alpha, a == Alt::Two ? Alt::Two : Alt::Upper));
    }
    if (solved) {
      out.steps.push_back(n_search([&](double x) { return m.power(x); }, n, s("power"), 2));
      out.answer = size_answer(n, power);
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
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
// The variance of the covariate: p(1 − p) for a binary one, σx² for a continuous one.
double cox_var(const Scenario& s) { return s.is("xtype", "binary") ? s("px") * (1 - s("px")) : s("sx") * s("sx"); }

// The mean of the test statistic per √event (Hsieh and Lavori): √((1 − R²) Var(X)) ln HR.
double cox_drift_per_event(double r2, double var, double hr) { return std::sqrt((1 - r2) * var) * std::log(hr); }

// The power with n subjects of whom a share pe has the event, the statistic's mean being drift_per_event per √event.
double cox_power(double n, double pe, double drift_per_event, double alpha, Alt a) {
  return z_power(std::sqrt(n * pe) * drift_per_event, alpha, a);
}

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
    // with no effect the power is alpha, which no sample size raises
    if (hr == 1 && s.solve == "n") throw CalcError("a hazard ratio of 1 has no effect to detect");
    const double var = cox_var(s);
    const double drift_per_event = cox_drift_per_event(r2, var, hr);
    auto P = [&](double n) { return cox_power(n, pe, drift_per_event, alpha, a); };
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const double n = r.num("n"), hr = r.num("hr"), pe = r.num("pevent"), r2 = r.num("r2"), alpha = s("alpha");
    const Alt a = alt(s);
    const bool solved = s.solve == "n";
    const double var = cox_var(s), drift_per_event = cox_drift_per_event(r2, var, hr);
    const double events = n * pe, ncp = std::sqrt(events) * drift_per_event;
    const double power = cox_power(n, pe, drift_per_event, alpha, a);
    const std::string PE = gsig(pe), R2 = gsig(r2);
    Derivation out;
    if (s.is("xtype", "binary"))
      out.steps.push_back({"Variance of X", "Var(X) = p (1 − p), p = P(X = 1)",
                           gsig(s("px")) + " × (1 − " + gsig(s("px")) + ") = " + dsig(var)});
    else
      out.steps.push_back({"Variance of X", "Var(X) = σx²", gsig(s("sx")) + "² = " + dsig(var)});
    if (solved) {
      // Hsieh and Lavori's closed form for the events, and the subjects they take
      out.find = find_size("N", s("power"));
      const double need = z_factor(alpha, s("power"), a) / (drift_per_event * drift_per_event);
      out.steps.push_back({"Events needed (Hsieh–Lavori)",
                           "D ≈ " + z_sum_formula(a) +
                               " / ((1 − R²) Var(X) (ln HR)²) and N ≈ D / P(Event), for comparison; the steps below "
                               "give the answer",
                           "D ≈ " + z_sum_numbers(alpha, a, s("power")) + " / ((1 − " + R2 + ") × " + dsig(var) +
                               " × (ln " + gv(hr) + ")²) = " + dv(need, 2) + ", N ≈ " + dv(need, 2) + " / " + PE +
                               " = " + dv(need / pe, 2)});
    } else {
      out.find = find_power("with N = " + gvi(n));
    }
    out.steps.push_back(
        {"Expected events", "D = N × P(Event)", size(n, solved) + " × " + PE + " = " + dv(events, 2)});
    out.steps.push_back({"Noncentrality (Hsieh–Lavori)", "λ = √D × √((1 − R²) Var(X)) × ln HR",
                         "√" + dv(events, 2) + " × √((1 − " + R2 + ") × " + dsig(var) + ") × ln " + gv(hr) + " = " +
                             dv(ncp)});
    out.steps.push_back(z_critical(alpha, a));
    out.steps.push_back(z_power_step(ncp, alpha, a));
    if (solved) {
      out.steps.push_back(
          n_search([&](double x) { return cox_power(x, pe, drift_per_event, alpha, a); }, n, s("power"), 1));
      out.answer = "N = " + dvi(n) + ", with " + dv(events, 2) + " expected events and power " + dvp(power);
    } else {
      out.answer = "power = " + dvp(power) + ", with " + dv(events, 2) + " expected events";
    }
    return out;
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
