// Shared steps for the derivations that show how a result came about. Each procedure plugs its own steps together
// (p.derivation) from these and its own quantities, computing every value with the same functions as its run, so
// that the steps and the answer cannot disagree.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "procedures/common.hpp"

namespace ssc::procs::derive {

using blocks::fint;
using blocks::fnum;
using Steps = std::vector<Step>;

/// A number as written in a formula: the minus sign of print, not a hyphen.
inline std::string num(double x, int decimals = 4) {
  std::string t = fnum(x, decimals);
  if (!t.empty() && t[0] == '-') t = "−" + t.substr(1);
  return t;
}
/// A number substituted into a formula: a negative one in parentheses, "0.5 − (−0.2)".
inline std::string sub(double x, int decimals = 4) {
  const std::string t = num(x, decimals);
  return x < 0 ? "(" + t + ")" : t;
}
/// A probability with the decimals it takes to tell it apart from the target it is compared with (0.79995 against
/// 0.8 would read 0.8 at four decimals).
inline std::string prob(double p, double against = NAN) {
  std::string t = fnum(p, 4);
  if (std::isfinite(against) && p != against)
    for (int d = 5; d <= 12 && t == fnum(against, d - 1); ++d) t = fnum(p, d);
  return t;
}

inline const char* level_symbol(Alt a) { return a == Alt::Two ? "1 − α/2" : "1 − α"; }

/// The critical value of a test whose statistic is standard normal under H0.
inline Step z_critical(double alpha, Alt a) {
  const double tail = tail_alpha(alpha, a);
  return {"Critical value", std::string("z = Φ⁻¹(") + level_symbol(a) + ")",
          "Φ⁻¹(1 − " + num(tail) + ") = " + num(stats::qnorm(tail, false))};
}

/// The critical value of a t statistic with `df` degrees of freedom (`df_formula` says how they arise).
inline Step t_critical(double alpha, double df, Alt a, const std::string& df_formula) {
  const double tail = tail_alpha(alpha, a);
  return {"Critical value",
          std::string("t = t(") + level_symbol(a) + "; ν), the t quantile with ν = " + df_formula +
              " degrees of freedom",
          "t(1 − " + num(tail) + "; " + num(df) + ") = " + num(stats::qt(tail, df, false))};
}

/// The power of a test whose statistic is N(λ, 1) under H1 (as z_power).
inline Step z_power_step(double ncp, double alpha, Alt a) {
  const double z = stats::qnorm(tail_alpha(alpha, a), false), p = z_power(ncp, alpha, a);
  if (a == Alt::Two)
    return {"Power", "Φ(λ − z) + Φ(−λ − z)",
            "Φ(" + num(ncp) + " − " + num(z) + ") + Φ(" + num(-ncp) + " − " + num(z) + ") = " +
                prob(stats::pnorm(ncp - z)) + " + " + prob(stats::pnorm(-ncp - z)) + " = " + prob(p)};
  if (a == Alt::Upper) return {"Power", "Φ(λ − z)", "Φ(" + num(ncp) + " − " + num(z) + ") = " + prob(p)};
  return {"Power", "Φ(−λ − z)", "Φ(" + num(-ncp) + " − " + num(z) + ") = " + prob(p)};
}

/// The power of a t-test (as t_power): T′ has the noncentral t distribution with ν degrees of freedom and
/// noncentrality λ.
inline Step t_power_step(double ncp, double df, double alpha, Alt a) {
  const double t = stats::qt(tail_alpha(alpha, a), df, false), p = t_power(ncp, df, alpha, a);
  const std::string law = ", T′ ~ noncentral t(ν, λ)";
  if (a == Alt::Two)
    return {"Power", "P(T′ > t) + P(T′ < −t)" + law,
            "P(T′ > " + num(t) + ") + P(T′ < −" + num(t) + ") = " + prob(stats::pnt(t, df, ncp, false)) + " + " +
                prob(stats::pnt(-t, df, ncp, true)) + " = " + prob(p)};
  if (a == Alt::Upper) return {"Power", "P(T′ > t)" + law, "P(T′ > " + num(t) + ") = " + prob(p)};
  return {"Power", "P(T′ < −t)" + law, "P(T′ < −" + num(t) + ") = " + prob(p)};
}

/// How a single sample size was found: the smallest admissible size (from `lo`, in steps of `step`) whose power
/// reaches the target, shown by the power one step below it and at it.
inline Step n_search(const std::function<double(double)>& power, double n, double target, double lo,
                     const std::string& sym = "N", double step = 1) {
  const double at = power(n);
  if (n - step < lo)
    return {"Sample size", sym + " = the smallest size whose power reaches the target",
            sym + " = " + fint(n) + " is the smallest size allowed, and its power is already " + prob(at, target) +
                " ≥ " + fnum(target)};
  return {"Sample size", sym + " = the smallest whole number whose power reaches the target",
          "power at " + sym + " = " + fint(n - step) + ": " + prob(power(n - step), target) + " < " + fnum(target) +
              " ≤ power at " + sym + " = " + fint(n) + ": " + prob(at, target)};
}

/// How two group sizes were found under the scenario's allocation rule (see blocks::smaller_two_group).
inline Step two_group_search(const Scenario& s, const std::function<double(double, double)>& power,
                             const blocks::TwoN& n, double target, double min_n) {
  const std::string& a = s.s("alloc");
  const std::string rule =
      a == "equal"   ? "N1 = N2, the smallest whole number whose power reaches the target"
      : a == "ratio" ? "N1 = the smallest whole number whose power reaches the target, with N2 = ⌈r × N1⌉ (r = " +
                           fnum(s("r")) + ")"
      : a == "n1"    ? "N2 = the smallest whole number whose power reaches the target, with N1 fixed"
      : a == "n2"    ? "N1 = the smallest whole number whose power reaches the target, with N2 fixed"
                     : "N = the smallest total whose power reaches the target, " + fnum(s("pct1")) +
                           "% of it (rounded) in group 1";
  const auto sizes = [](const blocks::TwoN& t) { return "N1 = " + fint(t.n1) + ", N2 = " + fint(t.n2); };
  const double at = power(n.n1, n.n2);
  const blocks::TwoN less = blocks::smaller_two_group(s, n);
  if (less.n1 < min_n || less.n2 < min_n)
    return {"Group sizes", rule,
            sizes(n) + " are the smallest sizes allowed, and their power is already " + prob(at, target) + " ≥ " +
                fnum(target)};
  return {"Group sizes", rule,
          "power at " + sizes(less) + ": " + prob(power(less.n1, less.n2), target) + " < " + fnum(target) +
              " ≤ power at " + sizes(n) + ": " + prob(at, target)};
}

/// How a detectable effect was found: the value at which the power equals the target.
inline Step effect_search(const std::string& sym, double value, double power_at, double target) {
  return {"Detectable " + sym, sym + " = the value at which the power equals the target",
          "power at " + sym + " = " + num(value) + ": " + prob(power_at, target) + " (target " + fnum(target) + ")"};
}

}  // namespace ssc::procs::derive
