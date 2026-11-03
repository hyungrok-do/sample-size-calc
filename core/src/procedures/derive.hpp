// Shared steps for the derivations that show how a result came about. Each procedure plugs its own derivation
// together (p.derivation) from these and its own quantities, computing every value with the same functions as its
// run, so that the steps and the answer cannot disagree.
//
// A derivation follows what was solved for. Solving for a sample size, it finds the smallest size whose power reaches
// the target and shows the power there and one size below; solving for power, it computes the power at the sizes
// given; solving for an effect, it works back from the target power to the effect. Its values are marked as given (an
// input, gv) or derived (computed, dv), which the page shows apart; the inputs themselves are listed by the runner.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "procedures/common.hpp"

namespace ssc::procs::derive {

using blocks::fint;
using blocks::fnum;
using Steps = std::vector<Step>;
using Derivation = ssc::Derivation;

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

// ------------------------------------------------------------------ given and derived values
// Marked text: \x01 (given) or \x02 (derived) before a value and \x03 after it; the runner turns it into plain text
// and segments (see set_marked in runner.cpp).
inline std::string given(const std::string& text) { return "\x01" + text + "\x03"; }
inline std::string derived(const std::string& text) { return "\x02" + text + "\x03"; }
/// A given number (an input), and as substituted (a negative one in parentheses).
inline std::string gv(double x, int decimals = 4) { return given(num(x, decimals)); }
inline std::string gvs(double x, int decimals = 4) { return x < 0 ? "(" + gv(x, decimals) + ")" : gv(x, decimals); }
inline std::string gvi(double n) { return given(fint(n)); }
/// A derived number (computed), as substituted, a whole one, and a probability told apart from `against`.
inline std::string dv(double x, int decimals = 4) { return derived(num(x, decimals)); }
inline std::string dvs(double x, int decimals = 4) { return x < 0 ? "(" + dv(x, decimals) + ")" : dv(x, decimals); }
inline std::string dvi(double n) { return derived(fint(n)); }
inline std::string dvp(double p, double against = NAN) { return derived(prob(p, against)); }
/// A size that is given when solving for power or an effect, and derived when it is what was solved for.
inline std::string size(double n, bool solved) { return solved ? dvi(n) : gvi(n); }

inline const char* level_symbol(Alt a) { return a == Alt::Two ? "1 − α/2" : "1 − α"; }
/// The level of the critical value with the given alpha in place: "1 − 0.05/2" or "1 − 0.05".
inline std::string level_numbers(double alpha, Alt a) { return "1 − " + gv(alpha) + (a == Alt::Two ? "/2" : ""); }

/// The critical value of a test whose statistic is standard normal under H0.
inline Step z_critical(double alpha, Alt a) {
  const double tail = tail_alpha(alpha, a);
  return {"Critical value", std::string("z = Φ⁻¹(") + level_symbol(a) + ")",
          "Φ⁻¹(" + level_numbers(alpha, a) + ") = " + dv(stats::qnorm(tail, false))};
}

/// The critical value of a t statistic with `df` degrees of freedom (`df_formula` says how they arise).
inline Step t_critical(double alpha, double df, Alt a, const std::string& df_formula) {
  const double tail = tail_alpha(alpha, a);
  return {"Critical value",
          std::string("t = t(") + level_symbol(a) + "; ν), the t quantile with ν = " + df_formula +
              " degrees of freedom",
          "t(" + level_numbers(alpha, a) + "; " + dv(df) + ") = " + dv(stats::qt(tail, df, false))};
}

/// The power of a test whose statistic is N(λ, 1) under H1 (as z_power).
inline Step z_power_step(double ncp, double alpha, Alt a) {
  const double z = stats::qnorm(tail_alpha(alpha, a), false), p = z_power(ncp, alpha, a);
  if (a == Alt::Two)
    return {"Power", "Φ(λ − z) + Φ(−λ − z)",
            "Φ(" + dv(ncp) + " − " + dv(z) + ") + Φ(" + dv(-ncp) + " − " + dv(z) + ") = " +
                dvp(stats::pnorm(ncp - z)) + " + " + dvp(stats::pnorm(-ncp - z)) + " = " + dvp(p)};
  if (a == Alt::Upper) return {"Power", "Φ(λ − z)", "Φ(" + dv(ncp) + " − " + dv(z) + ") = " + dvp(p)};
  return {"Power", "Φ(−λ − z)", "Φ(" + dv(-ncp) + " − " + dv(z) + ") = " + dvp(p)};
}

/// The power of a t-test (as t_power): T′ has the noncentral t distribution with ν degrees of freedom and
/// noncentrality λ.
inline Step t_power_step(double ncp, double df, double alpha, Alt a) {
  const double t = stats::qt(tail_alpha(alpha, a), df, false), p = t_power(ncp, df, alpha, a);
  const std::string law = ", T′ ~ noncentral t(ν, λ)";
  if (a == Alt::Two)
    return {"Power", "P(T′ > t) + P(T′ < −t)" + law,
            "P(T′ > " + dv(t) + ") + P(T′ < −" + dv(t) + ") = " + dvp(stats::pnt(t, df, ncp, false)) + " + " +
                dvp(stats::pnt(-t, df, ncp, true)) + " = " + dvp(p)};
  if (a == Alt::Upper) return {"Power", "P(T′ > t)" + law, "P(T′ > " + dv(t) + ") = " + dvp(p)};
  return {"Power", "P(T′ < −t)" + law, "P(T′ < −" + dv(t) + ") = " + dvp(p)};
}

/// Solving for an effect: the noncentrality at which a test whose statistic is N(λ, 1) has the target power, which
/// the row's effect (noncentrality `ncp`) has, to the search's precision.
inline Step z_ncp_needed(double ncp, double alpha, Alt a, double target) {
  const double z = stats::qnorm(tail_alpha(alpha, a), false);
  const std::string eq = a == Alt::Two ? "Φ(λ − z) + Φ(−λ − z)" : a == Alt::Upper ? "Φ(λ − z)" : "Φ(−λ − z)";
  const std::string at = a == Alt::Two   ? "Φ(" + dv(ncp) + " − " + dv(z) + ") + Φ(" + dv(-ncp) + " − " + dv(z) + ")"
                         : a == Alt::Upper ? "Φ(" + dv(ncp) + " − " + dv(z) + ")"
                                           : "Φ(" + dv(-ncp) + " − " + dv(z) + ")";
  return {"Noncentrality needed", "λ solves " + eq + " = 1 − β", at + " = " + gv(target) + " when λ = " + dv(ncp)};
}

/// Solving for an effect with a t-test: the noncentrality at which the power reaches the target.
inline Step t_ncp_needed(double ncp, double df, double alpha, Alt a, double target) {
  const double t = stats::qt(tail_alpha(alpha, a), df, false);
  const std::string eq = a == Alt::Two ? "P(T′ > t) + P(T′ < −t)" : a == Alt::Upper ? "P(T′ > t)" : "P(T′ < −t)";
  const std::string at = a == Alt::Two   ? "P(T′ > " + dv(t) + ") + P(T′ < −" + dv(t) + ")"
                         : a == Alt::Upper ? "P(T′ > " + dv(t) + ")"
                                           : "P(T′ < −" + dv(t) + ")";
  return {"Noncentrality needed", "λ solves " + eq + " = 1 − β, T′ ~ noncentral t(ν, λ)",
          at + " = " + gv(target) + " when λ = " + dv(ncp)};
}

/// Solving for a sample size: the size found is the smallest admissible one (from `lo`, in steps of `step`) whose
/// power reaches the target, shown by the power one step below it and at it.
inline Step n_search(const std::function<double(double)>& power, double n, double target, double lo,
                     const std::string& sym = "N", double step = 1) {
  const double at = power(n);
  if (n - step < lo)
    return {"Smallest size", sym + " = the smallest size whose power reaches 1 − β",
            sym + " = " + dvi(n) + " is the smallest size allowed, and its power is already " + dvp(at, target) +
                " ≥ " + gv(target)};
  return {"Smallest size", sym + " = the smallest whole number whose power reaches 1 − β",
          "power at " + sym + " = " + dvi(n - step) + ": " + dvp(power(n - step), target) + " < " + gv(target) +
              " ≤ power at " + sym + " = " + dvi(n) + ": " + dvp(at, target)};
}

/// The unknown of a two-group size solve under the scenario's allocation rule, in words: what `find` names.
inline std::string two_group_unknown(const Scenario& s) {
  const std::string& a = s.s("alloc");
  if (a == "equal") return "N1 = N2";
  if (a == "ratio") return "N1, with N2 = ⌈r × N1⌉ (r = " + gv(s("r")) + ")";
  if (a == "n1") return "N2, with N1 = " + gvi(s("n1")) + " fixed";
  if (a == "n2") return "N1, with N2 = " + gvi(s("n2")) + " fixed";
  return "N, with " + gv(s("pct1")) + "% of it (rounded) in group 1";
}
/// Two group sizes as written, given or derived.
inline std::string two_sizes(const blocks::TwoN& n, bool solved) {
  return "N1 = " + size(n.n1, solved) + ", N2 = " + size(n.n2, solved);
}

/// Solving for two group sizes under the scenario's allocation rule (see blocks::smaller_two_group): the power at
/// the next smaller design falls short of the target, and at the answer reaches it.
inline Step two_group_search(const Scenario& s, const std::function<double(double, double)>& power,
                             const blocks::TwoN& n, double target, double min_n) {
  const std::string& a = s.s("alloc");
  const std::string rule =
      a == "equal"   ? "N1 = N2, the smallest whole number whose power reaches 1 − β"
      : a == "ratio" ? "N1 = the smallest whole number whose power reaches 1 − β, with N2 = ⌈r × N1⌉"
      : a == "n1"    ? "N2 = the smallest whole number whose power reaches 1 − β, with N1 fixed"
      : a == "n2"    ? "N1 = the smallest whole number whose power reaches 1 − β, with N2 fixed"
                     : "N = the smallest total whose power reaches 1 − β, a fixed percentage of it in group 1";
  const double at = power(n.n1, n.n2);
  const blocks::TwoN less = blocks::smaller_two_group(s, n);
  if (less.n1 < min_n || less.n2 < min_n)
    return {"Smallest sizes", rule,
            two_sizes(n, true) + " are the smallest sizes allowed, and their power is already " + dvp(at, target) +
                " ≥ " + gv(target)};
  return {"Smallest sizes", rule,
          "power at " + two_sizes(less, true) + ": " + dvp(power(less.n1, less.n2), target) + " < " + gv(target) +
              " ≤ power at " + two_sizes(n, true) + ": " + dvp(at, target)};
}

/// Solving for an effect by search (where no noncentrality step applies): the value at which the power equals the
/// target.
inline Step effect_search(const std::string& sym, double value, double power_at, double target) {
  return {"Detectable " + sym, sym + " = the value at which the power equals 1 − β",
          "power at " + sym + " = " + dv(value) + ": " + dvp(power_at, target) + " = " + gv(target)};
}

/// The usual words of `find`: the smallest size reaching the target power, the power at given sizes, and the
/// effect detected with the target power.
inline std::string find_size(const std::string& unknown, double target) {
  return "the smallest " + unknown + " whose power is at least 1 − β = " + gv(target);
}
inline std::string find_power(const std::string& at) { return "the power " + at; }
inline std::string find_effect(const std::string& unknown, double target, const std::string& at) {
  return "the " + unknown + " the test detects with power 1 − β = " + gv(target) + ", " + at;
}

}  // namespace ssc::procs::derive
