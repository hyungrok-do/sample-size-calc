// Shared steps for the derivations of the designs for means: the F-tests of the analyses of variance, the two
// one-sided tests of equivalence, the group sizes that follow from an allocation and the normal approximations shown
// beside a size found, each computed with the same functions as the runs they explain.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "procedures/derive.hpp"
#include "procedures/integrals.hpp"

namespace ssc::procs::derive {

/// Terms of a sum or a list as a step writes them, joined by `sep`: more than six as the first five, "…" and the last
/// (a design may have hundreds of groups or measurement times).
inline std::string joined(const std::vector<std::string>& t, const std::string& sep) {
  std::string out;
  for (std::size_t i = 0; i < t.size(); ++i) {
    if (t.size() > 6 && i == 5) {
      out += sep + "…";
      i = t.size() - 1;
    }
    out += (i ? sep : "") + t[i];
  }
  return out;
}

/// A number below 1 with four significant digits rather than four decimals (0.03125, not 0.0312), for a small
/// quantity that a step multiplies by large ones.
inline std::string num_sig(double x) {
  const double ax = std::fabs(x);
  if (!(ax > 0) || ax >= 1) return num(x);
  return num(x, std::min(12, 3 - static_cast<int>(std::floor(std::log10(ax)))));
}

/// A whole exponent as written after its base: "" for 1, "²" for 2, "¹²" for 12.
inline std::string superscript(int k) {
  static const char* const digit[] = {"⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹"};
  if (k == 1) return "";
  std::string out;
  for (const char c : std::to_string(k)) out += c == '-' ? "⁻" : digit[c - '0'];
  return out;
}

/// Marked text substituted after an operator: in parentheses when its value is negative, "0.5 − (−0.2)".
inline std::string means_paren(const std::string& marked, double value) {
  return value < 0 ? "(" + marked + ")" : marked;
}

// ------------------------------------------------------------------ F-tests

/// The critical value of an F-test with ν1 and ν2 degrees of freedom (`df_formula` says how they arise).
inline Step f_critical(double alpha, double df1, double df2, const std::string& df_formula) {
  return {"Critical value", "F = F(1 − α; ν1, ν2), the F quantile with " + df_formula + " degrees of freedom",
          "F(1 − " + gv(alpha) + "; " + dv(df1) + ", " + dv(df2) + ") = " + dv(stats::qf(alpha, df1, df2, false))};
}

/// The power of an F-test: F′ has the noncentral F distribution with ν1 and ν2 degrees of freedom and
/// noncentrality λ.
inline Step f_power_step(double df1, double df2, double lambda, double alpha) {
  const double f = stats::qf(alpha, df1, df2, false);
  return {"Power", "P(F′ > F), F′ ~ noncentral F(ν1, ν2, λ)",
          "P(F′ > " + dv(f) + ") = " + dvp(stats::pnf(f, df1, df2, lambda, false))};
}

/// Solving for an effect with an F-test: the noncentrality at which its power reaches the target, which the row's
/// effect (noncentrality `lambda`) has, to the search's precision.
inline Step means_f_ncp_needed(double df1, double df2, double lambda, double alpha, double target) {
  const double f = stats::qf(alpha, df1, df2, false);
  return {"Noncentrality needed", "λ solves P(F′ > F) = 1 − β, F′ ~ noncentral F(ν1, ν2, λ)",
          "P(F′ > " + dv(f) + ") = " + gv(target) + " when λ = " + dv(lambda)};
}

// ------------------------------------------------------------------ equivalence

/// The names of the limits and the true value in the steps of a TOST ("EL", "EU", "δ"; "ln θL", "ln θU", "ln θ0").
struct TostSymbols {
  std::string lower, upper, value;
};

/// The limits and the true value of a TOST as its steps write them, marked: "−[0.5]", "[0.5]" and "[0]" for limits
/// ±E around a difference entered, the derived logarithms for a ratio's.
struct TostShown {
  std::string lower, upper, value;
};

/// The two one-sided tests of equivalence, as tost_power: the estimate d ~ N(δ, SE²) must lie more than t estimated
/// standard errors inside both limits, the standard error estimated with ν degrees of freedom (`df_formula`; with ν
/// infinite or above 1e7 it is known, and the tests are z-tests). The steps are the standardized distances to the
/// limits, the critical value and the power, exactly: an expectation over the ratio U of the estimated to the true
/// standard error.
inline Steps tost_steps(double el, double eu, double delta, double se, double df, double alpha,
                        const std::string& df_formula, const TostSymbols& sym, const TostShown& shown) {
  const double a = (eu - delta) / se, b = (el - delta) / se;
  const double p = tost_power(el, eu, delta, se, df, alpha);
  const std::string value = means_paren(shown.value, delta);
  Steps out;
  out.push_back({"Distances to the limits", "a = (" + sym.upper + " − " + sym.value + ") / SE, b = (" + sym.lower +
                                                " − " + sym.value + ") / SE",
                 "a = (" + shown.upper + " − " + value + ") / " + dv(se) + " = " + dv(a) + ", b = (" + shown.lower +
                     " − " + value + ") / " + dv(se) + " = " + dv(b)});
  if (!std::isfinite(df) || df > 1e7) {  // as tost_power: the standard error is known
    const double z = stats::qnorm(alpha, false);
    out.push_back(z_critical(alpha, Alt::Upper));
    out.push_back({"Power", "P(both tests reject) = Φ(a − z) − Φ(b + z)",
                   "Φ(" + dv(a) + " − " + dv(z) + ") − Φ(" + dv(b) + " + " + dv(z) + ") = " +
                       dvp(stats::pnorm(a - z)) + " − " + dvp(stats::pnorm(b + z)) + " = " + dvp(p)});
    return out;
  }
  const double t = stats::qt(alpha, df, false), umax = (eu - el) / (2 * t * se);
  out.push_back(t_critical(alpha, df, Alt::Upper, df_formula));
  out.push_back({"Power",
                 "P(both tests reject) = E[Φ(a − tU) − Φ(b + tU); U < (a − b)/(2t)], over U = √(χ²ν/ν), the "
                 "estimated SE over the true one",
                 "E[Φ(" + dv(a) + " − " + dv(t) + " U) − Φ(" + dv(b) + " + " + dv(t) + " U); U < " + dv(umax) +
                     "] = " + dvp(p)});
  return out;
}

/// Solving for the size of a TOST: the closed form of the z-tests when only the nearer limit counts (Chow, Shao &
/// Wang), for comparison beside the exact answer: `lhs` ≈ (z(1 − α) + z(1 − β))² v / g², g the distance from the
/// true value to the nearer limit, with z(1 − β/2) when the true value lies midway between the limits (as the run's
/// starting guess decides it). `v_formula` and `v_numbers` write the variance factor v, its value `v`.
inline Step means_tost_approx(double el, double eu, double delta, double alpha, double target,
                               const TostSymbols& sym, const TostShown& shown, const std::string& lhs,
                               const std::string& v_formula, const std::string& v_numbers, double v) {
  const bool mid = std::fabs(delta - (el + eu) / 2) < 1e-12;
  const double gap = std::min(delta - el, eu - delta);
  const double za = stats::qnorm(alpha, false), zb = stats::qnorm(mid ? 1 - (1 - target) / 2 : target);
  return {"Normal approximation",
          lhs + " ≈ (z(1 − α) + z(1 − β" + (mid ? "/2" : "") + "))² " + v_formula + " / g², g = min(" + sym.value +
              " − " + sym.lower + ", " + sym.upper + " − " + sym.value + ")" +
              (mid ? ", with β/2 as " + sym.value + " lies midway between the limits" : "") +
              ", for comparison; the steps below give the answer",
          "g = min(" + shown.value + " − " + means_paren(shown.lower, el) + ", " + shown.upper + " − " +
              means_paren(shown.value, delta) + ") = " + dv(gap) + "; (" + dv(za) + " + " + dv(zb) + ")² × " +
              v_numbers + " / " + dv(gap) + "² = " + dv((za + zb) * (za + zb) * v / (gap * gap))};
}

/// How an even total was found (two sequences or groups of equal size, `halves`): as n_search, in steps of 2.
inline Step even_n_search(const std::function<double(double)>& power, double n, double target, double lo,
                          const std::string& halves) {
  Step st = n_search(power, n, target, lo, "N", 2);
  st.formula = "N = the smallest even number whose power reaches the target, " + halves;
  return st;
}

// ------------------------------------------------------------------ two groups

/// The two group sizes as the steps write them, each marked as it arises: solved for (derived, but for a group whose
/// size was entered and held fixed), entered (given), or computed from the inputs (N2 = ⌈r × N1⌉, a total split by a
/// percentage: derived).
inline std::pair<std::string, std::string> means_group_sizes(const Scenario& s, const blocks::TwoN& n, bool solved) {
  const std::string& a = s.s("alloc");
  if (solved) return {a == "n1" ? gvi(n.n1) : dvi(n.n1), a == "n2" ? gvi(n.n2) : dvi(n.n2)};
  if (a == "pct") return {dvi(n.n1), dvi(n.n2)};
  if (a == "ratio") return {gvi(n.n1), dvi(n.n2)};
  return {gvi(n.n1), gvi(n.n2)};
}

/// The group sizes a power or an effect is computed at, in the words of `find`: as entered, or the inputs they
/// follow from (means_alloc_step derives them).
inline std::string means_sizes_at(const Scenario& s, const blocks::TwoN& n) {
  const std::string& a = s.s("alloc");
  if (a == "ratio") return "N1 = " + gvi(n.n1) + ", N2 = ⌈r × N1⌉ (r = " + gv(s("r")) + ")";
  if (a == "pct") return "N = " + gvi(s("ntotal")) + ", " + gv(s("pct1")) + "% of it (rounded) in group 1";
  return two_sizes(n, false);
}

/// The group sizes that follow from the inputs when they are not entered as they are (solving for power or an
/// effect): N2 from N1 and a ratio, or both from a total and a percentage, as two_group_n.
inline void means_alloc_step(Steps& out, const Scenario& s, const blocks::TwoN& n) {
  const std::string& a = s.s("alloc");
  if (a == "ratio")
    out.push_back({"Group 2 size", "N2 = ⌈r × N1⌉", "⌈" + gv(s("r")) + " × " + gvi(n.n1) + "⌉ = " + dvi(n.n2)});
  else if (a == "pct")
    out.push_back({"Group sizes", "N1 = N × p/100 rounded to a whole number, N2 = N − N1",
                   "N1 = round(" + gvi(s("ntotal")) + " × " + gv(s("pct1")) + "/100) = " + dvi(n.n1) + ", N2 = " +
                       gvi(s("ntotal")) + " − " + dvi(n.n1) + " = " + dvi(n.n2)});
}

/// A term of a normal approximation: its formula, its numbers (marked) and its value.
struct MeansTerm {
  std::string formula, numbers;
  double value;
};

/// Solving for two group sizes with equal groups or N2 = ⌈r × N1⌉, where group 2 follows group 1: the closed form a
/// normal statistic gives, for comparison beside the exact answer, N1 ≈ (z(1 − α/2) + z(1 − β))² v (1 + 1/r) / d² for
/// a difference d whose standard error is √(v (1/N1 + 1/N2)). With one group fixed, or a percentage of a total, no
/// closed form applies as simply, and none is shown.
inline void means_alloc_approx(Steps& out, const Scenario& s, Alt a, const MeansTerm& v, const MeansTerm& d) {
  const std::string& alloc = s.s("alloc");
  if (alloc != "equal" && alloc != "ratio") return;
  const bool equal = alloc == "equal";
  const double k = equal ? 1 : s("r");
  const double za = stats::qnorm(tail_alpha(s("alpha"), a), false), zb = stats::qnorm(s("power"));
  out.push_back({"Normal approximation",
                 std::string(equal ? "N1 = N2" : "N1") + " ≈ (z(" + level_symbol(a) + ") + z(1 − β))² " +
                     (equal ? "2" + v.formula : v.formula + "(1 + 1/r)") + " / " + d.formula +
                     "², for comparison; the steps below give the answer",
                 "(" + dv(za) + " + " + dv(zb) + ")² × " +
                     (equal ? "2 × " + v.numbers : v.numbers + " × (1 + 1/" + gv(k) + ")") + " / " + d.numbers +
                     "² = " + dv((za + zb) * (za + zb) * v.value * (1 + 1 / k) / (d.value * d.value))});
}

/// two_group_search, with the group size that was entered and held fixed (one group's size given, the other's
/// found) marked as given: the shared step marks every size it shows as found.
inline Step means_two_search(const Scenario& s, const std::function<double(double, double)>& power,
                             const blocks::TwoN& n, double target, double min_n) {
  Step st = two_group_search(s, power, n, target, min_n);
  const std::string& a = s.s("alloc");
  if (a == "n1" || a == "n2") {
    const std::string sym = a == "n1" ? "N1 = " : "N2 = ";
    const double fixed = a == "n1" ? n.n1 : n.n2;
    const std::string from = sym + dvi(fixed), to = sym + gvi(fixed);
    for (std::size_t at = st.numbers.find(from); at != std::string::npos; at = st.numbers.find(from, at + to.size()))
      st.numbers.replace(at, from.size(), to);
  }
  return st;
}

/// The answer of a two-group size solved for: the sizes (a fixed one given), their total and the power there.
inline std::string means_two_answer(const Scenario& s, const blocks::TwoN& n, double power) {
  const auto sizes = means_group_sizes(s, n, true);
  return "N1 = " + sizes.first + ", N2 = " + sizes.second + " (N = " + dvi(n.total()) + "), with power " + dvp(power);
}

}  // namespace ssc::procs::derive
