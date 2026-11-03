// Formatting and steps shared by the derivations of the survival, rate, variance, correlation and regression
// procedures, in the manner of derive.hpp. Their steps carry small quantities (hazards, variances of log ratios,
// effect sizes) that four decimals would blur.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <utility>

#include "procedures/apart.hpp"
#include "procedures/derive.hpp"

namespace ssc::procs::derive {

/// A number to four significant digits however small (a hazard of 0.0005776, not 0.0006), as num writes it.
inline std::string sig(double x) {
  const double ax = std::fabs(x);
  const int decimals = ax > 0 && ax < 1 ? std::max(4, 3 - static_cast<int>(std::floor(std::log10(ax)))) : 4;
  return num(x, decimals);
}
/// The same, given (an input) or derived (computed).
inline std::string gsig(double x) { return given(sig(x)); }
inline std::string dsig(double x) { return derived(sig(x)); }

/// Two numbers with the digits that tell them apart (fnum_apart) and print's minus sign: an effect next to its value
/// under H0, which four decimals can show as equal (a rate found next to a large λ0).
inline std::pair<std::string, std::string> apart(double a, double b) {
  const auto printed = [](std::string t) { return !t.empty() && t[0] == '-' ? "−" + t.substr(1) : t; };
  const auto [ta, tb] = fnum_apart(a, b);
  return {printed(ta), printed(tb)};
}

/// A step with `from` written as `to` in its formula and numbers: θ for the drift of a logrank statistic, where λ
/// is a hazard.
inline Step renamed(Step st, const std::string& from, const std::string& to) {
  for (std::string* t : {&st.formula, &st.numbers})
    for (std::size_t at = 0; (at = t->find(from, at)) != std::string::npos; at += to.size())
      t->replace(at, from.size(), to);
  return st;
}

/// The normal quantiles of a closed form for comparison, z(1 − α/2) (or z(1 − α)) and z(1 − β): in symbols,
/// "(z(1 − α/2) + z(1 − β))²", and with their values.
inline std::string z_sum_formula(Alt a) { return std::string("(z(") + level_symbol(a) + ") + z(1 − β))²"; }
inline double z_alpha(double alpha, Alt a) { return stats::qnorm(tail_alpha(alpha, a), false); }
inline std::string z_sum_numbers(double alpha, Alt a, double target) {
  return "(" + dv(z_alpha(alpha, a)) + " + " + dv(stats::qnorm(target)) + ")²";
}

/// Solving for an effect by search, as effect_search, with the condition the effect meets written out (`power`,
/// the power in terms of the effect) and its terms at the value found (`terms`, ending before the power itself).
inline Step effect_found(const std::string& sym, const std::string& value, const std::string& power,
                         const std::string& terms, double power_at, double target) {
  return {"Detectable " + sym, sym + " = the value at which the power, " + power + ", equals 1 − β",
          "power at " + sym + " = " + value + ": " + (terms.empty() ? "" : terms + " = ") + dvp(power_at, target) +
              " = " + gv(target)};
}

/// The answer to a sample size of one group, and to two (N1, N2 and their total), with the power there.
inline std::string size_answer(double n, double power) { return "N = " + dvi(n) + ", with power " + dvp(power); }
inline std::string sizes_answer(const blocks::TwoN& n, const std::string& with) {
  return two_sizes(n, true) + " (N = " + dvi(n.total()) + "), with " + with;
}

}  // namespace ssc::procs::derive
