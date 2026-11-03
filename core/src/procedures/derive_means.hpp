// Shared steps for the derivations of the designs for means: the F-tests of the analyses of variance and the two
// one-sided tests of equivalence, each computed with the same functions as the runs they explain.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
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

/// The critical value of an F-test with ν1 and ν2 degrees of freedom (`df_formula` says how they arise).
inline Step f_critical(double alpha, double df1, double df2, const std::string& df_formula) {
  return {"Critical value", "F = F(1 − α; ν1, ν2), the F quantile with " + df_formula + " degrees of freedom",
          "F(1 − " + num(alpha) + "; " + num(df1) + ", " + num(df2) + ") = " + num(stats::qf(alpha, df1, df2, false))};
}

/// The power of an F-test: F′ has the noncentral F distribution with ν1 and ν2 degrees of freedom and
/// noncentrality λ.
inline Step f_power_step(double df1, double df2, double lambda, double alpha) {
  const double f = stats::qf(alpha, df1, df2, false);
  return {"Power", "P(F′ > F), F′ ~ noncentral F(ν1, ν2, λ)",
          "P(F′ > " + num(f) + ") = " + prob(stats::pnf(f, df1, df2, lambda, false))};
}

/// The names of the limits and the true value in the steps of a TOST ("EL", "EU", "δ"; "ln θL", "ln θU", "ln θ0").
struct TostSymbols {
  std::string lower, upper, value;
};

/// The two one-sided tests of equivalence, as tost_power: the estimate d ~ N(δ, SE²) must lie more than t estimated
/// standard errors inside both limits, the standard error estimated with ν degrees of freedom (`df_formula`; with ν
/// infinite or above 1e7 it is known, and the tests are z-tests). The steps are the standardized distances to the
/// limits, the critical value and the power, exactly: an expectation over the ratio U of the estimated to the true
/// standard error.
inline Steps tost_steps(double el, double eu, double delta, double se, double df, double alpha,
                        const std::string& df_formula, const TostSymbols& sym) {
  const double a = (eu - delta) / se, b = (el - delta) / se;
  const double p = tost_power(el, eu, delta, se, df, alpha);
  Steps out;
  out.push_back({"Distances to the limits", "a = (" + sym.upper + " − " + sym.value + ") / SE, b = (" + sym.lower +
                                                " − " + sym.value + ") / SE",
                 "a = (" + num(eu) + " − " + sub(delta) + ") / " + num(se) + " = " + num(a) + ", b = (" + num(el) +
                     " − " + sub(delta) + ") / " + num(se) + " = " + num(b)});
  if (!std::isfinite(df) || df > 1e7) {  // as tost_power: the standard error is known
    const double z = stats::qnorm(alpha, false);
    out.push_back(z_critical(alpha, Alt::Upper));
    out.push_back({"Power", "P(both tests reject) = Φ(a − z) − Φ(b + z)",
                   "Φ(" + num(a) + " − " + num(z) + ") − Φ(" + num(b) + " + " + num(z) + ") = " +
                       prob(stats::pnorm(a - z)) + " − " + prob(stats::pnorm(b + z)) + " = " + prob(p)});
    return out;
  }
  const double t = stats::qt(alpha, df, false), umax = (eu - el) / (2 * t * se);
  out.push_back(t_critical(alpha, df, Alt::Upper, df_formula));
  out.push_back({"Power",
                 "P(both tests reject) = E[Φ(a − tU) − Φ(b + tU); U < (a − b)/(2t)], over U = √(χ²ν/ν), the "
                 "estimated SE over the true one",
                 "E[Φ(" + num(a) + " − " + num(t) + " U) − Φ(" + num(b) + " + " + num(t) + " U); U < " + num(umax) +
                     "] = " + prob(p)});
  return out;
}

/// How an even total was found (two sequences or groups of equal size, `halves`): as n_search, in steps of 2.
inline Step even_n_search(const std::function<double(double)>& power, double n, double target, double lo,
                          const std::string& halves) {
  Step st = n_search(power, n, target, lo, "N", 2);
  st.formula = "N = the smallest even number whose power reaches the target, " + halves;
  return st;
}

}  // namespace ssc::procs::derive
