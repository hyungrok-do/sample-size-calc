// Shared helpers for procedure implementations.
#pragma once

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "engine/blocks.hpp"
#include "engine/procedure.hpp"
#include "engine/solve.hpp"
#include "stats/stats.hpp"

namespace ssc::procs {

using blocks::Alt;

/// A share α/k of the level: α/2 for each tail of a two-sided test, α/k for each of k comparisons.  The share of
/// the smallest alphas (α/2 of α = 5e-324) rounds to 0, whose critical value is infinite, so that no study would
/// reach any power; it stays at the smallest positive double, the tail a one-sided test at that α has.
inline double alpha_share(double alpha, double k) {
  const double share = alpha / k;
  return share == 0 && alpha > 0 ? std::numeric_limits<double>::denorm_min() : share;
}

/// α/2, the level of each tail of a two-sided test.
inline double half_alpha(double alpha) { return alpha_share(alpha, 2); }

/// The level of the tail in which a test rejects: α/2 for each tail of a two-sided test, α for a one-sided one.
inline double tail_alpha(double alpha, Alt alt) { return alt == Alt::Two ? half_alpha(alpha) : alpha; }

/// Power of a test whose statistic is N(ncp, 1) under H1 (known variance).
inline double z_power(double ncp, double alpha, Alt alt) {
  using namespace stats;
  if (alt == Alt::Two) {
    const double z = qnorm(half_alpha(alpha), false);
    return pnorm(ncp - z) + pnorm(-ncp - z);
  }
  const double z = qnorm(alpha, false);
  return alt == Alt::Upper ? pnorm(ncp - z) : pnorm(-ncp - z);
}

/// Power of a t-test with `df` degrees of freedom and noncentrality `ncp`.
inline double t_power(double ncp, double df, double alpha, Alt alt) {
  using namespace stats;
  if (!(df > 0)) return NAN;
  if (alt == Alt::Two) {
    const double t = qt(half_alpha(alpha), df, false);
    return pnt(t, df, ncp, false) + pnt(-t, df, ncp, true);
  }
  const double t = qt(alpha, df, false);
  return alt == Alt::Upper ? pnt(t, df, ncp, false) : pnt(-t, df, ncp, true);
}

/// Normal-approximation sample size multiplier ((z_α' + z_β)²).
inline double z_factor(double alpha, double power, Alt alt) {
  using namespace stats;
  const double za = qnorm(tail_alpha(alpha, alt), false);
  const double zb = qnorm(power);
  return (za + zb) * (za + zb);
}

/// "two-sided" / "one-sided (upper)" etc. for summary statements.
inline std::string alt_words(Alt a) {
  return a == Alt::Two ? "two-sided" : "one-sided";
}

/// The power column(s) and the target value for a row when solving for N.
inline void put_power(const Scenario& s, Row& r, double power) {
  if (s.solve == "n") r.set("target_power", s("power"));
  r.set("power", std::isnan(power) ? power : std::fmin(1.0, std::fmax(0.0, power)));  // no 1 + 7e-16 from rounding
}

/// Phrase: "achieve 80.7% power" (+ target when solving for N).
inline std::string achieves(const Scenario& s, const Row& r) {
  std::string t = "achieve" + std::string(" ") + blocks::power_text(r.num("power")) + " power";
  if (s.solve == "n") t += " (target " + blocks::power_text(s("power")) + ")";
  return t;
}

inline std::string achieves_singular(const Scenario& s, const Row& r) {
  std::string t = "achieves " + blocks::power_text(r.num("power")) + " power";
  if (s.solve == "n") t += " (target " + blocks::power_text(s("power")) + ")";
  return t;
}

/// Commonly cited methodological references.
inline const char* REF_CHOW =
    "Chow, S.-C., Shao, J., Wang, H. & Lokhnygina, Y. (2018). Sample Size Calculations in Clinical Research, "
    "3rd ed. Chapman & Hall/CRC.";
/// "0, 1, 2, 4" for list inputs shown in result tables and statements.
inline std::string list_text(const std::vector<double>& v) {
  std::string out;
  for (std::size_t i = 0; i < v.size(); ++i) out += (i ? ", " : "") + blocks::fnum(v[i]);
  return out;
}

/// Text stored in a result cell (empty if absent or numeric).
inline std::string text_of(const Row& r, const std::string& key) {
  for (const auto& kv : r.cells)
    if (kv.first == key)
      if (const std::string* p = std::get_if<std::string>(&kv.second)) return *p;
  return {};
}

inline const char* REF_JULIOUS = "Julious, S. A. (2010). Sample Sizes for Clinical Trials. Chapman & Hall/CRC.";
inline const char* REF_MACHIN =
    "Machin, D., Campbell, M. J., Tan, S. B. & Tan, S. H. (2018). Sample Sizes for Clinical, Laboratory and "
    "Epidemiology Studies, 4th ed. Wiley-Blackwell.";
inline const char* REF_COHEN =
    "Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences, 2nd ed. Lawrence Erlbaum.";
inline const char* REF_LENTH =
    "Lenth, R. V. (1989). Algorithm AS 243: cumulative distribution function of the non-central t distribution. "
    "Applied Statistics 38, 185–189.";
inline const char* REF_BK =
    "Benton, D. & Krishnamoorthy, K. (2003). Computing discrete mixtures of continuous distributions: noncentral "
    "chisquare, noncentral t and the distribution of the square of the sample multiple correlation coefficient. "
    "Computational Statistics & Data Analysis 43, 249–267.";
inline const char* REF_FLEISS =
    "Fleiss, J. L., Levin, B. & Paik, M. C. (2003). Statistical Methods for Rates and Proportions, 3rd ed. Wiley.";

}  // namespace ssc::procs
