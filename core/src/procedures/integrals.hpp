// Expectations over the distribution of a sample-SD ratio and exact TOST power.
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "stats/numerics.hpp"
#include "stats/stats.hpp"

namespace ssc::procs {

/// E[ f(U) · 1{U < umax} ] where U = √(χ²_df / df) (the ratio S/σ of a sample
/// SD to the true SD).  Adaptive Gauss–Kronrod on pieces placed around the bulk
/// of U (mode ≈ 1, scale ≈ 1/√(2 df)).
template <class F>
double expect_sd_ratio(F&& f, double df, double umax = INFINITY, double eps_abs = 1e-15, double eps_rel = 1e-12) {
  using namespace stats;
  if (!(umax > 0)) return 0;
  const double lc = std::log(2.0) + (df / 2) * std::log(df / 2) - lgammafn(df / 2);
  auto g = [&](double u) {
    if (u <= 0) return df == 1 ? std::exp(lc) : 0.0;
    return std::exp(lc + (df - 1) * std::log(u) - df * u * u / 2);
  };
  const double s = 1 / std::sqrt(2 * df);
  std::vector<double> pts = {0.0};
  for (double z : {-14.0, -7.0, -3.5, -1.5, 0.0, 1.5, 3.5, 7.0, 14.0, 28.0, 45.0}) {
    const double p = 1 + z * s;
    if (p > pts.back() && p < umax) pts.push_back(p);
  }
  pts.push_back(std::min(umax, 1 + 60 * s + 1));
  double total = 0;
  for (std::size_t i = 0; i + 1 < pts.size(); ++i) {
    if (pts[i + 1] <= pts[i]) continue;
    total += numerics::integrate([&](double u) { return f(u) * g(u); }, pts[i], pts[i + 1], eps_abs, eps_rel, 25);
  }
  return total;
}

/// Exact power of the two one-sided tests (TOST) for equivalence.
/// The estimate d ~ N(delta, se²) and its standard error is estimated with
/// `df` degrees of freedom (df = ∞ → known variance, z-tests).  Equivalence is
/// concluded when  EL + t·SÊ < d < EU − t·SÊ  with t the 1 − α quantile.
inline double tost_power(double el, double eu, double delta, double se, double df, double alpha) {
  using namespace stats;
  if (!(eu > el) || !(se > 0)) return 0;
  const double a = (eu - delta) / se;
  const double b = (el - delta) / se;
  if (!std::isfinite(df) || df > 1e7) {
    const double z = qnorm(alpha, false);
    return std::max(0.0, pnorm(a - z) - pnorm(b + z));
  }
  const double t = qt(alpha, df, false);
  const double umax = (eu - el) / (2 * t * se);
  const double p = expect_sd_ratio([&](double u) { return std::max(0.0, pnorm(a - t * u) - pnorm(b + t * u)); }, df,
                                   umax);
  return clamp01(p);
}

}  // namespace ssc::procs
