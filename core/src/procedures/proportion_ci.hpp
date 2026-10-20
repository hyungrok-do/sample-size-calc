// Confidence limits for one binomial proportion, shared by the precision
// procedures for one proportion and for the difference of two proportions.
// The count x need not be an integer: planning calculations evaluate the
// interval at the anticipated count x = n·p.
#pragma once

#include <algorithm>
#include <cmath>
#include <string>

#include "stats/stats.hpp"

namespace ssc::procs {

struct Interval {
  double lo, hi;
};

inline Interval prop_ci(const std::string& method, double n, double x, double conf, bool two_sided) {
  using namespace stats;
  const double alpha = 1 - conf;
  const double a = two_sided ? alpha / 2 : alpha;
  const double z = qnorm(a, false);
  const double ph = x / n, qh = 1 - ph;
  if (method == "wald") {
    const double h = z * std::sqrt(ph * qh / n);
    return {std::max(0.0, ph - h), std::min(1.0, ph + h)};
  }
  if (method == "wilson") {
    const double den = 1 + z * z / n;
    const double c = (ph + z * z / (2 * n)) / den;
    const double h = z / den * std::sqrt(ph * qh / n + z * z / (4 * n * n));
    return {std::max(0.0, c - h), std::min(1.0, c + h)};
  }
  if (method == "wilsoncc") {
    const double z2 = z * z;
    double lo = 0, hi = 1;
    if (ph > 0) lo = (2 * n * ph + z2 - 1 - z * std::sqrt(std::max(0.0, z2 - 2 - 1 / n + 4 * ph * (n * qh + 1)))) / (2 * (n + z2));
    if (ph < 1) hi = (2 * n * ph + z2 + 1 + z * std::sqrt(std::max(0.0, z2 + 2 - 1 / n + 4 * ph * (n * qh - 1)))) / (2 * (n + z2));
    return {std::max(0.0, lo), std::min(1.0, hi)};
  }
  if (method == "ac") {
    const double nt = n + z * z;
    const double pt_ = (x + z * z / 2) / nt;
    const double h = z * std::sqrt(pt_ * (1 - pt_) / nt);
    return {std::max(0.0, pt_ - h), std::min(1.0, pt_ + h)};
  }
  if (method == "jeffreys") {
    return {x <= 0 ? 0.0 : qbeta(a, x + 0.5, n - x + 0.5), x >= n ? 1.0 : qbeta(a, x + 0.5, n - x + 0.5, false)};
  }
  // exact (Clopper–Pearson)
  return {x <= 0 ? 0.0 : qbeta(a, x, n - x + 1), x >= n ? 1.0 : qbeta(a, x + 1, n - x, false)};
}

}  // namespace ssc::procs
