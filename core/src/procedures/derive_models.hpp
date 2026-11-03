// Formatting shared by the derivations of the survival, rate, variance, correlation and regression procedures, whose
// steps carry small quantities (hazards, variances of log ratios, effect sizes) that four decimals would blur.
#pragma once

#include <algorithm>
#include <cmath>
#include <string>

#include "procedures/derive.hpp"

namespace ssc::procs::derive {

/// A number to four significant digits however small (a hazard of 0.0005776, not 0.0006), as num writes it.
inline std::string sig(double x) {
  const double ax = std::fabs(x);
  const int decimals = ax > 0 && ax < 1 ? std::max(4, 3 - static_cast<int>(std::floor(std::log10(ax)))) : 4;
  return num(x, decimals);
}

}  // namespace ssc::procs::derive
