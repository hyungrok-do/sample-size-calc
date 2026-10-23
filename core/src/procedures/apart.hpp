// Values that a statement shows side by side, with the digits it takes to tell them apart.
#pragma once

#include <cmath>
#include <initializer_list>
#include <string>
#include <utility>

#include "procedures/common.hpp"

namespace ssc::procs {

/// x as fnum_fine shows it (`extra` = 0), or with `extra` more digits: decimals where fnum shows a fixed number of
/// them, significant digits where it shows an exponent or a number of 1e6 or more.  Whole numbers are exact, and,
/// as with fnum_fine, a number that is not whole never shows as one (0.9999999 as 1, 1000000000.66 as 1000000001):
/// it takes the digits it needs to show a decimal point.
inline std::string fnum_digits(double x, int extra) {
  if (extra == 0 || !std::isfinite(x) || (x == std::floor(x) && std::fabs(x) < 1e15)) return fnum_fine(x);
  const double ax = std::fabs(x);
  std::string t;
  for (int d = 4 + extra; d <= 17; ++d) {
    t = ax >= 1e6 || ax < 1e-4 ? text::format_g(x, d) : blocks::fnum(x, d);
    if (t.find_first_of(".e") != std::string::npos) break;
  }
  return t;
}

/// The texts of a and b, two values a statement shows side by side (an effect and its value under H0), with the
/// digits it takes to tell them apart: fnum_fine's when the difference of the values shown is theirs to within 5%,
/// more when it is not (when they show both alike, say), or when either shows as one of the `limits` it is not
/// (ρ1 = −0.4999999956 when ρ0 = −0.49999999, not −0.5 for both and for an ICC's limit −1/(k − 1) = −0.5;
/// μD = 1000000000.66 when μD0 = 1000000000, not 1e+09).  Whole limits (a correlation's ±1) need no mention, as
/// fnum_digits never shows a number as a whole number it is not.
inline std::pair<std::string, std::string> fnum_apart(double a, double b, std::initializer_list<double> limits = {}) {
  auto value = [](const std::string& t) {
    double v;
    return text::parse_double(t, v) ? v : NAN;
  };
  auto at_limit = [&](double x, const std::string& t) {
    for (double l : limits)
      if (x != l && value(t) == l) return true;
    return false;
  };
  std::pair<std::string, std::string> t;
  for (int extra = 0; extra <= 13; ++extra) {
    t = {fnum_digits(a, extra), fnum_digits(b, extra)};
    const double shown = value(t.first) - value(t.second);
    if (std::fabs(shown - (a - b)) <= 0.05 * std::fabs(a - b) && !at_limit(a, t.first) && !at_limit(b, t.second))
      break;
  }
  return t;
}

}  // namespace ssc::procs
