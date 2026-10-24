// Values that a statement shows side by side, with the digits it takes to tell them apart.
#pragma once

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <string>
#include <utility>
#include <vector>

#include "procedures/common.hpp"

namespace ssc::procs {

/// Whether t shows x, a number that is not whole, as the whole number next to it (0.9999999 as 1, 999999999.9 as
/// 1e+09).
inline bool shows_whole(double x, const std::string& t) {
  double v;
  return x != std::floor(x) && text::parse_double(t, v) && v == std::floor(v) && std::fabs(v - x) < 1;
}

/// x as fnum_fine shows it (`extra` = 0), or with `extra` more digits: decimals where fnum shows a fixed number of
/// them, significant digits where it shows an exponent or a number of 1e6 or more.  Whole numbers are exact, and,
/// as with fnum_fine, a number that is not whole never shows as the whole number next to it (0.9999999 as 1,
/// 999999999.9 as 1e+09): it takes the digits it needs to show a decimal point.
inline std::string fnum_digits(double x, int extra) {
  if (extra == 0 || !std::isfinite(x) || (x == std::floor(x) && std::fabs(x) < 1e15)) return fnum_fine(x);
  const double ax = std::fabs(x);
  std::string t;
  for (int d = 4 + extra; d <= 17; ++d) {
    t = ax >= 1e6 || ax < 1e-4 ? text::format_g(x, d) : blocks::fnum(x, d);
    if (!shows_whole(x, t)) break;
  }
  return t;
}

/// x to the decimal place 10^−q (q < 0 for tens, hundreds, ...) in fnum's forms, fixed below 1e6, and with an
/// exponent and the significant digits that place takes at or above it and below 1e-4, and never with fewer digits
/// than fnum shows (four decimals, or four significant digits); a whole number below 1e15 exactly.  A number that is
/// not whole never shows as the whole number next to it (999999999.9, not 1e+09): it then takes the decimals it needs.
/// `place` is the place shown.
inline std::string fnum_place(double x, int q, int& place) {
  place = q;
  if (!std::isfinite(x) || x == 0 || (x == std::floor(x) && std::fabs(x) < 1e15)) return fnum_fine(x);
  const double ax = std::fabs(x);
  const int e = static_cast<int>(std::floor(std::log10(ax)));
  const bool exponent = ax >= 1e6 || ax < 1e-4;
  auto at = [&](int p) {
    return exponent ? text::format_g(x, std::min(p + e + 1, 17)) : blocks::fnum(x, std::min(p, 17));
  };
  place = exponent ? std::max(q, 3 - e) : std::max(q, 4);
  std::string t = at(place);
  while (shows_whole(x, t) && place < (exponent ? 16 - e : 17)) t = at(++place);
  return t;
}

/// The texts of values a statement shows together (an effect and its value under H0, an interval's bounds and its
/// estimate), with the digits it takes to tell each from every other: fnum_fine's when the differences of the values
/// shown are theirs to within 5%, and otherwise all to the decimal place that keeps them so (the same place for all,
/// whatever their size: μD = 1000000.00018 when μD0 = 999999.9999, not 999999.999900000053, which had decimals below
/// 1e6 and significant digits above it).  None shows as one of the `limits` it is not
/// (ρ1 = −0.4999999956 when ρ0 = −0.49999999, not −0.5 for both and for an ICC's limit −1/(k − 1) = −0.5), nor as
/// the whole number next to it (μD = 1000000000.66 when μD0 = 1000000000, not 1e+09).
inline std::vector<std::string> fnum_apart_list(const std::vector<double>& x,
                                                std::initializer_list<double> limits = {}) {
  auto value = [](const std::string& t) {
    double v;
    return text::parse_double(t, v) ? v : NAN;
  };
  const std::size_t k = x.size();
  std::vector<std::string> t(k);
  auto told = [&] {
    for (std::size_t i = 0; i < k; ++i) {
      if (shows_whole(x[i], t[i])) return false;
      for (double l : limits)
        if (x[i] != l && value(t[i]) == l) return false;
      for (std::size_t j = i + 1; j < k; ++j)
        if (!(std::fabs((value(t[i]) - value(t[j])) - (x[i] - x[j])) <= 0.05 * std::fabs(x[i] - x[j]))) return false;
    }
    return true;
  };
  for (std::size_t i = 0; i < k; ++i) t[i] = fnum_fine(x[i]);
  if (told()) return t;
  // from the place of the first digit of the smallest difference down
  double gap = HUGE_VAL;
  for (std::size_t i = 0; i < k; ++i)
    for (std::size_t j = i + 1; j < k; ++j)
      if (x[i] != x[j]) gap = std::fmin(gap, std::fabs(x[i] - x[j]));
  const int q0 = gap < HUGE_VAL ? -static_cast<int>(std::floor(std::log10(gap))) : 0;
  for (int q = q0 - 1; q <= q0 + 16; ++q) {
    // numbers of 1e6 or more all to the place that one of them takes to show a decimal point, so that they show
    // alike (-999999908.8 next to -999999999.7, not -9.9999991e+08)
    int large = q, p;
    for (double v : x)
      if (std::fabs(v) >= 1e6) {
        fnum_place(v, q, p);
        large = std::max(large, p);
      }
    for (std::size_t i = 0; i < k; ++i) t[i] = fnum_place(x[i], std::fabs(x[i]) >= 1e6 ? large : q, p);
    if (told()) break;
  }
  return t;
}

/// The texts of a and b, two values a statement shows side by side (an effect and its value under H0), told apart
/// as by fnum_apart_list.  Whole limits (a correlation's ±1) need no mention, as a number that is not whole never
/// shows as one.
inline std::pair<std::string, std::string> fnum_apart(double a, double b, std::initializer_list<double> limits = {}) {
  const auto texts = fnum_apart_list({a, b}, limits);
  const std::pair<std::string, std::string> t{texts[0], texts[1]};
  return t;
}

}  // namespace ssc::procs
