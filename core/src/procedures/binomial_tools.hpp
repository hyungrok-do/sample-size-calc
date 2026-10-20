// Helpers for exact (enumeration-based) power of tests on binomial counts.
#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "stats/stats.hpp"

namespace ssc::procs {

/// Binomial probabilities b(x; n, p) for x in [lo, hi] — the window outside of
/// which the total mass is below `tail` — computed from the mode outwards.
struct BinomWindow {
  int lo = 0, hi = -1;
  std::vector<double> pmf;  // pmf[x - lo]
  double at(int x) const { return (x < lo || x > hi) ? 0.0 : pmf[static_cast<std::size_t>(x - lo)]; }
};

inline BinomWindow binom_window(int n, double p, double tail = 1e-16) {
  using namespace stats;
  BinomWindow w;
  if (n <= 0) {
    w.lo = w.hi = 0;
    w.pmf = {1.0};
    return w;
  }
  if (p <= 0 || p >= 1) {
    w.lo = w.hi = (p <= 0 ? 0 : n);
    w.pmf = {1.0};
    return w;
  }
  const double q = 1 - p;
  int lo = 0, hi = n;
  if (n > 50) {
    lo = static_cast<int>(std::max(0.0, qbinom(tail, n, p) - 1));
    hi = static_cast<int>(std::min<double>(n, qbinom(1 - tail, n, p) + 1));
  }
  w.lo = lo;
  w.hi = hi;
  w.pmf.assign(static_cast<std::size_t>(hi - lo + 1), 0.0);
  int mode = static_cast<int>(std::floor((n + 1) * p));
  if (mode > hi) mode = hi;
  if (mode < lo) mode = lo;
  double v = dbinom_raw(mode, n, p, q);
  w.pmf[static_cast<std::size_t>(mode - lo)] = v;
  const double r = p / q;
  for (int x = mode; x < hi; ++x) {
    v *= static_cast<double>(n - x) / (x + 1) * r;
    w.pmf[static_cast<std::size_t>(x + 1 - lo)] = v;
  }
  v = w.pmf[static_cast<std::size_t>(mode - lo)];
  for (int x = mode; x > lo; --x) {
    v *= static_cast<double>(x) / (n - x + 1) / r;
    w.pmf[static_cast<std::size_t>(x - 1 - lo)] = v;
  }
  return w;
}

/// Smallest c with P(X ≥ c | n, p) ≤ a  (n + 1 if none).
inline int upper_critical(int n, double p, double a) {
  using namespace stats;
  int lo = -1, hi = n + 1;  // P(X ≥ hi) = 0 ≤ a ; find boundary
  // P(X ≥ c) = pbinom(c − 1, upper)
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (pbinom(mid - 1, n, p, false) <= a)
      hi = mid;
    else
      lo = mid;
  }
  return hi;
}

/// Largest c with P(X ≤ c | n, p) ≤ a  (−1 if none).
inline int lower_critical(int n, double p, double a) {
  using namespace stats;
  int lo = -1, hi = n + 1;  // invariant: P(X ≤ lo) ≤ a (lo = −1 trivially), P(X ≤ hi) > a
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (pbinom(mid, n, p, true) <= a)
      lo = mid;
    else
      hi = mid;
  }
  return lo;
}

/// P(X ≤ c) + P(X ≥ d) for X ~ Bin(n, p).
inline double two_tail_prob(int n, double p, int c_low, int c_high) {
  using namespace stats;
  double s = 0;
  if (c_low >= 0) s += pbinom(c_low, n, p, true);
  if (c_high <= n) s += pbinom(c_high - 1, n, p, false);
  return std::min(1.0, s);
}

/// Rough cost of an exact sample-size search for two binomial samples of about
/// n1 and n2: the cells of the joint window (binom_window keeps about 16.5
/// standard deviations of each) times the power evaluations of the saw-tooth
/// search in solve::min_n (a scan over 15% of n below the bisection result).
inline double exact_search_work(double n1, double n2, double p1, double p2) {
  auto width = [](double n, double p) { return 16.5 * std::sqrt(n * p * (1 - p)) + 3; };
  return width(n1, p1) * width(n2, p2) * (std::max(25.0, 0.15 * std::max(n1, n2)) + 40);
}

/// Above this work (a few seconds), exact searches use the normal approximation instead.
constexpr double EXACT_SEARCH_BUDGET = 2.5e8;

}  // namespace ssc::procs
