// Generic solvers used by the procedures.
//
//  * min_n():   smallest integer sample size with power ≥ target.  Monotone
//               power functions use galloping + bisection (O(log n) power
//               evaluations, started from a normal-approximation guess);
//               saw-toothed exact-test power uses bisection on the envelope
//               followed by a bounded upward scan for the first crossing
//               (over 15% of n, at most 10,000 sizes).
//  * solve_increasing(): continuous root of a monotone function (Brent),
//               with automatic bracket expansion.
#pragma once

#include <cmath>
#include <limits>
#include <string>

#include "engine/procedure.hpp"
#include "stats/numerics.hpp"
#include "util/numtext.hpp"

namespace ssc::solve {

struct NSearch {
  double lo = 2;            // smallest admissible n
  double hi = 1e8;          // largest admissible n
  double guess = NAN;       // starting point (e.g. normal approximation)
  double step = 1;          // grid spacing (2 for even totals, ...)
  bool monotone = true;     // power non-decreasing in n
  double scan_window = 0;   // non-monotone: width of the upward scan (0 = automatic)
  double period = 0;        // non-monotone: the saw-tooth's period in n, if known (the scan covers five of them)
  const char* goal = "the target power";   // what cannot be reached, for the message
};

// The smallest whole n ≥ from at which ok(n) holds, for a condition that stays true once true (such as group 2,
// the allocation ratio times n rounded up, having at least 2 subjects).  Doubles then halves, so that a tiny
// ratio does not count up one at a time; throws `none` when no n up to `cap` qualifies.
template <class F>
double first_whole(F&& ok, double from, double cap, const char* none) {
  if (from > cap) throw CalcError(none);
  if (ok(from)) return from;
  double lo = from, step = 1;  // ok(lo) is false
  while (!ok(lo + step)) {
    lo += step;
    step *= 2;
    if (lo > cap) throw CalcError(none);
  }
  double hi = lo + step;  // ok(hi) is true
  while (hi - lo > 1) {
    const double mid = std::floor(lo + (hi - lo) / 2);
    (ok(mid) ? hi : lo) = mid;
  }
  if (hi > cap) throw CalcError(none);
  return hi;
}

inline double snap(double n, const NSearch& o) {
  const double k = std::ceil((n - o.lo) / o.step - 1e-9);
  return o.lo + std::fmax(k, 0.0) * o.step;
}

template <class F>
double min_n(F&& power, double target, const NSearch& o, const RunContext* ctx = nullptr) {
  auto P = [&](double n) {
    if (ctx) ctx->check();
    return power(n);
  };
  const double lo0 = snap(o.lo, o);
  if (P(lo0) >= target) return lo0;
  // --- find a bracket lo < hi with P(lo) < target ≤ P(hi) -------------------
  double lo = lo0, hi = lo0;
  bool bracketed = false;
  if (std::isfinite(o.guess) && o.guess > lo0) {
    const double g = snap(std::fmin(o.guess, o.hi), o);
    if (P(g) >= target) {
      hi = g;
      bracketed = true;
      // walk down from the guess with growing steps to tighten `lo`
      double d = std::fmax(o.step, std::ceil(0.02 * g / o.step) * o.step);
      while (true) {
        const double c = snap(std::fmax(lo0, hi - d), o);
        if (c <= lo0) break;
        if (P(c) < target) {
          lo = c;
          break;
        }
        hi = c;
        d *= 2;
      }
    } else {
      lo = g;
    }
  }
  if (!bracketed) {
    double width = std::fmax(o.step, lo - lo0 + o.step);
    while (true) {
      if (lo >= o.hi)
        throw CalcError(std::string(o.goal) + " cannot be reached with a sample size up to " + text::format_g(o.hi, 10));
      hi = snap(std::fmin(o.hi, lo + width), o);
      if (P(hi) >= target) break;
      lo = hi;
      width *= 2;
    }
  }
  // --- bisection on the grid -------------------------------------------------
  while (hi - lo > o.step * 1.5) {
    const double mid = snap(lo + std::floor((hi - lo) / o.step / 2) * o.step, o);
    if (mid <= lo || mid >= hi) break;
    if (P(mid) >= target)
      hi = mid;
    else
      lo = mid;
  }
  if (o.monotone) return hi;
  // --- saw-tooth power: the first n reaching the target may lie below `hi` ---
  // (15% of n below; at most 10,000 sizes, which covers the saw-tooth even for
  // rare events until n is in the hundreds of thousands)
  // rare events make long teeth (period ≈ 1/p0 subjects for a binomial test): at least five periods
  const double window =
      o.scan_window > 0 ? o.scan_window
                        : std::fmax(25 * o.step, std::fmin(std::fmax(std::ceil(0.15 * hi), std::ceil(5 * o.period)),
                                                           10000 * o.step));
  const double start = snap(std::fmax(lo0, hi - window), o);
  for (double n = start; n < hi; n += o.step)
    if (P(n) >= target) return n;
  return hi;
}

/// Smallest n from which ok(n) holds for every larger n, as far as the next sizes show: for interval widths,
/// which dip under the target at a few degenerate small sizes (a bound on the wrong side of the estimate, or
/// one pressed against 0 or 1) before they settle into their decline.
template <class F>
double min_n_staying(F&& ok, NSearch o, const RunContext* ctx = nullptr) {
  o.monotone = true;
  auto P = [&](double n) { return ok(n) ? 1.0 : 0.0; };
  for (;;) {
    const double r = min_n(P, 0.5, o, ctx);
    const double last = std::fmin(o.hi, r + std::fmin(std::fmax(1000 * o.step, r), 10000 * o.step));
    double bad = -1;
    for (double n = snap(r + o.step, o); n <= last; n += o.step) {
      if (ctx) ctx->check();
      if (!ok(n)) bad = n;  // the last size in the window that falls short
    }
    if (bad < 0) return r;
    o.lo = bad + o.step;
  }
}

/// x0 + d, the value of an effect found as its distance d from x0, the value of no effect or of the other group
/// (μ1 = μ2 + d).  Numbers as large as x0 lie apart by a step that grows with them, and next to a large x0 the sum
/// keeps few of d's digits or none (μ2 = 1e16 keeps nothing of a difference of 0.9: μ1 rounds back to μ2), so that
/// the answer would stand for another effect than the one whose power is reported: it is refused unless the step
/// there is at most a millionth of d.  `what` and `from` name the two values (μ1 and μ2); `shift` says that only
/// their difference matters, so that they can be measured from a reference nearer to them.
inline double offset(double x0, double d, const std::string& what, const std::string& from, bool shift = false) {
  const double x = x0 + d, step = std::nextafter(std::fabs(x), HUGE_VAL) - std::fabs(x);
  if (!(std::fabs(d) >= 1e6 * step))
    throw CalcError(what + " cannot be told apart from " + from + " = " + text::format_roundtrip(x0) +
                    ": the difference" + (d != 0 ? " of " + text::format_g(std::fabs(d), 4) : std::string()) +
                    " it needs is lost to rounding at a value of that size" +
                    (shift ? "; measure the values from a nearby reference instead (" + from + " = 0, for example)"
                           : ""));
  return x;
}

/// Root of f(x) = target for f increasing on [a, ...); f(a) < target expected.
/// The upper end b is doubled (relative to a) until f(b) ≥ target or b_max; tol is relative to the root's
/// distance from a, the value of no effect, which `from` names (the answer is refused when it cannot be told apart
/// from a, as by offset()).
template <class F>
double solve_increasing(F&& f, double target, double a, double b, double b_max, double tol = 1e-10,
                        const char* what = "the effect size", const char* from = "the value of no effect") {
  b_max = std::fmin(b_max, std::numeric_limits<double>::max());  // 1e9·σ overflows for σ = 1e300
  b = std::fmin(b, b_max);
  // the search is for the distance t from x0 = a, so that tol is relative to the effect, not to the value (an ICC of
  // 0.9999996 against one of 0.999999), and does not depend on the units of the inputs: a tolerance taken from the
  // bracket (the SD) left the fourth digit wrong when the root lies far below it (the mean difference for an N of
  // 1e8), and one absolute below 1 left it wrong for an SD of 1e-8
  const double x0 = a;
  auto g = [&](double t) { return f(x0 + t) - target; };
  double ta = 0, tb = b - x0;
  const double t_max = b_max - x0;
  double fa = g(ta);
  if (fa >= 0)  // a is no effect: the power there is the significance level
    throw CalcError("the target power must exceed the power with no effect (the significance level)");
  double fb = g(tb);
  int guard = 0;
  while (fb < 0) {
    if (tb >= t_max || ++guard > 200)
      throw CalcError(std::string("the target power cannot be reached for any admissible value of ") + what);
    const double width = tb - ta;
    ta = tb;
    fa = fb;
    tb = std::fmin(t_max, ta + 2 * width + 1e-12 * std::fabs(ta) + 1e-300);
    fb = g(tb);
  }
  if (std::isnan(fa) || std::isnan(fb))  // the power is undefined there (an overflow at extreme inputs)
    throw CalcError(std::string("the power is undefined for some values of ") + what +
                    "; check for an extreme value among the inputs");
  return offset(x0, numerics::brent(g, ta, tb, fa, fb, 0, 300, tol), what, from);
}

/// Root of f(x) = target for f decreasing on [a, b] (e.g. power vs. a quantity
/// that must shrink); expands b downwards towards b_min.
template <class F>
double solve_decreasing(F&& f, double target, double a, double b, double b_min, double tol = 1e-10,
                        const char* what = "the effect size", const char* from = "the value of no effect") {
  const double t = solve_increasing([&](double x) { return f(a - x); }, target, 0, a - b, a - b_min, tol, what);
  return offset(a, -t, what, from);
}

}  // namespace ssc::solve
