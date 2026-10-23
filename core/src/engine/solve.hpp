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

/// The largest difference that two values inside the accepted inputs (at most 1e30 in size) can have: a search for
/// a difference reaches at least this far, so that it finds every answer that an input could state.
constexpr double BOUND_SPAN = 2e30;

/// How far a power at the value reported may lie from the target: within what the table shows (four decimals),
/// and within half the shortfall 1 − target, which statements show to more decimals near 100% (99.9999%).
inline double power_slack(double target) { return std::fmin(5e-5, 0.5 * (1 - target)); }

/// A power in a message, with three digits of its shortfall 1 − p near 1 (0.999998628, not 0.999999).
inline std::string power_digits(double p) {
  return text::format_g(p, p > 0.99 && p < 1 ? 3 + static_cast<int>(std::ceil(-std::log10(1 - p))) : 6);
}

/// x0 + d, the value of an effect found as its distance d from x0, the value of no effect or of the other group
/// (μ1 = μ2 + d).  Numbers as large as x0 lie apart by a step that grows with them, and next to a large x0 the sum
/// keeps few of d's digits or none (μ2 = 1e16 keeps nothing of a difference of 0.9: μ1 rounds back to μ2), so that
/// the answer would stand for another effect than the one whose power is reported: it is refused unless the
/// difference it keeps, x − x0, is d to the four significant digits that statements and tables show.  `what` and
/// `from` name the two values (μ1 and μ2); `shift` says that only their difference matters, so that they can be
/// measured from a reference nearer to them.  check_kept() is the test, for a value x found otherwise.
inline void check_kept(double x0, double d, double x, const std::string& what, const std::string& from,
                       bool shift = false) {
  const double kept = x - x0;
  if (!(std::fabs(kept - d) <= 5e-5 * std::fabs(d))) {
    const double step = std::nextafter(std::fabs(x), HUGE_VAL) - std::fabs(x);
    throw CalcError(what + " cannot be told apart from " + from + " = " + text::format_roundtrip(x0) +
                    ": the difference of " + text::format_g(std::fabs(d), 4) + " it needs " +
                    (kept == 0 ? "is lost to rounding" : "keeps fewer than four significant digits") +
                    " next to a value of that size, where numbers lie " + text::format_g(step, 3) + " apart" +
                    (shift ? "; measure the values from a nearby reference instead (" + from + " = 0, for example)"
                           : ""));
  }
}

inline double offset(double x0, double d, const std::string& what, const std::string& from, bool shift = false) {
  const double x = x0 + d;
  check_kept(x0, d, x, what, from, shift);
  return x;
}

/// offset(), with the power at the value it gives, power(x − x0), on the target to the precision shown: a
/// difference kept to four digits can still move a steep power off it.
template <class P>
double offset(double x0, double d, P&& power, double target, const std::string& what, const std::string& from,
              bool shift = false) {
  const double x = offset(x0, d, what, from, shift);
  if (x - x0 == d) return x;  // the power found for d
  const double pw = power(x - x0);
  if (!(std::fabs(pw - target) <= power_slack(target)))
    throw CalcError(what + " cannot be placed next to " + from + " = " + text::format_roundtrip(x0) +
                    " closely enough: the nearest value that can be stored, " + text::format_roundtrip(x) +
                    ", gives a power of " + power_digits(pw) + " for the target of " + power_digits(target) +
                    (shift ? "; measure the values from a nearby reference instead (" + from + " = 0, for example)"
                           : ""));
  return x;
}

/// The value at which the power reaches the target, for x, a root of power(x) = target found on a finer scale than
/// the doubles next to it (near a correlation of 1, those lie a hundredth of 1e-14 apart, and the power can change
/// by more than 1e-9 from one to the next): the first double, counted from the side of no effect up to `last` (the
/// last admissible value), at which the power no longer falls short (`up`: the power grows towards larger values;
/// x, the root converted to a value, may lie a few doubles to either side of it).  When the power there overshoots
/// the target by more than the precision shown, the double before it, if that one falls short by no more; the
/// answer is refused when even `last` falls short (with `unreachable`), or when neither double meets the target to
/// the precision shown.
template <class P>
double reach(P&& power, double x, double last, bool up, double target, const std::string& what,
             const std::string& unreachable) {
  const double away = up ? -HUGE_VAL : HUGE_VAL, slack = power_slack(target);
  double pw = power(x), before = NAN, pw_before = NAN;  // `before`: the double next to x, short of the target
  int steps = 0;
  for (; pw < target - 1e-9 && x != last && steps < 16; ++steps) {
    before = x;
    pw_before = pw;
    x = std::nextafter(x, last);
    pw = power(x);
  }
  for (; std::isnan(before) && pw > target + 1e-9 && steps < 16; ++steps) {
    const double y = std::nextafter(x, away), py = power(y);
    if (py < target - 1e-9) {
      before = y;
      pw_before = py;
    } else {
      x = y;
      pw = py;
    }
  }
  if (pw < target - 1e-9 && x == last) throw CalcError(unreachable);
  const bool crossed = pw >= target - 1e-9 && (pw <= target + 1e-9 || !std::isnan(before));
  if (!crossed) {  // the root, converted to the value, lay further off than the 16 doubles tried (it should be a few)
    if (std::fabs(pw - target) <= slack) return x;
    throw CalcError("the answer cannot be computed for these inputs; check for an extreme value");
  }
  if (pw - target <= slack) return x;
  if (target - pw_before <= slack) return before;
  throw CalcError("the power cannot be brought to the target at any " + what +
                  " that can be stored: the numbers there lie so far apart that it moves from " +
                  power_digits(pw_before) + " at " + text::format_roundtrip(before) + " to " + power_digits(pw) +
                  " at the next one, " + text::format_roundtrip(x));
}

namespace detail {

/// The distance t from a at which f(a + t) = target, for f increasing; see solve_increasing.
template <class F>
double root_from(F&& f, double target, double a, double b, double b_max, double tol, const std::string& unreachable,
                 const char* what) {
  // the search is for the distance t from x0 = a, so that tol is relative to the effect, not to the value (an ICC of
  // 0.9999996 against one of 0.999999), and does not depend on the units of the inputs: a tolerance taken from the
  // bracket (the SD) left the fourth digit wrong when the root lies far below it (the mean difference for an N of
  // 1e8), and one absolute below 1 left it wrong for an SD of 1e-8
  const double x0 = a;
  // never beyond b_max: x0 + (b_max − x0) can round past it, onto a limit where the power is not defined
  auto g = [&](double t) { return f(std::fmin(x0 + t, b_max)) - target; };
  double ta = 0, tb = b - x0;
  const double t_max = b_max - x0;
  double fa = g(ta);
  if (fa >= 0)  // a is no effect: the power there is the significance level
    throw CalcError("the target power must exceed the power with no effect (the significance level)");
  double fb = g(tb);
  int guard = 0;
  while (fb < 0) {
    if (tb >= t_max || ++guard > 200) throw CalcError(unreachable);
    const double width = tb - ta;
    ta = tb;
    fa = fb;
    tb = std::fmin(t_max, ta + 2 * width + 1e-12 * std::fabs(ta) + 1e-300);
    fb = g(tb);
  }
  if (std::isnan(fa) || std::isnan(fb))  // the power is undefined there (an overflow at extreme inputs)
    throw CalcError(std::string("the power is undefined for some values of ") + what +
                    "; check for an extreme value among the inputs");
  return numerics::brent(g, ta, tb, fa, fb, 0, 300, tol);
}

}  // namespace detail

/// Root of f(x) = target for f increasing on [a, ...); f(a) < target expected.
/// The upper end b is doubled (relative to a) until f(b) ≥ target or b_max; tol is relative to the root's
/// distance from a, the value of no effect, which `from` names (the answer is refused when it cannot be told apart
/// from a, as by offset(), and moved to the next double when the power there falls short, as by reach()).
template <class F>
double solve_increasing(F&& f, double target, double a, double b, double b_max, double tol = 1e-10,
                        const char* what = "the effect size", const char* from = "the value of no effect") {
  b_max = std::fmin(b_max, std::numeric_limits<double>::max());  // 1e9·σ overflows for σ = 1e300
  b = std::fmin(b, b_max);
  const std::string unreachable = std::string("the target power cannot be reached for any admissible value of ") + what;
  const double t = detail::root_from(f, target, a, b, b_max, tol, unreachable, what);
  const double x = offset(a, t, what, from);
  if (x - a == t) return x;  // not rounded (a = 0): the power is the one found
  const double y = reach(f, std::fmin(x, b_max), b_max, true, target, what, unreachable);
  check_kept(a, t, y, what, from);
  return y;
}

/// Root of f(x) = target for f decreasing on [a, b] (e.g. power vs. a quantity
/// that must shrink); expands b downwards towards b_min.
template <class F>
double solve_decreasing(F&& f, double target, double a, double b, double b_min, double tol = 1e-10,
                        const char* what = "the effect size", const char* from = "the value of no effect") {
  const std::string unreachable = std::string("the target power cannot be reached for any admissible value of ") + what;
  // as in solve_increasing, never below b_min
  const double t = detail::root_from([&](double x) { return f(std::fmax(a - x, b_min)); }, target, 0, a - b,
                                     a - b_min, tol, unreachable, what);
  const double x = offset(a, -t, what, from);
  if (x - a == -t) return x;
  const double y = reach(f, std::fmax(x, b_min), b_min, false, target, what, unreachable);
  check_kept(a, -t, y, what, from);
  return y;
}

}  // namespace ssc::solve
