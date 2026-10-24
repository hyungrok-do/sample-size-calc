// Monte Carlo power: run many simulated studies in parallel, summarize the
// rejections (with the Monte Carlo standard error and a 95% interval), and
// find the smallest sample size whose simulated power reaches a target.
//
// Studies are reduced in replicate order, so a result depends only on the
// seed, never on the number of threads.  The sample-size search relies on
// common random numbers (see random.hpp): each candidate size reuses the same
// streams, so simulated power changes smoothly from one size to the next.
#pragma once

#include <algorithm>
#include <atomic>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "engine/procedure.hpp"
#include "util/numtext.hpp"

namespace ssc::sim {

/// Outcome of one simulated study.
struct Outcome {
  bool reject = false;
  bool failed = false;   // the analysis failed (counted as not rejecting)
  bool design_failure = false;  // with failed: the data's design gave no estimate (not the arithmetic)
  double estimate = NAN;  // estimated effect
  double se = NAN;        // its standard error
  double df = NAN;        // degrees of freedom of the test
  bool df_fallback = false;  // the test's degrees of freedom could not be computed, and a fixed number stood in
  double extra = NAN;     // a design quantity, e.g. the number of events
  unsigned also = 0;      // other tests of the same data: bit i set when test i rejected
};

constexpr int ALSO_MAX = 16;  // other tests counted per study

/// Summary of a set of simulated studies.
struct Summary {
  long reps = 0;
  long rejections = 0;
  long failures = 0;
  long design_failures = 0;  // failures that the data's design caused (Outcome::design_failure)
  long df_fallbacks = 0;  // studies whose test used the fixed degrees of freedom in place of its own
  double power = NAN;  // rejection rate
  double se = NAN;     // Monte Carlo standard error √(p(1 − p)/R)
  double lo = NAN;     // 95% Wilson interval
  double hi = NAN;
  double mean_estimate = NAN;
  double sd_estimate = NAN;
  double mean_se = NAN;
  double mean_df = NAN;
  double mean_extra = NAN;
  long also[ALSO_MAX] = {};  // rejections of the other tests (bit i of Outcome::also)
  /// Rejection rate of other test i.
  double also_power(int i) const { return reps > 0 ? static_cast<double>(also[i]) / static_cast<double>(reps) : NAN; }
};

/// The summary `s` with its rejection count replaced by `rejections` (the rate, its standard error and
/// interval recomputed), for a power that is not one study's event, such as the smallest of several rates.
inline Summary with_rejections(Summary s, long rejections) {
  s.rejections = rejections;
  if (s.reps == 0) return s;
  const double r = static_cast<double>(s.reps), p = rejections / r;
  s.power = p;
  s.se = std::sqrt(p * (1 - p) / r);
  const double z = 1.959963984540054, z2 = z * z;
  const double center = (p + z2 / (2 * r)) / (1 + z2 / r);
  const double half = z * std::sqrt(p * (1 - p) / r + z2 / (4 * r * r)) / (1 + z2 / r);
  s.lo = std::fmax(0.0, center - half);
  s.hi = std::fmin(1.0, center + half);
  return s;
}

/// Summary of outcomes added one by one, in replicate order (so that the sums, and the result, do not
/// depend on how the replicates were spread over threads).  A study whose analysis failed counts as a failure and
/// as not rejecting; its other tests (Outcome::also) count as on any study, and the averages take the values it
/// holds, so that it leaves out of them only what its analysis did not give (NaN).
class Accumulator {
 public:
  void add(const Outcome& o) {
    ++s_.reps;
    if (o.failed) {
      ++s_.failures;
      if (o.design_failure) ++s_.design_failures;
    } else {
      if (o.reject) ++s_.rejections;
      if (o.df_fallback) ++s_.df_fallbacks;
    }
    for (int i = 0; i < ALSO_MAX; ++i)
      if (o.also >> i & 1u) ++s_.also[i];
    if (std::isfinite(o.estimate)) {
      ++ok_;
      sum_ += o.estimate;
      sum2_ += o.estimate * o.estimate;
      se_ += o.se;
    }
    if (std::isfinite(o.df)) {
      ++dfn_;
      df_ += o.df;
    }
    if (std::isfinite(o.extra)) {
      ++extran_;
      extra_ += o.extra;
    }
  }

  Summary finish() const {
    Summary s = s_;
    if (s.reps == 0) return s;
    if (ok_ > 0) {
      s.mean_estimate = sum_ / ok_;
      s.sd_estimate = ok_ > 1 ? std::sqrt(std::fmax(0.0, (sum2_ - sum_ * sum_ / ok_) / (ok_ - 1))) : NAN;
      s.mean_se = se_ / ok_;
    }
    if (dfn_ > 0) s.mean_df = df_ / dfn_;
    if (extran_ > 0) s.mean_extra = extra_ / extran_;
    return with_rejections(s, s.rejections);
  }

 private:
  Summary s_;
  double se_ = 0, sum_ = 0, sum2_ = 0, df_ = 0, extra_ = 0;
  long ok_ = 0, dfn_ = 0, extran_ = 0;
};

inline Summary summarize(const std::vector<Outcome>& out) {
  Accumulator acc;
  for (const auto& o : out) acc.add(o);
  return acc.finish();
}

/// Replicates simulated at a time: their outcomes are held until they are added to the summary, so the
/// memory a simulation needs does not grow with the number of replicates (16,384 outcomes, about 1 MB).
constexpr long BLOCK = 16384;

/// Simulate `reps` studies; study(r) gives replicate r's outcome and must draw
/// its random numbers only from streams of r.  Progress is reported from p0 to p1.
template <class F>
Summary simulate(long reps, F&& study, const RunContext& ctx, double p0 = 0, double p1 = 1) {
  Accumulator acc;
  std::vector<Outcome> out(static_cast<std::size_t>(std::min(std::max(0L, reps), BLOCK)));
  std::atomic<long> done{0};
  for (long start = 0; start < reps; start += BLOCK) {
    const long m = std::min(BLOCK, reps - start);
    parallel_for(
        static_cast<std::size_t>(m), ctx.threads,
        [&](std::size_t i) {
          out[i] = study(start + static_cast<long>(i));
          const long d = done.fetch_add(1, std::memory_order_relaxed) + 1;
          if ((d & 31) == 0 || d == reps) ctx.report(p0 + (p1 - p0) * static_cast<double>(d) / reps);
        },
        ctx.token);
    for (long i = 0; i < m; ++i) acc.add(out[static_cast<std::size_t>(i)]);
  }
  return acc.finish();
}

// ------------------------------------------------------------------ sample-size search

struct SizeSearch {
  double lo = 2;       // smallest admissible size
  double hi = 1e5;     // largest admissible size
  double step = 1;     // sizes lo, lo + step, ...
  double guess = NAN;  // starting point, e.g. an analytic approximation
  long reps = 1000;    // replicates of the final evaluations
  long explore_reps = 250;  // replicates while bracketing and bisecting
};

struct SizePoint {
  double n;
  Summary s;
};

struct SizeResult {
  double n = NAN;
  Summary s;
  std::vector<SizePoint> evaluated;  // every size simulated with the full replicates
};

/// Progress of a search whose number of evaluations is not known in advance:
/// each evaluation gets a share of the remaining bar.
class SearchProgress {
 public:
  explicit SearchProgress(double budget) : budget_(budget) {}
  /// The part of the bar for an evaluation with `reps` replicates.
  std::pair<double, double> next(double reps) {
    const double a = pos_;
    used_ += reps;
    double b = used_ / budget_;
    if (b > 0.95) b = a + (1 - a) / 2;  // past the budget: halve the remaining distance
    pos_ = b;
    return {a, b};
  }

 private:
  double budget_;
  double used_ = 0;
  double pos_ = 0;
};

/// Smallest size n (on the grid) whose simulated power reaches `target`.
/// power_at(n, reps, p0, p1) simulates `reps` studies of size n with the same
/// seed for every n.  Exploratory evaluations use o.explore_reps replicates; the
/// answer is then bracketed again and bisected with o.reps replicates, so that
/// it reaches the target and one step less does not (with the full replicates).
template <class F>
SizeResult search_size(F&& power_at, double target, const SizeSearch& o, const RunContext& ctx) {
  auto snap = [&](double n) {
    const double k = std::ceil((n - o.lo) / o.step - 1e-9);
    return std::fmin(o.hi, o.lo + std::fmax(k, 0.0) * o.step);
  };
  std::map<std::pair<double, long>, Summary> memo;
  SearchProgress prog(12.0 * o.explore_reps + 4.0 * o.reps);
  SizeResult res;
  auto P = [&](double n, long reps) -> const Summary& {
    auto key = std::make_pair(n, reps);
    auto it = memo.find(key);
    if (it != memo.end()) return it->second;
    ctx.check();
    const auto span = prog.next(static_cast<double>(reps));
    Summary s = power_at(n, reps, span.first, span.second);
    if (reps == o.reps) res.evaluated.push_back({n, s});
    return memo.emplace(key, s).first->second;
  };
  const long e = std::min(o.explore_reps, o.reps);
  auto reaches = [&](double n, long reps) { return P(n, reps).power >= target; };

  // --- bracket with exploratory replicates --------------------------------
  double lo = snap(o.lo), hi;
  const double g = snap(std::isfinite(o.guess) ? o.guess : o.lo);
  if (reaches(g, e)) {
    hi = g;
    double d = std::fmax(o.step, 0.1 * g);
    while (hi > lo) {
      const double c = snap(std::fmax(o.lo, hi - d));
      if (c >= hi) break;
      if (reaches(c, e)) {
        hi = c;
        d *= 2;
      } else {
        lo = c;
        break;
      }
    }
  } else {
    lo = g;
    double d = std::fmax(o.step, 0.25 * g);
    while (true) {
      if (lo >= o.hi)
        throw CalcError("the target power is not reached with up to " + text::format_g(o.hi, 10) +
                        " (simulated power " + text::format_g(100 * P(lo, e).power, 3) + "%)");
      const double c = snap(lo + d);
      if (reaches(c, e)) {
        hi = c;
        break;
      }
      lo = c;
      d *= 2;
    }
  }
  // --- bisection, until the bracket is as narrow as the first step of the confirmation below ---
  auto coarse = [&](double n) { return std::fmax(o.step, std::floor(0.02 * n / o.step) * o.step); };
  while (hi - lo > std::fmax(1.5 * o.step, coarse(hi))) {
    const double mid = snap(lo + std::floor((hi - lo) / o.step / 2) * o.step);
    if (mid <= lo || mid >= hi) break;
    (reaches(mid, e) ? hi : lo) = mid;
  }
  // --- confirm with the full replicates: bracket from the exploratory answer, then bisect ---
  // (with common random numbers the simulated power is close to monotone in the size)
  double yes = hi, no = NAN;  // yes reaches the target with the full replicates, no does not
  const double first = coarse(hi);  // the exploratory answer is close
  if (reaches(yes, o.reps)) {
    for (double d = first;; d *= 2) {
      const double c = snap(std::fmax(o.lo, yes - d));
      if (c >= yes) break;  // the smallest admissible size reaches it
      if (!reaches(c, o.reps)) {
        no = c;
        break;
      }
      yes = c;
    }
  } else {
    no = yes;
    for (double d = first;; d *= 2) {
      if (no >= o.hi) throw CalcError("the target power is not reached with up to " + text::format_g(o.hi, 10));
      const double c = snap(no + d);
      if (reaches(c, o.reps)) {
        yes = c;
        break;
      }
      no = c;
    }
  }
  while (std::isfinite(no) && yes - no > 1.5 * o.step) {
    const double mid = snap(no + std::floor((yes - no) / o.step / 2) * o.step);
    if (mid <= no || mid >= yes) break;
    (reaches(mid, o.reps) ? yes : no) = mid;
  }
  const double n = yes;
  res.n = n;
  res.s = P(n, o.reps);
  std::sort(res.evaluated.begin(), res.evaluated.end(), [](const SizePoint& a, const SizePoint& b) {
    return a.n < b.n;
  });
  ctx.report(1);
  return res;
}

}  // namespace ssc::sim
