// Shared numerical runtime for compiled study expressions and standalone C++
// output.
#pragma once
#include "sim/random.hpp"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace ssc_study {
using Vec = std::vector<double>;
inline double finite(double x) {
  if (!std::isfinite(x))
    throw std::runtime_error(
        "Non-finite arithmetic result; check the expression's domain.");
  return x;
}
inline int count(double x, int minimum = 0, int maximum = 100000) {
  if (!std::isfinite(x) || x != std::floor(x) || x < minimum || x > maximum)
    throw std::runtime_error("Count is not an integer in the supported range.");
  return static_cast<int>(x);
}
struct Context {
  ssc::sim::Rng rng{1};
  std::uint64_t operations = 0, limit = 100000000;
  std::function<void()> check;
  void tick(std::uint64_t n = 1) {
    operations += n;
    if (operations > limit)
      throw std::runtime_error("Execution budget exceeded. Reduce the search "
                               "range, sample size, or simulation replicates.");
    if (check && (operations % 4096 < n))
      check();
  }
};
inline double checkedSD(double x) {
  if (!(x >= 0))
    throw std::runtime_error("Normal SD must be nonnegative.");
  return finite(x);
}
inline double normalCDF(double x) {
  return 0.5 * std::erfc(-x / std::sqrt(2.0));
}
inline double normalSF(double x) { return 0.5 * std::erfc(x / std::sqrt(2.0)); }
inline double normalQuantile(double p) {
  if (!(p > 0 && p < 1))
    throw std::runtime_error("normalQuantile needs 0 < p < 1.");
  double lo = -40, hi = 40;
  for (int i = 0; i < 110; ++i) {
    double m = lo + (hi - lo) / 2;
    if (p <= 0.5 ? normalCDF(m) < p : normalSF(m) > 1 - p)
      lo = m;
    else
      hi = m;
  }
  return (lo + hi) / 2;
}
inline double normal(Context &c, double mu, double sd) {
  c.tick();
  return finite(mu + checkedSD(sd) * c.rng.normal());
}
inline double uniform(Context &c, double lo, double hi) {
  c.tick();
  if (!(lo < hi))
    throw std::runtime_error("uniform needs lower < upper.");
  return finite(lo + (hi - lo) * c.rng.uniform());
}
inline double bernoulli(Context &c, double p) {
  c.tick();
  if (!(p >= 0 && p <= 1))
    throw std::runtime_error("bernoulli needs a probability in [0,1].");
  return c.rng.uniform() < p ? 1 : 0;
}
inline double exponential(Context &c, double rate) {
  c.tick();
  if (!(rate > 0))
    throw std::runtime_error("exponential needs a positive rate.");
  return finite(c.rng.exponential() / rate);
}
inline double gamma(Context &c, double shape, double scale) {
  c.tick();
  if (!(shape > 0 && scale > 0))
    throw std::runtime_error("gamma needs positive shape and scale.");
  return finite(c.rng.gamma(shape) * scale);
}
inline double binomial(Context &c, double n, double p) {
  c.tick();
  int m = count(n);
  if (!(p >= 0 && p <= 1))
    throw std::runtime_error("binomial needs a probability in [0,1].");
  return double(c.rng.binomial(m, p));
}
inline double natural(double x) {
  count(x, 0, 10000000);
  return x;
}
inline double divide(double x, double y) {
  if (y == 0)
    throw std::runtime_error("Division by zero.");
  return finite(x / y);
}
inline double power(double x, double y) { return finite(std::pow(x, y)); }
inline double squareRoot(double x) {
  if (x < 0)
    throw std::runtime_error("sqrt needs a nonnegative argument.");
  return std::sqrt(x);
}
inline double logarithm(double x) {
  if (!(x > 0))
    throw std::runtime_error("log needs a positive argument.");
  return std::log(x);
}
inline double logistic(double x) {
  return x >= 0 ? 1 / (1 + std::exp(-x)) : std::exp(x) / (1 + std::exp(x));
}
inline double sum(const Vec &a) {
  double s = 0, correction = 0;
  for (double x : a) {
    double y = x - correction, t = s + y;
    correction = (t - s) - y;
    s = t;
  }
  return finite(s);
}
inline double mean(const Vec &a) {
  if (a.empty())
    throw std::runtime_error("mean needs a nonempty array.");
  return sum(a) / a.size();
}
inline double variance(const Vec &a) {
  if (a.size() < 2)
    throw std::runtime_error("variance needs at least two observations.");
  double m = mean(a), s = 0;
  for (double x : a)
    s += (x - m) * (x - m);
  return finite(s / (a.size() - 1));
}
inline double at(const Vec &a, double i) {
  int k = count(i);
  if (std::size_t(k) >= a.size())
    throw std::runtime_error(
        "Array index out of range (indices start at zero).");
  return a[k];
}
inline double dot(const Vec &a, const Vec &b) {
  if (a.size() != b.size())
    throw std::runtime_error("dot requires equal array lengths.");
  double s = 0;
  for (std::size_t i = 0; i < a.size(); ++i)
    s += a[i] * b[i];
  return finite(s);
}
inline double quantile(Vec a, double p) {
  if (a.empty() || !(p >= 0 && p <= 1))
    throw std::runtime_error(
        "quantile needs a nonempty array and a probability.");
  std::sort(a.begin(), a.end());
  double i = p * (a.size() - 1);
  std::size_t k = std::size_t(i);
  return k + 1 < a.size() ? a[k] + (i - k) * (a[k + 1] - a[k]) : a[k];
}
template <class F> Vec sample(Context &c, double n, F f) {
  int m = count(n, 1);
  Vec a;
  a.reserve(m);
  for (int i = 0; i < m; ++i) {
    c.tick();
    a.push_back(finite(f(double(i))));
  }
  return a;
}
template <class F> Vec map(Context &c, const Vec &a, F f) {
  Vec out;
  out.reserve(a.size());
  for (double x : a) {
    c.tick();
    out.push_back(finite(f(x)));
  }
  return out;
}
template <class F>
double integrate(Context &c, F f, double lo, double hi, double panels) {
  int n = count(panels, 2, 100000);
  if (n % 2)
    throw std::runtime_error(
        "integrate needs an even number of Simpson panels.");
  double h = (hi - lo) / n, s = finite(f(lo)) + finite(f(hi));
  for (int i = 1; i < n; ++i) {
    c.tick();
    s += (i % 2 ? 4 : 2) * finite(f(lo + i * h));
  }
  return finite(s * h / 3);
}
inline Vec normalSample(Context &c, double n, double mu, double sd) {
  return sample(c, n, [&](double) { return normal(c, mu, sd); });
}
inline Vec bernoulliSample(Context &c, double n, double p) {
  return sample(c, n, [&](double) { return bernoulli(c, p); });
}
inline Vec exponentialSample(Context &c, double n, double rate) {
  return sample(c, n, [&](double) { return exponential(c, rate); });
}
inline Vec lognormalSample(Context &c, double n, double mu, double sd) {
  return sample(c, n,
                [&](double) { return finite(std::exp(normal(c, mu, sd))); });
}
inline std::string quote(const std::string &s) {
  std::string o = "\"";
  for (unsigned char c : s) {
    if (c == '"' || c == '\\') {
      o += '\\';
      o += c;
    } else if (c == '\n')
      o += "\\n";
    else if (c == '\r')
      o += "\\r";
    else if (c == '\t')
      o += "\\t";
    else if (c < 32)
      o += ' ';
    else
      o += c;
  }
  return o + '"';
}
inline std::string number(double x) {
  if (!std::isfinite(x))
    return "null";
  std::ostringstream s;
  s << std::setprecision(17) << x;
  return s.str();
}
struct Report {
  std::map<std::string, std::string> fields;
  std::vector<std::string> points;
  void put(const std::string &k, const std::string &v) { fields[k] = quote(v); }
  void put(const std::string &k, const char *v) { put(k, std::string(v)); }
  void put(const std::string &k, double v) { fields[k] = number(v); }
  void point(double x, double estimate, double lo, double hi, int reps) {
    if (points.size() < 2000)
      points.push_back("{\"x\":" + number(x) +
                       ",\"estimate\":" + number(estimate) +
                       ",\"lower\":" + number(lo) + ",\"upper\":" + number(hi) +
                       ",\"replicates\":" + std::to_string(reps) + "}");
  }
  std::string json() const {
    std::string s = "{";
    for (const auto &kv : fields) {
      if (s.size() > 1)
        s += ',';
      s += quote(kv.first) + ':' + kv.second;
    }
    s += ",\"trace\":[";
    for (std::size_t i = 0; i < points.size(); ++i) {
      if (i)
        s += ',';
      s += points[i];
    }
    return s + "]}";
  }
};
struct Query {
  bool evaluate = false, integer = false, monteCarlo = false;
  std::string variable = "", comparison = "=",
              evaluation = "deterministic_expression";
  double lo = 0, hi = 1, step = 1, target = 0, tolerance = 1e-8;
  int replicates = 10000;
  std::uint64_t seed = 2026;
  double confidence = 0.95;
  std::vector<double> polynomial;
};
inline bool satisfies(double v, double target, const std::string &op) {
  if (op == ">=")
    return v >= target;
  if (op == ">")
    return v > target;
  if (op == "<=")
    return v <= target;
  if (op == "<")
    return v < target;
  return v == target;
}
using Function = std::function<double(double, Context &)>;
using Event = std::function<bool(double, Context &)>;
inline Report run(const Query &q, const Function &f, const Event &event,
                  Context &ctx) {
  Report r;
  r.put("kind", "program");
  r.put("evaluation", q.evaluation);
  r.put("variable", q.variable);
  r.put("target", q.target);
  r.put("comparison", q.comparison);
  r.put("numerics", "binary64_not_formally_certified");
  if (q.evaluate && !q.monteCarlo) {
    r.put("method", "expression_evaluation");
    r.put("status", "evaluated");
    r.put("value", finite(f(0, ctx)));
    return r;
  }
  if (q.monteCarlo) {
    // A fixed number of replicates at every grid point. Bonferroni + two-sided
    // Hoeffding bounds cover every candidate simultaneously under iid sampling.
    if (!q.evaluate && !q.integer)
      throw std::runtime_error("Simulation solving requires a finite grid: use "
                               "Nat or Grid and a step.");
    int candidates =
        q.evaluate
            ? 1
            : count(std::floor((q.hi - q.lo) / q.step + 1e-10) + 1, 1, 2000);
    if (std::uint64_t(candidates) * q.replicates > 10000000)
      throw std::runtime_error("Simulation exceeds ten million trial "
                               "evaluations; reduce the grid or replicates.");
    double radius = std::sqrt(std::log(2.0 * candidates / (1 - q.confidence)) /
                              (2 * q.replicates));
    double firstPossible = std::numeric_limits<double>::quiet_NaN(),
           firstEstimate = firstPossible;
    r.put("method", "monte_carlo");
    r.put("confidence", q.confidence);
    r.put("seed", double(q.seed));
    r.put("replicatesPerCandidate", double(q.replicates));
    r.put("intervalMethod", "Hoeffding bounds simultaneous over the specified "
                            "finite grid, conditional on iid simulation");
    r.put("scope", "Minimum on the supplied grid. Monte Carlo uncertainty "
                   "excludes model misspecification and numerical error.");
    for (int k = 0; k < candidates; ++k) {
      double x = q.evaluate ? 0 : q.lo + k * q.step;
      int hits = 0;
      for (int j = 0; j < q.replicates; ++j) {
        ctx.tick();
        ctx.rng = ssc::sim::Rng(q.seed, std::uint64_t(k), std::uint64_t(j));
        if (event(x, ctx))
          ++hits;
      }
      double p = double(hits) / q.replicates, lo = std::max(0.0, p - radius),
             hi = std::min(1.0, p + radius);
      r.point(x, p, lo, hi, q.replicates);
      r.put("value", p);
      r.put("lower", lo);
      r.put("upper", hi);
      if (q.evaluate) {
        r.put("status", "estimated");
        return r;
      }
      bool greater = q.comparison == ">=" || q.comparison == ">";
      bool good = satisfies(greater ? lo : hi, q.target, q.comparison);
      bool bad = !satisfies(greater ? hi : lo, q.target, q.comparison);
      if (!bad && !std::isfinite(firstPossible))
        firstPossible = x;
      if (satisfies(p, q.target, q.comparison) && !std::isfinite(firstEstimate))
        firstEstimate = x;
      if (good) {
        r.put("status", firstPossible == x ? "found" : "bounded");
        r.put("solutionLower", firstPossible);
        r.put("solutionUpper", x);
        r.put("estimatedSolution", firstEstimate);
        if (firstPossible == x)
          r.put("solution", x);
        r.put("evaluatedCandidates", double(k + 1));
        return r;
      }
    }
    r.put("status",
          std::isfinite(firstPossible) ? "unresolved" : "not_found_in_domain");
    r.put("solutionLower", firstPossible);
    r.put("estimatedSolution", firstEstimate);
    r.put("evaluatedCandidates", double(candidates));
    return r;
  }
  if (q.integer) {
    int candidates =
        count(std::floor((q.hi - q.lo) / q.step + 1e-10) + 1, 1, 10000000);
    r.put("method", "ascending_grid_search");
    r.put("scope", "First satisfying value on the supplied grid under "
                   "floating-point evaluation.");
    for (int k = 0; k < candidates; ++k) {
      ctx.tick();
      double x = q.lo + k * q.step, v = finite(f(x, ctx));
      if (q.evaluation == "analytic_gaussian_probability" &&
          std::abs(v - q.target) <= 1e-12 * std::max(1.0, std::abs(q.target)) &&
          q.comparison != "=") {
        r.put("status", "numerically_unresolved");
        r.put("candidate", x);
        r.put("value", v);
        return r;
      }
      if (satisfies(v, q.target, q.comparison)) {
        r.put("status", "found");
        r.put("solution", x);
        r.put("value", v);
        r.put("evaluatedCandidates", double(k + 1));
        return r;
      }
    }
    r.put("status", "not_found_in_domain");
    return r;
  }
  r.put("scope",
        "A root within the supplied interval; uniqueness is not assumed.");
  if (!q.polynomial.empty()) {
    double c = q.polynomial[0] - q.target,
           b = q.polynomial.size() > 1 ? q.polynomial[1] : 0,
           a = q.polynomial.size() > 2 ? q.polynomial[2] : 0;
    std::vector<double> roots;
    if (a == 0 && b == 0) {
      r.put("method", "algebraic_constant");
      r.put("status", c == 0 ? "all_values_satisfy" : "not_found_in_domain");
      return r;
    }
    if (a == 0)
      roots.push_back(-c / b);
    else {
      double scale = std::max({std::abs(a), std::abs(b), std::abs(c)});
      a /= scale;
      b /= scale;
      c /= scale;
      double disc = b * b - 4 * a * c;
      if (disc >= 0) {
        double t = -0.5 * (b + std::copysign(std::sqrt(disc), b));
        if (t == 0)
          roots.push_back(-b / (2 * a));
        else {
          roots.push_back(t / a);
          roots.push_back(c / t);
        }
      }
    }
    std::sort(roots.begin(), roots.end());
    for (double x : roots)
      if (std::isfinite(x) && x >= q.lo && x <= q.hi) {
        double value = finite(f(x, ctx)), residual = value - q.target;
        if (std::abs(residual) <=
            q.tolerance * std::max(1.0, std::abs(q.target))) {
          r.put("method", q.polynomial.size() > 2 ? "algebraic_quadratic"
                                                  : "algebraic_linear");
          r.put("status", "found");
          r.put("solution", x);
          r.put("value", value);
          r.put("residual", residual);
          return r;
        }
      }
    // If rounding invalidated an algebraic candidate, the bracketing path still
    // checks the original expression. No root is inferred from a failed
    // formula.
  }
  r.put("method", "bracketed_bisection");
  double lo = q.lo, hi = q.hi, fl = finite(f(lo, ctx)) - q.target,
         fh = finite(f(hi, ctx)) - q.target;
  auto done = [&](double x, double residual) {
    r.put("status", "found");
    r.put("solution", x);
    r.put("value", residual + q.target);
    r.put("residual", residual);
    r.put("bracketLower", lo);
    r.put("bracketUpper", hi);
    return r;
  };
  if (std::abs(fl) <= q.tolerance)
    return done(lo, fl);
  if (std::abs(fh) <= q.tolerance)
    return done(hi, fh);
  if ((fl < 0) == (fh < 0)) {
    r.put("status", "root_not_bracketed");
    r.put("detail", "Use a sign-changing interval. This does not establish "
                    "that no root exists.");
    return r;
  }
  for (int i = 0; i < 160; ++i) {
    ctx.tick();
    double m = lo + (hi - lo) / 2, fm = finite(f(m, ctx)) - q.target;
    if (std::abs(fm) <= q.tolerance)
      return done(m, fm);
    if (m == lo || m == hi)
      break;
    if ((fm < 0) == (fl < 0)) {
      lo = m;
      fl = fm;
    } else {
      hi = m;
      fh = fm;
    }
  }
  r.put("status", "numerically_unresolved");
  r.put("bracketLower", lo);
  r.put("bracketUpper", hi);
  r.put("detail", "The interval narrowed without a small residual; check "
                  "continuity and numerical conditioning.");
  return r;
}
} // namespace ssc_study
