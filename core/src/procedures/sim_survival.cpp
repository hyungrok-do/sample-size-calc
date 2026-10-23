// Simulation-based power for comparing two survival curves whose hazards are
// not proportional (sim.surv.nph): a treatment effect that starts late (as
// often with immunotherapy), wanes, or reverses.  The control group has
// exponential survival; the treatment group's hazard ratio changes from HR1 to
// HR2 at a time T* after entry.  Subjects enter uniformly during accrual and
// are followed until the analysis or loss to follow-up, as in surv.logrank.
//
// Each simulated study is analyzed with one of four tests: the logrank test, a
// Fleming–Harrington weighted logrank test, the MaxCombo test (the largest of
// four FH statistics, with a multivariate normal p-value) or the difference in
// restricted mean survival time.  The quick approximation used for live
// previews, power curves and the start of the sample-size search integrates
// the expected events over time (Schoenfeld's formula extended to a hazard
// ratio that changes); with proportional hazards and the logrank test it equals
// surv.logrank.
//
// Sign convention: every statistic is oriented like the log hazard ratio, so
// negative values favor the treatment (fewer events than expected in group 1,
// or a longer restricted mean survival time) and the one-sided alternative
// "HR < 1" rejects for small values.
#include "procedures/sim_survival.hpp"

#include <algorithm>

#include "procedures/sim_common.hpp"
#include "sim/random.hpp"
#include "stats/numerics.hpp"

namespace ssc::procs {
namespace {

/// x^e, exact and fast for the exponents 0 and 1 of the usual weights.
double power_of(double x, double e) { return e == 0 ? 1 : (e == 1 ? x : std::pow(x, e)); }

}  // namespace

// =====================================================================================
// Test statistics of one data set and the multivariate normal probability
// =====================================================================================
namespace nph {
namespace {

// Rank-1 lattice rules for the separated integrand: in one or two dimensions the Fibonacci lattice with
// 2584 points (optimal in two), in three the Korobov lattice with 2039 points and the generator
// (1, 653, 653² mod 2039), the best for the P2 criterion.  With the baker's (tent) transformation and a
// fixed shift their error is about 1e-5 (checked against mvtnorm and Fibonacci lattices of 800,000 points).
struct Lattice {
  long n;
  long z[3];
};
constexpr Lattice LATTICE_2D = {2584, {1, 1597, 0}};
constexpr Lattice LATTICE_3D = {2039, {1, 653, 653L * 653 % 2039}};
constexpr double LATTICE_SHIFT[3] = {0.41421356237309515, 0.7320508075688772, 0.2360679774997898};  // √2, √3, √5

}  // namespace

WeightedLogrank weighted_logrank(const std::vector<Obs>& sorted, long n1, long n2, int k, const double* rho,
                                 const double* gamma) {
  WeightedLogrank res;
  res.k = k;
  double at[2] = {static_cast<double>(n1), static_cast<double>(n2)};  // at risk
  double surv = 1;  // pooled Kaplan–Meier estimate just before the current time
  const std::size_t n = sorted.size();
  for (std::size_t i = 0; i < n;) {
    const double t = sorted[i].time;
    double d[2] = {0, 0}, c[2] = {0, 0};  // events and censorings at t (censored subjects are at risk at t)
    for (; i < n && sorted[i].time == t; ++i) (sorted[i].event ? d : c)[sorted[i].group] += 1;
    const double events = d[0] + d[1], total = at[0] + at[1];
    if (events > 0) {
      const double o = d[0] - events * at[0] / total;  // observed minus expected in group 1
      const double v = total > 1 ? events * at[0] * at[1] / (total * total) * (total - events) / (total - 1) : 0;
      double w[4];
      for (int j = 0; j < k; ++j) w[j] = power_of(surv, rho[j]) * power_of(1 - surv, gamma[j]);
      for (int j = 0; j < k; ++j) {
        res.u[j] += w[j] * o;
        for (int l = 0; l <= j; ++l) res.v[j * k + l] += w[j] * w[l] * v;
      }
      surv *= 1 - events / total;
    }
    at[0] -= d[0] + c[0];
    at[1] -= d[1] + c[1];
  }
  for (int j = 0; j < k; ++j)
    for (int l = 0; l < j; ++l) res.v[l * k + j] = res.v[j * k + l];
  return res;
}

Rmst restricted_mean(const std::vector<Obs>& sorted, long n1, long n2, double tau) {
  // one pass: Σ A_j² q_j = R² Σ q_j − 2R Σ q_j a_j + Σ q_j a_j², with a_j the area up to the j-th event time,
  // R the area up to tau and q_j = d_j / (n_j (n_j − d_j)) (zero when every subject at risk has the event)
  double at[2] = {static_cast<double>(n1), static_cast<double>(n2)}, surv[2] = {1, 1}, last[2] = {0, 0};
  double area[2] = {0, 0}, q0[2] = {0, 0}, q1[2] = {0, 0}, q2[2] = {0, 0};
  const std::size_t n = sorted.size();
  for (std::size_t i = 0; i < n && sorted[i].time <= tau;) {
    const double t = sorted[i].time;
    double d[2] = {0, 0}, c[2] = {0, 0};
    for (; i < n && sorted[i].time == t; ++i) (sorted[i].event ? d : c)[sorted[i].group] += 1;
    for (int g = 0; g < 2; ++g) {
      if (d[g] > 0) {
        area[g] += surv[g] * (t - last[g]);
        last[g] = t;
        if (at[g] > d[g]) {
          const double q = d[g] / (at[g] * (at[g] - d[g]));
          q0[g] += q;
          q1[g] += q * area[g];
          q2[g] += q * area[g] * area[g];
        }
        surv[g] *= 1 - d[g] / at[g];
      }
      at[g] -= d[g] + c[g];
    }
  }
  Rmst res;
  for (int g = 0; g < 2; ++g) {
    const double r = area[g] + surv[g] * (tau - last[g]);  // the curve stays flat after the last observation
    res.rmst[g] = r;
    res.var[g] = std::fmax(0.0, r * r * q0[g] - 2 * r * q1[g] + q2[g]);
  }
  return res;
}

double mvn_probability(int m, const double* a, const double* b, const double* r, int threads) {
  if (m == 1) return stats::clamp01(stats::pnorm(b[0]) - stats::pnorm(a[0]));
  // Cholesky factor; a variable that (nearly) depends linearly on the previous ones gets a zero pivot and
  // its limits become a condition on that combination
  double c[16] = {};
  for (int i = 0; i < m; ++i)
    for (int j = 0; j <= i; ++j) {
      double s = r[i * m + j];
      for (int l = 0; l < j; ++l) s -= c[i * 4 + l] * c[j * 4 + l];
      if (i == j)
        c[i * 4 + i] = s > 1e-10 ? std::sqrt(s) : 0;
      else
        c[i * 4 + j] = c[j * 4 + j] > 0 ? s / c[j * 4 + j] : 0;
    }
  // a dependent variable's limits restrict the last free variable it involves (MaxCombo's matrix always has
  // one: FH(1, 0) = FH(0, 0) − FH(0, 1) in the scores), which keeps the integrand smooth
  int fold[4] = {-1, -1, -1, -1}, last = -1, nfree = 0;
  for (int i = 0; i < m; ++i) {
    if (c[i * 4 + i] > 0) {
      last = i;
      ++nfree;
      continue;
    }
    for (int j = 0; j < i; ++j)
      if (c[j * 4 + j] > 0 && std::fabs(c[i * 4 + j]) > 1e-12) fold[i] = j;
    if (fold[i] < 0 && (a[i] > 0 || b[i] < 0)) return 0;  // a degenerate variable fixed at 0
  }
  // separation of variables (Genz 1992): X = C y with y standard normal, each free y_i drawn within the
  // limits left by the previous ones; the integrand over the unit cube is the product of their probabilities
  auto interval = [&](int i, const double* y, double& plo, double& phi, bool& upper) {
    double s = 0;
    for (int j = 0; j < i; ++j) s += c[i * 4 + j] * y[j];
    double lo = (a[i] - s) / c[i * 4 + i], hi = (b[i] - s) / c[i * 4 + i];
    for (int q = i + 1; q < m; ++q) {
      if (fold[q] != i) continue;
      double sq = 0;
      for (int j = 0; j < i; ++j) sq += c[q * 4 + j] * y[j];
      const double cq = c[q * 4 + i], u = (a[q] - sq) / cq, v = (b[q] - sq) / cq;
      lo = std::fmax(lo, std::fmin(u, v));
      hi = std::fmin(hi, std::fmax(u, v));
    }
    if (!(hi > lo)) return false;
    upper = lo > 0;  // then from the upper tail, without cancellation near 1
    plo = stats::pnorm(upper ? -hi : lo);
    phi = stats::pnorm(upper ? -lo : hi);
    return true;
  };
  double plo0 = 0, phi0 = 0;  // the first variable's interval is the same at every point
  bool upper0 = false;
  const double origin[4] = {0, 0, 0, 0};
  if (!interval(0, origin, plo0, phi0, upper0)) return 0;
  const Lattice& rule = nfree > 3 ? LATTICE_3D : LATTICE_2D;
  const long points = nfree > 1 ? rule.n : 1;  // with one free variable the integrand is constant
  std::vector<double> value(static_cast<std::size_t>(points));  // the points are shared between threads
  parallel_for(value.size(), threads, [&](std::size_t pt) {
    const long p = static_cast<long>(pt);
    double f = 1, y[4] = {0, 0, 0, 0};
    for (int i = 0, dim = 0; i <= last && f > 0; ++i) {
      if (!(c[i * 4 + i] > 0)) continue;
      double plo = plo0, phi = phi0;
      bool upper = upper0;
      if (i > 0 && !interval(i, y, plo, phi, upper)) {
        f = 0;
        break;
      }
      f *= phi - plo;
      if (i < last && f > 0) {
        double w = static_cast<double>(p * rule.z[dim] % rule.n) / static_cast<double>(rule.n) + LATTICE_SHIFT[dim];
        ++dim;
        w -= std::floor(w);
        w = 1 - std::fabs(2 * w - 1);  // baker's transformation
        const double x = stats::qnorm(std::fmin(std::fmax(plo + w * (phi - plo), 1e-300), 1 - 1e-16));
        y[i] = upper ? -x : x;
      }
    }
    value[pt] = f;
  });
  double total = 0;  // added in order: the same sum on any number of threads
  for (double f : value) total += f;
  // a NaN (from a NaN limit or correlation) stays NaN: std::fmin(1.0, x) would make it a probability of 1, the
  // co-primary endpoints' joint power 100% and MaxCombo's power 0
  return stats::clamp01(total / static_cast<double>(points));
}

double mvn_outside(int m, const double* a, const double* b, const double* r, int threads) {
  // X_0 outside its limits, then for each later X_j the probability that it is outside while X_0..X_(j−1) are
  // inside, with X_j first in the separation of variables
  double total = (b[0] < stats::INF ? stats::pnorm(-b[0]) : 0.0) + (a[0] > -stats::INF ? stats::pnorm(a[0]) : 0.0);
  for (int j = 1; j < m; ++j) {
    const int dim = j + 1;
    int order[4] = {j, 0, 1, 2};
    double rr[16], lo[4], hi[4];
    for (int u = 0; u < dim; ++u) {
      for (int v = 0; v < dim; ++v) rr[u * dim + v] = r[order[u] * m + order[v]];
      lo[u] = a[order[u]];
      hi[u] = b[order[u]];
    }
    if (b[j] < stats::INF) {
      lo[0] = b[j];
      hi[0] = stats::INF;
      total += mvn_probability(dim, lo, hi, rr, threads);
    }
    if (a[j] > -stats::INF) {
      lo[0] = -stats::INF;
      hi[0] = a[j];
      total += mvn_probability(dim, lo, hi, rr, threads);
    }
  }
  return stats::clamp01(total);
}

}  // namespace nph

namespace {

using namespace blocks;
using sim::Rng;
using simc::p_value;

enum class Test { Logrank, FH, MaxCombo, Rmst };

// FH(ρ, γ) weights of the MaxCombo test
constexpr double COMBO_RHO[4] = {0, 0, 1, 1};
constexpr double COMBO_GAMMA[4] = {0, 1, 0, 1};

// Below this probability of a MaxCombo statistic leaving its box (a level, a p-value or a power), it is computed
// as such (nph::mvn_outside) rather than as 1 − P(box), whose lattice error of about 1e-5 is then more than a
// thousandth of it: at 1e-4, up to a third (four statistics correlated 0.95), and below about 1e-16 all of it, so
// that the critical value was a single test's and the power of a study 0.
constexpr double COMBO_TAIL = 0.01;

// Largest total sample size N1 + N2 of the sample-size search.
constexpr double MAX_N = 2e5;

struct Design {
  double lam2 = 1;                  // control hazard ln 2 / median
  double hr1 = 1, hr2 = 1, tc = 0;  // treatment hazard ratio before and from tc on
  double A = 0, F = 0, eta = 0;     // accrual, follow-up after accrual, loss-to-follow-up hazard
  Test test = Test::Logrank;
  double rho = 0, gamma = 0;  // FH weights (logrank: 0, 0)
  double tau = 0;             // RMST truncation time

  /// Hazard of group g (0: treatment, 1: control) at time t after entry.
  double hazard(int g, double t) const { return g == 1 ? lam2 : lam2 * (t < tc ? hr1 : hr2); }
  double cumhaz(int g, double t) const {
    if (g == 1) return lam2 * t;
    return t < tc ? lam2 * hr1 * t : lam2 * (hr1 * tc + hr2 * (t - tc));
  }
  double surv(int g, double t) const { return std::exp(-cumhaz(g, t)); }
  /// Event time whose cumulative hazard is e (a unit exponential variate).
  double event_time(int g, double e) const {
    if (g == 1) return e / lam2;
    const double at_tc = lam2 * hr1 * tc;
    return e < at_tc ? e / (lam2 * hr1) : tc + (e - at_tc) / (lam2 * hr2);
  }
  /// ∫_{t0}^{t1} S_g(u) du / S_g(t0), in closed form (piecewise exponential) and without underflow.
  double residual(int g, double t0, double t1) const {
    if (!(t1 > t0)) return 0;
    if (g == 0 && t0 < tc && t1 > tc) {
      const double h = lam2 * hr1;
      return -std::expm1(-h * (tc - t0)) / h + std::exp(-h * (tc - t0)) * residual(g, tc, t1);
    }
    const double h = hazard(g, t0);
    return -std::expm1(-h * (t1 - t0)) / h;
  }
  /// G(t): probability that a subject is still under observation t after entry (not lost, not at the analysis).
  double observed(double t) const {
    const double kept = eta > 0 ? std::exp(-eta * t) : 1;
    if (t <= F) return kept;
    if (t >= A + F) return 0;
    return kept * (A + F - t) / A;
  }
};

Design nph_design(const Scenario& s) {
  Design d;
  d.lam2 = std::log(2.0) / s("med2");
  if (!std::isfinite(d.lam2)) throw CalcError("the median survival in group 2 is too short to compute with");
  d.hr1 = s("hr1");
  d.hr2 = s("hr2");
  d.tc = s("tchange");
  d.A = s("accrual");
  d.F = s("follow");
  d.eta = s("eta");
  if (!(d.A + d.F > 0)) throw CalcError("the accrual and follow-up times cannot both be zero");
  if (s.is("test", "fh")) {
    d.test = Test::FH;
    d.rho = s("rho");
    d.gamma = s("gamma");
  } else if (s.is("test", "maxcombo")) {
    d.test = Test::MaxCombo;
  } else if (s.is("test", "rmst")) {
    d.test = Test::Rmst;
    d.tau = s("tau");
    if (d.tau > d.A + d.F)
      throw CalcError("τ must not exceed the study duration (accrual + follow-up = " + fnum(d.A + d.F) + ")");
  }
  return d;
}

/// Calls f(t, w) at the nodes of a composite Gauss–Legendre rule on [0, end], split where the integrands
/// have kinks (at T* and at the end of follow-up of the last subjects to enter) and in pieces that double in
/// length from an eighth of the shortest half-life (of the events or of the loss to follow-up), so that fast
/// hazards are integrated as accurately as slow ones.
template <class Fn>
void over_time(const Design& d, double end, Fn&& f) {
  std::vector<double> cuts = {0, end};
  for (double c : {d.tc, d.F})
    if (c > 0 && c < end) cuts.push_back(c);
  const double fastest = d.lam2 * std::fmax(1.0, std::fmax(d.hr1, d.hr2)) + d.eta;
  const double first = std::log(2.0) / fastest / 8;
  if (first > 0 && std::isfinite(first))  // a doubling from 0 would never reach the end
    for (double t = first; t < end; t *= 2) cuts.push_back(t);
  std::sort(cuts.begin(), cuts.end());
  cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
  const stats::GLRule& gl = stats::gauss_legendre(8);
  const int parts = 4;
  for (std::size_t k = 0; k + 1 < cuts.size(); ++k) {
    const double h = (cuts[k + 1] - cuts[k]) / parts;
    for (int q = 0; q < parts; ++q) {
      const double mid = cuts[k] + (q + 0.5) * h;
      for (std::size_t i = 0; i < gl.x.size(); ++i) f(mid + 0.5 * h * gl.x[i], 0.5 * h * gl.w[i]);
    }
  }
}

/// Expected events per subject (allocation p = N1/N) over the study, their mean log hazard ratio, and the
/// drift E[Z]/√N of k weighted logrank statistics: √(p(1 − p)) ∫ w log HR dD / √(∫ w² dD), dD the expected
/// events and w the weight at the expected pooled survival p·S1 + (1 − p)·S2 (Schoenfeld's approximation
/// over fine time intervals).
struct Expected {
  double events = 0;
  double log_hr = 0;
  double drift[4] = {0, 0, 0, 0};
  double cross[16] = {};  // ∫ w_j w_l dD, row-major k × k: the statistics' covariance up to p(1 − p)
};

Expected expected(const Design& d, double p, int k, const double* rho, const double* gamma) {
  const double l1 = std::log(d.hr1), l2 = std::log(d.hr2);
  double events = 0, lhr = 0, num[4] = {0, 0, 0, 0};
  Expected e;
  over_time(d, d.A + d.F, [&](double t, double w) {
    const double s1 = d.surv(0, t), s2 = d.surv(1, t);
    const double dd = (p * d.hazard(0, t) * s1 + (1 - p) * d.hazard(1, t) * s2) * d.observed(t) * w;
    const double l = t < d.tc ? l1 : l2, pooled = p * s1 + (1 - p) * s2;
    events += dd;
    lhr += l * dd;
    double wt[4];
    for (int j = 0; j < k; ++j) {
      wt[j] = power_of(pooled, rho[j]) * power_of(1 - pooled, gamma[j]);
      num[j] += wt[j] * l * dd;
      for (int i = 0; i <= j; ++i) e.cross[j * 4 + i] += wt[j] * wt[i] * dd;
    }
  });
  e.events = events;
  e.log_hr = events > 0 ? lhr / events : 0;
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < j; ++i) e.cross[i * 4 + j] = e.cross[j * 4 + i];
    const double den = e.cross[j * 4 + j];
    e.drift[j] = den > 0 ? std::sqrt(p * (1 - p)) * num[j] / std::sqrt(den) : 0;
  }
  return e;
}

/// The RMST difference RMST2 − RMST1 up to τ (oriented like the log hazard ratio) and the asymptotic variance
/// per subject of each group's estimate, ∫_0^τ (∫_t^τ S(u) du)² λ(t) / (S(t) G(t)) dt.
struct RmstExpected {
  double diff = 0;
  double var[2] = {0, 0};
};

RmstExpected rmst_expected(const Design& d) {
  RmstExpected r;
  r.diff = d.residual(1, 0, d.tau) - d.residual(0, 0, d.tau);
  over_time(d, d.tau, [&](double t, double w) {
    const double g = d.observed(t);
    if (!(g > 0)) return;
    for (int k = 0; k < 2; ++k) {
      const double res = d.residual(k, t, d.tau);
      r.var[k] += res * res * d.surv(k, t) * d.hazard(k, t) / g * w;
    }
  });
  return r;
}

/// The asymptotic null distribution of MaxCombo's four statistics, Z ~ N(0, R) with R from the expected
/// events and weights, and the critical value c of their maximum: P(max T_j ≥ c) = α, T_j = |Z_j|
/// two-sided and ±Z_j one-sided (the same c for both directions).
struct ComboNull {
  double r[16] = {};
  double c = 0;
};

ComboNull combo_null(const Design& d, double p, double alpha, Alt a, int threads) {
  // the same design comes back for every point of a power curve: keep the last few
  struct Cached {
    double key[10];
    ComboNull value;
  };
  thread_local std::vector<Cached> cache;
  const double key[10] = {d.lam2, d.hr1, d.hr2, d.tc, d.A, d.F, d.eta, p, alpha, static_cast<double>(a)};
  for (const Cached& c : cache)
    if (std::equal(key, key + 10, c.key)) return c.value;
  const Expected e = expected(d, p, 4, COMBO_RHO, COMBO_GAMMA);
  ComboNull out;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      const double den = std::sqrt(e.cross[i * 5] * e.cross[j * 5]);
      out.r[i * 4 + j] = den > 0 ? std::clamp(e.cross[i * 4 + j] / den, -1.0, 1.0) : (i == j ? 1 : 0);
    }
  auto exceed = [&](double c) {
    double lo[4], hi[4];
    for (int j = 0; j < 4; ++j) {
      lo[j] = a == Alt::Two ? -c : -stats::INF;
      hi[j] = c;
    }
    return alpha < COMBO_TAIL ? nph::mvn_outside(4, lo, hi, out.r, threads)
                              : 1 - nph::mvn_probability(4, lo, hi, out.r, threads);
  };
  // between the critical values of a single test and of Bonferroni's correction for four
  const double q = tail_alpha(alpha, a);
  const double lo = stats::qnorm(q, false), hi = stats::qnorm(q / 4, false);
  const double flo = exceed(lo) - alpha, fhi = exceed(hi) - alpha;
  out.c = flo <= 0 ? lo : (fhi >= 0 ? hi : numerics::brent([&](double c) { return exceed(c) - alpha; }, lo, hi, flo,
                                                             fhi, 1e-8));
  Cached entry;
  std::copy(key, key + 10, entry.key);
  entry.value = out;
  if (cache.size() >= 8) cache.erase(cache.begin());
  cache.push_back(entry);
  return out;
}

double approx_power(const Design& d, const RmstExpected& rm, const ComboNull& combo, double n1, double n2,
                    double alpha, Alt a, int threads) {
  if (d.test == Test::Rmst) {
    const double se = std::sqrt(rm.var[0] / n1 + rm.var[1] / n2);
    return z_power(se > 0 ? rm.diff / se : 0, alpha, a);
  }
  const double root = std::sqrt(n1 + n2);
  if (d.test == Test::MaxCombo) {  // P(max T_j ≥ c) with Z ~ N(√N · drift, R)
    const Expected e = expected(d, n1 / (n1 + n2), 4, COMBO_RHO, COMBO_GAMMA);
    double lo[4], hi[4];
    for (int j = 0; j < 4; ++j) {
      const double mu = root * e.drift[j];
      lo[j] = a == Alt::Upper ? -stats::INF : -combo.c - mu;
      hi[j] = a == Alt::Lower ? stats::INF : combo.c - mu;
    }
    const double p = 1 - nph::mvn_probability(4, lo, hi, combo.r, threads);
    return p < COMBO_TAIL ? nph::mvn_outside(4, lo, hi, combo.r, threads) : p;
  }
  const Expected e = expected(d, n1 / (n1 + n2), 1, &d.rho, &d.gamma);
  return z_power(root * e.drift[0], alpha, a);
}

/// MaxCombo decision: the largest standardized statistic T = max T_i (T_i = |Z_i| two-sided, −Z_i for HR < 1,
/// Z_i for HR > 1) against its null distribution, with Z ~ N(0, R) and R estimated from the data.  The p-value
/// p = P(T ≥ t) is bracketed first by p_min ≤ p ≤ K·p_min (p_min = P(T_i ≥ t)), then by bivariate bounds: the
/// largest union of two events below, Hunter's (1976) spanning-tree bound above.  Only the studies left
/// undecided (about 1 in 50) need the multivariate normal probability.
bool maxcombo_rejects(const nph::WeightedLogrank& w, double alpha, Alt a) {
  int idx[4], k = 0;
  double z[4];
  for (int j = 0; j < 4; ++j)
    if (w.v[j * 5] > 0) {
      idx[k] = j;
      z[k++] = w.u[j] / std::sqrt(w.v[j * 5]);
    }
  if (k == 0) return false;  // no information (no events)
  double t = -stats::INF;
  for (int j = 0; j < k; ++j) t = std::fmax(t, a == Alt::Two ? std::fabs(z[j]) : (a == Alt::Lower ? -z[j] : z[j]));
  const double pmin = a == Alt::Two ? 2 * stats::pnorm(t, false) : stats::pnorm(t, false);
  if (pmin >= alpha) return false;
  if (k * pmin < alpha) return true;
  double r[16], both[16];
  for (int i = 0; i < k; ++i)
    for (int j = 0; j < k; ++j) {
      const double rij = w.v[idx[i] * 4 + idx[j]] / std::sqrt(w.v[idx[i] * 5] * w.v[idx[j] * 5]);
      r[i * k + j] = std::clamp(rij, -1.0, 1.0);  // a NaN stays NaN, not a correlation of 1
      if (i < j)  // P(T_i ≥ t, T_j ≥ t)
        both[i * k + j] = both[j * k + i] = a == Alt::Two ? 2 * (stats::bvn_upper(t, t, r[i * k + j]) +
                                                                  stats::bvn_upper(t, t, -r[i * k + j]))
                                                          : stats::bvn_upper(t, t, r[i * k + j]);
    }
  double lower = pmin, tree = 0;
  bool in[4] = {true, false, false, false};
  for (int i = 0; i < k; ++i)
    for (int j = i + 1; j < k; ++j) lower = std::fmax(lower, 2 * pmin - both[i * k + j]);
  for (int step = 1; step < k; ++step) {  // maximum spanning tree (Prim)
    double best = -1;
    int next = 0;
    for (int i = 0; i < k; ++i)
      for (int j = 0; j < k; ++j)
        if (in[i] && !in[j] && both[i * k + j] > best) {
          best = both[i * k + j];
          next = j;
        }
    in[next] = true;
    tree += best;
  }
  if (lower >= alpha) return false;
  if (k * pmin - tree < alpha) return true;
  double lo[4], hi[4];
  for (int i = 0; i < k; ++i) {
    lo[i] = a == Alt::Two ? -t : -stats::INF;
    hi[i] = t;
  }
  return (alpha < COMBO_TAIL ? nph::mvn_outside(k, lo, hi, r) : 1 - nph::mvn_probability(k, lo, hi, r)) < alpha;
}

/// Sorts observations with times in [0, end]: a counting pass into as many buckets as observations, then each
/// (usually tiny) bucket by comparison, so linear time for spread-out times and never worse than one sort.
void sort_by_time(std::vector<nph::Obs>& obs, double end) {
  const std::size_t n = obs.size();
  const double scale = static_cast<double>(n) / end;
  std::vector<std::uint32_t> first(n + 1, 0), slot(n);
  for (std::size_t i = 0; i < n; ++i) {
    slot[i] = static_cast<std::uint32_t>(std::min(n - 1, static_cast<std::size_t>(obs[i].time * scale)));
    ++first[slot[i] + 1];
  }
  for (std::size_t b = 0; b < n; ++b) first[b + 1] += first[b];
  std::vector<nph::Obs> out(n);
  for (std::size_t i = 0; i < n; ++i) out[first[slot[i]]++] = obs[i];  // first[b] ends as bucket b's end
  for (std::size_t b = 0, lo = 0; b < n; lo = first[b++])
    if (first[b] - lo > 1)
      std::sort(out.begin() + static_cast<std::ptrdiff_t>(lo), out.begin() + static_cast<std::ptrdiff_t>(first[b]),
                [](const nph::Obs& x, const nph::Obs& y) { return x.time < y.time; });
  obs.swap(out);
}

sim::Outcome nph_study(const Design& d, long n1, long n2, long rep, std::uint64_t seed, double alpha, Alt a) {
  std::vector<nph::Obs> obs(static_cast<std::size_t>(n1 + n2));
  long events = 0;
  std::size_t at = 0;
  for (int g = 0; g < 2; ++g) {
    for (long i = 0; i < (g == 0 ? n1 : n2); ++i) {
      Rng rng(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g + 1),
              static_cast<std::uint64_t>(i));
      const double entry = d.A * rng.uniform();
      const double t = d.event_time(g, rng.exponential());
      double c = d.A + d.F - entry;  // administrative censoring at the analysis
      if (d.eta > 0) c = std::fmin(c, rng.exponential() / d.eta);
      const bool event = t <= c;
      obs[at++] = {event ? t : c, g, event ? 1 : 0};
      events += event ? 1 : 0;
    }
  }
  sort_by_time(obs, d.A + d.F);
  sim::Outcome o;
  o.extra = static_cast<double>(events);
  double z = 0;
  if (d.test == Test::MaxCombo) {
    o.reject = maxcombo_rejects(nph::weighted_logrank(obs, n1, n2, 4, COMBO_RHO, COMBO_GAMMA), alpha, a);
    return o;
  }
  if (d.test == Test::Rmst) {
    const nph::Rmst r = nph::restricted_mean(obs, n1, n2, d.tau);
    o.estimate = r.rmst[0] - r.rmst[1];
    o.se = std::sqrt(r.var[0] + r.var[1]);
    if (o.se > 0) z = -o.estimate / o.se;
  } else {
    const nph::WeightedLogrank w = nph::weighted_logrank(obs, n1, n2, 1, &d.rho, &d.gamma);
    if (w.v[0] > 0) z = w.u[0] / std::sqrt(w.v[0]);
  }
  o.reject = p_value(z, stats::INF, a) < alpha;
  return o;
}

/// "to detect a delayed treatment effect: no difference for the first 3 time units, then a hazard ratio of 0.7".
std::string effect_words(const Row& r) {
  const double h1 = r.num("hr1"), h2 = r.num("hr2"), tc = r.num("tchange"), end = r.num("accrual") + r.num("follow");
  const std::string first = "for the first " + fnum(tc) + " time units";
  if (tc <= 0 || tc >= end || h1 == h2) {  // a single hazard ratio during the study
    const double h = tc <= 0 ? h2 : h1;
    if (h == 1) return "when there is no treatment effect (a hazard ratio of 1): this is the type I error rate";
    return "to detect a hazard ratio of " + fnum(h) + " (proportional hazards)";
  }
  if (h1 == 1)
    return "to detect a delayed treatment effect: no difference " + first + ", then a hazard ratio of " + fnum(h2);
  if (h2 == 1)
    return "to detect a waning treatment effect: a hazard ratio of " + fnum(h1) + " " + first + ", then no difference";
  return std::string("to detect ") + ((h1 < 1) != (h2 < 1) ? "crossing hazards" : "a changing treatment effect") +
         ": a hazard ratio of " + fnum(h1) + " " + first + ", then " + fnum(h2);
}

std::string test_words(const Scenario& s, const Row& r) {
  if (s.is("test", "fh"))
    return "Fleming–Harrington FH(" + fnum(r.num("rho")) + ", " + fnum(r.num("gamma")) + ") weighted logrank test";
  if (s.is("test", "maxcombo"))
    return "MaxCombo test (the largest of the FH(0, 0), FH(0, 1), FH(1, 0) and FH(1, 1) statistics)";
  if (s.is("test", "rmst"))
    return "test of the difference in restricted mean survival time up to τ = " + fnum(r.num("tau"));
  return "logrank test";
}

Procedure nph() {
  Procedure p;
  p.id = "sim.surv.nph";
  p.name = "Two Survival Curves with Non-Proportional Hazards, by Simulation";
  p.category = "Simulation";
  p.subcategory = "Survival";
  p.order = 1;
  p.simulation = true;
  p.summary =
      "Power and sample size for comparing two survival curves when the hazard ratio changes over time, such as "
      "a delayed treatment effect, a waning effect or crossing hazards, estimated by simulating whole studies "
      "with accrual, follow-up and loss to follow-up. Each study is analyzed with the logrank test, a "
      "Fleming–Harrington weighted logrank test, the MaxCombo test or the difference in restricted mean survival "
      "time.";
  p.keywords = {"simulation", "monte carlo", "non-proportional hazards", "nonproportional hazards", "nph",
                "delayed effect", "delayed treatment effect", "delayed separation", "immunotherapy",
                "immuno-oncology", "waning effect", "crossing hazards", "crossing survival curves",
                "weighted logrank", "fleming-harrington", "maxcombo", "max-combo", "combination test", "rmst",
                "restricted mean survival time", "piecewise exponential"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total Sample Size)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: HR ≠ 1)", "One-sided (H1: HR > 1)", "One-sided (H1: HR < 1)"),
      select("test", "Test",
             {{"logrank", "Logrank", {}},
              {"fh", "Fleming–Harrington weighted logrank FH(ρ, γ)", {}},
              {"maxcombo", "MaxCombo: largest of FH(0, 0), FH(0, 1), FH(1, 0), FH(1, 1)", {}},
              {"rmst", "Difference in restricted mean survival time (RMST)", {}}},
             "logrank")
          .grp(G_TEST)
          .hint("The analysis of each simulated study. Weighted logrank tests with γ > 0 emphasize late "
                "differences (delayed effects), ρ > 0 early ones (waning effects); MaxCombo takes the most "
                "significant of four weightings and adjusts for the choice."),
      num("rho", "ρ (FH Weight Exponent for Early Differences)", "0")
          .grp(G_TEST)
          .sym("ρ")
          .range(0, 10)
          .sugg({"0", "1"})
          .hint("Exponent ρ of the weight Ŝ(t−)^ρ (1 − Ŝ(t−))^γ, Ŝ the pooled Kaplan–Meier estimate: values above 0 "
                "give early differences more weight. FH(0, 0) is the logrank test.")
          .show(when(is("test", {"fh"}))),
      num("gamma", "γ (FH Weight Exponent for Late Differences)", "1")
          .grp(G_TEST)
          .sym("γ")
          .range(0, 10)
          .sugg({"1", "0 0.5 1"})
          .hint("Exponent γ of the weight (1 − Ŝ(t−))^γ: values above 0 give late differences more weight, as "
                "suits a delayed treatment effect.")
          .show(when(is("test", {"fh"}))),
      num("tau", "τ (Truncation Time of the RMST)", "18")
          .grp(G_TEST)
          .sym("τ")
          .above(0)
          .sugg({"18", "12 18 24"})
          .hint("The restricted mean survival time is the area under the survival curve up to τ, the mean "
                "survival time over the first τ time units. τ must not exceed accrual + follow-up.")
          .show(when(is("test", {"rmst"}))),
      power_param(),
      alpha_param(),
      num("n1", "N1 (Subjects in Group 1, Treatment)", "200")
          .grp(G_N)
          .sym("N1")
          .whole()
          .at_least(2)
          .sugg({"200", "100 to 400 by 50"})
          .hint("Number of subjects in the treatment group (group 1); the control group gets N1 × the allocation "
                "ratio, rounded up.")
          .show(when(not_solving({"n"}))),
      num("r", "Allocation Ratio N2/N1", "1")
          .grp(G_N)
          .sym("N2/N1")
          .range(0.01, 100)
          .sugg({"1", "0.5", "2"})
          .hint("Size of the control group relative to the treatment group (N2 = ratio × N1, rounded up); 1 gives "
                "equal groups."),
      num("hr1", "HR1 (Hazard Ratio before T*)", "1")
          .grp(G_EFFECT)
          .sym("HR1")
          .range(0.01, 100)
          .sugg({"1", "0.7", "1 0.8"})
          .hint("Hazard ratio of treatment to control from entry until T*; 1 for a delayed effect."),
      num("hr2", "HR2 (Hazard Ratio from T* on)", "0.7")
          .grp(G_EFFECT)
          .sym("HR2")
          .range(0.01, 100)
          .sugg({"0.7", "0.6 0.7 0.8", "1"})
          .hint("Hazard ratio of treatment to control after T*; 1 for an effect that wanes."),
      num("tchange", "T* (Time When the Hazard Ratio Changes)", "3")
          .grp(G_EFFECT)
          .sym("T*")
          .at_least(0)
          .sugg({"3", "0 3 6"})
          .hint("Time after entry at which the hazard ratio changes from HR1 to HR2. With HR1 = HR2, or T* = 0, "
                "the hazards are proportional."),
      num("med2", "Median Survival (Control)", "12")
          .grp(G_EFFECT)
          .sym("Median2")
          .above(0)
          .sugg({"12", "6 12 24"})
          .hint("Median survival time of the control group, whose survival is taken to be exponential, in the time "
                "unit of the study."),
      num("accrual", "Accrual Time (Uniform Entry)", "12")
          .grp("Study Duration")
          .sym("Accrual")
          .at_least(0)
          .sugg({"12", "24"})
          .hint("Length of the period over which subjects enter the study at an even rate, in the time unit of the "
                "study; 0 if all start together."),
      num("follow", "Additional Follow-Up after Accrual", "12")
          .grp("Study Duration")
          .sym("Follow-Up")
          .at_least(0)
          .sugg({"12", "6 12 24"})
          .hint("Time from the end of accrual to the final analysis; the last subject to enter is followed for "
                "this long."),
      num("eta", "Loss to Follow-Up Hazard Rate", "0")
          .grp("Study Duration")
          .sym("η")
          .at_least(0)
          .sugg({"0", "0.01"})
          .hint("Rate per unit of time at which subjects are lost to follow-up for reasons other than the end of "
                "the study (a constant hazard); 0 for none."),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double alpha = s("alpha"), ratio = s("r");
    const Design d = nph_design(s);
    const RmstExpected rm = d.test == Test::Rmst ? rmst_expected(d) : RmstExpected{};
    const ComboNull combo =
        d.test == Test::MaxCombo ? combo_null(d, 1 / (1 + ratio), alpha, a, ctx.threads) : ComboNull{};
    auto n2_of = [&](double n1) { return std::ceil(ratio * n1 - 1e-9); };
    auto approx = [&](double n1) { return approx_power(d, rm, combo, n1, n2_of(n1), alpha, a, ctx.threads); };
    if (s.solve == "n") {
      // the hazard ratios that apply during the study
      const bool before = d.tc > 0, after = d.tc < d.A + d.F;
      const double hmin = std::fmin(before ? d.hr1 : d.hr2, after ? d.hr2 : d.hr1);
      const double hmax = std::fmax(before ? d.hr1 : d.hr2, after ? d.hr2 : d.hr1);
      if (hmin == 1 && hmax == 1) throw CalcError("the hazard ratio is 1 throughout the study: there is no effect");
      if ((a == Alt::Upper && hmax <= 1) || (a == Alt::Lower && hmin >= 1))
        throw CalcError("the hazard ratios are on the side opposite to the one-sided alternative");
    }
    const double cap = std::fmax(2, std::floor(MAX_N / (1 + ratio)));
    const double lo = s.solve != "n" ? 2
                                     : solve::first_whole([&](double x) { return n2_of(x) >= 2; },
                                                          std::fmax(2, std::floor((1 + 1e-9) / ratio) + 1), cap,
                                                          "the allocation ratio is too small: group 2 would have "
                                                          "fewer than 2 subjects");
    const double hi = std::fmax(lo, cap);
    double n1;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = lo;
      o.hi = hi;
      n1 = solve::min_n(approx, s("power"), o, &ctx);
    } else {
      n1 = s("n1");
      if (n2_of(n1) < 2) throw CalcError("group 2 needs at least 2 subjects");
    }
    double events = NAN;
    if (ctx.approximate) {
      simc::put_approximation(s, r, approx(n1));
    } else {
      const long reps = static_cast<long>(s("reps"));
      const auto seed = static_cast<std::uint64_t>(s("seed"));
      auto power_at = [&](double m1, long R, double p0, double p1) {
        simc::check_study_size(m1 + n2_of(m1), "subjects");
        const long k1 = static_cast<long>(m1), k2 = static_cast<long>(n2_of(m1));
        return sim::simulate(R, [&](long rep) { return nph_study(d, k1, k2, rep, seed, alpha, a); }, ctx, p0, p1);
      };
      sim::Summary sm;
      if (s.solve == "n") {
        sim::SizeSearch o;
        o.lo = lo;
        o.hi = hi;
        o.guess = n1;
        o.reps = reps;
        o.explore_reps = simc::explore_reps(reps);
        const sim::SizeResult res = sim::search_size(power_at, s("power"), o, ctx);
        n1 = res.n;
        sm = res.s;
      } else {
        sm = power_at(n1, reps, 0, 1);
      }
      simc::put_summary(s, r, sm);
      r.set("approx_power", approx(n1));
      events = sm.mean_extra;
      if (d.test == Test::Rmst && std::isfinite(sm.mean_estimate)) r.set("mean_drmst", sm.mean_estimate);
    }
    const double n2 = n2_of(n1);
    const Expected e = expected(d, n1 / (n1 + n2), 0, nullptr, nullptr);
    put_two_group(r, {n1, n2});
    r.set("events", ctx.approximate ? (n1 + n2) * e.events : events);
    r.set("hr1", d.hr1);
    r.set("hr2", d.hr2);
    r.set("tchange", d.tc);
    r.set("ahr", std::exp(e.log_hr));
    r.set("med2", s("med2"));
    r.set("accrual", d.A);
    r.set("follow", d.F);
    r.set("eta", d.eta);
    if (d.test == Test::FH) {
      r.set("rho", d.rho);
      r.set("gamma", d.gamma);
    }
    if (d.test == Test::Rmst) {
      r.set("tau", d.tau);
      r.set("drmst", -rm.diff);
    }
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    c.push_back(col("n1", "N1", "int", "Subjects in group 1 (treatment)"));
    c.push_back(col("n2", "N2", "int", "Subjects in group 2 (control)"));
    c.push_back(col("n", "N", "int", "Total sample size N1 + N2"));
    for (auto x : {col("events", "Events", "num",
                       "Average number of events in the simulated studies (the expected number in an approximation)"),
                   col("hr1", "HR1", "num", "Hazard ratio (treatment/control) before T*"),
                   col("hr2", "HR2", "num", "Hazard ratio from T* on"),
                   col("tchange", "T*", "num", "Time after entry at which the hazard ratio changes"),
                   col("ahr", "AHR", "num", "Average hazard ratio: exp of the mean log hazard ratio over the events"),
                   col("med2", "Median2", "num", "Control median survival"),
                   col("accrual", "Accrual", "num"), col("follow", "Follow-Up", "num"),
                   col("eta", "η", "num", "Loss to follow-up hazard rate"),
                   col("rho", "ρ", "num", "FH weight exponent of Ŝ(t−)"),
                   col("gamma", "γ", "num", "FH weight exponent of 1 − Ŝ(t−)"),
                   col("tau", "τ", "num", "Truncation time of the RMST"),
                   col("drmst", "ΔRMST", "num", "Difference in restricted mean survival time up to τ (treatment − "
                       "control)"),
                   col("mean_drmst", "Mean ΔRMST", "num", "Average estimated RMST difference in the simulated studies"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns())
      if (x.key != "failures") c.push_back(x);  // every simulated study can be analyzed
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const bool approx = text_of(r, "method") == "approximation";
    std::string t = "Group sample sizes of " + fint(r.num("n1")) + " (treatment) and " + fint(r.num("n2")) +
                    " (control) " + simc::power_phrase(s, r) + " with a " + alt_words(alt(s)) + " " +
                    test_words(s, r) + " at α = " + fnum(r.num("alpha")) + " " + effect_words(r) + ".";
    if (s.is("test", "rmst")) {
      const double dr = r.num("drmst");
      t += " The restricted mean survival time up to τ is " + fnum(std::fabs(dr), 2) + " time units" +
           (dr >= 0 ? " longer" : " shorter") + " with the treatment.";
    }
    t += " Control survival is exponential with a median of " + fnum(r.num("med2")) + " time units; ";
    t += r.num("accrual") > 0 ? "subjects enter uniformly over " + fnum(r.num("accrual")) +
                                    " time units and are followed for " + fnum(r.num("follow")) + " more"
                              : "all subjects enter at the start and are followed for " + fnum(r.num("follow")) +
                                    " time units";
    if (r.num("eta") > 0) t += ", with loss to follow-up at a hazard rate of " + fnum(r.num("eta"));
    t += approx ? "; about " + fint(std::round(r.num("events"))) + " events are expected."
                : "; the simulated studies have " + fint(std::round(r.num("events"))) + " events on average.";
    return t + simc::method_note(r);
  };
  p.definitions = {
      {"N1, N2", "Number of subjects in group 1 (treatment) and group 2 (control); N = N1 + N2."},
      {"HR1, HR2, T*", "Hazard ratio of treatment to control before and from time T* after entry on."},
      {"Median2", "Median survival of the control group (exponential survival)."},
      {"Accrual, Follow-Up", "Length of uniform enrollment and of follow-up after enrollment ends; the analysis "
                             "is at Accrual + Follow-Up."},
      {"η", "Hazard rate of loss to follow-up (exponential, the same in both groups)."},
      {"FH(ρ, γ)", "Fleming–Harrington weighted logrank test with weights Ŝ(t−)^ρ (1 − Ŝ(t−))^γ, Ŝ the pooled "
                   "Kaplan–Meier estimate; FH(0, 0) is the logrank test."},
      {"MaxCombo", "Largest of the standardized FH(0, 0), FH(0, 1), FH(1, 0) and FH(1, 1) statistics, with its "
                   "p-value from their joint normal distribution."},
      {"RMST", "Restricted mean survival time: the area under the survival curve up to τ."},
      {"AHR", "Average hazard ratio: exp of the mean log hazard ratio, weighted by the expected events."},
      {"Events", "Average number of events in the simulated studies (the expected number in an approximation)."},
      {"Power", "Share of the simulated studies whose test was significant, with its 95% interval (Wilson) and "
                "Monte Carlo standard error."}};
  p.notes = {
      "Every simulated subject enters uniformly during accrual and is followed until the analysis at Accrual + "
      "Follow-Up or until lost to follow-up (exponential with rate η). Control event times are exponential with "
      "the given median; treatment event times have the hazard HR1·λ2 before T* and HR2·λ2 after it, and are "
      "drawn by inverting the cumulative hazard.",
      "All statistics are oriented like the log hazard ratio: negative values favor the treatment (fewer events "
      "than expected in group 1, or a longer restricted mean survival time), so the one-sided alternative "
      "HR < 1 means that the treatment is better.",
      "The weighted logrank statistics use the hypergeometric variance. The two-sided MaxCombo test takes the "
      "largest |Z| of its four components, the one-sided test the largest of the one-sided statistics; its "
      "p-value comes from the multivariate normal distribution with the correlation estimated from the data "
      "(Genz's method with a fixed lattice rule).",
      "The RMST of each group is the area under its Kaplan–Meier curve up to τ (flat after the last "
      "observation), with the variance Σ A_j² d_j / (n_j (n_j − d_j)), A_j the area from the j-th event time "
      "to τ; the difference is tested with a z statistic.",
      "The quick approximation (live previews, power curves and the start of the sample-size search) integrates "
      "the expected events over time: E[Z] ≈ √(p(1 − p)) Σ w_k log HR_k d_k / √(Σ w_k² d_k), with p = N1/N and "
      "w_k the weight at the expected pooled survival. With proportional hazards and the logrank test it is "
      "Schoenfeld's formula, as in Logrank Tests. For MaxCombo the four statistics are taken as jointly normal "
      "with these means and the correlation implied by the expected events, and the critical value of their "
      "maximum keeps the type I error at α; for RMST it uses the asymptotic variance of the Kaplan–Meier areas.",
      "Results depend only on the seed: each subject's data come from their own random stream, so a larger study "
      "contains the subjects of a smaller one, which keeps the sample-size search stable."};
  p.references = {"Fleming, T. R. & Harrington, D. P. (1991). Counting Processes and Survival Analysis. Wiley.",
                  "Lin, R. S., Lin, J., Roychoudhury, S. et al. (2020). Alternative analysis methods for time to "
                  "event endpoints under nonproportional hazards: a comparative analysis. Statistics in "
                  "Biopharmaceutical Research 12, 187–198.",
                  "Royston, P. & Parmar, M. K. B. (2013). Restricted mean survival time: an alternative to the "
                  "hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. BMC "
                  "Medical Research Methodology 13, 152.",
                  "Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival "
                  "distributions. Biometrika 68, 316–319.",
                  "Genz, A. (1992). Numerical computation of multivariate normal probabilities. Journal of "
                  "Computational and Graphical Statistics 1, 141–149.",
                  simc::REF_SIM};
  p.curve_param = "n1";
  p.curve_label = "N1 (Subjects in Group 1)";
  p.curve_min = 2;
  p.example_solve = "n";
  p.example = {{"power", "0.8"}, {"test", "maxcombo"}, {"hr1", "1"},      {"hr2", "0.6"},
               {"tchange", "3"}, {"med2", "12"},       {"accrual", "12"}, {"follow", "12"}, {"reps", "300"}};
  return p;
}
SSC_REGISTER(nph);

}  // namespace
}  // namespace ssc::procs
