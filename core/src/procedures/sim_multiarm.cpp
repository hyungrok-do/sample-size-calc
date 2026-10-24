// Simulation-based power for comparing several treatments with one control
// (sim.multiarm): K treatment arms and a control, normal outcomes with a
// common SD, each treatment compared with the control by a t-test on the
// pooled variance (N − K − 1 degrees of freedom), and the familywise error
// controlled by Dunnett's single-step test, Holm's or Hochberg's step-wise
// procedures or Bonferroni's correction.
//
// Power has three meanings with several hypotheses: finding at least one of
// the effective treatments (disjunctive power), finding all of them
// (conjunctive power), or finding each one (the smallest per-arm power).
// For the single-step tests (Dunnett, Bonferroni, none) all three are exact:
// given the control's mean and the pooled SD the treatment statistics are
// independent, so the probabilities are a double integral (over the common
// normal component and the SD's chi distribution).  The step-wise procedures
// have no closed form; their quick approximation is Bonferroni's exact power,
// a lower bound.  Each study is generated from its sufficient statistics: the
// arms' means and the pooled sum of squares.
#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

#include "procedures/sim_common.hpp"
#include "stats/numerics.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;
using sim::Rng;

constexpr int KMAX = 8;   // treatment arms
// bits of a study's outcome (sim::Outcome::also): arm i rejected (bits 0–7), then these events
constexpr int ANY_BIT = 8, ALL_BIT = 9, FALSE_BIT = 10;
constexpr int Z_NODES = 64, S_NODES = 48;
constexpr double Z_RANGE = 8.5;

enum class Adjust { Dunnett, Holm, Hochberg, Bonferroni, None };
enum class Goal { Any, All, Each };

struct Design {
  int k = 0;                       // treatment arms
  double mu0 = 0, sigma = 1;
  std::array<double, KMAX> mu{};   // treatment means
  double origin = 0;               // the arms are simulated around it: μ0 when the means are far from 0 next to σ
  Alt alt = Alt::Two;
  Adjust adjust = Adjust::Dunnett;
  Goal goal = Goal::Any;
  double alpha = 0.05, r0 = 1;
  /// Whether treatment i differs from the control in the direction of the alternative.
  bool effective(int i) const {
    const double d = mu[i] - mu0;
    if (alt == Alt::Two) return d != 0;
    return alt == Alt::Upper ? d > 0 : d < 0;
  }
  int effective_count() const {
    int e = 0;
    for (int i = 0; i < k; ++i) e += effective(i);
    return e;
  }
};

/// Sizes: n per treatment arm, n0 = ⌈r0·n⌉ in the control arm.
struct Sizes {
  double n = 0, n0 = 0;
  double df(int k) const { return k * n + n0 - k - 1; }
};

/// ln P(Z ≥ y), also where P(Z ≥ y) is below the smallest double.
double log_upper_normal(double y) {
  if (y < 30) return std::log(stats::pnorm(y, false));
  return -y * y / 2 - std::log(y) - stats::LN_SQRT_2PI + std::log1p(-1 / (y * y));
}

/// Nodes of the SD ratio S = √(χ²_ν/ν): Gauss–Legendre in t = log S, where the density is smooth and
/// unimodal, ∝ exp(g(t)) with g(t) = ν t − ν (e^{2t} − 1)/2, over the range where g > −40, in panels of at most 1.5
/// where it is wider than 3 (with few degrees of freedom: 22 with 2, 3 with 26), since P(T ≥ c | S) falls from about
/// 1 to 0 within a small part of it (with 48 nodes over the whole range the power was up to 3.5% off with 2 degrees
/// of freedom and 4e-6 with 18); the weights are normalized to sum to 1.  One node S = 1 when ν is huge.
struct SNodes {
  std::vector<double> s, w;
  explicit SNodes(double df) {
    if (df > 1e7) {
      s = {1};
      w = {1};
      return;
    }
    auto g = [&](double t) { return df * t - df * std::expm1(2 * t) / 2; };
    double lo = -1, hi = 1;
    while (g(lo) > -40) lo *= 2;
    while (g(hi) > -40) hi *= 2;
    a_ = numerics::brent([&](double t) { return g(t) + 40; }, lo, 0, 1e-12);
    b_ = numerics::brent([&](double t) { return g(t) + 40; }, 0, hi, 1e-12);
    const stats::GLRule& rule = stats::gauss_legendre(S_NODES);
    const int panels = b_ - a_ > 3 ? static_cast<int>(std::ceil((b_ - a_) / 1.5)) : 1;
    h_ = (b_ - a_) / panels;
    for (int q = 0; q < panels; ++q)
      for (std::size_t i = 0; i < rule.x.size(); ++i) {
        const double t = a_ + h_ * (q + (rule.x[i] + 1) / 2);
        s.push_back(std::exp(t));
        w.push_back(rule.w[i] * std::exp(g(t)));
        total_ += w.back();
      }
    for (double& x : w) x /= total_;
  }

  /// The same nodes and, where an arm with the standardized effect `low` (the arm least likely to reject) rejects at
  /// the critical value c mostly at a smaller S than they reach (a tiny alpha, more so with few degrees of freedom),
  /// panels there as well: over the range where the density times P(T ≥ c | S) is within e^−40 of its largest
  /// value, when more than e^−30 of that lies below the density's nodes.
  SNodes(double df, double c, double low) : SNodes(df) {
    if (df > 1e7) return;
    auto g = [&](double t) { return df * t - df * std::expm1(2 * t) / 2; };
    auto h = [&](double t) { return g(t) + log_upper_normal(c * std::exp(t) - low); };  // concave
    double a = -60, b = 5;  // its largest value, by golden section
    while (b - a > 1e-9) {
      const double m1 = b - 0.6180339887498949 * (b - a), m2 = a + 0.6180339887498949 * (b - a);
      if (h(m1) < h(m2))
        a = m1;
      else
        b = m2;
    }
    const double top = (a + b) / 2, level = h(top) - 40;
    if (top >= a_ && !(h(a_) > level + 10)) return;
    double lo = top - 1, hi = top + 1;
    while (h(lo) > level) lo -= 2 * (top - lo);
    while (h(hi) > level) hi += 2 * (hi - top);
    const double ta = numerics::brent([&](double t) { return h(t) - level; }, lo, top, 1e-12);
    const double tb = std::fmin(a_, numerics::brent([&](double t) { return h(t) - level; }, top, hi, 1e-12));
    if (!(tb > ta)) return;
    const stats::GLRule& rule = stats::gauss_legendre(S_NODES);
    const int panels = static_cast<int>(std::ceil((tb - ta) / 1.5));  // as narrow as the density's own
    const double width = (tb - ta) / panels;
    for (int q = 0; q < panels; ++q)
      for (std::size_t i = 0; i < rule.x.size(); ++i) {  // weighted as the density's panels
        const double t = ta + width * (q + (rule.x[i] + 1) / 2);
        s.push_back(std::exp(t));
        w.push_back(rule.w[i] * width / h_ * std::exp(g(t)) / total_);
      }
  }

 private:
  double a_ = 0, b_ = 0, h_ = 0, total_ = 0;  // the density's range in t, its panels' width, their weights' sum
};

using Probs = std::array<double, KMAX>;
constexpr int EVENTS_MAX = KMAX + 3;  // events integrated in one pass, and the rule's own total
using Events = std::array<double, EVENTS_MAX>;

/// For the single-step tests with critical value c on |T| (two-sided) or T (one-sided, after orienting the
/// statistics): given the common component z and the SD ratio s, arm i rejects with probability
/// P(i | z, s); the arms are then independent.  `events(probabilities, e)` sets the conditional
/// probabilities e[0], ..., e[m − 1] of m events, whose probabilities come out together from one pass over
/// the nodes.  Arms with the same effect (all of them under the null) have the same probabilities, which are
/// computed once.  The SD nodes are shared between threads, and their sums added in order.
template <class F>
Events single_step(const Design& d, const Sizes& n, double c, const SNodes& sn, int m, F&& events,
                   const RunContext& ctx) {
  const double lam = std::sqrt(n.n / (n.n + n.n0)), root = std::sqrt(1 - lam * lam);
  const double se = d.sigma * std::sqrt(1 / n.n + 1 / n.n0);
  std::array<double, KMAX> theta{};
  std::array<int, KMAX> same{};  // the first arm with arm i's effect
  const double sign = d.alt == Alt::Lower ? -1 : 1;
  for (int i = 0; i < d.k; ++i) {
    theta[i] = sign * (d.mu[i] - d.mu0) / se;
    same[i] = i;
    for (int j = 0; j < i && same[i] == i; ++j)
      if (theta[j] == theta[i]) same[i] = j;
  }
  const stats::GLRule& g = stats::gauss_legendre(Z_NODES);
  // an arm's rejection changes with z over a width of about root/λ: with a control arm under a third of a treatment
  // arm (root below 1/2), the nodes are split into two panels (with r0 = 0.1 and 8 arms, one panel gave the power
  // that all arms reject 0.7% low)
  const int panels = root >= 0.5 ? 1 : 2;
  // Given s, z given that an arm rejects by T ≥ c is normal around λ·E(X | X ≥ a), X standard normal and a = cs − θ,
  // with an SD of at most √(root² + λ²·Var(X | X ≥ a)) (by T ≤ −c the same with −z and a = cs + θ).  Where 8 of
  // these SDs reach past the nodes for an arm (two-sided, unless it rejects that way less than e^−40 as often as the
  // other), they go on in panels of the same width: with a tiny alpha the rejections lie beyond z = 8.5 (8 arms,
  // r0 = 0.1, α = 1e-20: an arm's power came out 4% of its value).  These panels add to the events but not to the
  // rule's own total, which they would change by less than its rounding.
  auto reach = [&](double a) {
    if (!(a < 38)) return 0.0;  // P(X ≥ a) is below the smallest normal double
    const double mills = std::exp(-a * a / 2 - stats::LN_SQRT_2PI - log_upper_normal(a));  // E(X | X ≥ a)
    const double var = std::clamp(1 + a * mills - mills * mills, 0.0, 1.0);
    return lam * mills + 8 * std::sqrt(root * root + lam * lam * var);
  };
  std::vector<Events> part(sn.s.size());
  parallel_for(sn.s.size(), ctx.threads, [&](std::size_t a) {
    const double cs = c * sn.s[a];
    Events inner{}, e{};
    e[m] = 1;  // the rule's own total, summed exactly as the events are
    Probs prob{};
    auto panel = [&](int q, int count) {  // panel q of [−8.5, 8.5] split into `panels`, also beyond it
      for (std::size_t b = 0; b < g.x.size(); ++b) {
        const double z = Z_RANGE * (g.x[b] + (2 * q + 1 - panels)) / panels;  // as before with one panel
        for (int i = 0; i < d.k; ++i) {
          if (same[i] != i) {
            prob[i] = prob[same[i]];
            continue;
          }
          const double mi = theta[i] + lam * z;
          prob[i] = d.alt == Alt::Two ? stats::pnorm((-cs - mi) / root) + stats::pnorm((mi - cs) / root)
                                      : stats::pnorm((mi - cs) / root);
        }
        const double weight = g.w[b] / panels * stats::dnorm(z);
        events(prob, e);
        for (int j = 0; j < count; ++j) inner[j] += weight * e[j];
      }
    };
    for (int q = 0; q < panels; ++q) panel(q, m + 1);
    double zhi = 0, zlo = 0;
    for (int i = 0; i < d.k; ++i) {
      if (same[i] != i) continue;
      const double above = cs - theta[i], below = cs + theta[i];
      const bool two = d.alt == Alt::Two;
      if (!two || log_upper_normal(above) > log_upper_normal(below) - 40) zhi = std::fmax(zhi, reach(above));
      if (two && log_upper_normal(below) > log_upper_normal(above) - 40) zlo = std::fmin(zlo, -reach(below));
    }
    for (int q = panels; Z_RANGE * (2 * q - panels) / panels < zhi; ++q) panel(q, m);
    for (int q = -1; Z_RANGE * (2 * q + 2 - panels) / panels > zlo; --q) panel(q, m);
    part[a] = inner;
  }, ctx.token);
  Events total{};
  for (std::size_t a = 0; a < sn.s.size(); ++a)
    for (int j = 0; j <= m; ++j) total[j] += sn.w[a] * Z_RANGE * part[a][j];
  // The weights add up to 1 only to rounding, which differs between processors (fused multiply-adds on Arm):
  // divided by their own total, an event that every node gives is exactly 1, and none goes beyond [0, 1].
  for (int j = 0; j < m; ++j) total[j] = stats::clamp01(total[j] / total[m]);
  total[m] = 0;
  return total;
}

/// Unadjusted critical value of one t-test at level `level`.
double t_critical(double level, double df, Alt a) {
  return stats::qt(tail_alpha(level, a), df, false);
}

/// The familywise error of a single-step test with critical value c, P(some arm rejects | every null true), from
/// the overlap of the arms' rejections rather than as 1 − P(no arm rejects), which single_step resolves only to
/// about 1e-16, and on its nodes only to 0.3% of α with a small control arm (k = 8, r0 = 0.1: 0.05014 for 0.05)
/// and 1e-4 with 3 degrees of freedom: k·P(T ≥ c) less E[k·p − (1 − (1 − p)^k)], with p = P(an arm rejects | z, s)
/// under the null and P(T ≥ c) (both tails if two-sided) exact.  The overlap is computed without cancellation, as
/// −expm1(k·log1p(−p)) for 1 − (1 − p)^k, and as a share of the nodes' own sum of p, so that their error cancels.
/// The nodes follow the rejections, also into the far tail where single_step's do not reach: in t = log s over the
/// range where the density times P(T ≥ c | s) is within e^−40 of its largest value, and in z out to λ·c·s (where
/// an arm's statistic reaches c) and 8.5 beyond.
double familywise_error(const Design& d, const Sizes& n, double c) {
  const double df = n.df(d.k);
  const double lam = std::sqrt(n.n / (n.n + n.n0)), root = std::sqrt(1 - lam * lam);
  const bool two = d.alt == Alt::Two;
  std::vector<double> s = {1}, w = {1};
  if (df <= 1e7) {
    auto g = [&](double t) { return df * t - df * std::expm1(2 * t) / 2; };
    auto h = [&](double t) { return g(t) + log_upper_normal(c * std::exp(t)); };  // concave
    double a = -60, b = 5;  // its largest value, by golden section
    while (b - a > 1e-9) {
      const double m1 = b - 0.6180339887498949 * (b - a), m2 = a + 0.6180339887498949 * (b - a);
      if (h(m1) < h(m2))
        a = m1;
      else
        b = m2;
    }
    const double top = (a + b) / 2, level = h(top) - 40;
    double lo = top - 1, hi = top + 1;
    while (h(lo) > level) lo -= 2 * (top - lo);
    while (h(hi) > level) hi += 2 * (hi - top);
    const double ta = numerics::brent([&](double t) { return h(t) - level; }, lo, top, 1e-12);
    const double tb = numerics::brent([&](double t) { return h(t) - level; }, top, hi, 1e-12);
    const stats::GLRule& rule = stats::gauss_legendre(S_NODES);
    s.clear();
    w.clear();
    for (std::size_t i = 0; i < rule.x.size(); ++i) {
      const double t = ta + (tb - ta) * (rule.x[i] + 1) / 2;
      s.push_back(std::exp(t));
      w.push_back(rule.w[i] * std::exp(g(t) - g(top)));  // only their ratios matter
    }
  }
  const stats::GLRule& rule = stats::gauss_legendre(Z_NODES);
  double sum_p = 0, overlap = 0;
  for (std::size_t a = 0; a < s.size(); ++a) {
    const double cs = c * s[a];
    // z from 0 when two-sided (the null's integrand is even), in panels as fine as single_step's nodes
    const double zlo = two ? 0 : -Z_RANGE, zhi = std::fmax(Z_RANGE, lam * cs + Z_RANGE);
    const int panels = static_cast<int>(std::ceil((zhi - zlo) / Z_RANGE - 1e-9));
    const double half = (zhi - zlo) / panels / 2;
    double inner_p = 0, inner_o = 0;
    for (int q = 0; q < panels; ++q)
      for (std::size_t b = 0; b < rule.x.size(); ++b) {
        const double z = zlo + half * (2 * q + 1 + rule.x[b]);
        double p = stats::pnorm((lam * z - cs) / root);
        if (two) p = stats::clamp01(p + stats::pnorm((-lam * z - cs) / root));
        const double weight = rule.w[b] * stats::dnorm(z);
        inner_p += weight * p;
        inner_o += weight * (d.k * p + std::expm1(d.k * std::log1p(-p)));
      }
    sum_p += w[a] * half * inner_p;
    overlap += w[a] * half * inner_o;
  }
  const double one = (two ? 2 : 1) * stats::pt(c, df, false);  // P(T ≥ c), or P(|T| ≥ c)
  return one * (d.k - (sum_p > 0 ? overlap / sum_p : 0));
}

/// Dunnett's critical value: P(some arm rejects | all nulls true) = α.
double dunnett_critical(const Design& d, const Sizes& n, const RunContext& ctx) {
  const double df = n.df(d.k);
  if (d.k == 1) return t_critical(d.alpha, df, d.alt);
  const double lo = t_critical(d.alpha, df, d.alt), hi = t_critical(alpha_share(d.alpha, d.k), df, d.alt) + 1e-9;
  auto f = [&](double c) {
    ctx.check();
    return std::log(familywise_error(d, n, c) / d.alpha);
  };
  const double flo = f(lo), fhi = f(hi);
  // within the nodes' error of Bonferroni's value (the arms' rejections hardly overlap so far out): that one, on
  // the safe side
  if (!(fhi < 0)) return hi;
  if (!(flo > 0)) return lo;
  return numerics::brent(f, lo, hi, flo, fhi, 1e-10);
}

double critical_value(const Design& d, const Sizes& n, const RunContext& ctx) {
  const double df = n.df(d.k);
  switch (d.adjust) {
    case Adjust::Dunnett: return dunnett_critical(d, n, ctx);
    case Adjust::None: return t_critical(d.alpha, df, d.alt);
    default: return t_critical(alpha_share(d.alpha, d.k), df, d.alt);  // Bonferroni (and its bound for Holm, Hochberg)
  }
}

/// The standardized effect of the arm least likely to reject: the smallest, or the smallest in size if two-sided.
double least_effect(const Design& d, const Sizes& n) {
  const double se = d.sigma * std::sqrt(1 / n.n + 1 / n.n0), sign = d.alt == Alt::Lower ? -1 : 1;
  double low = INFINITY;
  for (int i = 0; i < d.k; ++i) {
    const double theta = sign * (d.mu[i] - d.mu0) / se;
    low = std::fmin(low, d.alt == Alt::Two ? std::fabs(theta) : theta);
  }
  return low;
}

struct Exact {
  double any = 0, all = 0, each = 0;
  std::array<double, KMAX> arm{};
};

/// Exact power of a single-step test with critical value c (Dunnett, Bonferroni, none): that any effective
/// arm rejects, that all do, and that each arm does.
Exact exact_power(const Design& d, const Sizes& n, double c, const SNodes& sn, const RunContext& ctx) {
  Exact e;
  const int eff = d.effective_count();
  bool null = true;
  for (int i = 0; i < d.k; ++i) null = null && d.mu[i] == d.mu0;
  if (null) {
    // every null true: the familywise error, as Dunnett's value is found from it (exactly α), and each arm's exact
    // t tail
    e.any = familywise_error(d, n, c);
    for (int i = 0; i < d.k; ++i) e.arm[i] = (d.alt == Alt::Two ? 2 : 1) * stats::pt(c, n.df(d.k), false);
    return e;
  }
  const Events p = single_step(d, n, c, sn, 2 + d.k, [&](const Probs& q, Events& ev) {
    double none = 1, all = 1;
    for (int i = 0; i < d.k; ++i)
      if (eff == 0 || d.effective(i)) none *= 1 - q[i];
    for (int i = 0; i < d.k; ++i)
      if (d.effective(i)) all *= q[i];
    ev[0] = 1 - none;
    // without cancellation where none is near 1: 1 − none was 0 once every q_i was below about 1e-16, and below one
    // arm's own power before that (12.6% low at α = 1e-16 with three arms of tiny effects)
    if (none > 0.5) {
      double log_none = 0;
      for (int i = 0; i < d.k; ++i)
        if (eff == 0 || d.effective(i)) log_none += std::log1p(-q[i]);
      ev[0] = -std::expm1(log_none);
    }
    ev[1] = eff == 0 ? 0.0 : all;
    for (int i = 0; i < d.k; ++i) ev[2 + i] = q[i];
  }, ctx);
  e.any = p[0];
  e.all = p[1];
  e.each = INFINITY;
  for (int i = 0; i < d.k; ++i) {
    e.arm[i] = p[2 + i];
    if (d.effective(i)) e.each = std::fmin(e.each, e.arm[i]);
  }
  if (eff == 0) e.each = 0;
  return e;
}

double goal_power(const Design& d, const Exact& e) {
  return d.goal == Goal::Any ? e.any : (d.goal == Goal::All ? e.all : e.each);
}

/// One simulated study: the arms' means and the pooled SD from their sampling distributions.  The streams
/// are fixed per arm, and the pooled SD comes from the χ² quantile of a fixed uniform, so the simulated
/// statistics change smoothly with the sample size.
sim::Outcome multiarm_study(const Design& d, const Sizes& n, double c, long rep, std::uint64_t seed) {
  const double df = n.df(d.k);
  Rng r0(seed, static_cast<std::uint64_t>(rep), 0);
  const double ybar0 = (d.mu0 - d.origin) + d.sigma / std::sqrt(n.n0) * r0.normal();
  const double s = d.sigma * std::sqrt(stats::qchisq(r0.uniform(), df) / df);
  const double se = s * std::sqrt(1 / n.n + 1 / n.n0), sign = d.alt == Alt::Lower ? -1 : 1;
  std::array<double, KMAX> t{}, pv{};
  for (int i = 0; i < d.k; ++i) {
    Rng ri(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(i + 1));
    const double ybar = (d.mu[i] - d.origin) + d.sigma / std::sqrt(n.n) * ri.normal();
    t[i] = simc::t_stat(sign * (ybar - ybar0), se);
    pv[i] = d.alt == Alt::Two ? 2 * stats::pt(std::fabs(t[i]), df, false) : stats::pt(t[i], df, false);
  }
  unsigned rejected = 0;
  if (d.adjust == Adjust::Dunnett || d.adjust == Adjust::Bonferroni || d.adjust == Adjust::None) {
    for (int i = 0; i < d.k; ++i)
      if ((d.alt == Alt::Two ? std::fabs(t[i]) : t[i]) > c) rejected |= 1u << i;
  } else {
    std::array<int, KMAX> order{};
    for (int i = 0; i < d.k; ++i) order[i] = i;
    for (int i = 1; i < d.k; ++i)  // by p-value (insertion: a few arms)
      for (int j = i; j > 0 && pv[order[j]] < pv[order[j - 1]]; --j) std::swap(order[j], order[j - 1]);
    if (d.adjust == Adjust::Holm) {
      for (int j = 0; j < d.k; ++j) {
        if (pv[order[j]] > d.alpha / (d.k - j)) break;
        rejected |= 1u << order[j];
      }
    } else {  // Hochberg: the largest j with p(j) ≤ α/(K − j + 1), and every smaller p-value
      for (int j = d.k - 1; j >= 0; --j)
        if (pv[order[j]] <= d.alpha / (d.k - j)) {
          for (int m = 0; m <= j; ++m) rejected |= 1u << order[m];
          break;
        }
    }
  }
  sim::Outcome o;
  const int eff = d.effective_count();
  bool any = false, all = eff > 0, false_claim = false;
  for (int i = 0; i < d.k; ++i) {
    const bool rej = rejected >> i & 1u;
    if (d.effective(i)) {
      any = any || rej;
      all = all && rej;
    } else if (rej) {
      false_claim = true;
    }
  }
  if (eff == 0) any = rejected != 0;  // no effective treatment: the familywise error
  o.also = rejected | (any ? 1u << ANY_BIT : 0u) | (all ? 1u << ALL_BIT : 0u) | (false_claim ? 1u << FALSE_BIT : 0u);
  o.reject = d.goal == Goal::Any ? any : all;  // each: rebuilt from the per-arm counts
  return o;
}

Design design_of(const Scenario& s) {
  Design d;
  const std::vector<double>& mu = s.v("mu");
  if (mu.empty()) throw CalcError("enter the mean of at least one treatment arm");
  if (mu.size() > KMAX) throw CalcError("at most 8 treatment arms");
  d.k = static_cast<int>(mu.size());
  std::copy(mu.begin(), mu.end(), d.mu.begin());
  d.mu0 = s("mu0");
  d.sigma = s("sigma");
  // means far from 0 next to σ would swamp the variation of the simulated arms (μ0 = μi = 10¹⁶ with σ = 2 gave
  // no rejection at all): they are simulated around the control's mean instead (others exactly as given)
  double far = std::fabs(d.mu0);
  for (int i = 0; i < d.k; ++i) far = std::fmax(far, std::fabs(d.mu[i]));
  if (!(far < 1e6 * d.sigma)) d.origin = d.mu0;
  d.alt = alt(s);
  d.alpha = s("alpha");
  d.r0 = s("r0");
  const std::string& a = s.s("adjust");
  d.adjust = a == "holm" ? Adjust::Holm
             : a == "hochberg" ? Adjust::Hochberg
             : a == "bonferroni" ? Adjust::Bonferroni
             : a == "none" ? Adjust::None
                           : Adjust::Dunnett;
  const std::string& g = s.s("goal");
  d.goal = g == "all" ? Goal::All : (g == "each" ? Goal::Each : Goal::Any);
  return d;
}

const char* adjust_words(Adjust a) {
  switch (a) {
    case Adjust::Holm: return "Holm's step-down procedure";
    case Adjust::Hochberg: return "Hochberg's step-up procedure";
    case Adjust::Bonferroni: return "Bonferroni's correction";
    case Adjust::None: return "no adjustment for multiplicity";
    default: return "Dunnett's adjustment";
  }
}

Procedure multiarm() {
  Procedure p;
  p.id = "sim.multiarm";
  p.name = "Several Treatments against a Control with a Multiplicity Adjustment (Dunnett, Holm, Hochberg), by "
           "Simulation";
  p.category = "Simulation";
  p.subcategory = "Several Groups";
  p.order = 1;
  p.simulation = true;
  p.summary =
      "Subjects per arm (or power) for comparing each of several treatments (doses, regimens) with one control "
      "while controlling the familywise error with Dunnett's test, Holm's or Hochberg's procedure or "
      "Bonferroni's correction; power to find at least one, all or each of the effective treatments.";
  p.keywords = {"simulation", "monte carlo", "multi-arm", "multiarm", "several treatments", "multiple arms",
                "many-to-one", "dunnett", "holm", "hochberg", "bonferroni", "multiplicity", "multiple comparisons",
                "familywise error", "fwer", "dose finding", "doses against placebo", "each against control",
                "disjunctive power", "conjunctive power", "treatments versus control"};
  p.solve_for = {{"n", "Sample Size", "n", "n (Subjects per Treatment Arm)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: μi ≠ μ0)", "One-sided (H1: μi > μ0)", "One-sided (H1: μi < μ0)"),
      select("adjust", "Multiplicity Adjustment",
             {{"dunnett", "Dunnett's test (single-step)", {}},
              {"holm", "Holm's step-down procedure", {}},
              {"hochberg", "Hochberg's step-up procedure", {}},
              {"bonferroni", "Bonferroni's correction", {}},
              {"none", "None (each test at level α)", {}}},
             "dunnett")
          .grp(G_TEST)
          .hint("How the familywise error (the chance of any false claim) is kept at α. Dunnett's test uses the "
                "correlation of the comparisons, which share the control; Holm's and Hochberg's procedures "
                "work on the p-values and are uniformly more powerful than Bonferroni's correction."),
      select("goal", "Power to Find",
             {{"any", "At least one effective treatment (disjunctive)", {}},
              {"all", "All the effective treatments (conjunctive)", {}},
              {"each", "Each effective treatment (the smallest per-arm power)", {}}},
             "any")
          .grp(G_TEST)
          .hint("Which success the power refers to; the other two are reported too. Treatments whose mean "
                "equals the control's count as ineffective, and claiming them is a false positive."),
      power_param(),
      alpha_param("0.05", "Alpha (Familywise Error Rate)"),
      num("n", "n (Subjects per Treatment Arm)", "30")
          .grp(G_N)
          .sym("n")
          .whole()
          .at_least(2)
          .sugg({"30", "10 to 60 by 10"})
          .hint("Number of subjects in each treatment arm.")
          .show(when(not_solving({"n"}))),
      num("r0", "Control Allocation n0/n", "1")
          .grp(G_N)
          .sym("n0/n")
          .range(0.1, 10)
          .sugg({"1", "1.5", "1 1.7"})
          .hint("Size of the control arm relative to each treatment arm (n0 = ratio × n, rounded up). √K, the "
                "square root of the number of treatment arms, is close to optimal for Dunnett's test."),
      vec("mu", "μ1 … μK (Means of the Treatment Arms)", "0.6 0.6 0.3")
          .grp(G_EFFECT)
          .sym("μi")
          .hint("Expected mean outcome in each treatment arm, one value per arm (up to 8 arms)."),
      num("mu0", "μ0 (Mean of the Control Arm)", "0")
          .grp(G_EFFECT)
          .sym("μ0")
          .sugg({"0"})
          .hint("Expected mean outcome in the control arm."),
      num("sigma", "σ (Standard Deviation)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .sugg({"1"})
          .hint("Standard deviation of the outcome within each arm (common to all arms)."),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Design d = design_of(s);
    const int eff = d.effective_count();
    auto sizes = [&](double n) { return Sizes{n, std::ceil(d.r0 * n - 1e-9)}; };
    const bool single = d.adjust == Adjust::Dunnett || d.adjust == Adjust::Bonferroni || d.adjust == Adjust::None;
    // exact power of the single-step test (Bonferroni's for the step-wise procedures, a lower bound)
    auto exact_at = [&](double n, double* crit = nullptr) {
      const Sizes z = sizes(n);
      const double c = critical_value(d, z, ctx);
      const SNodes sn(z.df(d.k), c, least_effect(d, z));
      if (crit) *crit = c;
      return exact_power(d, z, c, sn, ctx);
    };
    auto approx = [&](double n) { return goal_power(d, exact_at(n)); };
    if (s.solve == "n" && eff == 0)
      throw CalcError(d.alt == Alt::Two ? "no treatment mean differs from the control's"
                                        : "no treatment mean differs from the control's in the direction of the "
                                          "one-sided alternative");
    const double lo = s.solve != "n" ? 2 : solve::first_whole([&](double x) { return sizes(x).n0 >= 2; }, 2, 1e6,
                                                           "the control arm's share is too small: it would have "
                                                           "fewer than 2 subjects");
    double n;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = lo;
      o.hi = 1e6;
      n = solve::min_n(approx, s("power"), o, &ctx);
    } else {
      n = s("n");
      if (sizes(n).n0 < 2) throw CalcError("the control arm needs at least 2 subjects");
    }
    double crit = NAN;
    const Exact ex = exact_at(n, &crit);
    std::array<double, KMAX> arm = ex.arm;
    double any = ex.any, all = ex.all, each = ex.each, fp = NAN;
    if (ctx.approximate) {
      simc::put_approximation(s, r, goal_power(d, ex));
    } else {
      const long reps = static_cast<long>(s("reps"));
      const auto seed = static_cast<std::uint64_t>(s("seed"));
      auto power_at = [&](double m, long R, double p0, double p1) {
        const Sizes z = sizes(m);
        const double c = critical_value(d, z, ctx);
        const sim::Summary sm =
            sim::simulate(R, [&](long rep) { return multiarm_study(d, z, c, rep, seed); }, ctx, p0, p1);
        if (d.goal != Goal::Each) return sm;
        long low = sm.reps;  // each: the smallest per-arm rejection count among the effective arms
        for (int i = 0; i < d.k; ++i)
          if (d.effective(i)) low = std::min(low, sm.also[i]);
        return sim::with_rejections(sm, low);
      };
      sim::Summary sm;
      if (s.solve == "n") {
        sim::SizeSearch o;
        o.lo = lo;
        o.hi = 1e6;
        o.guess = n;
        o.reps = reps;
        o.explore_reps = simc::explore_reps(reps);
        const sim::SizeResult res = sim::search_size(power_at, s("power"), o, ctx);
        n = res.n;
        sm = res.s;
      } else {
        sm = power_at(n, reps, 0, 1);
      }
      simc::put_summary(s, r, sm);
      r.set("approx_power", goal_power(d, exact_at(n, &crit)));
      // the three powers and the per-arm rates, from the same simulated studies
      any = sm.also_power(ANY_BIT);
      all = sm.also_power(ALL_BIT);
      each = eff == 0 ? 0 : INFINITY;
      for (int i = 0; i < d.k; ++i) {
        arm[i] = sm.also_power(i);
        if (d.effective(i)) each = std::fmin(each, arm[i]);
      }
      if (eff < d.k) fp = sm.also_power(FALSE_BIT);
    }
    const Sizes z = sizes(n);
    r.set("n", n);
    r.set("n0", z.n0);
    r.set("N", d.k * n + z.n0);
    r.set("k", static_cast<double>(d.k));
    r.set("effective", static_cast<double>(eff));
    r.set("power_any", any);
    r.set("power_all", all);
    r.set("power_each", each);
    if (std::isfinite(fp)) r.set("fwer", fp);
    if (single) r.set("crit", crit);
    r.set("df", z.df(d.k));
    r.set("mu", list_text(std::vector<double>(d.mu.begin(), d.mu.begin() + d.k)));
    r.set("mu0", d.mu0);
    r.set("sigma", d.sigma);
    r.set("alpha", d.alpha);
    Table tb;
    tb.title = "Power for Each Treatment";
    tb.columns = {{"arm", "Arm"}, {"mu", "Mean"}, {"diff", "Difference"}, {"es", "Effect Size"}, {"power", "Power"}};
    tb.formats = {"int", "num", "num", "num", "prob"};
    for (int i = 0; i < d.k; ++i)
      tb.rows.push_back({Cell(static_cast<double>(i + 1)), Cell(d.mu[i]), Cell(d.mu[i] - d.mu0),
                         Cell((d.mu[i] - d.mu0) / d.sigma), Cell(arm[i])});
    r.tables.push_back(std::move(tb));
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    for (auto x : {col("n", "n", "int", "Subjects per treatment arm"), col("n0", "n0", "int", "Subjects in the control arm"),
                   col("N", "N", "int", "Total number of subjects"), col("k", "K", "int", "Treatment arms"),
                   col("effective", "Effective", "int", "Treatment arms whose mean differs from the control's"),
                   col("power_any", "P(Any)", "prob", "Power to find at least one effective treatment"),
                   col("power_all", "P(All)", "prob", "Power to find every effective treatment"),
                   col("power_each", "P(Each)", "prob", "The smallest per-arm power among the effective treatments"),
                   col("fwer", "False Claims", "prob", "Probability of claiming at least one ineffective treatment"),
                   col("crit", "Critical t", "num", "Critical value of the single-step test"),
                   col("df", "df", "num", "Degrees of freedom of the pooled variance"),
                   col("mu", "μ1 … μK", "text", "Means of the treatment arms"), col("mu0", "μ0", "num", "Mean of the control arm"),
                   col("sigma", "σ", "num", "Standard deviation"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns()) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const Design d = design_of(s);
    const std::string goal = d.goal == Goal::All    ? "to find all " + fint(r.num("effective")) + " effective treatments"
                             : d.goal == Goal::Each ? "to find each effective treatment (the smallest per-arm power)"
                                                    : "to find at least one of the effective treatments";
    std::string t = "With " + fint(r.num("n")) + " subjects in each of the " + fint(r.num("k")) +
                    " treatment arms and " + fint(r.num("n0")) + " in the control arm (total " + fint(r.num("N")) +
                    "), the study would " + simc::power_phrase(s, r) + " " + goal + " (means " + text_of(r, "mu") +
                    " against " + fnum(r.num("mu0")) + " for the control; σ = " + fnum(r.num("sigma")) +
                    "), comparing each treatment with the control by a " + alt_words(d.alt) + " t-test with " +
                    adjust_words(d.adjust);
    if (r.has("crit") && std::isfinite(r.num("crit"))) t += " (critical value " + fnum(r.num("crit"), 3) + ")";
    t += " at a familywise α of " + fnum(r.num("alpha")) + ".";
    if (d.goal != Goal::Any) t += " The power to find at least one is " + fpct(r.num("power_any")) + ".";
    if (d.goal != Goal::All && r.num("effective") > 1) t += " The power to find all is " + fpct(r.num("power_all")) + ".";
    return t + simc::method_note(r);
  };
  p.definitions = {{"n, n0", "Subjects in each treatment arm and in the control arm."},
                   {"P(Any)", "Disjunctive power: the probability of finding at least one effective treatment."},
                   {"P(All)", "Conjunctive power: the probability of finding every effective treatment."},
                   {"P(Each)", "The smallest of the effective treatments' own powers."},
                   {"False Claims", "Probability of declaring at least one ineffective treatment better than the "
                                    "control (at most α)."},
                   {"Power", "Share of the simulated studies that met the chosen goal, with its 95% interval "
                             "(Wilson) and Monte Carlo standard error."}};
  p.notes = {"Each treatment is compared with the control by a t statistic on the pooled variance of all arms, "
             "with N − K − 1 degrees of freedom; the comparisons are correlated because they share the control.",
             "For Dunnett's test, Bonferroni's correction and unadjusted tests the power is also computed exactly "
             "(the analytic value); for Holm's and Hochberg's procedures the analytic value is Bonferroni's power, "
             "a lower bound.",
             "With no effective treatment, the power to find at least one is the familywise error rate.",
             "Hochberg's procedure assumes positively dependent tests, which comparisons with a shared control "
             "are."};
  p.references = {"Dunnett, C. W. (1955). A multiple comparison procedure for comparing several treatments with a "
                  "control. Journal of the American Statistical Association 50, 1096–1121.",
                  "Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal "
                  "of Statistics 6, 65–70.",
                  "Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. "
                  "Biometrika 75, 800–802.",
                  "Horn, M. & Vollandt, R. (1995). Multiple Tests und Auswahlverfahren. Gustav Fischer (sample sizes "
                  "for many-to-one comparisons).",
                  simc::REF_SIM};
  p.curve_param = "n";
  p.curve_label = "n (Subjects per Treatment Arm)";
  p.example_solve = "n";
  p.example = {{"power", "0.8"}, {"mu", "0.6 0.6 0.3"}, {"mu0", "0"}, {"sigma", "1"}, {"reps", "1000"}};
  return p;
}
SSC_REGISTER(multiarm);

}  // namespace
}  // namespace ssc::procs
