// Simulation-based power for comparing two means when the data are not
// normal (sim.means.robust): skewed (log-normal, gamma), heavy-tailed (t),
// contaminated by outliers, or light-tailed (uniform), with possibly unequal
// SDs.  Each group's outcomes are μ + σ·Z for a standardized variate Z (mean
// 0, SD 1) of the chosen shape, so μ1 − μ2 is the difference in means.
//
// Every simulated study is analyzed with the pooled t-test, Welch's t-test
// and the Mann–Whitney (Wilcoxon rank-sum) test; the chosen test decides the
// power and the sample size, and the other two are reported for the same
// studies.  The quick approximation is the normal-theory power of the t-tests
// and Noether's formula for the Mann–Whitney test with P(X1 > X2) integrated
// numerically for the chosen shape.
#include <algorithm>
#include <cstdint>
#include <utility>

#include "procedures/sim_common.hpp"
#include "stats/numerics.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;
using sim::Rng;

constexpr double SQRT12 = 3.4641016151377545870;

/// A standardized variate Z (mean 0, SD 1) of one shape.
struct Shape {
  enum Kind { Normal, Lognormal, Gamma, T, Contaminated, Uniform } kind = Normal;
  double sign = 1;                  // −1 mirrors a skewed shape (skewed to the left)
  double s = 0, lm = 1, lsd = 1;    // log-normal exp(sN): SD of the log, mean and SD
  double a = 0;                     // gamma shape
  double df = 0, tscale = 1;        // t: degrees of freedom and √((ν − 2)/ν)
  double eps = 0, k = 1, c = 1;     // contaminated: share of outliers, their SD multiple, the overall SD
  double skew = 0, exkurt = 0;      // skewness and excess kurtosis
  std::string words;                // "a log-normal distribution (skewness 2)"

  double draw(Rng& rng) const {
    switch (kind) {
      case Lognormal: return sign * (std::exp(s * rng.normal()) - lm) / lsd;
      case Gamma: return sign * (rng.gamma(a) - a) / std::sqrt(a);
      case T: return tscale * rng.normal() / std::sqrt(rng.chisq(df) / df);
      case Contaminated: {
        const bool outlier = rng.uniform() < eps;
        return (outlier ? k : 1) * rng.normal() / c;
      }
      case Uniform: return (rng.uniform() - 0.5) * SQRT12;
      default: return rng.normal();
    }
  }

  /// P(Z ≤ z) of the right-skewed (unmirrored) version.
  double cdf_right(double z) const {
    if (kind == Lognormal) {
      const double x = lm + lsd * z;
      return x <= 0 ? 0 : stats::pnorm(std::log(x) / s);
    }
    const double x = a + std::sqrt(a) * z;  // gamma
    return x <= 0 ? 0 : stats::pgamma_raw(x, a);
  }

  double cdf(double z) const {
    switch (kind) {
      case Lognormal:
      case Gamma: return sign > 0 ? cdf_right(z) : 1 - cdf_right(-z);
      case T: return stats::pt(z / tscale, df);
      case Contaminated: return (1 - eps) * stats::pnorm(z * c) + eps * stats::pnorm(z * c / k);
      case Uniform: return std::clamp(z / SQRT12 + 0.5, 0.0, 1.0);
      default: return stats::pnorm(z);
    }
  }

  double quantile_right(double u) const {
    if (kind == Lognormal) return (std::exp(s * stats::qnorm(u)) - lm) / lsd;
    return (stats::qgamma(u, a) - a) / std::sqrt(a);
  }

  double quantile(double u) const {
    switch (kind) {
      case Lognormal:
      case Gamma: return sign > 0 ? quantile_right(u) : -quantile_right(1 - u);
      case T: return tscale * stats::qt(u, df);
      case Uniform: return (u - 0.5) * SQRT12;
      case Contaminated: {
        double lo = -1, hi = 1;
        while (cdf(lo) > u) lo *= 2;
        while (cdf(hi) < u) hi *= 2;
        return numerics::brent([&](double z) { return cdf(z) - u; }, lo, hi, 1e-13);
      }
      default: return stats::qnorm(u);
    }
  }
};

/// SD of the log of a log-normal variate with skewness γ > 0: (w + 2)√(w − 1) = γ with w = exp(s²).
double lognormal_sdlog(double gamma) {
  const double w1 = numerics::brent([&](double v) { return (v + 3) * std::sqrt(v) - gamma; }, 0, 1e6, 1e-15);
  return std::sqrt(std::log1p(w1));
}

Shape shape_of(const Scenario& s) {
  Shape z;
  const std::string& d = s.s("dist");
  auto skew_words = [](double g) {
    return "skewness " + fnum(g, 3) + (g < 0 ? ", to the left" : "");
  };
  if (d == "lognormal" || d == "gamma") {
    const double g = s("skew");
    // below 1e-6 the skewed shapes' quantiles are rounding noise around the normal ones (and the integral of
    // P(X1 > X2) over them never meets its tolerance)
    if (std::fabs(g) < 1e-6) {
      z.words = "a normal distribution";
      return z;
    }
    z.sign = g < 0 ? -1 : 1;
    z.skew = g;
    if (d == "lognormal") {
      z.kind = Shape::Lognormal;
      z.s = lognormal_sdlog(std::fabs(g));
      const double w = std::exp(z.s * z.s);
      z.lm = std::exp(z.s * z.s / 2);
      z.lsd = std::sqrt((w - 1) * w);
      z.exkurt = w * w * w * w + 2 * w * w * w + 3 * w * w - 6;
      z.words = "a log-normal distribution (" + skew_words(g) + ")";
    } else {
      z.kind = Shape::Gamma;
      z.a = 4 / (g * g);
      z.exkurt = 6 / z.a;
      z.words = "a gamma distribution (" + skew_words(g) + ")";
    }
  } else if (d == "t") {
    z.kind = Shape::T;
    z.df = s("df");
    z.tscale = std::sqrt((z.df - 2) / z.df);
    z.exkurt = z.df > 4 ? 6 / (z.df - 4) : INFINITY;
    z.words = "a t distribution with " + fnum(z.df) + " degrees of freedom (heavy tails)";
  } else if (d == "contaminated") {
    z.kind = Shape::Contaminated;
    z.eps = s("eps");
    z.k = s("kout");
    const double var = (1 - z.eps) + z.eps * z.k * z.k;
    // past about 1.3e154 the outliers' variance passes the largest double: every value would be standardized to 0
    // (k² itself, as a tiny ε keeps the mixed variance finite while k² overflows)
    if (!std::isfinite(z.k * z.k))
      throw CalcError("the outliers' SD multiple k is too large to simulate (its square passes the largest number)");
    z.c = std::sqrt(var);
    z.exkurt = 3 * ((1 - z.eps) + z.eps * z.k * z.k * z.k * z.k) / (var * var) - 3;
    if (!std::isfinite(z.exkurt)) {  // k^4 passes the largest double where k^2 does not
      const double w = z.k * z.k / var;
      z.exkurt = 3 * ((1 - z.eps) / var / var + z.eps * w * w) - 3;
    }
    z.words = "a normal distribution with " + fpct(z.eps, 0) + " outliers (SD × " + fnum(z.k) + ")";
  } else if (d == "uniform") {
    z.kind = Shape::Uniform;
    z.exkurt = -1.2;
    z.words = "a uniform distribution (light tails)";
  } else {
    z.words = "a normal distribution";
  }
  return z;
}

struct Design {
  Shape shape;
  double mu1 = 0, mu2 = 0, sd1 = 1, sd2 = 1;
};

/// P(X1 > X2) for X_g = μ_g + σ_g Z_g: the mean over Z1's quantiles of P(Z2 < (μ1 − μ2 + σ1 Z1)/σ2).
double prob_greater(const Design& d) {
  const double delta = d.mu1 - d.mu2;
  return numerics::integrate(
      [&](double u) { return d.shape.cdf((delta + d.sd1 * d.shape.quantile(u)) / d.sd2); }, 0, 1, 1e-10, 1e-9,
      20);  // at most about 2 million panels, whatever the integrand
}

enum Test { POOLED = 0, WELCH = 1, RANKS = 2 };

sim::Outcome robust_study(const Design& d, long n1, long n2, long rep, std::uint64_t seed, Test test, double alpha,
                          Alt a) {
  // each group's subjects are drawn in order from one stream, so a larger study contains a smaller one
  std::vector<std::pair<double, int>> all;  // freed with the study
  all.reserve(static_cast<std::size_t>(n1 + n2));
  double sum[2] = {0, 0}, ss[2] = {0, 0};
  const long n[2] = {n1, n2};
  sim::Outcome o;
  for (int g = 0; g < 2; ++g) {
    Rng rng(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g + 1));
    const double mu = g == 0 ? d.mu1 : d.mu2, sd = g == 0 ? d.sd1 : d.sd2;
    // Welford's update keeps the variance accurate for any location
    double mean = 0, m2 = 0;
    for (long i = 0; i < n[g]; ++i) {
      const double x = mu + sd * d.shape.draw(rng);
      // a value past the largest double (an outlier SD multiple near 1e308 gives inf / inf) leaves no study to
      // analyze: the ranks of a NaN never end, as its ties are never equal
      if (!std::isfinite(x)) {
        o.failed = true;
        return o;
      }
      all.emplace_back(x, g);
      const double delta = x - mean;
      mean += delta / static_cast<double>(i + 1);
      m2 += delta * (x - mean);
    }
    sum[g] = mean;
    ss[g] = m2;
  }
  const double a1 = static_cast<double>(n1), a2 = static_cast<double>(n2);
  const double diff = sum[0] - sum[1], v1 = ss[0] / (a1 - 1), v2 = ss[1] / (a2 - 1);
  // pooled t
  const double dfp = a1 + a2 - 2, sp2 = (ss[0] + ss[1]) / dfp, sep = std::sqrt(sp2 * (1 / a1 + 1 / a2));
  const bool rej_pooled = sep > 0 && simc::p_value(simc::t_stat(diff, sep), dfp, a) < alpha;
  // Welch
  const double q1 = v1 / a1, q2 = v2 / a2, sew = std::sqrt(q1 + q2);
  const double dfw = (q1 + q2) * (q1 + q2) / (q1 * q1 / (a1 - 1) + q2 * q2 / (a2 - 1));
  const bool rej_welch = sew > 0 && simc::p_value(simc::t_stat(diff, sew), dfw, a) < alpha;
  // Mann–Whitney with midranks
  std::sort(all.begin(), all.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
  double r1 = 0, ties = 0;
  for (std::size_t i = 0, j = 0; i < all.size(); i = j) {
    long in1 = 0;
    while (j < all.size() && all[j].first == all[i].first) in1 += all[j++].second == 0;
    const double t = static_cast<double>(j - i), rank = (static_cast<double>(i + j) + 1) / 2;
    r1 += rank * static_cast<double>(in1);
    ties += t * t * t - t;
  }
  const bool rej_ranks = simc::rank_sum_p(r1 - a1 * (a1 + 1) / 2, a1, a2, ties, a) < alpha;
  o.also = (rej_pooled ? 1u : 0u) | (rej_welch ? 2u : 0u) | (rej_ranks ? 4u : 0u);
  o.reject = test == POOLED ? rej_pooled : (test == WELCH ? rej_welch : rej_ranks);
  o.estimate = diff;
  o.se = test == POOLED ? sep : sew;
  if (test != RANKS) o.df = test == POOLED ? dfp : dfw;
  return o;
}

Procedure robust() {
  Procedure p;
  p.id = "sim.means.robust";
  p.name = "Two Means with Skewed or Heavy-Tailed Data (t, Welch, Mann–Whitney), by Simulation";
  p.category = "Simulation";
  p.subcategory = "Means";
  p.order = 1;
  p.simulation = true;
  p.summary =
      "Sample size (or power) for comparing two groups when the outcome is skewed, heavy-tailed or has outliers, "
      "estimated by simulating studies from the chosen distribution and analyzing each with the pooled t-test, "
      "Welch's t-test and the Mann–Whitney test, so the tests can be compared on the same data.";
  p.keywords = {"simulation", "monte carlo", "non-normal", "nonnormal", "skewed", "skewness", "log-normal",
                "lognormal", "gamma", "heavy tails", "heavy-tailed", "outliers", "contaminated", "robust",
                "t-test", "welch", "unequal variances", "mann-whitney", "wilcoxon rank-sum", "rank sum",
                "nonparametric", "distribution", "kurtosis", "costs", "length of stay"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total Sample Size)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: μ1 ≠ μ2)", "One-sided (H1: μ1 > μ2)", "One-sided (H1: μ1 < μ2)"),
      select("test", "Test",
             {{"welch", "Welch's t-test (unequal variances)", {}},
              {"t", "Pooled t-test (equal variances)", {}},
              {"mw", "Mann–Whitney (Wilcoxon rank-sum) test", {}}},
             "welch")
          .grp(G_TEST)
          .hint("The test that decides the power and the sample size; the other two are run on the same simulated "
                "studies and reported alongside. The Mann–Whitney test compares P(X1 > X2) with 1/2, not the "
                "means."),
      power_param(),
      alpha_param(),
      num("n1", "N1 (Subjects in Group 1)", "30")
          .grp(G_N)
          .sym("N1")
          .whole()
          .at_least(2)
          .sugg({"30", "10 to 100 by 10"})
          .hint("Number of subjects in group 1; group 2 gets N1 × the allocation ratio, rounded up.")
          .show(when(not_solving({"n"}))),
      num("r", "Allocation Ratio N2/N1", "1")
          .grp(G_N)
          .sym("N2/N1")
          .above(0)
          .sugg({"1", "2"})
          .hint("Size of group 2 relative to group 1 (N2 = ratio × N1, rounded up); 1 gives equal groups."),
      num("mu1", "μ1 (Mean of Group 1)", "1")
          .grp(G_EFFECT)
          .sym("μ1")
          .sugg({"1", "0.5 1"})
          .hint("Mean outcome in group 1 (e.g. the treatment) under the alternative."),
      num("mu2", "μ2 (Mean of Group 2)", "0")
          .grp(G_EFFECT)
          .sym("μ2")
          .sugg({"0"})
          .hint("Mean outcome in group 2 (e.g. the control); μ1 − μ2 is the difference to detect."),
      num("sd1", "σ1 (SD of Group 1)", "2")
          .grp(G_EFFECT)
          .sym("σ1")
          .above(0)
          .sugg({"2", "1 2"})
          .hint("Standard deviation of the outcome in group 1."),
      num("sd2", "σ2 (SD of Group 2)", "2")
          .grp(G_EFFECT)
          .sym("σ2")
          .above(0)
          .sugg({"2", "1 2"})
          .hint("Standard deviation of the outcome in group 2; unequal SDs affect the pooled t-test most."),
      select("dist", "Distribution of the Outcome",
             {{"lognormal", "Log-normal (skewed)", {}},
              {"gamma", "Gamma (skewed)", {}},
              {"t", "t distribution (heavy tails)", {}},
              {"contaminated", "Normal with outliers (contaminated normal)", {}},
              {"uniform", "Uniform (light tails)", {}},
              {"normal", "Normal", {}}},
             "lognormal")
          .grp("Distribution")
          .hint("The shape of the outcome's distribution in both groups; each group's outcomes are its mean plus "
                "its SD times a variate of this shape with mean 0 and SD 1."),
      num("skew", "γ (Skewness)", "2")
          .grp("Distribution")
          .sym("γ")
          .range(-20, 20)
          .sugg({"2", "1 2 4"})
          .hint("Skewness of the outcome: 0 is symmetric, 1 is moderate and 2 or more strong (costs, lengths of "
                "stay, counts of cells); negative values skew to the left.")
          .show(when(is("dist", {"lognormal", "gamma"}))),
      num("df", "ν (Degrees of Freedom of the t Distribution)", "4")
          .grp("Distribution")
          .sym("ν")
          .above(2)
          .sugg({"4", "3 5 10"})
          .hint("Smaller values give heavier tails (more extreme values); above 30 the t distribution is almost "
                "normal. It must exceed 2 for the SD to exist.")
          .show(when(is("dist", {"t"}))),
      num("eps", "ε (Share of Outliers)", "0.1")
          .grp("Distribution")
          .sym("ε")
          .range(0, 1, false, true)
          .sugg({"0.1", "0.05 0.1 0.2"})
          .hint("Proportion of subjects whose outcome comes from a wider normal distribution.")
          .show(when(is("dist", {"contaminated"}))),
      num("kout", "k (SD of the Outliers, × the Others' SD)", "5")
          .grp("Distribution")
          .sym("k")
          .at_least(1)
          .sugg({"5", "3 5 10"})
          .hint("How many times wider the outliers' distribution is than the others'; the overall SD stays σ.")
          .show(when(is("dist", {"contaminated"}))),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double alpha = s("alpha"), ratio = s("r");
    Design d;
    d.shape = shape_of(s);
    d.mu1 = s("mu1");
    d.mu2 = s("mu2");
    d.sd1 = s("sd1");
    d.sd2 = s("sd2");
    // the three tests are unchanged by a common change of location and scale: extreme values are simulated on the
    // scale of the larger SD, where their squares cannot overflow, and around μ2 (group 1 around μ1 − μ2), so
    // that means far from 0 next to the SDs do not swamp the variation (μ1 = μ2 = 10¹⁶ with SDs of 2 gave the
    // t-tests a type I error of 35%, and μ1 − μ2 = 2 there a power of 76% for 97%); others exactly as given
    Design w = d;
    double unit = 1;
    const double sd_min = std::fmin(d.sd1, d.sd2), mu_max = std::fmax(std::fabs(d.mu1), std::fabs(d.mu2));
    if (!(std::fmax(d.sd1, d.sd2) < 1e100 && sd_min > 1e-100 && mu_max < 1e100 && mu_max < 1e6 * sd_min)) {
      unit = std::fmax(d.sd1, d.sd2);
      w.mu1 = (d.mu1 - d.mu2) / unit;
      w.mu2 = 0;
      w.sd1 = d.sd1 / unit;
      w.sd2 = d.sd2 / unit;
    }
    const Test test = s.is("test", "t") ? POOLED : (s.is("test", "mw") ? RANKS : WELCH);
    const double delta = d.mu1 - d.mu2, pgt = prob_greater(w);
    auto n2_of = [&](double n1) { return std::ceil(ratio * n1 - 1e-9); };
    auto approx = [&](double n1) {
      const double n2 = n2_of(n1);
      if (test == RANKS) return z_power((pgt - 0.5) / std::sqrt((n1 + n2 + 1) / (12 * n1 * n2)), alpha, a);
      const double q1 = w.sd1 * w.sd1 / n1, q2 = w.sd2 * w.sd2 / n2, ncp = (w.mu1 - w.mu2) / std::sqrt(q1 + q2);
      const double df = test == POOLED ? n1 + n2 - 2 : (q1 + q2) * (q1 + q2) / (q1 * q1 / (n1 - 1) + q2 * q2 / (n2 - 1));
      return t_power(ncp, df, alpha, a);
    };
    if (s.solve == "n") {
      const double effect = test == RANKS ? pgt - 0.5 : delta;
      if (std::fabs(effect) < 1e-12)
        throw CalcError(test == RANKS ? "P(X1 > X2) is 1/2; the Mann–Whitney test has no power beyond alpha"
                                      : "μ1 equals μ2; power cannot exceed alpha");
      if ((a == Alt::Upper && effect < 0) || (a == Alt::Lower && effect > 0))
        throw CalcError(test == RANKS ? "P(X1 > X2) − 1/2 has the sign opposite to the one-sided alternative"
                                      : "μ1 − μ2 has the sign opposite to the one-sided alternative");
    }
    const double lo = s.solve != "n" ? 2 : solve::first_whole([&](double x) { return n2_of(x) >= 2; }, 2, 1e5,
                                                           "the allocation ratio is too small: group 2 would have "
                                                           "fewer than 2 subjects");
    double n1;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = lo;
      o.hi = 1e5;
      n1 = solve::min_n(approx, s("power"), o, &ctx);
    } else {
      n1 = s("n1");
      if (n2_of(n1) < 2) throw CalcError("group 2 needs at least 2 subjects");
    }
    if (ctx.approximate) {
      simc::put_approximation(s, r, approx(n1));
    } else {
      const long reps = static_cast<long>(s("reps"));
      const auto seed = static_cast<std::uint64_t>(s("seed"));
      auto power_at = [&](double m1, long R, double p0, double p1) {
        simc::check_study_size(m1 + n2_of(m1), "subjects");
        const long k1 = static_cast<long>(m1), k2 = static_cast<long>(n2_of(m1));
        return sim::simulate(R, [&](long rep) { return robust_study(w, k1, k2, rep, seed, test, alpha, a); }, ctx,
                             p0, p1);
      };
      sim::Summary sm;
      if (s.solve == "n") {
        sim::SizeSearch o;
        o.lo = lo;
        o.hi = 1e5;
        o.guess = n1;
        o.reps = reps;
        o.explore_reps = simc::explore_reps(reps);
        const sim::SizeResult res = sim::search_size(power_at, s("power"), o, ctx);
        n1 = res.n;
        sm = res.s;
      } else {
        sm = power_at(n1, reps, 0, 1);
      }
      simc::put_summary(s, r, sm, simc::Failures::Numerical);
      r.set("approx_power", approx(n1));
      r.set("power_t", sm.also_power(POOLED));
      r.set("power_welch", sm.also_power(WELCH));
      r.set("power_mw", sm.also_power(RANKS));
      if (std::isfinite(sm.mean_estimate)) r.set("mean_delta", sm.mean_estimate * unit);
    }
    put_two_group(r, {n1, n2_of(n1)});
    r.set("mu1", d.mu1);
    r.set("mu2", d.mu2);
    r.set("delta", delta);
    r.set("sd1", d.sd1);
    r.set("sd2", d.sd2);
    r.set("dist", d.shape.words);
    r.set("skewness", d.shape.skew);
    r.set("exkurt", d.shape.exkurt);
    r.set("pgt", pgt);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("mu1", "μ1", "num", "Mean of group 1"), col("mu2", "μ2", "num", "Mean of group 2"),
                   col("delta", "μ1 − μ2", "num", "Difference in means"),
                   col("mean_delta", "Mean Diff.", "num", "Average observed difference over the simulated studies"),
                   col("sd1", "σ1", "num", "SD of group 1"), col("sd2", "σ2", "num", "SD of group 2"),
                   col("skewness", "Skewness", "num", "Skewness of the outcome's distribution"),
                   col("exkurt", "Excess Kurtosis", "num", "Excess kurtosis of the outcome's distribution (0 normal)"),
                   col("pgt", "P(X1 > X2)", "prob", "Probability that a subject of group 1 has the larger outcome"),
                   col("power_t", "Power (t)", "prob", "Simulated power of the pooled t-test on the same studies"),
                   col("power_welch", "Power (Welch)", "prob", "Simulated power of Welch's t-test on the same studies"),
                   col("power_mw", "Power (M–W)", "prob", "Simulated power of the Mann–Whitney test on the same studies"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns()) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string test = s.is("test", "t")    ? "pooled t-test"
                             : s.is("test", "mw") ? "Mann–Whitney test"
                                                  : "Welch t-test";
    std::string t = "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " (total " +
                    fint(r.num("n")) + ") " + simc::power_phrase(s, r) + " to detect a difference of " +
                    fnum(r.num("delta")) + " between the means (μ1 = " + fnum(r.num("mu1")) + ", μ2 = " +
                    fnum(r.num("mu2")) + "; SDs " + fnum(r.num("sd1")) + " and " + fnum(r.num("sd2")) +
                    ") when the outcome follows " + text_of(r, "dist") + ", using a " + alt_words(alt(s)) + " " +
                    test + " at α = " + fnum(r.num("alpha")) + ".";
    if (r.has("power_t") && std::isfinite(r.num("power_t")))
      t += " On the same simulated studies the pooled t-test had a power of " + fpct(r.num("power_t")) +
           ", Welch's t-test " + fpct(r.num("power_welch")) + " and the Mann–Whitney test " +
           fpct(r.num("power_mw")) + ".";
    return t + simc::method_note(r);
  };
  p.definitions = {{"N1, N2", "Number of subjects in groups 1 and 2; N = N1 + N2."},
                   {"μ1, μ2", "Means of the outcome in the two groups."},
                   {"σ1, σ2", "Standard deviations of the outcome in the two groups."},
                   {"Skewness", "Asymmetry of the distribution: 0 for a symmetric one, positive for a long right "
                                "tail."},
                   {"P(X1 > X2)", "Probability that a random subject of group 1 has a larger outcome than one of "
                                  "group 2: the quantity the Mann–Whitney test compares with 1/2."},
                   {"Power", "Share of the simulated studies whose test was significant, with its 95% interval "
                             "(Wilson) and Monte Carlo standard error."}};
  p.notes = {"Each group's outcomes are its mean plus its SD times a variate of the chosen shape standardized to "
             "mean 0 and SD 1, so the groups differ in location (and in scale when the SDs differ).",
             "The Mann–Whitney test asks whether P(X1 > X2) = 1/2. With skewed data and unequal SDs this "
             "probability differs from 1/2 even when the means are equal, so the test can reject without a "
             "difference in means. Set μ1 = μ2 to see each test's type I error for your design.",
             "With unequal SDs and unequal group sizes the pooled t-test's type I error drifts away from α; "
             "Welch's test keeps it close.",
             "The Mann–Whitney p-value uses the normal approximation with the tie and continuity corrections (as "
             "R's wilcox.test with exact = FALSE).",
             "The quick approximation is the normal-theory power of the t-tests and Noether's formula for the "
             "Mann–Whitney test, with P(X1 > X2) integrated numerically for the chosen distribution."};
  p.references = {"Noether, G. E. (1987). Sample size determination for some common nonparametric tests. Journal "
                  "of the American Statistical Association 82, 645–647.",
                  "Fagerland, M. W. & Sandvik, L. (2009). The Wilcoxon–Mann–Whitney test under scrutiny. Statistics "
                  "in Medicine 28, 1487–1497.",
                  "Welch, B. L. (1947). The generalization of 'Student's' problem when several different population "
                  "variances are involved. Biometrika 34, 28–35.",
                  simc::REF_SIM};
  p.curve_param = "n1";
  p.curve_label = "N1 (Subjects in Group 1)";
  p.example_solve = "n";
  p.example = {{"power", "0.8"}, {"mu1", "1"}, {"mu2", "0"}, {"sd1", "2"}, {"sd2", "2"}, {"dist", "lognormal"},
               {"skew", "2"}, {"reps", "500"}};
  return p;
}
SSC_REGISTER(robust);

}  // namespace
}  // namespace ssc::procs
