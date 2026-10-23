// Simulation-based power for linear mixed models:
//  * sim.lmm.longitudinal: two groups measured repeatedly, compared on their
//    mean slopes (group × time interaction), with random intercepts (and
//    slopes) and monotone dropout;
//  * sim.lmm.cluster: cluster-randomized comparison of two means with a random
//    cluster intercept and cluster sizes that may vary.
//
// Each simulated study is analyzed as it would be in practice: a REML fit of
// the mixed model (sim/lmm.hpp) and a t test of the effect with Satterthwaite's
// (or a fixed) degrees of freedom.  The quick approximation used for live
// previews and power curves is the generalized least squares variance of the
// effect at the assumed variances; for complete balanced designs it equals the
// closed-form procedures means.two.slopes and means.two.cluster.
#include <cstdint>
#include <utility>

#include "procedures/sim_common.hpp"
#include "sim/lmm.hpp"
#include "sim/random.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;
using sim::Rng;
namespace lmm = sim::lmm;
using simc::df_method;
using simc::p_value;

/// The sentence on the simulated studies whose Satterthwaite degrees of freedom could not be computed, so that
/// their test used `fixed` degrees of freedom instead.
std::string fallback_note(const Row& r, const std::string& fixed) {
  if (!r.has("df_fallbacks")) return "";
  return " In " + fint(r.num("df_fallbacks")) + " of the " + fint(r.num("reps")) +
         " simulated studies, Satterthwaite's degrees of freedom could not be computed, and the test used " + fixed +
         ".";
}

// =====================================================================================
// Longitudinal: difference between two mean slopes
// =====================================================================================
// The analysis does not depend on the time origin and unit or on the outcome's unit (the estimates are only
// rescaled), so each study is simulated and fitted with time measured from the visit nearest to time 0 (from 0
// when the visits span it) in units of the farthest visit, t' = (t − t0)/τ, and the outcome in units of its
// largest SD: a follow-up of 10⁶ time units, visits 10⁶ time units after the origin, or an intercept SD of 10¹⁶⁰
// is then no different from the defaults.  With random slopes, t0 moves to where the estimates of a group's mean
// level and mean slope are uncorrelated when at the visit nearest 0 they are nearly collinear (slopes far more
// variable than the intercepts at time 0, visits away from it).  When the visits span 0, or all precede it, or t0
// has moved, the rows of the subjects' Z (with random slopes) measure time from the first visit instead,
// t'' = (t − t1)/τ, where Z'Z is far from singular whatever the subject's visits: from t0, two early visits close
// together and far from it (times −1000001, −1000000 and 0) would have made Z'Z look singular and lost the
// information in their difference.
struct LongDesign {
  std::vector<double> times;
  bool slopes = true;  // random intercepts and slopes (otherwise intercepts only)
  double sd0 = 1, sd1 = 0, rho = 0, sde = 1;
  double hazard = 0;  // probability of dropping out after each visit
  double l00 = 1, l10 = 0, l11 = 0;  // Cholesky factor of D
  double t0 = 0, t1 = 0, tau = 1, unit = 1;  // time origins of the analysis and of Z's rows, time unit, outcome unit
  // in those units: the residual SD and a factor M of the random effects' covariance, so that the intercept at
  // t0 and the slope per τ are (m00 z0 + m01 z1, m10 z0 + m11 z1) for standard normal z
  double se1 = 1, m00 = 1, m01 = 0, m10 = 0, m11 = 0;
  int q() const { return slopes ? 2 : 1; }
  double time(std::size_t j) const { return (times[j] - t0) / tau; }
  double z_time(std::size_t j) const { return (times[j] - t1) / tau; }
};

std::vector<double> visit_times(const Scenario& s) {
  if (s.is("timing", "custom")) {
    std::vector<double> t = s.v("times");
    if (t.size() < 2) throw CalcError("enter at least two measurement times");
    std::sort(t.begin(), t.end());
    if (t.front() == t.back()) throw CalcError("the measurement times must not all be equal");
    if (t.size() > 50) throw CalcError("at most 50 measurement times");
    return t;
  }
  const int m = static_cast<int>(s("m"));
  if (m > 50) throw CalcError("at most 50 measurements per subject");
  std::vector<double> t(static_cast<std::size_t>(m));
  for (int j = 0; j < m; ++j) t[static_cast<std::size_t>(j)] = s("duration") * j / (m - 1);
  return t;
}

/// Probability that a subject is measured at exactly the first k visits (k = 1..K).
std::vector<double> visit_probs(const LongDesign& d) {
  const int K = static_cast<int>(d.times.size());
  std::vector<double> p(static_cast<std::size_t>(K));
  for (int k = 1; k <= K; ++k)
    p[static_cast<std::size_t>(k - 1)] = std::pow(1 - d.hazard, k - 1) * (k < K ? d.hazard : 1);
  return p;
}

/// Set the time origin t0 of the analysis, with the time unit τ (the distance to the farthest visit), the outcome
/// unit (the largest SD) and, in those units, the residual SD and M.
void set_origin(LongDesign& d, double t0) {
  d.t0 = t0;
  d.tau = std::fmax(std::fabs(d.times.front() - d.t0), std::fabs(d.times.back() - d.t0));
  // b0 + b1 t = (b0 + b1 t0) + (b1 τ) t'
  const double intercept0 = d.l00 + d.t0 * d.l10, intercept1 = d.t0 * d.l11;
  d.unit = std::fmax(std::fmax(std::hypot(intercept0, intercept1), d.sde), d.sd1 * d.tau);
  d.se1 = d.sde / d.unit;
  d.m00 = intercept0 / d.unit, d.m01 = intercept1 / d.unit;
  d.m10 = d.tau * d.l10 / d.unit, d.m11 = d.tau * d.l11 / d.unit;
}

/// With random slopes: how far (in units of τ) to move the time origin so that a group's estimated mean level and
/// mean slope are uncorrelated at the assumed variances, and (share) 1 − their squared correlation at t0.  Far from
/// the time where the subjects' lines cross with the least spread, both estimates follow the slope's: with slopes
/// 10⁷ times as variable as the intercepts at time 0 and visits at 1, 2 and 4, X'V⁻¹X at the first visit is singular
/// to within 10⁻¹³ and its inverse has lost all but three digits.  From each pattern of visits (the first k), with
/// mean t̄ and spread S = Σ (t − t̄)², E = σ² (Z'Z)⁻¹ + G is the covariance of a subject's least-squares level and
/// slope, |E| = σ⁴/(k S) + σ² (v/S + s²/k) + |M|² (v and s² the variances of the level at t̄ and of the slope),
/// a sum of positive terms, and the information is Σ_k p_k E⁻¹ (a subject seen at one time tells only the level
/// there).  The visits' mean and spread are taken from the visit nearest t0 (from t0 itself when the visits span
/// it): from a t0 far beyond them, the visits' times lose their spacing (lines that cross 10¹⁸ time units before
/// visits at 0, 1 and 2 put all three at 1 in units of τ, and the origin went back to 0).
double steady_shift(const LongDesign& d, double& share) {
  const std::vector<double> prob = visit_probs(d);
  const double s2 = d.se1 * d.se1, slope2 = d.m10 * d.m10 + d.m11 * d.m11, det_m = d.m00 * d.m11 - d.m01 * d.m10;
  const double g_aa = d.m00 * d.m00 + d.m01 * d.m01, g_as = d.m00 * d.m10 + d.m01 * d.m11;
  const double ref = std::clamp(d.t0, d.times.front(), d.times.back()), offset = (ref - d.t0) / d.tau;
  double maa = 0, mas = 0, mss = 0;
  for (std::size_t k = 1; k <= d.times.size(); ++k) {
    const double p = prob[k - 1], dk = static_cast<double>(k);
    double mean = 0, ss = 0;
    for (std::size_t j = 0; j < k; ++j) mean += (d.times[j] - ref) / d.tau / dk;
    for (std::size_t j = 0; j < k; ++j) {
      const double t = (d.times[j] - ref) / d.tau;
      ss += (t - mean) * (t - mean);
    }
    mean = offset + mean;  // from t0
    const double la = d.m00 + mean * d.m10, lb = d.m01 + mean * d.m11, v = la * la + lb * lb;
    if (ss > 0) {
      const double det = s2 * s2 / (dk * ss) + s2 * (v / ss + slope2 / dk) + det_m * det_m;
      const double e_aa = s2 * (1 / dk + mean * mean / ss) + g_aa, e_as = g_as - s2 * mean / ss;
      const double e_ss = s2 / ss + slope2;
      maa += p * e_ss / det, mas -= p * e_as / det, mss += p * e_aa / det;
    } else {
      const double w = p / (s2 / dk + v);
      maa += w, mas += w * mean, mss += w * mean * mean;
    }
  }
  share = 1 - mas / maa * (mas / mss);
  return mas / maa;
}

LongDesign long_design(const Scenario& s) {
  LongDesign d;
  d.times = visit_times(s);
  d.slopes = !s.is("re", "intercepts");
  d.sd0 = s("sd0");
  d.sde = s("sde");
  if (d.slopes) {
    d.sd1 = s("sd1");
    d.rho = s("rho");
  }
  d.l00 = d.sd0;
  d.l10 = d.rho * d.sd1;
  d.l11 = d.sd1 * std::sqrt(1 - d.rho * d.rho);
  const double dropout = s("dropout");
  const double k = static_cast<double>(d.times.size());
  d.hazard = dropout > 0 ? 1 - std::pow(1 - dropout, 1 / (k - 1)) : 0;
  set_origin(d, std::clamp(0.0, d.times.front(), d.times.back()));
  d.t1 = d.slopes ? d.times.front() : d.t0;  // a random intercept alone is the same at every time
  if (d.slopes) {
    // a new origin where the estimates at t0 have lost four digits or more, from whichever of t0 and the time
    // c = −ρσ0/σ1, where the lines cross with the least spread, it lies nearer to: from the other, the sum of a time
    // and a long shift would lose its place (visits at 10²⁹ and 10³⁰, and lines that cross at 0 within the residual
    // SD, put t0 + shift 10¹³ time units off)
    double share;
    const double shift = steady_shift(d, share);
    if (std::isfinite(shift) && share < 1e-4) {
      double t = d.t0 + shift * d.tau;
      const double c = d.sd1 > 0 ? -d.rho * d.sd0 / d.sd1 : t;
      if (std::fabs(t - c) < std::fabs(t - d.t0)) {
        LongDesign from_c = d;
        set_origin(from_c, c);
        t = c + steady_shift(from_c, share) * from_c.tau;
      }
      set_origin(d, t);
    }
  }
  return d;
}

/// Λ (lower triangular, G = ΛΛ' = MM'/σe² in the analysis units): the starting value of the fit, and G at the
/// assumed variances.  M is triangular when the time origin is 0; otherwise |M| = τ l00 l11 / unit² is taken from
/// the factors of D, not from a difference.
void long_lambda(const LongDesign& d, double* theta) {
  if (d.m01 == 0) {
    theta[0] = d.m00 / d.se1, theta[1] = d.m10 / d.se1, theta[2] = d.m11 / d.se1;
    return;
  }
  const double a = std::hypot(d.m00, d.m01), det = d.tau * d.l00 / d.unit * (d.l11 / d.unit);
  theta[0] = a / d.se1;
  theta[1] = (d.m10 * d.m00 + d.m11 * d.m01) / a / d.se1;
  theta[2] = det / a / d.se1;
}

/// Rows of X (intercept, group, t', group × t') and Z (1, t'') of a subject of group g (0 is group 1) at every
/// visit.
void long_rows(const LongDesign& d, int g, std::vector<double>& x, std::vector<double>& z) {
  const std::size_t K = d.times.size();
  const int q = d.q();
  const double ind = g == 0 ? 1 : 0;
  x.resize(4 * K);
  z.resize(static_cast<std::size_t>(q) * K);
  for (std::size_t j = 0; j < K; ++j) {
    const double t = d.time(j);
    double* xr = &x[4 * j];
    xr[0] = 1, xr[1] = ind, xr[2] = t, xr[3] = ind * t;
    z[static_cast<std::size_t>(q) * j] = 1;
    if (q == 2) z[2 * j + 1] = d.z_time(j);
  }
}

/// Designs of the patterns: index = group × K + (visits − 1); group 0 is group 1 (indicator 1).  With each, X's
/// rows within a subject (the rows of (I − P)X), which depend only on the pattern: formed once, they leave a
/// subject's sums little more work than its errors.
struct LongPatterns {
  lmm::Data data;
  std::vector<std::vector<double>> within;
};

LongPatterns long_patterns(const LongDesign& d) {
  const int K = static_cast<int>(d.times.size()), q = d.q();
  LongPatterns lp;
  lmm::Data& data = lp.data;
  data.p = 4;
  data.q = q;
  data.two_rows = q == 2 && K == 2;
  data.patterns.resize(static_cast<std::size_t>(2 * K));
  lp.within.resize(data.patterns.size());
  // with random slopes and the visits from t1 ≠ t0, (1, t') = (1, t'') S, and X = (1, group, t', group × t') =
  // Z S C exactly; G at the assumed variances decides whether a subject's nearly coinciding visits tell anything
  // apart.  Otherwise Z = (1, t'), and the design comes from X and Z as it always has
  const double s[4] = {1, (d.t1 - d.t0) / d.tau, 0, 1};
  double theta[3], g[4];
  long_lambda(d, theta);
  g[0] = theta[0] * theta[0], g[1] = g[2] = theta[0] * theta[1], g[3] = theta[1] * theta[1] + theta[2] * theta[2];
  std::vector<double> x, z;
  for (int grp = 0; grp < 2; ++grp) {
    long_rows(d, grp, x, z);
    const double ind = grp == 0 ? 1 : 0, c[8] = {1, ind, 0, 0, 0, 0, 1, ind};
    for (int k = 1; k <= K; ++k) {
      const auto i = static_cast<std::size_t>(grp * K + k - 1);
      lmm::Pattern& pt = data.patterns[i];
      lp.within[i].resize(static_cast<std::size_t>(4 * k));
      if (q == 2 && d.t1 != d.t0) {
        pt.set_spanned(k, 4, q, z.data(), c, s, g);  // nothing of X within a subject
      } else {
        pt.set_design(k, 4, q, x.data(), z.data());
        pt.within_rows(4, q, x.data(), z.data(), lp.within[i].data());
      }
    }
  }
  return lp;
}

/// Standard error of the slope difference by generalized least squares at the assumed variances.
double long_se(const LongDesign& d, double n1, double n2) {
  lmm::Data data = long_patterns(d).data;
  const std::vector<double> prob = visit_probs(d);
  const std::size_t K = d.times.size();
  for (std::size_t k = 0; k < K; ++k) {
    data.patterns[k].count = n1 * prob[k];
    data.patterns[K + k].count = n2 * prob[k];
  }
  double theta[3], g[4], cov[lmm::PMAX * lmm::PMAX];
  long_lambda(d, theta);
  if (d.slopes) {
    g[0] = theta[0] * theta[0];
    g[1] = g[2] = theta[0] * theta[1];
    g[3] = theta[1] * theta[1] + theta[2] * theta[2];
  } else {
    g[0] = theta[0] * theta[0];
  }
  // variance components more than 10¹⁵⁴ times apart square to infinity, whatever the units (|G| from the factor:
  // from G's entries it is the small difference of large numbers when the visits are far from time 0)
  const double det_g = d.slopes ? theta[0] * theta[2] * (theta[0] * theta[2]) : g[0];
  const char* extreme = "the answer cannot be computed for these inputs; check for an extreme value";
  if (!std::isfinite(det_g) || !std::isfinite(g[0]) || (d.slopes && !std::isfinite(g[3]))) throw CalcError(extreme);
  if (!lmm::beta_covariance(data, d.se1 * d.se1, g, cov, det_g)) {
    // the slopes can be estimated as soon as some subjects are measured at two times
    if ((n1 + n2) * (1 - prob[0]) < 1)
      throw CalcError("the slopes cannot be estimated with this design (too few subjects measured twice)");
    throw CalcError(extreme);
  }
  return std::sqrt(cov[3 * 4 + 3]) * (d.unit / d.tau);
}

/// Whether a simulated study's subjects leave the model without an estimate, whatever the arithmetic: a group none of
/// whose subjects was measured twice tells nothing of its slope, and with random slopes and no subject measured three
/// times (two visits, or heavy dropout), the subjects seen twice, less one for each group's means, must be two or
/// more, or the covariance of the two measurements is not determined (with one subject left over, the REML criterion
/// falls without bound as it turns singular).  With heavy dropout and few subjects this is common, and it is the
/// design's doing, not an extreme input's: such a study counts as not rejecting.
bool long_unidentified(const lmm::Data& data, int K) {
  double twice[2] = {0, 0}, thrice = 0;
  for (int g = 0; g < 2; ++g)
    for (int k = 2; k <= K; ++k) {
      const double count = data.patterns[static_cast<std::size_t>(g * K + k - 1)].count;
      twice[g] += count;
      if (k > 2) thrice += count;
    }
  if (twice[0] == 0 || twice[1] == 0) return true;
  return data.q == 2 && thrice == 0 && (twice[0] - 1) + (twice[1] - 1) < 2;
}

/// One simulated study.  The analysis estimates the group means and slopes freely, so adding them to the
/// outcome changes the estimate of δ by exactly δ and nothing else: the subjects are simulated without them
/// (in the analysis units, only their random effects and errors), and δ is added to the estimate, which keeps
/// a δ of 10¹⁶⁰ out of the sums of squares.
sim::Outcome long_study(const LongDesign& d, const LongPatterns& designs, long n1, long n2, double delta, long rep,
                        std::uint64_t seed, int dfm, double alpha, Alt a) {
  lmm::Data data = designs.data;
  const int K = static_cast<int>(d.times.size()), q = d.q();
  std::vector<double> x, z, e(static_cast<std::size_t>(K));
  for (int g = 0; g < 2; ++g) {
    long_rows(d, g, x, z);
    const long n = g == 0 ? n1 : n2;
    for (long u = 0; u < n; ++u) {
      Rng rng(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g + 1),
              static_cast<std::uint64_t>(u));
      const double z0 = rng.normal(), z1 = d.slopes ? rng.normal() : 0;
      const double b[2] = {d.m00 * z0 + d.m01 * z1, d.m10 * z0 + d.m11 * z1};  // the intercept at t0, the slope per τ
      int k = 1;
      while (k < K && rng.uniform() >= d.hazard) ++k;
      const auto i = static_cast<std::size_t>(g * K + k - 1);
      const double* rx = designs.within[i].data();
      lmm::Pattern& pt = data.patterns[i];
      if (q == 2) {
        for (int j = 0; j < k; ++j) e[static_cast<std::size_t>(j)] = d.se1 * rng.normal();
        pt.add_effects<4, 2>(4, rx, z.data(), b, e.data());
        continue;
      }
      // a random intercept alone (Z = 1): what is left of X within a subject is (0, 0, t' − t̄', group × (t' − t̄')),
      // so the subject's sums need only the errors' mean, their deviations from it and the time's column; these are
      // Pattern::add_effects' sums to the bit, with one pass over the errors fewer
      double sum = 0, sx = 0, syy = 0;
      for (int j = 0; j < k; ++j) {
        e[static_cast<std::size_t>(j)] = d.se1 * rng.normal();
        sum += e[static_cast<std::size_t>(j)];
      }
      const double mean = pt.zl[0] * sum;  // (Z'Z)⁻¹ Z'e
      if (k > 1) {
        for (int j = 0; j < k; ++j) {
          const double v = e[static_cast<std::size_t>(j)] - mean;
          sx += rx[4 * j + 2] * v;
          syy += v * v;
        }
      }
      const double xy[4] = {0, 0, sx, g == 0 ? sx : 0}, level = b[0] + mean;
      pt.add_sums(4, 1, xy, &level, syy);
    }
  }
  sim::Outcome o;
  if (long_unidentified(data, K)) {
    o.failed = o.design_failure = true;
    return o;
  }
  double theta[3];
  long_lambda(d, theta);
  const lmm::Fit f = lmm::fit(data, theta);
  if (!f.ok) {
    o.failed = true;
    return o;
  }
  const double c[4] = {0, 0, 0, 1};
  const lmm::Wald w = lmm::wald(data, f, c, dfm, static_cast<double>(n1 + n2 - 2));
  const double estimate = w.estimate + delta * (d.tau / d.unit);
  o.reject = p_value(simc::t_stat(estimate, w.se), w.df, a) < alpha;
  o.estimate = estimate * (d.unit / d.tau);
  o.se = w.se * (d.unit / d.tau);
  o.df = w.df;
  o.df_fallback = w.fallback;
  return o;
}

Procedure longitudinal() {
  Procedure p;
  p.id = "sim.lmm.longitudinal";
  p.name = "Two Mean Slopes with a Mixed Model, by Simulation";
  p.category = "Simulation";
  p.subcategory = "Longitudinal Studies";
  p.order = 1;
  p.simulation = true;
  p.summary =
      "Power and sample size for comparing the mean rate of change of two groups measured repeatedly, estimated "
      "by simulating whole studies and analyzing each with a linear mixed model (random intercepts and slopes, "
      "REML, Satterthwaite t-test). Handles dropout between visits and any measurement schedule.";
  p.keywords = {"simulation", "monte carlo", "mixed model", "lmm", "linear mixed model", "random slope",
                "random intercept", "longitudinal", "slope", "rate of change", "growth curve", "repeated measures",
                "dropout", "attrition", "missing data", "satterthwaite", "reml", "lme4", "lmer"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total Sample Size)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: β1 ≠ β2)", "One-sided (H1: β1 > β2)", "One-sided (H1: β1 < β2)"),
      simc::df_param("t with N − 2 df (subjects)"),
      power_param(),
      alpha_param(),
      num("n1", "N1 (Subjects in Group 1)", "50")
          .grp(G_N)
          .sym("N1")
          .whole()
          .at_least(3)
          .sugg({"50", "20 to 100 by 20"})
          .hint("Number of subjects in group 1; group 2 gets N1 × the allocation ratio, rounded up.")
          .show(when(not_solving({"n"}))),
      num("r", "Allocation Ratio N2/N1", "1")
          .grp(G_N)
          .sym("N2/N1")
          .above(0)
          .sugg({"1", "2"})
          .hint("Size of group 2 relative to group 1 (N2 = ratio × N1, rounded up); 1 gives equal groups."),
      num("dropout", "Dropout by the Last Visit", "0")
          .grp(G_N)
          .sym("Dropout")
          .range(0, 1, false, true)
          .sugg({"0", "0.1", "0.2", "0 0.1 0.2"})
          .hint("Proportion of subjects who have left the study by the last visit. Subjects leave at a constant "
                "rate between visits, completely at random, and their earlier measurements are analyzed."),
      num("delta", "δ (Difference in Mean Slopes, β1 − β2)", "0.5")
          .grp(G_EFFECT)
          .sym("δ")
          .sugg({"0.5", "0.25 0.5"})
          .hint("Difference between the groups' mean changes per time unit."),
      select("re", "Random Effects",
             {{"slopes", "Random intercepts and slopes", {}}, {"intercepts", "Random intercepts only", {}}},
             "slopes")
          .grp(simc::G_MODEL)
          .hint("The variation between subjects, used both to simulate the data and in the analysis model."),
      num("sd0", "σ0 (SD of the Subjects' Intercepts)", "1")
          .grp(simc::G_MODEL)
          .sym("σ0")
          .at_least(0)
          .sugg({"1", "0.5 1 2"})
          .hint("Between-subject standard deviation of the true values at time 0."),
      num("sd1", "σ1 (SD of the Subjects' Slopes)", "0.5")
          .grp(simc::G_MODEL)
          .sym("σ1")
          .at_least(0)
          .sugg({"0.5", "0.25 0.5 1"})
          .hint("Between-subject standard deviation of the true slopes within a group.")
          .show(when(is("re", {"slopes"}))),
      num("rho", "ρ (Correlation of Intercepts and Slopes)", "0")
          .grp(simc::G_MODEL)
          .sym("ρ")
          .range(-1, 1, true, true)
          .sugg({"0", "-0.3 0 0.3"})
          .hint("Correlation between a subject's true value at time 0 and true slope; 0 when the starting level "
                "says nothing about the rate of change.")
          .show(when(is("re", {"slopes"}))),
      num("sde", "σe (Residual SD)", "1")
          .grp(simc::G_MODEL)
          .sym("σe")
          .above(0)
          .sugg({"1", "0.5 1 2"})
          .hint("Within-subject standard deviation of the measurements around the subject's line."),
      select("timing", "Spacing of the Measurements",
             {{"equal", "Equally spaced", {}}, {"custom", "Custom times", {}}}, "equal")
          .grp("Measurement Schedule")
          .hint("Plan M measurements equally spaced from time 0 to D, or list the measurement times yourself."),
      num("m", "M (Measurements per Subject)", "5")
          .grp("Measurement Schedule")
          .sym("M")
          .at_least(2)
          .whole()
          .sugg({"5", "3 5 7"})
          .hint("How many measurements are planned for each subject, counting the first (baseline) "
                "measurement.")
          .show(when(is("timing", {"equal"}))),
      num("duration", "D (Time from First to Last Measurement)", "4")
          .grp("Measurement Schedule")
          .sym("D")
          .above(0)
          .sugg({"4", "12"})
          .hint("Length of follow-up from the first to the last measurement, in the time unit of the slopes "
                "(e.g. years if δ is a change per year).")
          .show(when(is("timing", {"equal"}))),
      vec("times", "Measurement Times (e.g. 0 1 2 4 8)", "0 1 2 4 8")
          .grp("Measurement Schedule")
          .sym("Times")
          .hint("The planned times of the measurements, in the time unit of the slopes.")
          .show(when(is("timing", {"custom"}))),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double alpha = s("alpha"), ratio = s("r"), delta = s("delta");
    const LongDesign d = long_design(s);
    const int dfm = df_method(s);
    auto n2_of = [&](double n1) { return std::ceil(ratio * n1 - 1e-9); };
    auto approx = [&](double n1) {
      const double n2 = n2_of(n1), ncp = delta / long_se(d, n1, n2);
      return dfm == 2 ? z_power(ncp, alpha, a) : t_power(ncp, n1 + n2 - 2, alpha, a);
    };
    if (s.solve == "n") {
      if (delta == 0) throw CalcError("δ is zero; power cannot exceed alpha");
      if ((a == Alt::Upper && delta < 0) || (a == Alt::Lower && delta > 0))
        throw CalcError("δ has the sign opposite to the one-sided alternative");
    }
    const double lo = s.solve != "n" ? 3 : solve::first_whole([&](double x) { return n2_of(x) >= 3; }, 3, 1e5,
                                                           "the allocation ratio is too small: group 2 would have "
                                                           "fewer than 3 subjects");
    double n1;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = lo;
      o.hi = 1e5;
      n1 = solve::min_n(approx, s("power"), o, &ctx);
    } else {
      n1 = s("n1");
      if (n2_of(n1) < 3) throw CalcError("group 2 needs at least 3 subjects");
    }
    if (ctx.approximate) {
      simc::put_approximation(s, r, approx(n1));
    } else {
      const long reps = static_cast<long>(s("reps"));
      const auto seed = static_cast<std::uint64_t>(s("seed"));
      const LongPatterns designs = long_patterns(d);
      auto power_at = [&](double m1, long R, double p0, double p1) {
        simc::check_study_size((m1 + n2_of(m1)) * static_cast<double>(d.times.size()), "measurements");
        const long k1 = static_cast<long>(m1), k2 = static_cast<long>(n2_of(m1));
        return sim::simulate(
            R, [&](long rep) { return long_study(d, designs, k1, k2, delta, rep, seed, dfm, alpha, a); }, ctx, p0,
            p1);
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
        approx(n1);  // first: where it cannot be computed, neither can the simulated power
        sm = power_at(n1, reps, 0, 1);
      }
      // studies whose subjects leave the model without an estimate count as not rejecting, as a design that
      // often gives too few subjects measured twice deserves; only failures of the arithmetic in more than half of
      // the studies withhold the power
      const bool numerical = 2 * (sm.failures - sm.design_failures) > sm.reps;
      simc::put_summary(s, r, sm, numerical ? simc::Failures::Numerical : simc::Failures::Design,
                        sm.design_failures == sm.reps ? "each had a group with too few subjects measured twice" : "");
      r.set("approx_power", approx(n1));
      if (std::isfinite(sm.mean_estimate)) r.set("mean_delta", sm.mean_estimate);
      if (sm.df_fallbacks > 0) r.set("df_fallbacks", static_cast<double>(sm.df_fallbacks));
    }
    put_two_group(r, {n1, n2_of(n1)});
    r.set("delta", delta);
    r.set("sd0", d.sd0);
    if (d.slopes) {
      r.set("sd1", d.sd1);
      r.set("rho", d.rho);
    }
    r.set("sde", d.sde);
    r.set("m", static_cast<double>(d.times.size()));
    r.set("times", list_text(d.times));
    r.set("dropout", s("dropout"));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("delta", "δ", "num", "Difference in mean slopes β1 − β2"),
                   col("mean_delta", "Mean δ̂", "num", "Average estimated slope difference over the simulated studies"),
                   col("sd0", "σ0", "num", "SD of the subjects' intercepts"),
                   col("sd1", "σ1", "num", "SD of the subjects' slopes"),
                   col("rho", "ρ", "num", "Correlation of intercepts and slopes"),
                   col("sde", "σe", "num", "Residual SD"),
                   col("times", "Times", "text", "Measurement times"),
                   col("dropout", "Dropout", "num", "Proportion lost by the last visit"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns()) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string model =
        s.is("re", "intercepts") ? "random intercepts" : "random intercepts and slopes";
    const std::string test =
        s.is("df", "z") ? "z-test" : (s.is("df", "fixed") ? "t-test (N − 2 df)" : "Satterthwaite t-test");
    std::string t = "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " (total " +
                    fint(r.num("n")) + ") " + simc::power_phrase(s, r) + " to detect a difference of " +
                    fnum(r.num("delta")) + " between the groups' mean slopes, analyzing each study with a linear " +
                    "mixed model with " + model + " and a " + alt_words(alt(s)) + " " + test +
                    " at α = " + fnum(r.num("alpha")) + ". Subjects are measured at times " + text_of(r, "times");
    t += r.num("dropout") > 0 ? " and " + fpct(r.num("dropout"), 0) + " drop out by the last visit." : ".";
    return t + simc::method_note(r) + fallback_note(r, "N − 2");
  };
  p.definitions = {{"N1, N2", "Number of subjects in groups 1 and 2; N = N1 + N2."},
                   {"δ", "Difference β1 − β2 between the groups' mean slopes (change per time unit)."},
                   {"σ0, σ1, ρ", "SDs of the subjects' true intercepts and slopes, and their correlation."},
                   {"σe", "Residual standard deviation of a measurement around the subject's line."},
                   {"Dropout", "Proportion of subjects who have left by the last visit (monotone, completely at "
                               "random)."},
                   {"Power", "Share of the simulated studies whose test of δ = 0 was significant, with its 95% "
                             "interval (Wilson) and Monte Carlo standard error."}};
  p.notes = {"Every simulated study is generated from the model y = β0 + β·group + (β_t + δ·group)·t + b0 + b1·t + "
             "e and analyzed with the same model by REML; the effect δ is tested with a t statistic.",
             "With complete, balanced data the generalized least squares estimate of δ is the difference of the "
             "subjects' mean least-squares slopes, so with many subjects the simulated power agrees with Two Mean "
             "Slopes in a Longitudinal Design. With few subjects the REML estimate of the slope variance is often "
             "zero, which makes the test conservative (type I error about 4% at 20 subjects per group, as with "
             "lme4), so the simulated power is a few points lower.",
             "Results depend only on the seed: each subject's data come from their own random stream, so a "
             "larger study contains the subjects of a smaller one, which keeps the sample-size search stable."};
  p.references = {"Laird, N. M. & Ware, J. H. (1982). Random-effects models for longitudinal data. Biometrics 38, "
                  "963–974.",
                  "Satterthwaite, F. E. (1946). An approximate distribution of estimates of variance components. "
                  "Biometrics Bulletin 2, 110–114.",
                  "Kuznetsova, A., Brockhoff, P. B. & Christensen, R. H. B. (2017). lmerTest package: tests in "
                  "linear mixed effects models. Journal of Statistical Software 82(13).",
                  simc::REF_SIM, simc::REF_ARNOLD};
  p.curve_param = "n1";
  p.curve_label = "N1 (Subjects in Group 1)";
  p.curve_min = 3;
  p.example_solve = "n";
  p.example = {{"power", "0.8"}, {"delta", "0.5"}, {"sd1", "0.5"}, {"dropout", "0.1"}, {"reps", "200"}};
  return p;
}
SSC_REGISTER(longitudinal);

// =====================================================================================
// Cluster-randomized: difference between two means
// =====================================================================================
struct ClusterDesign {
  simc::ClusterSizes sizes;
  double sb = 0, se = 1;  // SDs of the cluster effects and residuals in units of the total SD σ
  double unit = 1;        // σ: the outcome's unit in the simulated studies and in the approximation
  double cv() const { return sizes.cv; }
};

ClusterDesign cluster_design(const Scenario& s) {
  const double icc = s("icc");
  return {simc::ClusterSizes(s("m"), s("cv")), std::sqrt(icc), std::sqrt(1 - icc), s("sigma")};
}

/// Mean of a cluster's weight m/(σe² + m σb²) (in units of 1/σ²) and of its size, over the size distribution.
std::pair<double, double> expected_weight(const ClusterDesign& d) {
  double w = 0, m = 0;
  const std::vector<double> grid = d.sizes.grid();
  for (double size : grid) {
    w += size / (d.se * d.se + size * d.sb * d.sb);
    m += size;
  }
  return {w / static_cast<double>(grid.size()), m / static_cast<double>(grid.size())};
}

/// An arm's sizes in the order they were first drawn, as indices into its sizes: an open-addressing table, 4 bytes a
/// slot and at most half full (a map took 64 bytes a size, and sizes that vary widely give almost every
/// cluster one of its own).
class SizeIndex {
 public:
  explicit SizeIndex(lmm::Arm& arm) : arm_(arm), slots_(64, EMPTY) {}
  /// The arm's sums for clusters of size m, added (empty) for the first one.
  lmm::Arm::Size& operator[](long m) {
    std::size_t i = slot(m);
    if (slots_[i] == EMPTY) {
      if (2 * (arm_.sizes.size() + 1) > slots_.size()) {
        grow();
        i = slot(m);
      }
      slots_[i] = static_cast<std::uint32_t>(arm_.sizes.size());
      arm_.sizes.push_back({static_cast<double>(m)});
    }
    return arm_.sizes[slots_[i]];
  }

 private:
  static constexpr std::uint32_t EMPTY = 0xffffffff;
  lmm::Arm& arm_;
  std::vector<std::uint32_t> slots_;
  /// The slot that holds size m, or the empty one where it goes (linear probing from a Fibonacci hash).
  std::size_t slot(long m) const {
    const std::size_t mask = slots_.size() - 1;
    std::size_t i = static_cast<std::size_t>(static_cast<std::uint64_t>(m) * 0x9e3779b97f4a7c15ULL >> 32) & mask;
    while (slots_[i] != EMPTY && arm_.sizes[slots_[i]].n != static_cast<double>(m)) i = (i + 1) & mask;
    return i;
  }
  void grow() {
    std::vector<std::uint32_t>(2 * slots_.size(), EMPTY).swap(slots_);
    for (std::size_t k = 0; k < arm_.sizes.size(); ++k)
      slots_[slot(static_cast<long>(arm_.sizes[k].n))] = static_cast<std::uint32_t>(k);
  }
};

/// One simulated study.  The analysis estimates both means freely, so adding them to the outcome changes the
/// estimated difference by exactly μ1 − μ2 and nothing else: the clusters are simulated around 0, in units of
/// σ, and the difference is added to the estimate (sums of squares of outcomes far from 0, or of 10¹⁶⁰ σ, would
/// lose the digits the variance estimates need, or overflow).  A study of a million clusters takes about a second,
/// so the token is checked within it: as the clusters are drawn, and as each evaluation of the REML criterion begins,
/// in the fit and in Satterthwaite's degrees of freedom (some 20 ms apart with a million sizes).
sim::Outcome cluster_study(const ClusterDesign& d, long k1, long k2, double diff, long rep, std::uint64_t seed,
                           int dfm, double alpha, Alt a, const RunContext& ctx) {
  lmm::Data data;
  data.p = 2;
  data.q = 1;
  data.arms.resize(2);
  data.poll = [](const void* c) { static_cast<const RunContext*>(c)->check(); };
  data.poll_arg = &ctx;
  for (int g = 0; g < 2; ++g) {
    lmm::Arm& arm = data.arms[static_cast<std::size_t>(g)];
    arm.x[0] = 1, arm.x[1] = g == 0 ? 1 : 0;  // X = (1, group), the same on every row of a cluster
    SizeIndex index(arm);
    for (long c = 0; c < (g == 0 ? k1 : k2); ++c) {
      if (c % 4096 == 0) ctx.check();
      Rng rng(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g + 1),
              static_cast<std::uint64_t>(c));
      const long m = d.sizes.draw(rng);
      const double dm = static_cast<double>(m);
      // the cluster's mean and within-cluster sum of squares (normal theory)
      const double mean = d.sb * rng.normal() + d.se / std::sqrt(dm) * rng.normal();
      const double ss = m > 1 ? d.se * d.se * rng.chisq(dm - 1) : 0;
      lmm::Arm::Size& sz = index[m];
      sz.count += 1;
      sz.b += mean;
      sz.bb += mean * mean;
      sz.yy += ss;
    }
  }
  sim::Outcome o;
  const double theta = d.sb / d.se;
  const lmm::Fit f = lmm::fit(data, &theta);
  if (!f.ok) {
    o.failed = true;
    return o;
  }
  const double c[2] = {0, 1};
  const lmm::Wald w = lmm::wald(data, f, c, dfm, static_cast<double>(k1 + k2 - 2));
  const double estimate = w.estimate + diff / d.unit;
  o.reject = p_value(simc::t_stat(estimate, w.se), w.df, a) < alpha;
  o.estimate = estimate * d.unit;
  o.se = w.se * d.unit;
  o.df = w.df;
  o.df_fallback = w.fallback;
  return o;
}

Procedure cluster() {
  Procedure p;
  p.id = "sim.lmm.cluster";
  p.name = "Two Means in a Cluster-Randomized Design with a Mixed Model, by Simulation";
  p.category = "Simulation";
  p.subcategory = "Cluster-Randomized Trials";
  p.order = 1;
  p.simulation = true;
  p.summary =
      "Number of clusters per arm (or power) for comparing two means when whole clusters are randomized, "
      "estimated by simulating studies and analyzing each with a linear mixed model (random cluster intercept, "
      "REML, Satterthwaite t-test). Cluster sizes may vary, which a design effect formula only approximates.";
  p.keywords = {"simulation", "monte carlo", "mixed model", "lmm", "cluster", "cluster randomized",
                "group randomized", "icc", "intraclass", "unequal cluster sizes", "variable cluster size",
                "design effect", "community trial", "satterthwaite", "reml", "lme4", "lmer"};
  p.solve_for = {{"k", "Number of Clusters", "k1", "Clusters in Group 1"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: μ1 ≠ μ2)", "One-sided (H1: μ1 > μ2)", "One-sided (H1: μ1 < μ2)"),
      simc::df_param("t with k1 + k2 − 2 df (clusters)"),
      power_param(),
      alpha_param(),
      num("k1", "k1 (Clusters in Group 1)", "10")
          .grp(G_N)
          .sym("k1")
          .whole()
          .at_least(2)
          .sugg({"10", "5 to 30 by 5"})
          .hint("Number of clusters (such as clinics, schools or villages) randomized to group 1.")
          .show(when(not_solving({"k"}))),
      num("kr", "Cluster Ratio k2/k1", "1")
          .grp(G_N)
          .sym("k2/k1")
          .above(0)
          .sugg({"1", "2"})
          .hint("Clusters in group 2 per cluster in group 1; k2 = ratio × k1, rounded up. Use 1 for equal "
                "numbers of clusters."),
      simc::cluster_size_param(),
      simc::cluster_cv_param(),
      num("icc", "ρ (Intracluster Correlation)", "0.05")
          .grp(G_EFFECT)
          .sym("ICC")
          .range(0, 1, false, true)
          .sugg({"0.01", "0.05", "0.01 0.05 0.1"})
          .hint("The share of the outcome's variance that lies between clusters, from 0 to 1; values from 0.01 to "
                "0.1 are common. Take it from earlier studies of similar clusters."),
      num("mu1", "μ1 (Mean of Group 1)", "0.5")
          .grp(G_EFFECT)
          .sym("μ1")
          .sugg({"0.5", "0.25 0.5"})
          .hint("Mean outcome in group 1 (e.g. the intervention) that you expect under the alternative."),
      num("mu2", "μ2 (Mean of Group 2)", "0")
          .grp(G_EFFECT)
          .sym("μ2")
          .sugg({"0"})
          .hint("Mean outcome in group 2 (e.g. the control); the difference μ1 − μ2 is the effect to detect."),
      num("sigma", "σ (Total SD of a Subject's Outcome)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .sugg({"1"})
          .hint("Standard deviation of one subject's outcome within a group, including the variation between "
                "clusters; the ICC sets how much of it lies between clusters."),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double alpha = s("alpha"), kr = s("kr"), mu1 = s("mu1"), mu2 = s("mu2"), sigma = s("sigma");
    const double diff = mu1 - mu2;
    const ClusterDesign d = cluster_design(s);
    const auto ew = expected_weight(d);
    const int dfm = df_method(s);
    auto k2_of = [&](double k1) { return std::ceil(kr * k1 - 1e-9); };
    auto approx = [&](double k1) {
      const double k2 = k2_of(k1), se = std::sqrt((1 / k1 + 1 / k2) / ew.first), ncp = diff / sigma / se;
      return dfm == 2 ? z_power(ncp, alpha, a) : t_power(ncp, k1 + k2 - 2, alpha, a);
    };
    if (s.solve == "k") {
      if (diff == 0) throw CalcError("μ1 equals μ2");
      if ((a == Alt::Upper && diff < 0) || (a == Alt::Lower && diff > 0))
        throw CalcError("μ1 − μ2 has the sign opposite to the one-sided alternative");
    }
    const double lo = s.solve != "k" ? 2 : solve::first_whole([&](double x) { return k2_of(x) >= 2; }, 2, 1e5,
                                                           "the cluster ratio is too small: group 2 would have "
                                                           "fewer than 2 clusters");
    double k1;
    if (s.solve == "k") {
      solve::NSearch o;
      o.lo = lo;
      o.hi = 1e5;
      k1 = solve::min_n(approx, s("power"), o, &ctx);
    } else {
      k1 = s("k1");
      if (k2_of(k1) < 2) throw CalcError("group 2 needs at least 2 clusters");
    }
    if (ctx.approximate) {
      simc::put_approximation(s, r, approx(k1));
    } else {
      const long reps = static_cast<long>(s("reps"));
      const auto seed = static_cast<std::uint64_t>(s("seed"));
      auto power_at = [&](double c1, long R, double p0, double p1) {
        simc::check_study_size(c1 + k2_of(c1), "clusters");
        const long n1 = static_cast<long>(c1), n2 = static_cast<long>(k2_of(c1));
        return sim::simulate(
            R, [&](long rep) { return cluster_study(d, n1, n2, diff, rep, seed, dfm, alpha, a, ctx); }, ctx, p0,
            p1);
      };
      sim::Summary sm;
      if (s.solve == "k") {
        sim::SizeSearch o;
        o.lo = lo;
        o.hi = 1e5;
        o.guess = k1;
        o.reps = reps;
        o.explore_reps = simc::explore_reps(reps);
        const sim::SizeResult res = sim::search_size(power_at, s("power"), o, ctx);
        k1 = res.n;
        sm = res.s;
      } else {
        sm = power_at(k1, reps, 0, 1);
      }
      simc::put_summary(s, r, sm, simc::Failures::Numerical);
      r.set("approx_power", approx(k1));
      if (std::isfinite(sm.mean_estimate)) r.set("mean_delta", sm.mean_estimate);
      if (sm.df_fallbacks > 0) r.set("df_fallbacks", static_cast<double>(sm.df_fallbacks));
    }
    const double k2 = k2_of(k1);
    // the average size as entered (rounding the drawn sizes moves their mean only in the third decimal)
    const double mbar = d.cv() > 0 ? s("m") : ew.second;
    r.set("k1", k1);
    r.set("k2", k2);
    r.set("m", mbar);
    r.set("cv", d.cv());
    r.set("n1", k1 * mbar);
    r.set("n2", k2 * mbar);
    r.set("n", (k1 + k2) * mbar);
    r.set("mu1", mu1);
    r.set("mu2", mu2);
    r.set("delta", diff);
    r.set("sigma", sigma);
    r.set("icc", s("icc"));
    r.set("de", ew.second / ew.first);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    for (auto x : {col("k1", "k1", "int", "Clusters in group 1"), col("k2", "k2", "int", "Clusters in group 2"),
                   col("m", "m̄", "num", "Average cluster size"), col("cv", "CV", "num", "CV of the cluster sizes"),
                   col("n", "N", "num", "Expected number of subjects"), col("mu1", "μ1", "num"),
                   col("mu2", "μ2", "num"),
                   col("mean_delta", "Mean δ̂", "num", "Average estimated difference over the simulated studies"),
                   col("sigma", "σ", "num"), col("icc", "ICC", "num"),
                   col("de", "Design Effect", "num", "Variance inflation relative to individual randomization"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns()) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string test =
        s.is("df", "z") ? "z-test" : (s.is("df", "fixed") ? "t-test (k1 + k2 − 2 df)" : "Satterthwaite t-test");
    const std::string sizes = r.num("cv") > 0 ? "clusters of " + fnum(r.num("m")) + " subjects on average (CV " +
                                                     fnum(r.num("cv")) + ")"
                                               : "clusters of " + fnum(r.num("m")) + " subjects";
    return fint(r.num("k1")) + " and " + fint(r.num("k2")) + " " + sizes + " " + simc::power_phrase(s, r) +
           " to detect a difference of " + fnum(r.num("delta")) + " between the means (σ = " +
           fnum(r.num("sigma")) + ", ICC = " + fnum(r.num("icc")) + "), analyzing each study with a linear mixed " +
           "model with a random cluster intercept and a " + alt_words(alt(s)) + " " + test + " at α = " +
           fnum(r.num("alpha")) + "." + simc::method_note(r) + fallback_note(r, "k1 + k2 − 2");
  };
  p.definitions = {{"k1, k2", "Number of clusters randomized to each group."},
                   {"m̄, CV", "Average cluster size and coefficient of variation of the sizes."},
                   {"ICC", "Intracluster correlation coefficient: the share of the outcome's variance that lies "
                           "between clusters."},
                   {"Design Effect", "Variance of the estimated difference relative to individual randomization "
                                     "of the same number of subjects."},
                   {"Power", "Share of the simulated studies whose test of μ1 = μ2 was significant, with its 95% "
                             "interval (Wilson) and Monte Carlo standard error."}};
  p.notes = {"Each simulated cluster has a size drawn from the size distribution, a normal cluster effect and "
             "normal subject-level outcomes; the cluster means and within-cluster sums of squares are drawn "
             "directly from their exact distributions.",
             "With equal cluster sizes the analysis is equivalent to a t-test on the cluster means unless the "
             "estimated ICC is zero, so with many clusters the simulated power agrees with Two Means in a "
             "Cluster-Randomized Design. With few clusters a zero estimate is common and makes the test "
             "conservative (with 4 clusters per group and K − 2 degrees of freedom, a type I error of about 3% and "
             "power 0.60 against 0.62).",
             "The quick approximation uses the expected inverse-variance weight of a cluster over the size "
             "distribution. It is close to the simulated power when the sizes vary moderately, and above it when "
             "they vary widely and clusters are few (0.43 against 0.39 with a size CV of 1.2 and 6 clusters per "
             "group)."};
  p.references = {"Donner, A. & Klar, N. (2000). Design and Analysis of Cluster Randomization Trials in Health "
                  "Research. Arnold.",
                  "Eldridge, S. M., Ashby, D. & Kerry, S. (2006). Sample size for cluster randomized trials: effect "
                  "of coefficient of variation of cluster size and analysis method. International Journal of "
                  "Epidemiology 35, 1292–1300.",
                  "Kuznetsova, A., Brockhoff, P. B. & Christensen, R. H. B. (2017). lmerTest package: tests in "
                  "linear mixed effects models. Journal of Statistical Software 82(13).",
                  simc::REF_SIM};
  p.curve_param = "k1";
  p.curve_label = "k1 (Clusters in Group 1)";
  p.example_solve = "k";
  p.example = {{"power", "0.8"}, {"m", "20"}, {"cv", "0.6"}, {"icc", "0.05"}, {"mu1", "0.4"}, {"reps", "500"}};
  return p;
}
SSC_REGISTER(cluster);

}  // namespace
}  // namespace ssc::procs
