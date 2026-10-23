// Simulation-based power for a parallel-group trial with repeated measures
// and dropout, analyzed with a mixed model for repeated measures (sim.mmrm):
// the usual primary analysis of a continuous outcome in confirmatory trials
// (change from baseline at the last visit, with a mean for every group and
// visit and an unstructured covariance matrix).
//
// Patients are measured at T visits with a compound-symmetric or AR(1)
// correlation; each group loses patients at a constant rate between visits,
// completely at random, and patients who leave keep their earlier
// measurements.  Each simulated trial is fitted by REML (sim_mmrm.hpp) and the
// difference at the last visit is tested with a t statistic; a t-test on the
// completers is reported alongside.  The quick approximation is the expected
// generalized least squares variance under the true covariance and the
// expected dropout pattern (Lu, Luo & Chen 2008).
#include "procedures/sim_mmrm.hpp"

#include <algorithm>
#include <cstdint>

#include "procedures/sim_common.hpp"

namespace ssc::procs {

namespace mmrm {
namespace {

/// Cholesky factor of the leading n × n block of a (row-major, stride TMAX); false if not positive definite.
bool cholesky(const double* a, int n, double* l) {
  for (int j = 0; j < n; ++j) {
    double s = a[j * TMAX + j];
    for (int k = 0; k < j; ++k) s -= l[j * TMAX + k] * l[j * TMAX + k];
    if (!(s > 0)) return false;
    l[j * TMAX + j] = std::sqrt(s);
    for (int i = j + 1; i < n; ++i) {
      double t = a[i * TMAX + j];
      for (int k = 0; k < j; ++k) t -= l[i * TMAX + k] * l[j * TMAX + k];
      l[i * TMAX + j] = t / l[j * TMAX + j];
    }
  }
  return true;
}

/// Solve (L L') x = b in place, L from cholesky().
void chol_solve(const double* l, int n, double* b) {
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < i; ++k) b[i] -= l[i * TMAX + k] * b[k];
    b[i] /= l[i * TMAX + i];
  }
  for (int i = n - 1; i >= 0; --i) {
    for (int k = i + 1; k < n; ++k) b[i] -= l[k * TMAX + i] * b[k];
    b[i] /= l[i * TMAX + i];
  }
}

}  // namespace

void Stats::add(int g, int k, const double* y) {
  Pattern& p = pattern[g][k - 1];
  p.count += 1;
  for (int a = 0; a < k; ++a) {
    p.sum[a] += y[a];
    for (int b = 0; b <= a; ++b) p.cross[a * TMAX + b] += y[a] * y[b];
  }
}

Fit fit(const Stats& s, bool reml, const double* contrast) {
  Fit f;
  const int T = s.visits;
  auto fail = [&](Failure why) {
    f.failure = why;
    return f;
  };
  if (T < 1 || T > TMAX) return fail(Failure::Empty);
  // --- the regression of each visit on the group and the earlier visits, over the patients observed at it
  // z = (group 1, group 2, y1, …, yt); C = Σ z z'
  constexpr int Z = TMAX + 2;
  std::array<double, TMAX * TMAX> sigma{};
  std::array<std::array<double, TMAX>, TMAX> beta{};  // beta[t][u]: coefficient of visit u in visit t's regression
  std::array<double, TMAX> resid{};
  for (int t = 0; t < T; ++t) {
    std::array<double, Z * Z> c{};
    double n = 0, ng[2] = {0, 0};
    for (int g = 0; g < 2; ++g)
      for (int k = t + 1; k <= T; ++k) {
        const Stats::Pattern& p = s.pattern[g][k - 1];
        if (p.count == 0) continue;
        n += p.count;
        ng[g] += p.count;
        c[g * Z + g] += p.count;
        for (int a = 0; a <= t; ++a) {
          c[(2 + a) * Z + g] += p.sum[a];
          for (int b = 0; b <= a; ++b) c[(2 + a) * Z + 2 + b] += p.cross[a * TMAX + b];
        }
      }
    const int q = t + 2;  // predictors: two group means and t earlier visits
    if (ng[0] == 0 || ng[1] == 0) return fail(Failure::Empty);
    if (n - q < 1) return fail(Failure::Few);
    std::array<double, Z * Z> a{};
    std::array<double, Z> rhs{};
    for (int i = 0; i < q; ++i) {
      for (int j = 0; j < q; ++j) a[i * Z + j] = i >= j ? c[i * Z + j] : c[j * Z + i];
      rhs[i] = c[(2 + t) * Z + i];  // i < 2 + t
    }
    // solve by Cholesky (stride Z)
    std::array<double, Z * Z> l{};
    for (int j = 0; j < q; ++j) {
      double v = a[j * Z + j];
      for (int k = 0; k < j; ++k) v -= l[j * Z + k] * l[j * Z + k];
      if (!(v > 1e-12 * std::fmax(1.0, a[j * Z + j]))) return fail(Failure::Arithmetic);
      l[j * Z + j] = std::sqrt(v);
      for (int i = j + 1; i < q; ++i) {
        double w = a[i * Z + j];
        for (int k = 0; k < j; ++k) w -= l[i * Z + k] * l[j * Z + k];
        l[i * Z + j] = w / l[j * Z + j];
      }
    }
    std::array<double, Z> coef = rhs;
    for (int i = 0; i < q; ++i) {
      for (int k = 0; k < i; ++k) coef[i] -= l[i * Z + k] * coef[k];
      coef[i] /= l[i * Z + i];
    }
    for (int i = q - 1; i >= 0; --i) {
      for (int k = i + 1; k < q; ++k) coef[i] -= l[k * Z + i] * coef[k];
      coef[i] /= l[i * Z + i];
    }
    double rss = c[(2 + t) * Z + 2 + t];
    for (int i = 0; i < q; ++i) rss -= coef[i] * rhs[i];
    const double denom = reml ? n - 2 : n;
    if (!(denom > 0)) return fail(Failure::Few);
    resid[t] = std::fmax(rss, 0.0) / denom;
    for (int u = 0; u < t; ++u) beta[t][u] = coef[2 + u];
  }
  // --- the covariance matrix from the regressions: Σ_ts = Σ_u β_tu Σ_us (s < t), Σ_tt = σ²_t + Σ_u β_tu Σ_ut
  for (int t = 0; t < T; ++t) {
    for (int sv = 0; sv < t; ++sv) {
      double v = 0;
      for (int u = 0; u < t; ++u) v += beta[t][u] * sigma[u * TMAX + sv];
      sigma[t * TMAX + sv] = sigma[sv * TMAX + t] = v;
    }
    double v = resid[t];
    for (int u = 0; u < t; ++u) v += beta[t][u] * sigma[u * TMAX + t];
    sigma[t * TMAX + t] = v;
  }
  // --- generalized least squares for each group's means: the leading blocks of Σ⁻¹ come from M = L⁻¹
  std::array<double, TMAX * TMAX> lch{}, m{};
  if (!cholesky(sigma.data(), T, lch.data())) return fail(Failure::Arithmetic);
  for (int j = 0; j < T; ++j) {  // M = L⁻¹, lower triangular
    m[j * TMAX + j] = 1 / lch[j * TMAX + j];
    for (int i = j + 1; i < T; ++i) {
      double v = 0;
      for (int k = j; k < i; ++k) v -= lch[i * TMAX + k] * m[k * TMAX + j];
      m[i * TMAX + j] = v / lch[i * TMAX + i];
    }
  }
  // the combination of the visits' means compared (the last visit: c = e_T), and its variance c' info⁻¹ c by group
  std::array<double, TMAX> c{};
  if (contrast)
    std::copy(contrast, contrast + T, c.begin());
  else
    c[T - 1] = 1;
  double var_diff[2] = {0, 0};
  for (int g = 0; g < 2; ++g) {
    std::array<double, TMAX * TMAX> info{};
    std::array<double, TMAX> rhs{};
    for (int k = 1; k <= T; ++k) {
      const Stats::Pattern& p = s.pattern[g][k - 1];
      if (p.count == 0) continue;
      if (k == T) f.observed_last += p.count;
      // Σ_kk⁻¹[a][b] = Σ_{c<k} M[c][a] M[c][b]
      for (int a = 0; a < k; ++a)
        for (int b = 0; b <= a; ++b) {
          double v = 0;
          for (int cc = a; cc < k; ++cc) v += m[cc * TMAX + a] * m[cc * TMAX + b];
          info[a * TMAX + b] += p.count * v;
          if (a != b) info[b * TMAX + a] += p.count * v;
          rhs[a] += v * p.sum[b];
          if (a != b) rhs[b] += v * p.sum[a];
        }
    }
    std::array<double, TMAX * TMAX> li{};
    if (!cholesky(info.data(), T, li.data())) return fail(Failure::Arithmetic);
    std::array<double, TMAX> mu = rhs;
    chol_solve(li.data(), T, mu.data());
    for (int t = 0; t < T; ++t) f.mu[g][t] = mu[t];
    std::array<double, TMAX> e = c;
    chol_solve(li.data(), T, e.data());
    for (int t = 0; t < T; ++t) var_diff[g] += c[t] * e[t];
  }
  for (int t = 0; t < T; ++t) f.delta += c[t] * (f.mu[0][t] - f.mu[1][t]);
  f.se = std::sqrt(var_diff[0] + var_diff[1]);
  f.sigma = sigma;
  f.ok = f.se > 0;
  if (!f.ok) f.failure = Failure::Arithmetic;
  return f;
}

}  // namespace mmrm

namespace {

using namespace blocks;
using sim::Rng;
using mmrm::TMAX;

// The analyses do not depend on the outcome's unit or on the means of the groups at the visits (the MMRM
// estimates every one of them), so the trials are simulated with mean 0 in units of σ, and δ/σ is added to the
// estimated difference: a δ of 10¹⁰ σ, or a σ of 10¹⁶⁰, would otherwise lose every digit of the residual
// variances to the sums of squares, or overflow.
struct Design {
  int visits = 5;
  double delta = 0.5, sd = 1, rho = 0.6;
  bool ar1 = false;
  double drop[2] = {0, 0};  // proportion lost by the last visit, by group
  std::array<double, TMAX * TMAX> sigma{}, chol{};  // the correlation matrix (the covariance in units of σ²)
  /// Probability that a patient of group g is still observed at visit t (1 … T).
  double retained(int g, int t) const { return std::pow(1 - drop[g], static_cast<double>(t) / visits); }
};

Design design_of(const Scenario& s) {
  Design d;
  d.visits = static_cast<int>(s("visits"));
  d.delta = s("delta");
  d.sd = s("sd");
  d.rho = s("rho");
  d.ar1 = s.is("corr", "ar1");
  d.drop[0] = s("drop1");
  d.drop[1] = s("drop2");
  const int T = d.visits;
  for (int a = 0; a < T; ++a)
    for (int b = 0; b < T; ++b) {
      d.sigma[a * TMAX + b] = a == b ? 1 : (d.ar1 ? std::pow(d.rho, std::abs(a - b)) : d.rho);
    }
  // Cholesky factor of the true covariance (positive definite for −1/(T − 1) < ρ < 1, or |ρ| < 1 for AR(1))
  for (int j = 0; j < T; ++j) {
    double v = d.sigma[j * TMAX + j];
    for (int k = 0; k < j; ++k) v -= d.chol[j * TMAX + k] * d.chol[j * TMAX + k];
    if (!(v > 0))
      throw CalcError("the correlation " + fnum(d.rho) + " is not possible with " + std::to_string(T) +
                      " visits (exchangeable correlations must exceed −1/(T − 1))");
    d.chol[j * TMAX + j] = std::sqrt(v);
    for (int i = j + 1; i < T; ++i) {
      double w = d.sigma[i * TMAX + j];
      for (int k = 0; k < j; ++k) w -= d.chol[i * TMAX + k] * d.chol[j * TMAX + k];
      d.chol[i * TMAX + j] = w / d.chol[j * TMAX + j];
    }
  }
  return d;
}

/// Expected variance of the difference at the last visit with n1 and n2 patients, in units of σ²: the inverse of
/// the expected GLS information Σ_k n π_k E_k' Σ_kk⁻¹ E_k under the true covariance and dropout pattern.  In the
/// errors e = L⁻¹y (Σ = LL') the visits are independent with variance 1, and the information on each visit's mean
/// is n times the probability of being observed there; the last visit's mean is the combination Σ_t L_Tt e_t, so
/// its variance is Σ_t L_Tt² / (n P(observed at visit t)).  (Inverting the information of the visits themselves
/// lost every digit once the correlation was within about 10⁻¹⁵ of 1.)
double expected_variance(const Design& d, double n1, double n2) {
  const int T = d.visits;
  double v = 0;
  for (int g = 0; g < 2; ++g) {
    const double n = g == 0 ? n1 : n2;
    for (int t = 0; t < T; ++t) {
      const double l = d.chol[(T - 1) * TMAX + t];
      v += l * l / (n * d.retained(g, t + 1));
    }
  }
  return v;
}

enum Analysis { MMRM = 0, COMPLETERS = 1 };

// Bits of sim::Outcome::also: each analysis's rejection (bit MMRM or COMPLETERS), then why an analysis failed.
constexpr unsigned EMPTY = 2, FEW = 3, ARITHMETIC = 4, NO_COMPLETERS = 5, FLAT_COMPLETERS = 6;

const std::vector<simc::Cause>& causes_of(Analysis an) {
  using simc::Failures;
  static const std::vector<simc::Cause> mmrm = {
      {EMPTY, Failures::Design, "a group had no patients left at a visit"},
      {FEW, Failures::Design,
       "too few patients were left at a visit to estimate its covariances with the earlier visits"},
      {ARITHMETIC, Failures::Numerical, "the arithmetic of the MMRM fit broke down"}};
  static const std::vector<simc::Cause> completers = {
      {NO_COMPLETERS, Failures::Design, "fewer than 2 patients of a group completed the study"},
      {FLAT_COMPLETERS, Failures::Numerical, "the completers' outcomes had no spread"}};
  return an == MMRM ? mmrm : completers;
}

sim::Outcome mmrm_study(const Design& d, long n1, long n2, long rep, std::uint64_t seed, Analysis an, bool z,
                        double alpha, Alt a) {
  const int T = d.visits;
  // The MMRM is fitted to the errors e = L⁻¹y, whose visits are independent with variance 1, rather than to the
  // outcomes y = Le: its fit is unchanged when each visit is replaced by a combination of itself and the earlier
  // visits (L⁻¹ is lower triangular, so a patient's observed visits stay observed, and the regression of each visit
  // on the earlier ones changes with them), and the difference at the last visit is that of the combination given
  // by L's last row.  With visits correlated within about 10⁻¹² of 1, their differences in y kept too few digits
  // in the sums of squares, and the fits failed.
  mmrm::Stats st;
  st.visits = T;
  double csum[2] = {0, 0}, css[2] = {0, 0}, ccount[2] = {0, 0};
  const long n[2] = {n1, n2};
  for (int g = 0; g < 2; ++g) {
    Rng rng(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g + 1));
    const double ldrop = d.drop[g] > 0 ? std::log1p(-d.drop[g]) : 0;
    for (long i = 0; i < n[g]; ++i) {
      double e[TMAX];
      for (int t = 0; t < T; ++t) e[t] = rng.normal();
      const double u = rng.uniform();
      // the last visit observed: retained through visit t when u ≤ (1 − dropout)^(t/T)
      int k = T;
      if (d.drop[g] > 0) {  // compared before conversion: a tiny dropout rate puts it far beyond any int
        const double last = std::floor(T * std::log(u) / ldrop + 1e-12);
        k = last >= T ? T : static_cast<int>(last);
      }
      if (k >= 1) st.add(g, k, e);
      if (k == T) {
        double y = 0;  // the outcome at the last visit
        for (int c = 0; c < T; ++c) y += d.chol[(T - 1) * TMAX + c] * e[c];
        csum[g] += y;
        css[g] += y * y;
        ccount[g] += 1;
      }
    }
  }
  sim::Outcome o;
  const double shift = d.delta / d.sd;  // the difference at the last visit, left out of the simulated outcomes
  // completers: pooled t-test at the last visit
  bool comp = false;
  unsigned comp_failed = NO_COMPLETERS;
  if (ccount[0] >= 2 && ccount[1] >= 2) {
    const double m1 = csum[0] / ccount[0], m2 = csum[1] / ccount[1];
    const double ss = css[0] - ccount[0] * m1 * m1 + css[1] - ccount[1] * m2 * m2, df = ccount[0] + ccount[1] - 2;
    const double se = std::sqrt(ss / df * (1 / ccount[0] + 1 / ccount[1]));
    comp_failed = se > 0 ? 0u : FLAT_COMPLETERS;
    comp = se > 0 && simc::p_value(simc::t_stat(m1 - m2 + shift, se), z ? INFINITY : df, a) < alpha;
  }
  const mmrm::Fit f = mmrm::fit(st, true, &d.chol[(T - 1) * TMAX]);
  bool mm = false;
  unsigned mm_failed = 0;
  if (f.ok) {
    o.df = z ? INFINITY : f.observed_last - 2;
    mm = o.df > 0 && simc::p_value(simc::t_stat(f.delta + shift, f.se), o.df, a) < alpha;
    o.estimate = (f.delta + shift) * d.sd;
    o.se = f.se * d.sd;
  } else {
    mm_failed = f.failure == mmrm::Failure::Empty ? EMPTY : (f.failure == mmrm::Failure::Few ? FEW : ARITHMETIC);
  }
  o.also = (mm ? 1u : 0u) | (comp ? 2u : 0u) | (mm_failed ? 1u << mm_failed : 0u) |
           (comp_failed ? 1u << comp_failed : 0u);
  if (an == MMRM ? !f.ok : comp_failed != 0) {
    o.failed = true;
  } else {
    o.reject = an == MMRM ? mm : comp;
  }
  o.extra = ccount[0] + ccount[1];
  return o;
}

Procedure mmrm_proc() {
  Procedure p;
  p.id = "sim.mmrm";
  p.name = "Two Means at the Last Visit with Dropout: Mixed Model for Repeated Measures (MMRM), by Simulation";
  p.category = "Simulation";
  p.subcategory = "Longitudinal Studies";
  p.order = 2;
  p.simulation = true;
  p.summary =
      "Sample size (or power) for a parallel-group trial whose primary analysis is a mixed model for repeated "
      "measures (MMRM: a mean for every group and visit, unstructured covariance) of the difference at the last "
      "visit, with patients dropping out over time; by simulating whole trials and fitting each by REML.";
  p.keywords = {"simulation", "monte carlo", "mmrm", "mixed model for repeated measures", "repeated measures",
                "longitudinal", "dropout", "attrition", "missing data", "missing at random", "unstructured covariance",
                "change from baseline", "last visit", "visits", "confirmatory trial", "completers", "reml",
                "kenward-roger", "estimand"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total Sample Size)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: δ ≠ 0)", "One-sided (H1: δ > 0)", "One-sided (H1: δ < 0)"),
      select("analysis", "Analysis",
             {{"mmrm", "MMRM: all visits, unstructured covariance (REML)", {}},
              {"completers", "t-test on the patients who complete the study", {}}},
             "mmrm")
          .grp(G_TEST)
          .hint("The analysis that decides the power and the sample size; the other is run on the same simulated "
                "trials and reported alongside. The MMRM uses the earlier visits of the patients who leave."),
      select("df", "Degrees of Freedom",
             {{"last", "t with the patients observed at the last visit − 2", {}},
              {"z", "None (normal approximation)", {}}},
             "last")
          .grp(G_TEST)
          .hint("Reference distribution of the test statistic at the last visit."),
      power_param(),
      alpha_param(),
      num("n1", "N1 (Patients Randomized to Group 1)", "100")
          .grp(G_N)
          .sym("N1")
          .whole()
          .at_least(3)
          .sugg({"100", "50 to 200 by 50"})
          .hint("Number of patients randomized to group 1, before dropout; group 2 gets N1 × the allocation ratio, "
                "rounded up.")
          .show(when(not_solving({"n"}))),
      num("r", "Allocation Ratio N2/N1", "1")
          .grp(G_N)
          .sym("N2/N1")
          .above(0)
          .sugg({"1", "2"})
          .hint("Size of group 2 relative to group 1; 1 gives equal groups."),
      num("drop1", "Dropout by the Last Visit, Group 1", "0.2")
          .grp(G_N)
          .sym("Dropout 1")
          .range(0, 1, false, true)
          .sugg({"0.2", "0.1 0.2 0.3"})
          .hint("Proportion of group 1 no longer observed at the last visit. Patients leave at a constant rate "
                "between visits, completely at random, and their earlier visits are analyzed."),
      num("drop2", "Dropout by the Last Visit, Group 2", "0.2")
          .grp(G_N)
          .sym("Dropout 2")
          .range(0, 1, false, true)
          .sugg({"0.2", "0.1 0.2 0.3"})
          .hint("Proportion of group 2 no longer observed at the last visit."),
      num("delta", "δ (Difference in Means at the Last Visit)", "0.4")
          .grp(G_EFFECT)
          .sym("δ")
          .sugg({"0.4", "0.3 0.4 0.5"})
          .hint("Mean of group 1 minus group 2 at the last visit (e.g. in the change from baseline); the difference "
                "grows linearly from the first visit."),
      num("sd", "σ (SD at Each Visit)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .sugg({"1"})
          .hint("Standard deviation of the outcome at each visit within a group."),
      num("visits", "T (Visits after Baseline)", "5")
          .grp(simc::G_MODEL)
          .sym("T")
          .whole()
          .range(1, TMAX)
          .single()
          .sugg({"5", "4", "8"})
          .hint("Number of scheduled visits after randomization, up to 12; the last one is the primary time point."),
      select("corr", "Correlation between Visits",
             {{"cs", "Exchangeable (the same for any two visits)", {}},
              {"ar1", "AR(1) (decreasing with the distance between visits)", {}}},
             "cs")
          .grp(simc::G_MODEL)
          .hint("How measurements of the same patient are correlated; the analysis estimates an unstructured "
                "covariance whatever the truth."),
      num("rho", "ρ (Correlation between Visits)", "0.6")
          .grp(simc::G_MODEL)
          .sym("ρ")
          .range(-1, 1, true, true)
          .sugg({"0.6", "0.4 0.6 0.8"})
          .hint("Correlation between two measurements of a patient (adjacent visits for AR(1)). Higher correlation "
                "lets the earlier visits of patients who leave say more about the last one."),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double alpha = s("alpha"), ratio = s("r");
    const Design d = design_of(s);
    const Analysis an = s.is("analysis", "completers") ? COMPLETERS : MMRM;
    const bool z = s.is("df", "z");
    auto n2_of = [&](double n1) { return std::ceil(ratio * n1 - 1e-9); };
    auto approx = [&](double n1) {
      const double n2 = n2_of(n1);
      const double keep1 = d.retained(0, d.visits), keep2 = d.retained(1, d.visits);
      const double var = an == MMRM ? expected_variance(d, n1, n2) : 1 / (n1 * keep1) + 1 / (n2 * keep2);
      const double df = n1 * keep1 + n2 * keep2 - 2, ncp = d.delta / d.sd / std::sqrt(var);
      return z || !(df > 0) ? z_power(ncp, alpha, a) : t_power(ncp, df, alpha, a);
    };
    if (s.solve == "n") {
      if (d.delta == 0) throw CalcError("δ is zero; power cannot exceed alpha");
      if ((a == Alt::Upper && d.delta < 0) || (a == Alt::Lower && d.delta > 0))
        throw CalcError("δ has the sign opposite to the one-sided alternative");
    }
    // the unstructured covariance regresses each visit on the group and the earlier visits: the last visit needs
    // at least T + 2 patients, so fewer expected there would fail most fits
    const double need = d.visits + 2;
    auto at_last = [&](double n1) { return n1 * d.retained(0, d.visits) + n2_of(n1) * d.retained(1, d.visits); };
    double lo = s.solve != "n" ? 3 : solve::first_whole([&](double x) { return n2_of(x) >= 3; }, 3, 1e5,
                                                     "the allocation ratio is too small: group 2 would have "
                                                     "fewer than 3 patients");
    if (s.solve == "n" && an == MMRM)
      lo = solve::first_whole([&](double x) { return at_last(x) >= need; }, lo, 1e5,
                              "too few patients would remain at the last visit for the MMRM fit");
    double n1;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = lo;
      o.hi = 1e5;
      n1 = solve::min_n(approx, s("power"), o, &ctx);
    } else {
      n1 = s("n1");
      if (n2_of(n1) < 3) throw CalcError("group 2 needs at least 3 patients");
      if (an == MMRM && at_last(n1) < need)
        throw CalcError("with " + std::to_string(d.visits) + " visits the MMRM needs at least " + fint(need) +
                        " patients at the last visit, but about " + fnum(at_last(n1), 3) +
                        " remain after dropout; randomize more patients");
    }
    if (ctx.approximate) {
      simc::put_approximation(s, r, approx(n1));
    } else {
      const long reps = static_cast<long>(s("reps"));
      const auto seed = static_cast<std::uint64_t>(s("seed"));
      auto power_at = [&](double m1, long R, double p0, double p1) {
        simc::check_study_size((m1 + n2_of(m1)) * d.visits, "measurements");
        const long k1 = static_cast<long>(m1), k2 = static_cast<long>(n2_of(m1));
        return simc::simulate_tests(
            R, [&](long rep) { return mmrm_study(d, k1, k2, rep, seed, an, z, alpha, a); }, ctx, p0, p1);
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
      simc::put_summary(s, r, sm, causes_of(an));
      r.set("approx_power", approx(n1));
      r.set("power_mmrm", simc::other_power(sm, MMRM, causes_of(MMRM)));
      r.set("power_completers", simc::other_power(sm, COMPLETERS, causes_of(COMPLETERS)));
      if (std::isfinite(sm.mean_estimate)) r.set("mean_delta", sm.mean_estimate);
      if (std::isfinite(sm.mean_extra)) r.set("completers", sm.mean_extra);
    }
    put_two_group(r, {n1, n2_of(n1)});
    if (!r.has("completers")) r.set("completers", n1 * d.retained(0, d.visits) + n2_of(n1) * d.retained(1, d.visits));
    r.set("delta", d.delta);
    r.set("sd", d.sd);
    r.set("rho", d.rho);
    r.set("visits", static_cast<double>(d.visits));
    r.set("drop1", d.drop[0]);
    r.set("drop2", d.drop[1]);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("completers", "Completers", "num", "Average number of patients observed at the last visit"),
                   col("delta", "δ", "num", "Difference in means at the last visit"),
                   col("mean_delta", "Mean δ̂", "num", "Average estimated difference over the simulated trials"),
                   col("sd", "σ", "num", "SD at each visit"), col("rho", "ρ", "num", "Correlation between visits"),
                   col("visits", "T", "int", "Visits after baseline"),
                   col("drop1", "Dropout 1", "num", "Proportion of group 1 lost by the last visit"),
                   col("drop2", "Dropout 2", "num", "Proportion of group 2 lost by the last visit"),
                   col("power_mmrm", "Power (MMRM)", "prob", "Simulated power of the MMRM analysis"),
                   col("power_completers", "Power (Completers)", "prob", "Simulated power of the completers' t-test"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns()) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string corr = s.is("corr", "ar1") ? "AR(1) correlation " : "exchangeable correlation ";
    const std::string analysis = s.is("analysis", "completers")
                                     ? "a t-test on the patients who complete the study"
                                     : "a mixed model for repeated measures with an unstructured covariance (REML)";
    std::string drop = r.num("drop1") == r.num("drop2")
                           ? fpct(r.num("drop1"), 0) + " of each group"
                           : fpct(r.num("drop1"), 0) + " of group 1 and " + fpct(r.num("drop2"), 0) + " of group 2";
    std::string t = "Randomizing " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " patients (total " +
                    fint(r.num("n")) + ") would " + simc::power_phrase(s, r) + " to detect a difference of " +
                    fnum(r.num("delta")) + " at the last of " + fint(r.num("visits")) + " visits (σ = " +
                    fnum(r.num("sd")) + ", " + corr + fnum(r.num("rho")) + "), with " + drop +
                    " lost by the last visit, analyzing with " + analysis + " and a " + alt_words(alt(s)) +
                    " test at α = " + fnum(r.num("alpha")) + ".";
    t += simc::others_sentence(r, "trials", {{"MMRM", "power_mmrm"}, {"completers' t-test", "power_completers"}});
    return t + simc::method_note(r);
  };
  p.definitions = {{"N1, N2", "Number of patients randomized to each group; N = N1 + N2."},
                   {"δ", "Difference between the groups' means at the last visit."},
                   {"Dropout", "Proportion of a group no longer observed at the last visit (monotone, completely at "
                               "random, at a constant rate between visits)."},
                   {"Completers", "Average number of patients observed at the last visit."},
                   {"Power", "Share of the simulated trials whose test was significant, with its 95% interval "
                             "(Wilson) and Monte Carlo standard error."}};
  p.notes = {"The MMRM has a mean for every group and visit and an unstructured covariance matrix common to the "
             "groups; with monotone dropout its REML fit has a closed form (each visit regressed on the earlier "
             "ones), and the difference at the last visit is the generalized least squares estimate.",
             "Degrees of freedom: the patients observed at the last visit minus 2, close to the Kenward–Roger value "
             "for these designs; software such as the mmrm R package or SAS PROC MIXED gives similar power.",
             "Dropout is completely at random. Under missingness at random the MMRM stays unbiased while the "
             "completers' analysis does not; the power differs little.",
             "A trial whose MMRM cannot be fitted (a group with no patients left at a visit, or too few to estimate "
             "a visit's covariances with the earlier visits) counts as a failed analysis of the MMRM (not "
             "rejecting), and one with fewer than 2 completers in a group as a failed analysis of the completers' "
             "t-test.",
             "The quick approximation is the expected generalized least squares variance under the true covariance "
             "and the expected dropout pattern (Lu, Luo & Chen 2008)."};
  p.references = {"Mallinckrodt, C. H., Clark, W. S. & David, S. R. (2001). Accounting for dropout bias using "
                  "mixed-effects models. Journal of Biopharmaceutical Statistics 11, 9–21.",
                  "Lu, K., Luo, X. & Chen, P.-Y. (2008). Sample size estimation for repeated measures analysis in "
                  "randomized clinical trials with missing data. International Journal of Biostatistics 4(1), 9.",
                  "Anderson, T. W. (1957). Maximum likelihood estimates for a multivariate normal distribution when "
                  "some observations are missing. Journal of the American Statistical Association 52, 200–203.",
                  simc::REF_SIM};
  p.curve_param = "n1";
  p.curve_label = "N1 (Patients Randomized to Group 1)";
  p.curve_min = 3;
  p.example_solve = "n";
  p.example = {{"power", "0.9"}, {"delta", "0.4"}, {"drop1", "0.2"}, {"drop2", "0.2"}, {"reps", "500"}};
  return p;
}
SSC_REGISTER(mmrm_proc);

}  // namespace
}  // namespace ssc::procs
