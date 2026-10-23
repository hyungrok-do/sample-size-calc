// Simulation-based power for a cross-sectional stepped-wedge trial analyzed
// with a linear mixed model (sim.sw): S groups of g clusters switch from
// control to intervention at successive steps, over S + 1 periods with m new
// subjects per cluster and period.
//
// The data follow the model of Hooper et al. (2016): period effects, the
// treatment effect, a random cluster effect and a random cluster-period
// effect (the cluster autocorrelation, CAC, is the correlation of a cluster's
// period means; CAC = 1 is the model of Hussey & Hughes 2007).  Each study is
// fitted by REML (sim/lmm.hpp) with a random cluster intercept on the
// subjects' outcomes (Hussey & Hughes), or with random cluster and
// cluster-period intercepts, fitted on the cluster-period means (equivalent
// with equal cluster-period sizes when the cluster-period variance is
// positive); the treatment effect is tested with a Satterthwaite t-test.
// The quick approximation is the generalized least squares variance under
// the true model, which for CAC = 1 is the variance of Hussey & Hughes.
#include <algorithm>
#include <cstdint>

#include "procedures/sim_common.hpp"
#include "sim/lmm.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;
using sim::Rng;
namespace lmm = sim::lmm;

constexpr int SMAX = 10;  // steps (S + 1 periods, S + 2 fixed effects ≤ PMAX)

struct Design {
  int steps = 4;
  double m = 20;
  double delta = 0.3;
  double sa = 0, sc = 0, se = 1;  // SDs of the cluster effect, the cluster-period effect and a subject's error
                                  // in units of the total SD σ
  double unit = 1;                // σ: the outcome's unit in the simulated trials and in the approximation
  bool hooper = false;            // analysis with a cluster-period effect (on the cluster-period means)
  int periods() const { return steps + 1; }
  int p() const { return steps + 2; }
};

/// Fixed-effect row of period j for a cluster switching at step s (1 … S): intercept, period j ≥ 1, treatment.
void design_row(const Design& d, int s, int j, double* x) {
  std::fill(x, x + d.p(), 0.0);
  x[0] = 1;
  if (j >= 1) x[j] = 1;
  x[d.p() - 1] = j >= s ? 1 : 0;
}

/// Row of period j for a cluster switching at step s minus the cluster's mean row (every period has as many
/// rows): what is left of it within the cluster, exactly 0 for the intercept, which a large cluster variance
/// would otherwise weigh against digits.  Each period is 1/T of the mean row, and the treatment (T − s)/T.
void within_row(const Design& d, int s, int j, double* dev) {
  const int T = d.periods(), p = d.p();
  design_row(d, s, j, dev);
  dev[0] = 0;
  for (int a = 1; a < T; ++a) dev[a] -= 1.0 / T;
  dev[p - 1] -= static_cast<double>(T - s) / T;
}

/// One pattern per sequence: the design of a cluster switching at step s.  Subject level (m rows per period,
/// a random cluster intercept) or cluster-period means (one row per period).
lmm::Data patterns(const Design& d, double clusters_per_step) {
  lmm::Data data;
  data.p = d.p();
  data.q = 1;
  const int T = d.periods(), p = d.p();
  const double w = d.hooper ? 1 : d.m;  // rows per period
  for (int s = 1; s <= d.steps; ++s) {
    lmm::Pattern pt;
    pt.count = clusters_per_step;
    pt.n = std::llround(w * T);
    double x[lmm::PMAX], dev[lmm::PMAX];
    for (int j = 0; j < T; ++j) {
      design_row(d, s, j, x);
      for (int a = 0; a < p; ++a) pt.zx[a] += w * x[a];
      pt.zz[0] += w;
    }
    pt.set_coefficients(p, 1);
    // X'(I − P)X: the rows' deviations from the cluster's mean row (every period has w rows)
    for (int j = 0; j < T; ++j) {
      within_row(d, s, j, dev);
      for (int a = 0; a < p; ++a)
        for (int b = 0; b < p; ++b) pt.xx[a * p + b] += w * dev[a] * dev[b];
    }
    data.patterns.push_back(pt);
  }
  return data;
}

/// Variance of the treatment estimate under the true model, in units of σ² (generalized least squares on the
/// cluster-period means: residual variance σc² + σe²/m, random cluster intercept σa²).
double gls_variance(const Design& d, double g) {
  Design means = d;
  means.hooper = true;
  const lmm::Data data = patterns(means, g);
  const double psi = d.sc * d.sc + d.se * d.se / d.m;
  const double rel[1] = {d.sa * d.sa / psi};
  double cov[lmm::PMAX * lmm::PMAX];
  if (!lmm::beta_covariance(data, psi, rel, cov)) return INFINITY;
  return cov[(d.p() - 1) * d.p() + d.p() - 1];
}

/// One simulated trial.  `means` and `subjects` are the patterns of the cluster-period means and of the
/// subjects; the analysis with cluster-period effects fits the means, and when that fit puts the
/// cluster-period variance below zero (the means vary less than their subjects alone would make them) it
/// refits on the subjects without it, as REML constrained to a non-negative variance (lme4) does.
///
/// The period effects and the treatment effect are estimated freely, so adding them to the outcome changes
/// the estimate of δ by exactly δ and nothing else: the trial is simulated without them, in units of σ, and δ/σ
/// is added to the estimate.
sim::Outcome sw_study(const Design& d, const lmm::Data& means, const lmm::Data& subjects, long g, long rep,
                      std::uint64_t seed, int dfm, double alpha, Alt a) {
  lmm::Data subj = subjects, mean = means;
  const int T = d.periods(), p = d.p();
  const double dm = d.m;
  double within = 0;  // the subjects' sum of squares around their cluster-period means
  for (int s = 1; s <= d.steps; ++s) {
    lmm::Pattern& ps = subj.patterns[static_cast<std::size_t>(s - 1)];  // add_sums counts the clusters
    lmm::Pattern& pm = mean.patterns[static_cast<std::size_t>(s - 1)];
    double dev_x[SMAX + 1][lmm::PMAX];  // the periods' rows within the cluster, the same for all its clusters
    for (int j = 0; j < T; ++j) within_row(d, s, j, dev_x[j]);
    for (long c = 0; c < g; ++c) {
      Rng rng(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(s), static_cast<std::uint64_t>(c));
      const double cluster = d.sa * rng.normal();
      // the cluster-period means are the cluster effect plus r_j; the deviations of r_j from their mean (and the
      // subjects' sums of squares around the means) are what lies within the cluster
      double r[SMAX + 1], rsum = 0, ss = 0;
      for (int j = 0; j < T; ++j) {
        r[j] = d.sc * rng.normal() + d.se / std::sqrt(dm) * rng.normal();
        ss += dm > 1 ? d.se * d.se * rng.chisq(dm - 1) : 0;
        rsum += r[j];
      }
      within += ss;
      double xs[lmm::PMAX] = {}, xm[lmm::PMAX] = {}, ym = 0;
      for (int j = 0; j < T; ++j) {
        const double* x = dev_x[j];
        const double dev = r[j] - rsum / T;
        for (int k = 0; k < p; ++k) {
          xs[k] += dm * dev * x[k];
          xm[k] += dev * x[k];
        }
        ym += dev * dev;
      }
      const double level = cluster + rsum / T;  // the cluster's mean
      ps.add_sums(p, 1, xs, &level, ss + dm * ym);
      if (d.hooper) pm.add_sums(p, 1, xm, &level, ym);
    }
  }
  sim::Outcome o;
  const lmm::Data* data = &subj;
  double theta = std::fmax(d.sa / d.se, 0.05);
  lmm::Fit f;
  if (d.hooper) {
    const double psi = d.sc * d.sc + d.se * d.se / dm;
    const double start = std::fmax(d.sa / std::sqrt(psi), 0.05);
    f = lmm::fit(mean, &start);
    const double cells = static_cast<double>(g * d.steps * T);
    const bool boundary = dm > 1 && f.ok && f.sigma2 < within / (cells * (dm - 1)) / dm;
    if (f.ok && !boundary) {
      data = &mean;
    } else {
      f = lmm::fit(subj, &theta);
    }
  } else {
    f = lmm::fit(subj, &theta);
  }
  if (!f.ok) {
    o.failed = true;
    return o;
  }
  double c[lmm::PMAX] = {};
  c[p - 1] = 1;
  const double clusters = static_cast<double>(g * d.steps);
  const lmm::Wald w = lmm::wald(*data, f, c, dfm, clusters - 2);
  const double estimate = w.estimate + d.delta / d.unit;
  o.reject = simc::p_value(simc::t_stat(estimate, w.se), w.df, a) < alpha;
  o.estimate = estimate * d.unit;
  o.se = w.se * d.unit;
  o.df = w.df;
  return o;
}

Design design_of(const Scenario& s) {
  Design d;
  d.steps = static_cast<int>(s("steps"));
  if (d.steps > SMAX) throw CalcError("at most 10 steps");
  d.m = s("m");
  d.delta = s("delta");
  const double sigma = s("sigma"), icc = s("icc"), cac = s("cac");
  d.sa = std::sqrt(icc * cac);
  d.sc = std::sqrt(icc * (1 - cac));
  d.se = std::sqrt(1 - icc);
  d.unit = sigma;
  d.hooper = s.is("analysis", "hooper");
  return d;
}

Procedure stepped_wedge_sim() {
  Procedure p;
  p.id = "sim.sw";
  p.name = "Stepped-Wedge Trial with a Mixed Model (Cluster-Period Effects), by Simulation";
  p.category = "Simulation";
  p.subcategory = "Cluster-Randomized Trials";
  p.order = 3;
  p.simulation = true;
  p.summary =
      "Clusters per step (or power) for a cross-sectional stepped-wedge trial, estimated by simulating whole "
      "trials and analyzing each with a linear mixed model with period effects (random cluster intercept, or "
      "random cluster and cluster-period intercepts; REML, Satterthwaite t-test). Handles a cluster "
      "autocorrelation below 1 and few clusters, where formulas are least reliable.";
  p.keywords = {"simulation", "monte carlo", "stepped wedge", "stepped-wedge", "sw-crt", "cluster randomized",
                "staggered rollout", "hussey hughes", "hooper", "cluster autocorrelation", "cac",
                "cluster-period", "mixed model", "lmm", "period effects", "satterthwaite", "reml", "icc",
                "implementation trial"};
  p.solve_for = {{"g", "Clusters per Step", "g", "Clusters per Step"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: δ ≠ 0)", "One-sided (H1: δ > 0)", "One-sided (H1: δ < 0)"),
      select("analysis", "Analysis Model",
             {{"hh", "Random cluster intercept (Hussey & Hughes)", {}},
              {"hooper", "Random cluster and cluster-period intercepts (Hooper et al.)", {}}},
             "hh")
          .grp(G_TEST)
          .hint("The mixed model fitted to each simulated trial. With a cluster autocorrelation below 1 the "
                "random-intercept model's test rejects too often; the model with cluster-period effects keeps "
                "the type I error near α, though with fewer than about 8 clusters both reject somewhat too often."),
      simc::df_param("t with K − 2 df (clusters)"),
      power_param(),
      alpha_param(),
      num("steps", "S (Number of Steps)", "4")
          .grp("Design")
          .sym("Steps")
          .whole()
          .range(2, SMAX)
          .single()
          .sugg({"4", "3", "6"})
          .hint("Number of times a new group of clusters switches from control to intervention; there are S + 1 "
                "measurement periods (2 to 10 steps: with one, the effect could not be told apart from the period)."),
      num("g", "g (Clusters Switching at Each Step)", "2")
          .grp("Design")
          .sym("Clusters/Step")
          .whole()
          .at_least(1)
          .sugg({"2", "1 2 3"})
          .hint("Number of clusters that switch at each step; the trial has K = S × g clusters.")
          .show(when(not_solving({"g"}))),
      num("m", "m (Subjects per Cluster per Period)", "20")
          .grp("Design")
          .sym("m")
          .whole()
          .at_least(1)
          .sugg({"20", "10 20 50"})
          .hint("Number of subjects measured in each cluster in each period (different subjects in every "
                "period)."),
      num("delta", "δ (Treatment Effect)", "0.3")
          .grp(G_EFFECT)
          .sym("δ")
          .sugg({"0.3", "0.2 0.3"})
          .hint("Difference between the mean outcome under intervention and under control."),
      num("sigma", "σ (Total Standard Deviation)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .sugg({"1"})
          .hint("Standard deviation of one subject's outcome, including the variation between clusters and "
                "periods."),
      num("icc", "ICC (Intracluster Correlation)", "0.05")
          .grp(G_EFFECT)
          .sym("ICC")
          .range(0, 1, false, true)
          .sugg({"0.01", "0.05", "0.1"})
          .hint("Correlation of two subjects of the same cluster in the same period: the share of the variance "
                "between clusters (and cluster-periods)."),
      num("cac", "CAC (Cluster Autocorrelation)", "1")
          .grp(G_EFFECT)
          .sym("CAC")
          .range(0, 1)
          .sugg({"1", "0.8", "0.5 0.8 1"})
          .hint("Correlation between a cluster's true means in two periods: 1 when the cluster effect is constant "
                "over time (Hussey & Hughes); values from 0.5 to 0.9 are common in practice."),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double alpha = s("alpha");
    const Design d = design_of(s);
    const int dfm = simc::df_method(s);
    auto approx = [&](double g) {
      const double v = gls_variance(d, g), clusters = g * d.steps;
      const double ncp = d.delta / d.unit / std::sqrt(v);
      // Satterthwaite's degrees of freedom for the treatment effect are large (it is estimated mostly within
      // clusters), so the normal distribution stands in for them
      return dfm != 1 || clusters <= 2 ? z_power(ncp, alpha, a) : t_power(ncp, clusters - 2, alpha, a);
    };
    if (s.solve == "g") {
      if (d.delta == 0) throw CalcError("δ is zero; power cannot exceed alpha");
      if ((a == Alt::Upper && d.delta < 0) || (a == Alt::Lower && d.delta > 0))
        throw CalcError("δ has the sign opposite to the one-sided alternative");
    }
    const double lo = d.steps >= 3 ? 1 : std::ceil(3.0 / d.steps);  // at least 3 clusters for K − 2 df
    double g;
    if (s.solve == "g") {
      solve::NSearch o;
      o.lo = lo;
      o.hi = 1e4;
      g = solve::min_n(approx, s("power"), o, &ctx);
    } else {
      g = s("g");
      if (g * d.steps < 3) throw CalcError("the trial needs at least 3 clusters");
    }
    if (ctx.approximate) {
      simc::put_approximation(s, r, approx(g));
    } else {
      const long reps = static_cast<long>(s("reps"));
      const auto seed = static_cast<std::uint64_t>(s("seed"));
      Design by_subject = d, by_mean = d;
      by_subject.hooper = false;
      by_mean.hooper = true;
      const lmm::Data subjects = patterns(by_subject, 0), means = patterns(by_mean, 0);
      auto power_at = [&](double c, long R, double p0, double p1) {
        simc::check_study_size(c * d.steps * d.periods(), "cluster-periods");
        const long k = static_cast<long>(c);
        return sim::simulate(
            R, [&](long rep) { return sw_study(d, means, subjects, k, rep, seed, dfm, alpha, a); }, ctx, p0, p1);
      };
      sim::Summary sm;
      if (s.solve == "g") {
        sim::SizeSearch o;
        o.lo = lo;
        o.hi = 1e4;
        o.guess = g;
        o.reps = reps;
        o.explore_reps = simc::explore_reps(reps);
        const sim::SizeResult res = sim::search_size(power_at, s("power"), o, ctx);
        g = res.n;
        sm = res.s;
      } else {
        sm = power_at(g, reps, 0, 1);
      }
      simc::put_summary(s, r, sm, simc::Failures::Numerical);
      r.set("approx_power", approx(g));
      if (std::isfinite(sm.mean_estimate)) r.set("mean_delta", sm.mean_estimate);
    }
    r.set("g", g);
    r.set("steps", static_cast<double>(d.steps));
    r.set("clusters", g * d.steps);
    r.set("periods", static_cast<double>(d.periods()));
    r.set("m", d.m);
    r.set("N", g * d.steps * d.periods() * d.m);
    r.set("delta", d.delta);
    r.set("icc", s("icc"));
    r.set("cac", s("cac"));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    for (auto x : {col("g", "Clusters/Step", "int", "Clusters switching at each step"),
                   col("clusters", "K", "int", "Clusters in the trial"), col("steps", "Steps", "int", "Steps"),
                   col("periods", "Periods", "int", "Measurement periods"),
                   col("m", "m", "int", "Subjects per cluster and period"),
                   col("N", "N", "int", "Observations in the trial"), col("delta", "δ", "num", "Treatment effect"),
                   col("mean_delta", "Mean δ̂", "num", "Average estimated treatment effect"),
                   col("icc", "ICC", "num", "Intracluster correlation"),
                   col("cac", "CAC", "num", "Cluster autocorrelation"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns()) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string model = s.is("analysis", "hooper")
                                  ? "random cluster and cluster-period intercepts"
                                  : "a random cluster intercept";
    const std::string test =
        s.is("df", "z") ? "z-test" : (s.is("df", "fixed") ? "t-test (K − 2 df)" : "Satterthwaite t-test");
    return "A stepped-wedge trial with " + fint(r.num("clusters")) + " clusters (" + fint(r.num("g")) +
           " switching to the intervention at each of " + fint(r.num("steps")) + " steps, " + fint(r.num("periods")) +
           " periods) and " + fint(r.num("m")) + " subjects per cluster and period (" + fint(r.num("N")) +
           " observations) would " + simc::power_phrase(s, r) + " to detect a treatment effect of " +
           fnum(r.num("delta")) + " (ICC " + fnum(r.num("icc")) + ", cluster autocorrelation " + fnum(r.num("cac")) +
           "), analyzing each trial with a linear mixed model with period effects and " + model + " and a " +
           alt_words(alt(s)) + " " + test + " at α = " + fnum(r.num("alpha")) + "." + simc::method_note(r);
  };
  p.definitions = {{"S, g, K", "Steps, clusters switching at each step, and clusters in the trial (K = S × g)."},
                   {"m", "Subjects measured in each cluster in each period (cross-sectional: new subjects every "
                         "period)."},
                   {"ICC", "Correlation of two subjects of the same cluster in the same period."},
                   {"CAC", "Correlation of a cluster's true means in two different periods (1: a constant cluster "
                           "effect)."},
                   {"Power", "Share of the simulated trials whose test of the treatment effect was significant, "
                             "with its 95% interval (Wilson) and Monte Carlo standard error."}};
  p.notes = {"Every simulated trial has period effects, the treatment effect, a random cluster effect and, when the "
             "cluster autocorrelation is below 1, a random cluster-period effect (Hooper et al. 2016).",
             "The model with cluster-period effects is fitted on the cluster-period means, which with equal "
             "cluster-period sizes gives the same estimate and test as the subject-level model whenever the "
             "cluster-period variance is estimated as positive; when it would be negative, the trial is refitted "
             "without it, as lme4's REML (which keeps variances non-negative) does.",
             "With CAC = 1 and the random-intercept analysis the quick approximation is Hussey & Hughes's variance, "
             "as in Stepped-Wedge Cluster-Randomized Designs; set CAC below 1 to see how much power the "
             "decay of the cluster effect costs.",
             "With few clusters, Satterthwaite's degrees of freedom keep the type I error closer to α than the "
             "normal approximation, but with fewer than about 8 clusters even a correctly specified model rejects "
             "somewhat too often (6–7% at α = 0.05, as lmerTest's tests do)."};
  p.references = {"Hussey, M. A. & Hughes, J. P. (2007). Design and analysis of stepped wedge cluster randomized "
                  "trials. Contemporary Clinical Trials 28, 182–191.",
                  "Hooper, R., Teerenstra, S., de Hoop, E. & Eldridge, S. (2016). Sample size calculation for "
                  "stepped wedge and other longitudinal cluster randomised trials. Statistics in Medicine 35, "
                  "4718–4728.",
                  "Baio, G. et al. (2015). Sample size calculation for a stepped wedge trial. Trials 16, 354.",
                  simc::REF_SIM};
  p.curve_param = "g";
  p.curve_label = "g (Clusters Switching at Each Step)";
  p.curve_min = 1;
  p.example_solve = "g";
  p.example = {{"power", "0.8"}, {"steps", "4"}, {"m", "20"}, {"delta", "0.3"}, {"icc", "0.05"}, {"cac", "0.8"},
               {"analysis", "hooper"}, {"reps", "300"}};
  return p;
}
SSC_REGISTER(stepped_wedge_sim);

}  // namespace
}  // namespace ssc::procs
