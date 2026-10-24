// Simulation-based power for a cluster-randomized comparison of two
// proportions analyzed with a logistic mixed model (sim.glmm.cluster).
//
// Each simulated cluster has a size from the size distribution, a normal
// random intercept on the logit scale and a binomial number of successes.
// Each study is analyzed by maximum likelihood with adaptive Gauss–Hermite
// quadrature (sim/glmm.hpp) and a Wald test of the group effect, with k1 + k2
// − 2 degrees of freedom by default, which keeps the type I error near α
// with few clusters (Li & Redden 2015).  The proportions may be given as
// population-averaged values or for a typical cluster, and the intracluster
// correlation on the proportion or the latent logistic scale; the quick
// approximation compares the population-averaged proportions with cluster
// weights averaged over the size distribution.
#include <algorithm>
#include <cstdint>

#include "procedures/sim_common.hpp"
#include "sim/glmm.hpp"
#include "sim/random.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;
using sim::Rng;
namespace glmm = sim::glmm;

constexpr double PI2_3 = 3.2898681336964528729;  // π²/3, the variance of the standard logistic distribution
constexpr double MAX_SIGMA = 6;                  // largest SD of the random intercept (latent ICC 0.916)

double logit(double p) { return std::log(p / (1 - p)); }
double expit(double t) { return 1 / (1 + std::exp(-t)); }

struct GlmmDesign {
  simc::ClusterSizes sizes{20, 0};
  double sigma = 0;   // SD of the random intercept (logit scale)
  double b0 = 0;      // cluster-specific intercept of group 2
  double b1 = 0;      // cluster-specific log odds ratio
  double pm1 = 0, pm2 = 0;      // population-averaged proportions
  double icc1 = 0, icc2 = 0;    // proportion-scale ICCs
  double icc_latent = 0;
};

GlmmDesign glmm_design(const Scenario& s) {
  GlmmDesign d;
  d.sizes = simc::ClusterSizes(s("m"), s("cv"));
  const double p1 = s("p1"), p2 = s("p2"), icc = s("icc");
  const bool marginal = !s.is("pscale", "conditional");
  auto beta_of = [&](double p, double sigma) { return marginal ? glmm::beta_for_marginal(p, sigma) : logit(p); };
  if (s.is("iccscale", "latent")) {
    d.sigma = std::sqrt(PI2_3 * icc / (1 - icc));
  } else {
    try {
      d.sigma = glmm::sigma_for_icc(icc, [&](double sigma) { return beta_of(p2, sigma); });
    } catch (const std::domain_error&) {
      throw CalcError("no logistic random-intercept model has this ICC with these proportions");
    }
  }
  // beyond this SD clusters are almost all-or-none, and the Wald test breaks down: at an SD of 8 it rejected
  // a true null hypothesis 7.5% of the time with 155 clusters per arm, and 13% of the fits failed; at 15
  // it rejected none of 40 studies of 300 clusters per arm with proportions of 0.3 and 0.2
  if (d.sigma > MAX_SIGMA) {
    const double latent = d.sigma * d.sigma / (d.sigma * d.sigma + PI2_3);
    throw CalcError("the ICC implies a between-cluster SD of " + text::format_g(d.sigma, 3) +
                    " on the logit scale (latent ICC " + text::format_g(latent, 3) + "); with an SD above 6 (latent "
                    "ICC 0.92) clusters are almost all-or-none, and the mixed model's Wald test breaks down");
  }
  d.b0 = beta_of(p2, d.sigma);
  d.b1 = beta_of(p1, d.sigma) - d.b0;
  d.pm1 = glmm::marginal_mean(d.b0 + d.b1, d.sigma);
  d.pm2 = glmm::marginal_mean(d.b0, d.sigma);
  d.icc1 = glmm::between_variance(d.b0 + d.b1, d.sigma) / (d.pm1 * (1 - d.pm1));
  d.icc2 = glmm::between_variance(d.b0, d.sigma) / (d.pm2 * (1 - d.pm2));
  d.icc_latent = d.sigma * d.sigma / (d.sigma * d.sigma + PI2_3);
  return d;
}

// Bits of sim::Outcome::also: why a study's fit failed, and whether its clusters all had one subject, so that it
// was analyzed by logistic regression (glmm::fit holds σ at 0).
constexpr unsigned SEPARATED = 0, RAN_OFF = 1, ARITHMETIC = 2, NO_SE = 3, CONCORDANT = 4, ONE_SUBJECT = 5;

const std::vector<simc::Cause>& causes() {
  using simc::Failures;
  static const std::vector<simc::Cause> c = {
      {SEPARATED, Failures::Design,
       "no subject of a group had the outcome, or every one had it, so the odds ratio had no finite estimate"},
      {CONCORDANT, Failures::Design,
       "the subjects of every cluster all had the outcome or all lacked it, so the between-cluster SD had no finite "
       "estimate"},
      {RAN_OFF, Failures::Design, "the estimates grew without bound"},
      {ARITHMETIC, Failures::Numerical, "the fit did not converge"},
      {NO_SE, Failures::Numerical,
       "the fit gave no standard error (its information matrix was not positive definite)"}};
  return c;
}

/// Whether no subject of a group had the outcome, or every one had it.  The likelihood then grows as the log odds
/// ratio runs off to ±∞, whatever σ is, so the odds ratio has no finite estimate, and the Wald statistic, whose
/// standard error grows faster, goes to 0.  Newton's method stops where the likelihood has grown flat, short of the
/// run-off guard or not, and may count as converged there: a study's data, not its fit, tell.
bool separated(const std::vector<glmm::Cluster>& data) {
  bool none[2] = {true, true}, all[2] = {true, true};
  for (const glmm::Cluster& c : data) {
    const int g = c.x == 1 ? 0 : 1;
    none[g] = none[g] && c.s == 0;
    all[g] = all[g] && c.s == c.m;
  }
  return none[0] || all[0] || none[1] || all[1];
}

/// Whether the subjects of every cluster agree (all with the outcome or all without it) and a cluster has more than
/// one.  For given population-averaged proportions such a cluster is likelier the larger σ is, so the likelihood
/// grows as σ runs off to ∞.  β1 grows with σ, and so does its standard error, so a fit that stops at a large σ
/// still tests (and may reject), but one that fails has met the design, not the arithmetic.
bool concordant(const std::vector<glmm::Cluster>& data) {
  bool agree = true, larger = false;
  for (const glmm::Cluster& c : data) {
    agree = agree && (c.s == 0 || c.s == c.m);
    larger = larger || c.m > 1;
  }
  return agree && larger;
}

/// Expected weight of a cluster, size / (p(1 − p)(1 + (size − 1)ρ)), over the size distribution.
double expected_weight(const GlmmDesign& d, double p, double icc) {
  double w = 0;
  const std::vector<double> grid = d.sizes.grid();
  for (double m : grid) w += m / (p * (1 - p) * (1 + (m - 1) * icc));
  return w / static_cast<double>(grid.size());
}

sim::Outcome glmm_study(const GlmmDesign& d, long k1, long k2, long rep, std::uint64_t seed, bool z, double alpha,
                        Alt a) {
  // clusters with the same (group, size, successes) enter the likelihood once, with their count:
  // sort the triples packed into 64 bits (group, then size and successes below 2³¹ each)
  std::vector<std::uint64_t> keys;  // freed with the study
  keys.reserve(static_cast<std::size_t>(k1 + k2));
  for (int g = 0; g < 2; ++g)
    for (long c = 0; c < (g == 0 ? k1 : k2); ++c) {
      Rng rng(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g + 1),
              static_cast<std::uint64_t>(c));
      const long m = d.sizes.draw(rng);
      const double eta = d.b0 + (g == 0 ? d.b1 : 0) + d.sigma * rng.normal();
      const long s = rng.binomial(m, expit(eta));
      keys.push_back(static_cast<std::uint64_t>(g) << 62 | static_cast<std::uint64_t>(m) << 31 |
                     static_cast<std::uint64_t>(s));
    }
  std::sort(keys.begin(), keys.end());
  std::vector<glmm::Cluster> data;
  for (std::size_t i = 0, j = 0; i < keys.size(); i = j) {
    while (j < keys.size() && keys[j] == keys[i]) ++j;
    glmm::Cluster c;
    c.x = keys[i] >> 62 == 0 ? 1 : 0;
    c.m = static_cast<double>(keys[i] >> 31 & 0x7fffffff);
    c.s = static_cast<double>(keys[i] & 0x7fffffff);
    c.count = static_cast<double>(j - i);
    data.push_back(c);
  }
  const unsigned one_subject =
      std::all_of(data.begin(), data.end(), [](const glmm::Cluster& c) { return c.m == 1; }) ? 1u << ONE_SUBJECT : 0u;
  if (separated(data)) {  // not fitted
    sim::Outcome o;
    o.failed = true;
    o.also = 1u << SEPARATED | one_subject;
    return o;
  }
  sim::Outcome o;
  const glmm::Fit f = glmm::fit(data, d.b0, d.b1, d.sigma);
  if (!f.ok || !(f.se_b1 > 0)) {
    o.failed = true;
    const glmm::Failure why = f.ok ? glmm::Failure::Information : f.failure;
    o.also = 1u << (concordant(data)                    ? CONCORDANT
                    : why == glmm::Failure::Separation  ? RAN_OFF
                    : why == glmm::Failure::Information ? NO_SE
                                                        : ARITHMETIC) |
             one_subject;
    return o;
  }
  o.also = one_subject;
  o.estimate = f.b1;
  o.se = f.se_b1;
  o.df = z ? INFINITY : static_cast<double>(k1 + k2 - 2);
  o.reject = simc::p_value(f.b1 / f.se_b1, o.df, a) < alpha;
  return o;
}

Procedure glmm_cluster() {
  Procedure p;
  p.id = "sim.glmm.cluster";
  p.name = "Two Proportions in a Cluster-Randomized Design with a Logistic Mixed Model, by Simulation";
  p.category = "Simulation";
  p.subcategory = "Cluster-Randomized Trials";
  p.order = 2;
  p.simulation = true;
  p.summary =
      "Number of clusters per arm (or power) for comparing two proportions when whole clusters are randomized, "
      "estimated by simulating studies and analyzing each with a logistic mixed model (random cluster intercept, "
      "adaptive Gauss–Hermite quadrature, Wald test). Cluster sizes may vary.";
  p.keywords = {"simulation", "monte carlo", "glmm", "generalized linear mixed model", "logistic mixed model",
                "mixed effects logistic regression", "random intercept", "cluster", "cluster randomized",
                "group randomized", "binary outcome", "yes no outcome", "proportions", "icc", "intraclass",
                "unequal cluster sizes", "glmer", "lme4", "adaptive quadrature"};
  p.solve_for = {{"k", "Number of Clusters", "k1", "Clusters in Group 1"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: P1 ≠ P2)", "One-sided (H1: P1 > P2)", "One-sided (H1: P1 < P2)"),
      select("test", "Test",
             {{"t", "Wald t-test (k1 + k2 − 2 df)", {}}, {"z", "Wald z-test (normal approximation)", {}}}, "t")
          .grp(G_TEST)
          .hint("With few clusters the Wald z-test rejects too often; the t reference distribution with k1 + "
                "k2 − 2 degrees of freedom keeps the type I error close to α."),
      power_param(),
      alpha_param(),
      num("k1", "k1 (Clusters in Group 1)", "15")
          .grp(G_N)
          .sym("k1")
          .whole()
          .at_least(2)
          .sugg({"15", "10 to 40 by 10"})
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
      num("p1", "P1 (Proportion in Group 1)", "0.3")
          .grp(G_EFFECT)
          .sym("P1")
          .range(0, 1, true, true)
          .sugg({"0.3", "0.25 0.3"})
          .hint("The share of subjects in group 1 (e.g. the intervention) with the outcome, between 0 and 1; "
                "whether over all clusters or in a typical cluster is chosen below."),
      num("p2", "P2 (Proportion in Group 2)", "0.2")
          .grp(G_EFFECT)
          .sym("P2")
          .range(0, 1, true, true)
          .sugg({"0.2"})
          .hint("The share of subjects in group 2 (e.g. the control) with the outcome, between 0 and 1, in the "
                "same sense as P1."),
      select("pscale", "The Proportions Are",
             {{"marginal", "Population-averaged (over all clusters)", {}},
              {"conditional", "For a typical cluster (random effect 0)", {}}},
             "marginal")
          .grp(G_EFFECT)
          .hint("Population-averaged proportions are what a group's overall share would be; the mixed model's "
                "odds ratio refers to a typical cluster and is larger."),
      num("icc", "ρ (Intracluster Correlation)", "0.05")
          .grp(G_EFFECT)
          .sym("ICC")
          .range(0, 1, false, true)
          .sugg({"0.01", "0.05", "0.01 0.05 0.1"})
          .hint("How alike the outcomes of subjects in the same cluster are, from 0 (not at all) to 1, on the "
                "scale chosen below; values from 0.01 to 0.1 are common on the proportion scale."),
      select("iccscale", "ICC Scale",
             {{"proportion", "Proportions: correlation of two outcomes in a cluster (group 2)", {}},
              {"latent", "Latent logistic scale: σ²/(σ² + π²/3)", {}}},
             "proportion")
          .grp(G_EFFECT)
          .hint("The proportion scale is the usual one for planning and for design effects; mixed-model "
                "software reports the latent scale."),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double alpha = s("alpha"), kr = s("kr"), p1 = s("p1"), p2 = s("p2");
    const GlmmDesign d = glmm_design(s);
    const bool z = s.is("test", "z");
    const double w1 = expected_weight(d, d.pm1, d.icc1), w2 = expected_weight(d, d.pm2, d.icc2);
    auto k2_of = [&](double k1) { return std::ceil(kr * k1 - 1e-9); };
    auto approx = [&](double k1) {
      const double se = std::sqrt(1 / (k1 * w1) + 1 / (k2_of(k1) * w2)), ncp = (d.pm1 - d.pm2) / se;
      return z ? z_power(ncp, alpha, a) : t_power(ncp, k1 + k2_of(k1) - 2, alpha, a);
    };
    if (s.solve == "k") {
      if (p1 == p2) throw CalcError("P1 equals P2");
      if ((a == Alt::Upper && p1 < p2) || (a == Alt::Lower && p1 > p2))
        throw CalcError("P1 − P2 has the sign opposite to the one-sided alternative");
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
      auto power_at = [&](double c1, long R, double p0, double p1_) {
        simc::check_study_size(c1 + k2_of(c1), "clusters");
        const long n1 = static_cast<long>(c1), n2 = static_cast<long>(k2_of(c1));
        return simc::simulate_tests(
            R, [&](long rep) { return glmm_study(d, n1, n2, rep, seed, z, alpha, a); }, ctx, p0, p1_);
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
      simc::put_summary(s, r, sm, causes());
      r.set("approx_power", approx(k1));
      if (std::isfinite(sm.mean_estimate)) r.set("mean_logor", sm.mean_estimate);
      if (sm.also[ONE_SUBJECT] > 0) r.set("logistic_studies", static_cast<double>(sm.also[ONE_SUBJECT]));
    }
    // the average size as entered (rounding the drawn sizes moves their mean only in the third decimal)
    const double msize = d.sizes.cv > 0 ? s("m") : d.sizes.grid().front();
    const double k2 = k2_of(k1);
    r.set("k1", k1);
    r.set("k2", k2);
    r.set("m", msize);
    r.set("cv", d.sizes.cv);
    r.set("n1", k1 * msize);
    r.set("n2", k2 * msize);
    r.set("n", (k1 + k2) * msize);
    r.set("p1", p1);
    r.set("p2", p2);
    r.set("pm1", d.pm1);
    r.set("pm2", d.pm2);
    r.set("or", std::exp(d.b1));
    r.set("icc", s("icc"));
    r.set("icc2", d.icc2);
    r.set("icc_latent", d.icc_latent);
    r.set("sigma_b", d.sigma);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    for (auto x : {col("k1", "k1", "int", "Clusters in group 1"), col("k2", "k2", "int", "Clusters in group 2"),
                   col("m", "m̄", "num", "Average cluster size"), col("cv", "CV", "num", "CV of the cluster sizes"),
                   col("n", "N", "num", "Expected number of subjects"),
                   col("pm1", "P1", "prob", "Population-averaged proportion in group 1"),
                   col("pm2", "P2", "prob", "Population-averaged proportion in group 2"),
                   col("or", "OR (cluster)", "num", "Odds ratio for a typical cluster (the model's group effect)"),
                   col("mean_logor", "Mean log OR̂", "num", "Average estimated log odds ratio"),
                   col("icc2", "ICC", "num", "Intracluster correlation of the outcomes in group 2 (proportion scale)"),
                   col("icc_latent", "ICC (latent)", "num", "Intracluster correlation on the latent logistic scale"),
                   col("sigma_b", "σb", "num", "SD of the random intercept (logit scale)"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns()) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string subjects = fnum(r.num("m")) == "1" ? " subject" : " subjects";
    const std::string sizes = r.num("cv") > 0 ? "clusters of " + fnum(r.num("m")) + subjects + " on average (CV " +
                                                     fnum(r.num("cv")) + ")"
                                               : "clusters of " + fnum(r.num("m")) + subjects;
    const std::string test = s.is("test", "z") ? "Wald z-test" : "Wald t-test with k1 + k2 − 2 degrees of freedom";
    // studies whose clusters all had one subject are analyzed by logistic regression, whose odds ratio is the
    // population-averaged one: every study with equal sizes of 1, otherwise those that drew no larger cluster
    const double logistic = r.has("logistic_studies") ? r.num("logistic_studies") : 0;
    const bool all_logistic = (r.num("cv") == 0 && r.num("m") == 1) || (logistic > 0 && logistic == r.num("reps"));
    const double pm1 = r.num("pm1"), pm2 = r.num("pm2");
    const double averaged = std::exp(std::log(pm1 / (1 - pm1)) - std::log(pm2 / (1 - pm2)));
    const std::string logistic_or = "the odds ratio of the population-averaged proportions" +
                                    (std::isfinite(averaged) ? " (" + fnum(averaged, 3) + ")" : std::string());
    const std::string analysis = all_logistic ? "by logistic regression (clusters of one subject leave the random "
                                                "intercept's SD unestimable, so it is held at 0) and a "
                                              : "with a logistic mixed model with a random cluster intercept and a ";
    std::string t = fint(r.num("k1")) + " and " + fint(r.num("k2")) + " " + sizes + " " + simc::power_phrase(s, r) +
                    " to detect population-averaged proportions of " + fnum(pm1) + " and " + fnum(pm2) +
                    " (an odds ratio of " + fnum(r.num("or"), 3) + " for a typical cluster; ICC " +
                    fnum(r.num("icc2"), 3) + " on the proportion scale, " + fnum(r.num("icc_latent"), 3) +
                    " on the latent scale), analyzing each study " + analysis + alt_words(alt(s)) + " " + test +
                    " at α = " + fnum(r.num("alpha")) + ".";
    if (all_logistic)
      t += " Logistic regression estimates " + logistic_or + ".";
    else if (logistic > 0)
      t += " In " + simc::grouped(logistic) + " of the " + simc::grouped(r.num("reps")) + " simulated studies every "
           "cluster had one subject, which leaves the random intercept's SD unestimable; " +
           (logistic == 1 ? "it was" : "those were") + " analyzed by logistic regression, with the SD held at 0, which "
           "estimates " + logistic_or + ".";
    return t + simc::method_note(r);
  };
  p.definitions = {{"k1, k2", "Number of clusters randomized to each group."},
                   {"P1, P2", "Population-averaged proportions: the overall share with the outcome in each group."},
                   {"OR (cluster)", "Odds ratio between the groups for a typical cluster (random intercept 0), the "
                                    "effect the mixed model estimates. A study whose clusters all have one subject "
                                    "is analyzed by logistic regression, which estimates the odds ratio of the "
                                    "population-averaged proportions, nearer 1 when the ICC is above 0."},
                   {"ICC", "Correlation between the outcomes of two subjects in the same cluster (proportion "
                           "scale), or σb²/(σb² + π²/3) on the latent logistic scale."},
                   {"Power", "Share of the simulated studies whose test of OR = 1 was significant, with its 95% "
                             "interval (Wilson) and Monte Carlo standard error."}};
  p.notes = {"Each cluster's number of successes is binomial given its random intercept; the model is fitted by "
             "maximum likelihood with 21-point adaptive Gauss–Hermite quadrature (as glmer with nAGQ = 21).",
             "With few clusters, maximum likelihood underestimates the between-cluster variance, so Wald tests "
             "reject somewhat more often than α (about 5.5% at α = 5% with 15 clusters of 30 per arm, with the t "
             "reference distribution). Set P1 = P2 to see the type I error of your design.",
             "The SD of the random intercept may be at most 6 on the logit scale (a latent ICC of 0.92): beyond "
             "it, clusters are almost all-or-none and the Wald test breaks down.",
             "With clusters of one subject the SD of the random intercept cannot be told apart from the groups' "
             "proportions, so it is held at 0 and each study is analyzed by ordinary logistic regression, which "
             "estimates the odds ratio of the population-averaged proportions; with unequal sizes, so is each study "
             "that draws no larger cluster.",
             "Studies whose fit fails count as not rejecting and are reported as failed fits, with the reason: no "
             "subject of a group with the outcome, or every one with it (no finite odds ratio; such a study is not "
             "fitted), the subjects of every cluster all alike (no finite between-cluster SD), or a fit that did not "
             "converge. No simulated power is given when no study could be fitted, or when fits did not converge in "
             "more than half of them.",
             "The quick approximation compares the population-averaged proportions with each cluster weighted by "
             "size / (p(1 − p)(1 + (size − 1)ρ)), averaged over the size distribution."};
  p.references = {"Li, P. & Redden, D. T. (2015). Comparing denominator degrees of freedom approximations for the "
                  "generalized linear mixed model in analyzing binary outcome in small sample cluster-randomized "
                  "trials. BMC Medical Research Methodology 15, 38.",
                  "Eldridge, S. M., Ashby, D. & Kerry, S. (2006). Sample size for cluster randomized trials: effect "
                  "of coefficient of variation of cluster size and analysis method. International Journal of "
                  "Epidemiology 35, 1292–1300.",
                  "Pinheiro, J. C. & Bates, D. M. (1995). Approximations to the log-likelihood function in the "
                  "nonlinear mixed-effects model. Journal of Computational and Graphical Statistics 4, 12–35.",
                  simc::REF_SIM};
  p.curve_param = "k1";
  p.curve_label = "k1 (Clusters in Group 1)";
  p.example_solve = "k";
  p.example = {{"power", "0.8"}, {"m", "30"}, {"cv", "0.4"}, {"p1", "0.3"}, {"p2", "0.2"}, {"icc", "0.05"},
               {"reps", "300"}};
  return p;
}
SSC_REGISTER(glmm_cluster);

}  // namespace
}  // namespace ssc::procs
