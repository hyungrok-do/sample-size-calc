// Shared parts of the simulation procedures: inputs (number of simulated
// studies, seed, degrees of freedom), test decisions, result cells, columns
// and wording.
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "procedures/common.hpp"
#include "sim/random.hpp"
#include "sim/simulate.hpp"

namespace ssc::procs::simc {

using blocks::Alt;

inline const char* G_SIM = "Simulation";
inline const char* G_MODEL = "Model";

/// The most units (subjects, measurements, or clusters where a cluster is simulated through its sums) that one
/// simulated study may have.  Stop takes effect between studies, so one study must stay short: a few tenths of
/// a second at this size.
constexpr double STUDY_MAX = 1e6;

/// A whole number with its digits grouped by commas (1,000,000).
inline std::string grouped(double x) {
  std::string s = blocks::fint(x);
  for (int i = static_cast<int>(s.size()) - 3; i > (s[0] == '-' ? 1 : 0); i -= 3)
    s.insert(static_cast<std::size_t>(i), ",");
  return s;
}

/// Refuses to simulate studies of more than `max` `units` each.  Check a size before converting it to `long`,
/// which has only 32 bits on Windows and in WebAssembly.
inline void check_study_size(double size, const char* units, double max = STUDY_MAX) {
  if (size > max)
    throw CalcError("a simulated study is limited to " + grouped(max) + " " + units + " (this one has " +
                    grouped(size) + "); the approximation shown before simulating applies to larger designs");
}

inline ParamDef reps_param(const std::string& def = "1000") {
  return num("reps", "Simulated Studies", def)
      .grp(G_SIM)
      .sym("Simulations")
      .whole()
      .range(100, 1e6)
      .single()
      .sugg({"1000", "2000", "5000"})
      .hint("Number of studies to simulate. The Monte Carlo standard error of a power estimate p is "
            "√(p(1 − p)/R): about 1.3 percentage points for 1,000 studies at 80% power.");
}

inline ParamDef seed_param() {
  return num("seed", "Random Seed", "2026")
      .grp(G_SIM)
      .sym("Seed")
      .whole()
      .range(0, 1e9)
      .single()
      .hint("The same seed gives the same simulated studies, so every result can be reproduced exactly.");
}

/// Degrees of freedom of the t test: Satterthwaite, a fixed design-based value, or the normal limit.
inline ParamDef df_param(const std::string& fixed_label) {
  return select("df", "Degrees of Freedom",
                {{"satterthwaite", "Satterthwaite (as in lmerTest and SAS)", {}},
                 {"fixed", fixed_label, {}},
                 {"z", "None (normal approximation)", {}}},
                "satterthwaite")
      .grp(blocks::G_TEST)
      .hint("Satterthwaite's approximation estimates the degrees of freedom of each simulated study's t test; "
            "the normal approximation is anticonservative with few subjects or clusters.");
}

inline int df_method(const Scenario& s) { return s.is("df", "z") ? 2 : (s.is("df", "fixed") ? 1 : 0); }

/// p-value of a t (df finite) or z statistic for the alternative.
inline double p_value(double t, double df, Alt a) {
  using namespace stats;
  if (!std::isfinite(t)) return 1;
  if (std::isinf(df)) {
    if (a == Alt::Two) return 2 * pnorm(-std::fabs(t));
    return a == Alt::Upper ? pnorm(t, false) : pnorm(t, true);
  }
  if (a == Alt::Two) return 2 * pt(std::fabs(t), df, false);
  return a == Alt::Upper ? pt(t, df, false) : pt(t, df, true);
}

/// t statistic estimate / se (0 without a positive se); one too large to represent is still a rejection, so it is
/// kept at the largest finite value rather than infinity, which p_value treats as no statistic.
inline double t_stat(double estimate, double se) {
  if (!(se > 0)) return 0;
  const double t = estimate / se;
  return std::isinf(t) ? std::copysign(std::numeric_limits<double>::max(), t) : t;
}

/// p-value of the Mann–Whitney test with the normal approximation, the tie correction and a continuity
/// correction (as R's wilcox.test with exact = FALSE); `w` the rank sum of group 1 minus n1(n1 + 1)/2.
inline double rank_sum_p(double w, double n1, double n2, double ties, Alt a) {
  const double n = n1 + n2, mean = n1 * n2 / 2;
  const double var = n1 * n2 / 12 * ((n + 1) - ties / (n * (n - 1)));
  if (!(var > 0)) return 1;
  const double d = w - mean;
  const double corr = a == Alt::Two ? (d > 0 ? 0.5 : (d < 0 ? -0.5 : 0)) : (a == Alt::Upper ? 0.5 : -0.5);
  const double z = (d - corr) / std::sqrt(var);
  if (a == Alt::Two) return 2 * std::fmin(stats::pnorm(z), stats::pnorm(z, false));
  return a == Alt::Upper ? stats::pnorm(z, false) : stats::pnorm(z);
}

/// Replicates of the exploratory evaluations of a sample-size search.
inline long explore_reps(long reps) { return std::min(reps, std::max(100L, reps / 4)); }

/// Cells of a simulated power estimate.  Studies whose analysis failed count as not rejecting, which is a fair
/// account of a few of them; when they are most of the studies, the share that rejected says more about the
/// analysis than about the design, and no power is given.
inline void put_summary(const Scenario& s, Row& r, const sim::Summary& sm) {
  if (2 * sm.failures > sm.reps)
    throw CalcError("the simulated power cannot be computed: " + grouped(static_cast<double>(sm.failures)) +
                    " of the " + grouped(static_cast<double>(sm.reps)) +
                    " simulated studies could not be analyzed; check for an extreme value");
  if (s.solve != "power") r.set("target_power", s("power"));
  r.set("power", sm.power);
  r.set("power_se", sm.se);
  r.set("power_lo", sm.lo);
  r.set("power_hi", sm.hi);
  r.set("reps", static_cast<double>(sm.reps));
  r.set("seed", s("seed"));
  r.set("failures", static_cast<double>(sm.failures));
  if (std::isfinite(sm.mean_df)) r.set("mean_df", sm.mean_df);
  r.set("method", "simulation");
}

/// Cells of an analytic approximation (the quick preview of a simulation procedure).
inline void put_approximation(const Scenario& s, Row& r, double power) {
  if (s.solve != "power") r.set("target_power", s("power"));
  r.set("power", power);
  r.set("method", "approximation");
}

/// Power columns of a simulation procedure (the Monte Carlo columns stay empty, and are
/// dropped, in an approximation).
inline std::vector<Column> power_columns(const Scenario& s) {
  using blocks::col;
  std::vector<Column> c;
  if (s.solve != "power") c.push_back(col("target_power", "Target Power", "prob"));
  c.push_back(col("power", s.solve != "power" ? "Actual Power" : "Power", "prob",
                  "Simulated power: the share of simulated studies whose test rejected the null hypothesis"));
  c.push_back(col("power_lo", "Lower 95%", "prob", "Lower end of the 95% interval for the power (Wilson)"));
  c.push_back(col("power_hi", "Upper 95%", "prob", "Upper end of the 95% interval for the power (Wilson)"));
  c.push_back(col("power_se", "MC SE", "prob", "Monte Carlo standard error of the power estimate"));
  c.push_back(col("approx_power", "Approx. Power", "prob", "Analytic approximation of the power"));
  return c;
}

inline std::vector<Column> run_columns() {
  using blocks::col;
  return {col("mean_df", "Mean df", "num", "Average degrees of freedom of the simulated studies' tests"),
          col("reps", "Studies", "int", "Number of simulated studies"),
          col("failures", "Failed Fits", "int", "Simulated studies whose model could not be fitted (not rejecting)"),
          col("seed", "Seed", "int", "Random seed")};
}

/// "achieve a simulated power of 81.2% (95% CI 78.7% to 83.5%, 1,000 simulated studies; target 80.0%)".
inline std::string power_phrase(const Scenario& s, const Row& r) {
  using blocks::fint;
  using blocks::fpct;
  const bool approx = text_of(r, "method") == "approximation";
  std::string t = approx ? "achieve an approximate power of " + fpct(r.num("power"))
                         : "achieve a simulated power of " + fpct(r.num("power")) + " (95% CI " +
                               fpct(r.num("power_lo")) + " to " + fpct(r.num("power_hi")) + ", " +
                               fint(r.num("reps")) + " simulated studies";
  if (s.solve != "power") t += std::string(approx ? " (" : "; ") + "target " + fpct(s("power")) + ")";
  else if (!approx) t += ")";
  return t;
}

/// Closing sentences: what the approximation is, or how the simulation compares with it.
inline std::string method_note(const Row& r) {
  using blocks::fpct;
  if (text_of(r, "method") == "approximation")
    return " This is a quick analytic approximation; run the simulation to confirm it.";
  std::string t;
  if (r.has("approx_power") && std::isfinite(r.num("approx_power")))
    t += " The analytic approximation gives " + fpct(r.num("approx_power")) + ".";
  const double failed = r.num("failures");
  if (failed > 0)
    t += " " + blocks::fint(failed) + " simulated " + (failed == 1 ? "study" : "studies") +
         " could not be analyzed and count as not rejecting.";
  return t;
}

/// Cluster sizes: all equal (CV 0) or gamma-distributed with mean m̄ and coefficient of variation CV,
/// rounded to whole subjects (at least 1).
struct ClusterSizes {
  double mean = 20, cv = 0, shape = 0, scale = 0;
  ClusterSizes(double m, double c) : mean(m), cv(c) {
    if (cv > 0) {
      shape = 1 / (cv * cv);
      scale = mean / shape;
    }
  }
  /// A size drawn from the distribution (capped at 10⁹, which keeps it within a 32-bit long).
  long draw(sim::Rng& rng) const {
    if (cv <= 0) return std::max(1L, std::lround(mean));
    return std::max(1L, std::lround(std::fmin(scale * rng.gamma(shape), 1e9)));
  }
  /// Equally likely sizes spanning the distribution (512 quantiles), for expected values over it.
  std::vector<double> grid() const {
    if (cv <= 0) return {std::max(1.0, std::round(mean))};
    std::vector<double> g(512);
    for (std::size_t i = 0; i < g.size(); ++i)
      g[i] = std::max(1.0, std::round(scale * stats::qgamma((i + 0.5) / static_cast<double>(g.size()), shape)));
    return g;
  }
};

inline ParamDef cluster_size_param() {
  return num("m", "m̄ (Average Cluster Size)", "20")
      .grp(blocks::G_N)
      .sym("m̄")
      .range(1, 1e7)
      .sugg({"20", "10 20 50"})
      .hint("Average number of subjects per cluster (rounded to whole subjects when the sizes are equal).");
}

inline ParamDef cluster_cv_param() {
  return num("cv", "CV of the Cluster Sizes", "0")
      .grp(blocks::G_N)
      .sym("CV")
      .range(0, 5)
      .sugg({"0", "0.4", "0 0.4 0.8"})
      .hint("Coefficient of variation (SD/mean) of the cluster sizes; 0 when every cluster has the same size. "
            "Sizes are drawn from a gamma distribution with this mean and CV, rounded to whole subjects (at "
            "least 1).");
}

inline const char* REF_SIM =
    "Landau, S. & Stahl, D. (2013). Sample size and power calculations for medical studies by simulation when "
    "closed form expressions are not available. Statistical Methods in Medical Research 22, 324–345.";
inline const char* REF_ARNOLD =
    "Arnold, B. F., Hogan, D. R., Colford, J. M. & Hubbard, A. E. (2011). Simulation methods to estimate design "
    "power: an overview for applied research. BMC Medical Research Methodology 11, 94.";

}  // namespace ssc::procs::simc
