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
/// a second at this size.  (The cluster mixed models check for Stop within a study too, which there can take
/// minutes: 10⁶ clusters of 10⁷ subjects, half of them with an outcome as rare as 1e-9, whose binomial counts by
/// inversion take some 0.3 ms each.)
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

/// What a simulated study whose analysis failed says about the design.
enum class Failures {
  Design,     // the design gives that study no answer (no wins or no losses, every cluster of a group with the same
              // outcome, too few patients left at a visit), so that not rejecting is a fair account of it
  Numerical,  // the arithmetic broke down (a likelihood past the largest number), which says nothing of the design
};

/// One reason why a simulated study's analysis fails: the bit of sim::Outcome::also that the study sets for it (a
/// failed analysis sets the bit of one reason), what it says about the design, and its words, a clause that follows
/// "could not be analyzed (…)" or "in 150, …" ("a group had no patients left at a visit").
struct Cause {
  int bit;
  Failures kind;
  std::string what;
};

/// sim::simulate, for studies that run other tests on the same data, or tell why an analysis failed, in the bits of
/// sim::Outcome::also: the name says what such studies rely on, that those bits count on every study, also where the
/// analysis that decides the power failed, as sim::Accumulator counts them.
template <class F>
sim::Summary simulate_tests(long reps, F&& study, const RunContext& ctx, double p0 = 0, double p1 = 1) {
  return sim::simulate(reps, study, ctx, p0, p1);
}

/// The reasons of the failed studies: "a group had no patients left at a visit", or "in 150, …; in 13, …" when
/// there were several.
inline std::string causes_text(const sim::Summary& sm, const std::vector<Cause>& causes) {
  std::vector<const Cause*> seen;
  for (const Cause& c : causes)
    if (sm.also[c.bit] > 0) seen.push_back(&c);
  if (seen.size() == 1) return seen[0]->what;
  std::string t;
  for (const Cause* c : seen)
    t += (t.empty() ? "in " : "; in ") + grouped(static_cast<double>(sm.also[c->bit])) + ", " + c->what;
  return t;
}

/// Studies whose analysis broke down in its arithmetic, of the reasons `causes`.
inline long numerical_failures(const sim::Summary& sm, const std::vector<Cause>& causes) {
  long n = 0;
  for (const Cause& c : causes)
    if (c.kind == Failures::Numerical) n += sm.also[c.bit];
  return n;
}

/// Studies whose analysis failed, of the reasons `causes` (a failed analysis sets the bit of one reason).
inline long failures_of(const sim::Summary& sm, const std::vector<Cause>& causes) {
  long n = 0;
  for (const Cause& c : causes) n += sm.also[c.bit];
  return n;
}

/// Simulated power of another test run on the same studies (bit `bit` of sim::Outcome::also), whose analysis fails
/// for the reasons `causes` (such studies count as not rejecting); NaN when it failed in every study, or when the
/// arithmetic broke down in more than half of them, as put_summary refuses.
inline double other_power(const sim::Summary& sm, int bit, const std::vector<Cause>& causes) {
  const bool none = sm.reps > 0 && failures_of(sm, causes) == sm.reps;
  return none || 2 * numerical_failures(sm, causes) > sm.reps ? NAN : sm.also_power(bit);
}

/// Cells of a test run on the same studies: its power under `key` (other_power) and, where its analysis failed,
/// the number of those studies and their reasons under key + "_failures" and key + "_why", which others_sentence
/// reads.  The test that decides the power (`deciding`) leaves those to put_summary's failures and failures_why.
inline void put_other(Row& r, const sim::Summary& sm, const std::string& key, int bit,
                      const std::vector<Cause>& causes, bool deciding) {
  r.set(key, other_power(sm, bit, causes));
  const long failed = failures_of(sm, causes);
  if (deciding || failed == 0) return;
  r.set(key + "_failures", static_cast<double>(failed));
  r.set(key + "_why", causes_text(sm, causes));
}

/// Cells of a simulated power estimate.  Studies whose analysis failed count as not rejecting.  No power is given
/// when none could be analyzed, or when failures are numerical breakdowns and more than half of the studies: the
/// share that rejected would then say more about the arithmetic than about the design.  `why` says what makes a
/// study fail, for the message when none could be analyzed.
inline void put_summary(const Scenario& s, Row& r, const sim::Summary& sm, Failures kind = Failures::Design,
                        const std::string& why = "") {
  const std::string extreme = kind == Failures::Numerical ? "; check for an extreme value" : "";
  if (sm.reps > 0 && sm.failures == sm.reps)
    throw CalcError("the simulated power cannot be computed: none of the " +
                    grouped(static_cast<double>(sm.reps)) + " simulated studies could be analyzed" +
                    (why.empty() ? "" : " (" + why + ")") + extreme);
  if (kind == Failures::Numerical && 2 * sm.failures > sm.reps)
    throw CalcError("the simulated power cannot be computed: " + grouped(static_cast<double>(sm.failures)) +
                    " of the " + grouped(static_cast<double>(sm.reps)) +
                    " simulated studies could not be analyzed" + extreme);
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

/// put_summary for studies that tell why their analysis failed (simulate_tests and the bits of `causes`): design
/// failures and numerical ones are counted apart, the rule above applies to the numerical ones, and the messages and
/// the statement say what made the studies fail.  `studies` names them in the messages, as the statement does
/// ("trials").
inline void put_summary(const Scenario& s, Row& r, const sim::Summary& sm, const std::vector<Cause>& causes,
                        const std::string& studies = "studies") {
  const bool numerical = 2 * numerical_failures(sm, causes) > sm.reps;
  const std::string why = causes_text(sm, causes), extreme = numerical ? "; check for an extreme value" : "";
  if (sm.reps > 0 && sm.failures == sm.reps)
    throw CalcError("the simulated power cannot be computed: none of the " +
                    grouped(static_cast<double>(sm.reps)) + " simulated " + studies + " could be analyzed (" + why +
                    ")" + extreme);
  if (numerical)
    throw CalcError("the simulated power cannot be computed: " + grouped(static_cast<double>(sm.failures)) +
                    " of the " + grouped(static_cast<double>(sm.reps)) + " simulated " + studies +
                    " could not be analyzed (" + why + ")" + extreme);
  put_summary(s, r, sm);
  if (sm.failures > 0) r.set("failures_why", why);
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
          col("failures", "Failed Fits", "int",
              "Simulated studies whose analysis failed, as when its model could not be fitted (not rejecting)"),
          col("seed", "Seed", "int", "Random seed")};
}

/// "achieve a simulated power of 81.2% (95% CI 78.7% to 83.5%, 1,000 simulated studies; target 80.0%)", with the
/// studies named as the rest of the statement names them ("trials").
inline std::string power_phrase(const Scenario& s, const Row& r, const std::string& studies = "studies") {
  using blocks::fint;
  using blocks::fpct;
  const bool approx = text_of(r, "method") == "approximation";
  std::string t = approx ? "achieve an approximate power of " + fpct(r.num("power"))
                         : "achieve a simulated power of " + fpct(r.num("power")) + " (95% CI " +
                               fpct(r.num("power_lo")) + " to " + fpct(r.num("power_hi")) + ", " +
                               fint(r.num("reps")) + " simulated " + studies;
  if (s.solve != "power") t += std::string(approx ? " (" : "; ") + "target " + fpct(s("power")) + ")";
  else if (!approx) t += ")";
  return t;
}

/// " On the same simulated trials the MMRM had a power of 80.5% and the completers' t-test 71.5%.": the powers of
/// the tests named in `tests` (a name and the key of its power each), run on the same `studies`, then, for each test
/// whose analysis failed in some of them (put_other), in how many and why, and why a power that is not given
/// (other_power) is not.  Empty in an approximation.
inline std::string others_sentence(const Row& r, const std::string& studies,
                                   const std::vector<std::pair<std::string, std::string>>& tests) {
  std::vector<std::string> given;
  std::string failed_text;
  const std::string all = grouped(r.num("reps")) + " simulated " + studies;
  for (const auto& t : tests) {
    if (!r.has(t.second)) return "";
    const double failed = r.has(t.second + "_failures") ? r.num(t.second + "_failures") : 0;
    const std::string why = " (" + text_of(r, t.second + "_why") + ")";
    if (std::isfinite(r.num(t.second))) {
      given.push_back(t.first + (given.empty() ? " had a power of " : " ") + blocks::fpct(r.num(t.second)));
      if (failed > 0)
        failed_text += " " + grouped(failed) + " of the " + all + " could not be analyzed by the " + t.first +
                       (failed == 1 ? " and counts" : " and count") + " as not rejecting" + why + ".";
    } else if (failed >= r.num("reps")) {
      failed_text += " None of the " + all + " could be analyzed by the " + t.first + why +
                     ", so its power is not given.";
    } else {
      failed_text += " " + grouped(failed) + " of the " + all + " could not be analyzed by the " + t.first + why +
                     ", and as its arithmetic broke down in more than half of them, its power is not given.";
    }
  }
  std::string s;
  for (std::size_t i = 0; i < given.size(); ++i)
    s += (i == 0 ? " On the same simulated " + studies + " the " : (i + 1 == given.size() ? " and the " : ", the ")) +
         given[i];
  if (!given.empty()) s += ".";
  return s + failed_text;
}

/// Closing sentences: what the approximation is, or how the simulation compares with it.  `studies` names the
/// simulated studies as others_sentence does ("trials"; the singular drops the s).
inline std::string method_note(const Row& r, const std::string& studies = "studies") {
  using blocks::fpct;
  if (text_of(r, "method") == "approximation")
    return " This is a quick analytic approximation; run the simulation to confirm it.";
  std::string t;
  if (r.has("approx_power") && std::isfinite(r.num("approx_power")))
    t += " The analytic approximation gives " + fpct(r.num("approx_power")) + ".";
  const double failed = r.num("failures");
  const std::string why = text_of(r, "failures_why");
  const std::string one = studies == "studies" ? "study" : studies.substr(0, studies.size() - 1);
  if (failed > 0)
    t += " " + grouped(failed) + " simulated " + (failed == 1 ? one : studies) + " could not be analyzed and " +
         (failed == 1 ? "counts" : "count") + " as not rejecting" + (why.empty() ? "" : " (" + why + ")") + ".";
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
  /// Whether every size is the mean: at a CV of 0, and at one whose square underflows to 0 (below about 1e-154),
  /// whose shape 1/CV² is infinite (every draw would be 10⁹); the sizes are the mean to the nearest subject well
  /// before that.
  bool equal() const { return cv <= 0 || std::isinf(shape); }
  /// A size drawn from the distribution (capped at 10⁹, which keeps it within a 32-bit long).
  long draw(sim::Rng& rng) const {
    if (equal()) return std::max(1L, std::lround(mean));
    return std::max(1L, std::lround(std::fmin(scale * rng.gamma(shape), 1e9)));
  }
  /// Equally likely sizes spanning the distribution (512 quantiles), for expected values over it.
  std::vector<double> grid() const {
    if (equal()) return {std::max(1.0, std::round(mean))};
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
