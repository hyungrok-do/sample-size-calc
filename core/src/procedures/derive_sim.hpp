// Shared steps of the simulation procedures' derivations: the studies that rejected H0, the power as their share, its
// Monte Carlo standard error and Wilson interval, and how a sample size was found by simulation.  They read the
// row's own values: a derivation runs for up to 1,000 rows and never simulates again.
#pragma once

#include <cmath>
#include <string>

#include "procedures/derive.hpp"
#include "procedures/sim_common.hpp"

namespace ssc::procs::derive {

/// Whether the row holds the analytic approximation of a simulation procedure rather than its simulated power.
inline bool approximated(const Row& r) { return text_of(r, "method") == "approximation"; }

/// k, the studies counted in the simulated power p = k / R, which is p·R to the nearest whole number.
inline double sim_count(const Row& r) { return std::round(r.num("power") * r.num("reps")); }

/// "two-sided, α = 0.05": the alternative and the level of each simulated study's test.
inline std::string sim_level(Alt a, double alpha) {
  return std::string(a == Alt::Two ? "two-sided" : "one-sided") + ", α = " + fnum(alpha);
}

/// The test of each study of a mixed model with the degrees of freedom of simc::df_param: the model's `fit`, then
/// t = `estimate` / SE with Satterthwaite's ν (the row's mean over the studies, and the studies whose own could not
/// be computed and took the fixed ones instead), with the fixed ν (`fixed`, such as "k1 + k2 − 2", which is
/// `fixed_numbers`), or z without degrees of freedom.  `studies` names the simulated studies ("trials").
inline Step sim_mixed_test(const Scenario& s, const Row& r, const std::string& fit, const std::string& estimate,
                           const std::string& fixed, const std::string& fixed_numbers, double fixed_df,
                           const std::string& studies = "studies") {
  const int dfm = simc::df_method(s);
  const std::string ratio = estimate + " / SE(" + estimate + ")";
  const std::string test = dfm == 2   ? "z = " + ratio + ", normal"
                           : dfm == 1 ? "t = " + ratio + " with ν = " + fixed
                                      : "t = " + ratio + " with Satterthwaite's ν from the study's fit (" + fixed +
                                            " where it cannot be computed)";
  std::string numbers = sim_level(blocks::alt(s), s("alpha"));
  if (dfm == 1) numbers += "; ν = " + fixed_numbers + " = " + fint(fixed_df);
  if (dfm == 0 && r.has("df_fallbacks"))
    numbers += "; " + simc::grouped(r.num("df_fallbacks")) + " " + studies + " took ν = " + fint(fixed_df);
  if (dfm == 0 && r.has("mean_df")) numbers += "; mean ν over the " + studies + " = " + num(r.num("mean_df"), 2);
  const std::string one = studies == "studies" ? "study" : studies.substr(0, studies.size() - 1);
  return {"Analysis of each " + one, fit + "; " + test + "; H0 rejected when p < α", numbers};
}

/// The studies counted: k of the R simulated with the seed (`counted` says which: "the studies whose test rejects
/// H0"), with those whose analysis failed, which count as not rejecting, and why.  `count` replaces "k = 812" where k
/// is not a plain count.
inline Step sim_count_step(const Row& r, const std::string& what, const std::string& counted,
                           const std::string& count = "") {
  std::string numbers = "R = " + simc::grouped(r.num("reps")) + " (seed " + fint(r.num("seed")) + ")";
  const double failed = r.num("failures");
  if (failed > 0) {
    const std::string why = text_of(r, "failures_why");
    numbers += ", " + simc::grouped(failed) + " not analyzed" + (why.empty() ? "" : " (" + why + ")") +
               " and counted as not rejecting";
  }
  return {what, "k = " + counted + " among the R simulated, each from its own random streams of the seed",
          numbers + ": " + (count.empty() ? "k = " + simc::grouped(sim_count(r)) : count)};
}

/// The simulated power p = k / R, its Monte Carlo standard error and its 95% Wilson interval, as
/// sim::with_rejections computes them (`target`: the power the row's sample size was to reach).
inline Steps sim_power_steps(const Row& r, double target = NAN) {
  const double p = r.num("power"), reps = r.num("reps");
  const double z = 1.959963984540054, z2 = z * z;  // as sim::with_rejections
  const double shift = z2 / (2 * reps), spread = std::sqrt(p * (1 - p) / reps + z2 / (4 * reps * reps));
  Steps out;
  out.push_back({"Simulated power", "p = k / R",
                 simc::grouped(sim_count(r)) + " / " + simc::grouped(reps) + " = " + prob(p, target)});
  out.push_back({"Monte Carlo standard error", "SE = √(p(1 − p) / R)",
                 "√(" + prob(p) + " × " + prob(1 - p) + " / " + simc::grouped(reps) + ") = " +
                     num(r.num("power_se"))});
  out.push_back({"95% interval (Wilson)", "(p + z²/2R ± z √(p(1 − p)/R + z²/4R²)) / (1 + z²/R), z = Φ⁻¹(0.975)",
                 "(" + prob(p) + " + " + num(shift) + " ± " + num(z, 2) + " × " + num(spread) + ") / " +
                     num(1 + z2 / reps) + " = " + prob(r.num("power_lo")) + " to " + prob(r.num("power_hi"))});
  return out;
}

/// How a sample size `sym` = n was found by simulation (sim::search_size): from the approximation's answer `start`
/// (with its approximate power), simulated with fewer studies a size (simc::explore_reps) and then with all R, until
/// the answer reaches the target and one size less does not, both with all R studies of the same seed; or the
/// smallest size allowed (`lo`), which already reaches it.  `studies` names them ("trials").
inline Step sim_size_search(const Row& r, const std::string& what, const std::string& sym, double n, double lo,
                            double start, double start_power, double target, const std::string& studies = "studies") {
  const double p = r.num("power"), reps = r.num("reps");
  if (n - 1 < lo)
    return {what, sym + " = the smallest size whose simulated power reaches the target",
            sym + " = " + fint(n) + " is the smallest size allowed, and its simulated power is already " +
                prob(p, target) + " ≥ " + fnum(target)};
  const double explore = static_cast<double>(simc::explore_reps(static_cast<long>(reps)));
  const std::string counts = explore < reps ? "R′ = " + simc::grouped(explore) + " " + studies +
                                                  " a size, then all R = " + simc::grouped(reps)
                                            : "R = " + simc::grouped(reps) + " " + studies + " a size";
  return {what,
          sym + " = the size at which the simulated power with all R " + studies +
              " reaches the target and one size less does not, searched from the approximation's answer " + sym +
              "₀ with R′ = min(R, max(100, ⌊R/4⌋)) " + studies + " a size first",
          sym + "₀ = " + fint(start) + " (approximate power " + prob(start_power, target) + "); " + counts +
              ": power at " + sym + " = " + fint(n - 1) + " < " + fnum(target) + " ≤ power at " + sym + " = " +
              fint(n) + ": " + prob(p, target)};
}

}  // namespace ssc::procs::derive
