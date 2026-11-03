// Shared steps and frames of the simulation procedures' derivations.  Solving for a size, a simulated row finds the
// smallest size whose simulated power reaches the target: its steps show how sim::search_size went from the
// approximation's answer, the studies that rejected H0 at the answer, the power as their share with its Monte Carlo
// standard error and Wilson interval, and what the search established one size below.  Solving for power, they show
// the studies counted at the sizes given and the power.  An approximation is framed like an analytic procedure.
// They read the row's own values: a derivation runs for up to 1,000 rows and never simulates again.
#pragma once

#include <cctype>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "procedures/derive.hpp"
#include "procedures/sim_common.hpp"

namespace ssc::procs::derive {

/// Whether the row holds the analytic approximation of a simulation procedure rather than its simulated power.
inline bool approximated(const Row& r) { return text_of(r, "method") == "approximation"; }

/// k, the studies counted in the simulated power p = k / R, which is p·R to the nearest whole number.
inline double sim_count(const Row& r) { return std::round(r.num("power") * r.num("reps")); }

/// A count with its digits grouped (1,000), given or derived.
inline std::string gvg(double x) { return given(simc::grouped(x)); }
inline std::string dvg(double x) { return derived(simc::grouped(x)); }

/// Numbers as a list ("0.2, 0.25, 0.3"), each given or derived.
inline std::string gv_list(const std::vector<double>& v) {
  std::string t;
  for (double x : v) t += (t.empty() ? "" : ", ") + gv(x);
  return t;
}
inline std::string dv_list(const std::vector<double>& v) {
  std::string t;
  for (double x : v) t += (t.empty() ? "" : ", ") + dv(x);
  return t;
}

/// "two-sided, α = 0.05": the alternative and the level of each simulated study's test.
inline std::string sim_level(Alt a, double alpha) {
  return std::string(a == Alt::Two ? "two-sided" : "one-sided") + ", α = " + gv(alpha);
}

/// Two groups of sizes N1 and N2 = ⌈r N1⌉ as substituted: "N1 = 30, N2 = ⌈1 × 30⌉ = 30", N1 given or derived.
inline std::string ratio_sizes(const std::string& s1, const std::string& s2, double n1, double ratio, double n2,
                               bool solved) {
  return s1 + " = " + size(n1, solved) + ", " + s2 + " = ⌈" + gv(ratio) + " × " + size(n1, solved) + "⌉ = " +
         dvi(n2);
}

// ------------------------------------------------------------------ what is found, and the answer

/// Solving for a size by simulation: the smallest `unknown` whose simulated power (`which`: " to find …", where the
/// power has several meanings) reaches the target.
inline std::string find_sim_size(const std::string& unknown, double target, const std::string& which = "") {
  return "the smallest " + unknown + " whose simulated power" + which + " is at least 1 − β = " + gv(target);
}
/// Solving for power by simulation: the simulated power `at` the sizes given.
inline std::string find_sim_power(const std::string& at, const std::string& which = "") {
  return "the simulated power" + which + " " + at;
}
/// The same, by the analytic approximation shown before simulating.
inline std::string find_approx_size(const std::string& unknown, double target, const std::string& which = "") {
  return "the smallest " + unknown + " whose power" + which + " is at least 1 − β = " + gv(target) +
         ", by the analytic approximation (nothing simulated)";
}
inline std::string find_approx_power(const std::string& at, const std::string& which = "") {
  return "the power" + which + " " + at + ", by the analytic approximation (nothing simulated)";
}

/// The answer of a simulated row: the sizes found (`sizes`, solving for one) with the simulated power and its
/// interval, or the simulated power alone.
inline std::string sim_answer(const Row& r, const std::string& sizes, bool solved) {
  const std::string p = dvp(r.num("power")) + " (95% CI " + dvp(r.num("power_lo")) + " to " +
                        dvp(r.num("power_hi")) + ")";
  return solved ? sizes + ", with simulated power " + p : "simulated power = " + p;
}
/// The answer of an approximation.
inline std::string approx_answer(const Row& r, const std::string& sizes, bool solved) {
  const std::string p = dvp(r.num("power"));
  return solved ? sizes + ", with approximate power " + p : "approximate power = " + p;
}

// ------------------------------------------------------------------ the steps of a simulated row

/// The test of each study of a mixed model with the degrees of freedom of simc::df_param: the model's `fit`, then
/// t = `estimate` / SE with Satterthwaite's ν (the row's mean over the studies, and the studies whose own could not
/// be computed and took the fixed ones instead), with the fixed ν (`fixed`, such as "k1 + k2 − 2", which is
/// `fixed_numbers` with the values marked), or z without degrees of freedom.  `studies` names the simulated studies
/// ("trials").
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
  if (dfm == 1) numbers += "; ν = " + fixed_numbers + " = " + dvi(fixed_df);
  if (dfm == 0 && r.has("df_fallbacks"))
    numbers += "; " + dvg(r.num("df_fallbacks")) + " " + studies + " took ν = " + fixed_numbers + " = " +
               dvi(fixed_df);
  if (dfm == 0 && r.has("mean_df")) numbers += "; mean ν over the " + studies + " = " + dv(r.num("mean_df"), 2);
  const std::string one = studies == "studies" ? "study" : studies.substr(0, studies.size() - 1);
  return {"Analysis of each " + one, fit + "; " + test + "; H0 rejected when p < α", numbers};
}

/// The reasons of the failed studies (simc::causes_text, "in 150, …; in 13, …"), their counts marked as derived.
inline std::string causes_marked(const std::string& why) {
  std::string out;
  std::size_t i = 0;
  while (i < why.size()) {
    const bool lead = why.compare(i, 3, "in ") == 0 && (i == 0 || (i >= 2 && why.compare(i - 2, 2, "; ") == 0));
    std::size_t j = i + 3;  // a grouped count ends at the comma before a space
    while (lead && j < why.size() && (std::isdigit(static_cast<unsigned char>(why[j])) ||
                                      (why[j] == ',' && j + 1 < why.size() && why[j + 1] != ' ')))
      ++j;
    if (lead && j > i + 3) {
      out += "in " + derived(why.substr(i + 3, j - i - 3));
      i = j;
    } else {
      out += why[i++];
    }
  }
  return out;
}

/// The studies counted: k of the R simulated with the seed (`counted` says which: "the studies whose test rejects
/// H0"), `at` the sizes where they were counted ("N1 = 17, N2 = 17"; empty for the sizes given), with those whose
/// analysis failed, which count as not rejecting, and why.  `count` replaces "k = 812" where k is not a plain count.
inline Step sim_count_step(const Row& r, const std::string& what, const std::string& counted, const std::string& at,
                           const std::string& count = "") {
  std::string numbers = (at.empty() ? "" : "at " + at + ", ") + "R = " + gvg(r.num("reps")) + " (seed " +
                        gvi(r.num("seed")) + ")";
  const double failed = r.num("failures");
  if (failed > 0) {
    const std::string why = causes_marked(text_of(r, "failures_why"));
    numbers += ", " + dvg(failed) + " not analyzed" + (why.empty() ? "" : " (" + why + ")") +
               " and counted as not rejecting";
  }
  return {what, "k = " + counted + " among the R simulated, each from its own random streams of the seed",
          numbers + ": " + (count.empty() ? "k = " + dvg(sim_count(r)) : count)};
}

/// The simulated power p = k / R, its Monte Carlo standard error and its 95% Wilson interval, as
/// sim::with_rejections computes them (`target`: the power the row's sample size was to reach).
inline Steps sim_power_steps(const Row& r, double target = NAN) {
  const double p = r.num("power"), reps = r.num("reps");
  const double z = 1.959963984540054, z2 = z * z;  // as sim::with_rejections
  const double shift = z2 / (2 * reps), spread = std::sqrt(p * (1 - p) / reps + z2 / (4 * reps * reps));
  Steps out;
  out.push_back({"Simulated power", "p = k / R", dvg(sim_count(r)) + " / " + gvg(reps) + " = " + dvp(p, target)});
  out.push_back({"Monte Carlo standard error", "SE = √(p(1 − p) / R)",
                 "√(" + dvp(p) + " × " + dvp(1 - p) + " / " + gvg(reps) + ") = " + dv(r.num("power_se"))});
  out.push_back({"95% interval (Wilson)", "(p + z²/2R ± z √(p(1 − p)/R + z²/4R²)) / (1 + z²/R), z = Φ⁻¹(0.975)",
                 "(" + dvp(p) + " + " + dv(shift) + " ± " + num(z, 2) + " × " + dv(spread) + ") / " +
                     dv(1 + z2 / reps) + " = " + dvp(r.num("power_lo")) + " to " + dvp(r.num("power_hi"))});
  return out;
}

/// A size solved for by simulation (sim::search_size), as its derivation shows it.
struct SimSize {
  std::string sym;           // the size's symbol: "N1", "k1", "g", "n"
  double n = NAN;            // the size found
  double lo = NAN;           // the smallest size allowed
  double start = NAN;        // the approximation's smallest size reaching the target, where the search began
  double start_power = NAN;  // its approximate power
  double target = NAN;       // 1 − β
  /// The design at a size, as substituted ("N1 = 16, N2 = 16"); the size alone when unset.
  std::function<std::string(double)> design;
  std::string at(double x) const { return design ? design(x) : sym + " = " + dvi(x); }
};

/// How sim::search_size went: from the approximation's answer, sizes simulated with R′ studies each, bracketing the
/// target and halving the bracket; then with all R, bracketing again from there and halving until a size reaches
/// the target and the one below does not.  Every size is simulated with the same seed's streams.
inline Steps sim_search_steps(const Row& r, const SimSize& z, const std::string& studies = "studies") {
  const double reps = r.num("reps");
  const double explore = static_cast<double>(simc::explore_reps(static_cast<long>(reps)));
  const std::string s0 = z.sym + "₀";
  Steps out;
  out.push_back({"Search from the approximation",
                 s0 + " = the approximation's smallest " + z.sym + " whose power reaches 1 − β; sizes simulated "
                     "with R′ = min(R, max(100, ⌊R/4⌋)) " + studies + " each, from " + s0 + " down in steps of 10% "
                     "of it (up in steps of 25% if it falls short), each step twice the last, until the target is "
                     "bracketed, then the bracket halved until it spans at most 2% of its upper end",
                 s0 + " = " + dvi(z.start) + ", approximate power " + dvp(z.start_power, z.target) +
                     "; R′ = min(" + gvg(reps) + ", max(100, ⌊" + gvg(reps) + "/4⌋)) = " + dvg(explore)});
  out.push_back({"Confirmation with all R " + studies,
                 "from the upper end of that bracket, sizes simulated with all R " + studies + " each: down in steps "
                     "of 2% of it (up if it falls short), each step twice the last, until the target is bracketed, "
                     "then the bracket halved until a size reaches 1 − β and the one below does not; every size with "
                     "the same random streams of the seed",
                 "R = " + gvg(reps) + " " + studies + " a size, seed " + gvi(r.num("seed"))});
  return out;
}

/// What the search established below the answer: with all R studies the size one less falls short of the target
/// (its power is not kept, and is not simulated again), or the answer is the smallest size allowed.
inline Step sim_below_step(const Row& r, const SimSize& z, const std::string& studies = "studies") {
  const double p = r.num("power");
  if (z.n - 1 < z.lo)
    return {"Smallest size", z.sym + " = the smallest size whose simulated power reaches 1 − β",
            z.sym + " = " + dvi(z.n) + " is the smallest size allowed, and its simulated power is already " +
                dvp(p, z.target) + " ≥ " + gv(z.target)};
  return {"Smallest size",
          z.sym + " = the size whose simulated power with all R " + studies + " reaches 1 − β while that of " +
              z.sym + " − 1 does not",
          "power at " + z.at(z.n - 1) + " < " + gv(z.target) + " ≤ power at " + z.at(z.n) + ": " + dvp(p, z.target) +
              ", both with R = " + gvg(r.num("reps")) + " " + studies + " of seed " + gvi(r.num("seed"))};
}

/// The steps of a simulated row after its model and the analysis of each study: solving for a size (`size` set),
/// how the search went, the studies counted at the answer (`count`), the power with its standard error and
/// interval, and what the search established below the answer; solving for power, the studies counted, the power,
/// its standard error and interval.
inline void sim_result_steps(Steps& out, const Row& r, const Step& count, const SimSize* size,
                             const std::string& studies = "studies") {
  if (size)
    for (auto& st : sim_search_steps(r, *size, studies)) out.push_back(st);
  out.push_back(count);
  for (auto& st : sim_power_steps(r, size ? size->target : NAN)) out.push_back(st);
  if (size) out.push_back(sim_below_step(r, *size, studies));
}

}  // namespace ssc::procs::derive
