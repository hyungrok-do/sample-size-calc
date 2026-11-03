// Steps shared by the derivations of the agreement, diagnostic, nonparametric, confidence-interval and
// group-sequential procedures, in the manner of derive.hpp: searches whose power is saw-toothed in N, group-size
// searches that reach a target of another quantity than the power, and the closed forms shown for comparison.
#pragma once

#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "procedures/derive.hpp"

namespace ssc::procs::derive {

/// Derived numbers listed in a formula: "0.3333, 0.6667, 1", each marked.
inline std::string dvlist(const std::vector<double>& v, int decimals = 4) {
  std::string t;
  for (std::size_t i = 0; i < v.size(); ++i) t += (i ? ", " : "") + dv(v[i], decimals);
  return t;
}

/// Two group sizes as substituted: when solved for, derived but for a group fixed by the allocation; otherwise as
/// entered, or derived from the entries under the allocation rule (N2 = ⌈r × N1⌉, or a percentage of the total).
inline std::pair<std::string, std::string> group_sizes(const Scenario& s, const blocks::TwoN& n, bool solved) {
  const std::string& a = s.s("alloc");
  if (solved) return {a == "n1" ? gvi(n.n1) : dvi(n.n1), a == "n2" ? gvi(n.n2) : dvi(n.n2)};
  if (a == "pct") return {dvi(n.n1), dvi(n.n2)};
  return {gvi(n.n1), a == "ratio" ? dvi(n.n2) : gvi(n.n2)};
}
inline std::string group_sizes_text(const Scenario& s, const blocks::TwoN& n, bool solved) {
  const auto [n1, n2] = group_sizes(s, n, solved);
  return "N1 = " + n1 + ", N2 = " + n2;
}

/// The words that end a closed form shown beside the exact answer.
inline const std::string FOR_COMPARISON = ", for comparison; the steps below give the answer";

/// (z(1 − α/2) + z(1 − β)) as written in a closed form, and its numbers.
inline std::string z_sum_symbol(Alt a) { return std::string("(z(") + level_symbol(a) + ") + z(1 − β))"; }
inline std::string z_sum_numbers(double alpha, double target, Alt a) {
  return "(" + dv(stats::qnorm(tail_alpha(alpha, a), false)) + " + " + dv(stats::qnorm(target)) + ")";
}

/// How a single sample size was found where the power of an exact test is saw-toothed in N: the search's scan
/// upward from below its bracket returns the first size that reaches the target. Shown by the power one step below
/// and at the answer; when the size below reaches the target too, the answer is where that scan began.
inline Step saw_search(const std::function<double(double)>& power, double n, double target, double lo,
                       const std::string& sym = "N", double step = 1) {
  const std::string rule =
      sym + " = the smallest whole number whose power reaches 1 − β (exact power is not monotone in " + sym +
      ": an upward scan finds the first size that reaches it)";
  const double at = power(n);
  if (n - step < lo)
    return {"Smallest size", rule,
            sym + " = " + dvi(n) + " is the smallest size allowed, and its power is already " + dvp(at, target) +
                " ≥ " + gv(target)};
  const double below = power(n - step);
  if (below < target)
    return {"Smallest size", rule,
            "power at " + sym + " = " + dvi(n - step) + ": " + dvp(below, target) + " < " + gv(target) +
                " ≤ power at " + sym + " = " + dvi(n) + ": " + dvp(at, target)};
  return {"Smallest size", rule,
          "power at " + sym + " = " + dvi(n) + ": " + dvp(at, target) + " ≥ " + gv(target) +
              ", the first size of the scan (power at " + sym + " = " + dvi(n - step) + ", below it: " +
              dvp(below, target) + ")"};
}

/// How two group sizes were found under the scenario's allocation rule when the search reaches a target of another
/// quantity than the power (`quantity`, such as "drift", the target `target` written `target_sym`, itself
/// computed): the value one design below the answer and at it (as derive::two_group_search).
inline Step two_group_search_of(const Scenario& s, const std::function<double(double, double)>& f,
                                const blocks::TwoN& n, double target, double min_n, const std::string& quantity,
                                const std::string& target_sym) {
  const std::string& a = s.s("alloc");
  const std::string reach = "the smallest whole number whose " + quantity + " reaches " + target_sym;
  const std::string rule =
      a == "equal"   ? "N1 = N2, " + reach
      : a == "ratio" ? "N1 = " + reach + ", with N2 = ⌈r × N1⌉"
      : a == "n1"    ? "N2 = " + reach + ", with N1 fixed"
      : a == "n2"    ? "N1 = " + reach + ", with N2 fixed"
                     : "N = the smallest total whose " + quantity + " reaches " + target_sym +
                           ", a fixed percentage of it in group 1";
  const double at = f(n.n1, n.n2);
  const blocks::TwoN less = blocks::smaller_two_group(s, n);
  if (less.n1 < min_n || less.n2 < min_n)
    return {"Smallest sizes", rule,
            group_sizes_text(s, n, true) + " are the smallest sizes allowed, and their " + quantity +
                " is already " + dvp(at, target) + " ≥ " + dv(target)};
  return {"Smallest sizes", rule,
          quantity + " at " + group_sizes_text(s, less, true) + ": " + dvp(f(less.n1, less.n2), target) + " < " +
              dv(target) + " ≤ " + quantity + " at " + group_sizes_text(s, n, true) + ": " + dvp(at, target)};
}

}  // namespace ssc::procs::derive
