// Steps shared by the derivations of the agreement, diagnostic, nonparametric and group-sequential procedures, in
// the manner of derive.hpp: searches whose power is saw-toothed in N, and group-size searches in units other than
// the power.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "procedures/derive.hpp"

namespace ssc::procs::derive {

/// Numbers listed in a formula: "0.3333, 0.6667, 1".
inline std::string nums(const std::vector<double>& v, int decimals = 4) {
  std::string t;
  for (std::size_t i = 0; i < v.size(); ++i) t += (i ? ", " : "") + num(v[i], decimals);
  return t;
}

/// How a single sample size was found where the power of an exact test is saw-toothed in N: the search's scan
/// upward from below its bracket returns the first size that reaches the target. Shown by the power one step below
/// and at the answer; when the size below reaches the target too, the answer is where that scan began.
inline Step saw_search(const std::function<double(double)>& power, double n, double target, double lo,
                       const std::string& sym = "N", double step = 1) {
  const std::string rule =
      sym + " = the smallest whole number whose power reaches the target (exact power is not monotone in " + sym +
      ": an upward scan finds the first size that reaches it)";
  const double at = power(n);
  if (n - step < lo)
    return {"Sample size", rule,
            sym + " = " + fint(n) + " is the smallest size allowed, and its power is already " + prob(at, target) +
                " ≥ " + fnum(target)};
  const double below = power(n - step);
  if (below < target)
    return {"Sample size", rule,
            "power at " + sym + " = " + fint(n - step) + ": " + prob(below, target) + " < " + fnum(target) +
                " ≤ power at " + sym + " = " + fint(n) + ": " + prob(at, target)};
  return {"Sample size", rule,
          "power at " + sym + " = " + fint(n) + ": " + prob(at, target) + " ≥ " + fnum(target) +
              ", the first size of the scan (power at " + sym + " = " + fint(n - step) + ", below it: " +
              prob(below, target) + ")"};
}

/// How two group sizes were found under the scenario's allocation rule when the search reaches a target of another
/// quantity than the power (`quantity`, such as "drift", with the target's symbol `target_sym`): the value one design
/// below the answer and at it (as derive::two_group_search).
inline Step two_group_search_of(const Scenario& s, const std::function<double(double, double)>& f,
                                const blocks::TwoN& n, double target, double min_n, const std::string& quantity,
                                const std::string& target_sym) {
  const std::string& a = s.s("alloc");
  const std::string reach = "the smallest whole number whose " + quantity + " reaches " + target_sym;
  const std::string rule =
      a == "equal"   ? "N1 = N2, " + reach
      : a == "ratio" ? "N1 = " + reach + ", with N2 = ⌈r × N1⌉ (r = " + fnum(s("r")) + ")"
      : a == "n1"    ? "N2 = " + reach + ", with N1 fixed"
      : a == "n2"    ? "N1 = " + reach + ", with N2 fixed"
                     : "N = the smallest total whose " + quantity + " reaches " + target_sym + ", " +
                           fnum(s("pct1")) + "% of it (rounded) in group 1";
  const auto sizes = [](const blocks::TwoN& t) { return "N1 = " + fint(t.n1) + ", N2 = " + fint(t.n2); };
  const double at = f(n.n1, n.n2);
  const blocks::TwoN less = blocks::smaller_two_group(s, n);
  if (less.n1 < min_n || less.n2 < min_n)
    return {"Group sizes", rule,
            sizes(n) + " are the smallest sizes allowed, and their " + quantity + " is already " + prob(at, target) +
                " ≥ " + num(target)};
  return {"Group sizes", rule,
          quantity + " at " + sizes(less) + ": " + prob(f(less.n1, less.n2), target) + " < " + num(target) + " ≤ " +
              quantity + " at " + sizes(n) + ": " + prob(at, target)};
}

}  // namespace ssc::procs::derive
