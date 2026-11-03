// Shared steps for the derivations of the proportion procedures: normal tail powers, binomial windows, and the
// size searches whose power is costly to compute again (an enumeration, taken from the row at the answer) or
// saw-toothed in the size.
#pragma once

#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "procedures/apart.hpp"
#include "procedures/binomial_tools.hpp"
#include "procedures/derive.hpp"

namespace ssc::procs::derive {

/// Whether the row carries the note (the normal approximation used where an enumeration was asked for, say).
inline bool has_note(const Row& r, const std::string& note) {
  for (const auto& n : r.notes)
    if (n == note) return true;
  return false;
}

/// A number below 1 in size, such as a standard error, with four significant digits (0.003569, not the 0.0036 of
/// four decimals, from which the steps after it would not give the answer back), in a formula's minus sign; num()
/// from 1 on.
inline std::string sig(double x) {
  if (!std::isfinite(x) || x == 0 || std::fabs(x) >= 1) return num(x);
  const std::string t = text::format_g(x, 4);
  return !t.empty() && t[0] == '-' ? "−" + t.substr(1) : t;
}

/// The power of a test whose numerator is normal with mean m and SD s1 and that rejects beyond z × s0 (less a
/// continuity correction c), one tail at a time: Φ((m − c − z s0)/s1) above and Φ((−m − c − z s0)/s1) below.
/// `up` and `lo` are the tails as the run computed them and `power` what it made of them; an empty `c_sym` leaves the
/// correction out.
inline Step tails_power(Alt a, const std::string& m_sym, const std::string& c_sym, const std::string& s0_sym,
                        const std::string& s1_sym, double m, double c, double z, double s0, double s1, double up,
                        double lo, double power) {
  const std::string cs = c_sym.empty() ? "" : " − " + c_sym, cn = c_sym.empty() ? "" : " − " + sig(c);
  const std::string rest = cs + " − z " + s0_sym + ")/" + s1_sym + ")";
  const std::string upper = "Φ((" + m_sym + rest, lower = "Φ((−" + m_sym + rest;
  auto at = [&](double mean) {
    return "Φ((" + sig(mean) + cn + " − " + sub(z) + " × " + sig(s0) + ")/" + sig(s1) + ")";
  };
  if (a == Alt::Two)
    return {"Power", upper + " + " + lower,
            at(m) + " + " + at(-m) + " = " + prob(up) + " + " + prob(lo) + " = " + prob(power)};
  if (a == Alt::Upper) return {"Power", upper, at(m) + " = " + prob(power)};
  return {"Power", lower, at(-m) + " = " + prob(power)};
}

/// "a − b" for a difference shown as `shown`: with that many decimals where it has more than four (0.1 − 0.05001 for
/// 0.04999), and otherwise with the digits it takes to tell the two apart (a one-sided bound next to its estimate:
/// 0.05 − 0.0499981, not 0.05 − 0.05), in a formula's minus signs.
inline std::string difference(double a, double b, const std::string& shown) {
  auto sign = [](std::string t) { return !t.empty() && t[0] == '-' ? "−" + t.substr(1) : t; };
  auto both = [&](const std::string& ta, const std::string& tb) {
    return sign(ta) + " − " + (b < 0 ? "(" + sign(tb) + ")" : sign(tb));
  };
  const std::size_t dot = shown.find('.');
  const int decimals = dot == std::string::npos ? 0 : static_cast<int>(shown.size() - dot - 1);
  if (shown.find('e') == std::string::npos && decimals > 4) return both(fnum(a, decimals), fnum(b, decimals));
  const auto [ta, tb] = fnum_apart(a, b);
  return both(ta, tb);
}

/// The counts an enumeration covers for one binomial (those carrying all but 1e-16 of its probability): "9–48"; q,
/// when given, is 1 − p (as binom_window_pq takes it).
inline std::string window_text(int n, double p, double q = NAN) {
  const BinomWindow w = binom_window_pq(n, p, std::isnan(q) ? 1 - p : q);
  return std::to_string(w.lo) + "–" + std::to_string(w.hi);
}

/// The exact power of a test on two binomial counts, groups a and b ("1" and "2"), as the probability of the tables
/// it rejects (`over`: "over the tables with |Z| > z").
inline std::string table_sum(const std::string& a, const std::string& b, const std::string& over) {
  return "Σ b(x" + a + "; N" + a + ", P" + a + ") b(x" + b + "; N" + b + ", P" + b + ") " + over +
         ", b(x; N, P) the binomial probability, over the counts that carry all but 1e-16 of each distribution";
}

/// The search's rule for an exact power, which is saw-toothed in the size, after the rule it shares with a smooth one.
inline std::string sawtooth_rule(const std::string& rule) {
  return rule + "; the exact power is saw-toothed in the size, so the search takes the first size that reaches the "
                "target in an upward scan from below where a bisection finds it reached (some larger sizes fall short "
                "again)";
}

/// How P1 follows from the effect given with P2 (as p1_from computes it); nothing when P1 itself was given.
inline void p1_step(const Scenario& s, double p1, double p2, Steps& out) {
  if (s.is("es", "diff"))
    out.push_back({"Group 1 proportion", "P1 = P2 + D1", num(p2) + " + " + sub(s("d1")) + " = " + num(p1)});
  else if (s.is("es", "ratio"))
    out.push_back({"Group 1 proportion", "P1 = R1 × P2", num(s("r1")) + " × " + num(p2) + " = " + num(p1)});
  else if (s.is("es", "or")) {
    const double o = s("or1") * p2 / (1 - p2);
    out.push_back({"Group 1 proportion", "P1 = O/(1 + O), O = OR1 × P2/(1 − P2), the odds of group 1",
                   "O = " + num(s("or1")) + " × " + num(p2) + "/" + num(1 - p2) + " = " + num(o) + ", P1 = " + num(o) +
                       "/" + num(1 + o) + " = " + num(p1)});
  }
}

/// How a size was found: the power one step below the answer and at it (`at`, the row's), the designs named by
/// `below_text` and `at_text` ("N = 52"); an empty `below_text`: the answer is the smallest design allowed.  Where a
/// saw-toothed search began its scan at the answer, the size below it can reach the target too, and is said so.
inline Step size_step(const std::string& what, const std::string& rule, const std::string& at_text, double at,
                      const std::string& below_text, double below, double target) {
  if (below_text.empty())
    return {what, rule,
            "at " + at_text + ", the smallest allowed, the power is already " + prob(at, target) + " ≥ " +
                fnum(target)};
  if (below >= target)
    return {what, rule,
            "power at " + below_text + ": " + prob(below, target) + " ≥ " + fnum(target) +
                " as well, but the search's scan began above it; power at " + at_text + ": " + prob(at, target)};
  return {what, rule,
          "power at " + below_text + ": " + prob(below, target) + " < " + fnum(target) + " ≤ power at " + at_text +
              ": " + prob(at, target)};
}

/// The rule of a two-group search under the scenario's allocation (as derive::two_group_search words it).
inline std::string two_group_rule(const Scenario& s) {
  const std::string& a = s.s("alloc");
  return a == "equal"   ? "N1 = N2, the smallest whole number whose power reaches the target"
         : a == "ratio" ? "N1 = the smallest whole number whose power reaches the target, with N2 = ⌈r × N1⌉ (r = " +
                              fnum(s("r")) + ")"
         : a == "n1"    ? "N2 = the smallest whole number whose power reaches the target, with N1 fixed"
         : a == "n2"    ? "N1 = the smallest whole number whose power reaches the target, with N2 fixed"
                        : "N = the smallest total whose power reaches the target, " + fnum(s("pct1")) +
                              "% of it (rounded) in group 1";
}

/// derive::two_group_search for a power that is costly (`at`, the power at the answer, is the row's; only the
/// smaller design's is computed) or saw-toothed (`sawtooth`).
inline Step two_group_step(const Scenario& s, const std::function<double(double, double)>& power,
                           const blocks::TwoN& n, double at, double target, double min_n, bool sawtooth) {
  const blocks::TwoN less = blocks::smaller_two_group(s, n);
  const auto sizes = [](const blocks::TwoN& t) { return "N1 = " + fint(t.n1) + ", N2 = " + fint(t.n2); };
  const std::string rule = sawtooth ? sawtooth_rule(two_group_rule(s)) : two_group_rule(s);
  if (less.n1 < min_n || less.n2 < min_n) return size_step("Group sizes", rule, sizes(n), at, "", NAN, target);
  return size_step("Group sizes", rule, sizes(n), at, sizes(less), power(less.n1, less.n2), target);
}

}  // namespace ssc::procs::derive
