// Shared steps for the derivations of the proportion procedures: normal tail powers, binomial windows, the closed
// forms shown for comparison, and the size searches whose power is costly to compute again (an enumeration, taken
// from the row at the answer) or saw-toothed in the size. Values are marked given or derived (see derive.hpp).
#pragma once

#include <cmath>
#include <functional>
#include <string>
#include <utility>
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
/// A derived number as sig() writes it.
inline std::string dsig(double x) { return derived(sig(x)); }
/// A number marked given (an input) or derived, by `typed`.
inline std::string mv(double x, bool typed) { return typed ? gv(x) : dv(x); }
/// The level of the tail in which the test rejects: α/2, computed, or α itself.
inline std::string tail_level(double alpha, Alt a) { return a == Alt::Two ? dv(tail_alpha(alpha, a)) : gv(alpha); }

/// The power of a test whose numerator is normal with mean m and SD s1 and that rejects beyond z × s0 (less a
/// continuity correction c), one tail at a time: Φ((m − c − z s0)/s1) above and Φ((−m − c − z s0)/s1) below.
/// `up` and `lo` are the tails as the run computed them and `power` what it made of them; an empty `c_sym` leaves the
/// correction out. Every value is derived.
inline Step tails_power(Alt a, const std::string& m_sym, const std::string& c_sym, const std::string& s0_sym,
                        const std::string& s1_sym, double m, double c, double z, double s0, double s1, double up,
                        double lo, double power) {
  const std::string cs = c_sym.empty() ? "" : " − " + c_sym, cn = c_sym.empty() ? "" : " − " + dsig(c);
  const std::string rest = cs + " − z " + s0_sym + ")/" + s1_sym + ")";
  const std::string upper = "Φ((" + m_sym + rest, lower = "Φ((−" + m_sym + rest;
  auto at = [&](double mean) {
    return "Φ((" + dsig(mean) + cn + " − " + dvs(z) + " × " + dsig(s0) + ")/" + dsig(s1) + ")";
  };
  if (a == Alt::Two)
    return {"Power", upper + " + " + lower,
            at(m) + " + " + at(-m) + " = " + dvp(up) + " + " + dvp(lo) + " = " + dvp(power)};
  if (a == Alt::Upper) return {"Power", upper, at(m) + " = " + dvp(power)};
  return {"Power", lower, at(-m) + " = " + dvp(power)};
}

/// "a − b" for a difference shown as `shown`: with that many decimals where it has more than four (0.1 − 0.05001 for
/// 0.04999), and otherwise with the digits it takes to tell the two apart (a one-sided bound next to its estimate:
/// 0.05 − 0.0499981, not 0.05 − 0.05), in a formula's minus signs; each marked given or derived.
inline std::string difference(double a, double b, const std::string& shown, bool a_given = false,
                              bool b_given = false) {
  auto sign = [](std::string t) { return !t.empty() && t[0] == '-' ? "−" + t.substr(1) : t; };
  auto both = [&](const std::string& ta, const std::string& tb) {
    const std::string mb = b_given ? given(sign(tb)) : derived(sign(tb));
    return (a_given ? given(sign(ta)) : derived(sign(ta))) + " − " + (b < 0 ? "(" + mb + ")" : mb);
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
  return dvi(w.lo) + "–" + dvi(w.hi);
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

/// Whether P1 was typed as it is (not derived from a difference, ratio or odds ratio with P2, nor solved for).
inline bool p1_typed(const Scenario& s) { return s.solve != "p1" && s.is("es", "props"); }

/// How P1 follows from the effect given with P2 (as p1_from computes it); nothing when P1 itself was given.
inline void p1_step(const Scenario& s, double p1, double p2, Steps& out) {
  if (s.is("es", "diff"))
    out.push_back({"Group 1 proportion", "P1 = P2 + D1", gv(p2) + " + " + gvs(s("d1")) + " = " + dv(p1)});
  else if (s.is("es", "ratio"))
    out.push_back({"Group 1 proportion", "P1 = R1 × P2", gv(s("r1")) + " × " + gv(p2) + " = " + dv(p1)});
  else if (s.is("es", "or")) {
    const double o = s("or1") * p2 / (1 - p2);
    out.push_back({"Group 1 proportion", "P1 = O/(1 + O), O = OR1 × P2/(1 − P2), the odds of group 1",
                   "O = " + gv(s("or1")) + " × " + gv(p2) + "/" + dv(1 - p2) + " = " + dv(o) + ", P1 = " + dv(o) +
                       "/" + dv(1 + o) + " = " + dv(p1)});
  }
}

/// How a size was found: the power one step below the answer and at it (`at`, the row's), the designs named by
/// `below_text` and `at_text` (marked: "N = {52}"); an empty `below_text`: the answer is the smallest design allowed.
/// Where a saw-toothed search began its scan at the answer, the size below it can reach the target too, and is said
/// so.
inline Step size_step(const std::string& what, const std::string& rule, const std::string& at_text, double at,
                      const std::string& below_text, double below, double target) {
  if (below_text.empty())
    return {what, rule,
            "at " + at_text + ", the smallest allowed, the power is already " + dvp(at, target) + " ≥ " + gv(target)};
  if (below >= target)
    return {what, rule,
            "power at " + below_text + ": " + dvp(below, target) + " ≥ " + gv(target) +
                " as well, but the search's scan began above it; power at " + at_text + ": " + dvp(at, target)};
  return {what, rule,
          "power at " + below_text + ": " + dvp(below, target) + " < " + gv(target) + " ≤ power at " + at_text + ": " +
              dvp(at, target)};
}

/// The two group sizes, each marked by the allocation: typed (N1 of equal, ratio and separate groups, N2 of equal and
/// separate ones, the fixed group of a size solve) or computed (the sizes solved for, N2 = ⌈r × N1⌉, a percentage's
/// split of N).
inline std::pair<std::string, std::string> group_sizes(const Scenario& s, const blocks::TwoN& n) {
  const std::string& a = s.s("alloc");
  const bool solved = s.solve == "n";
  const bool typed1 = solved ? a == "n1" : (a == "equal" || a == "ratio" || a == "sep");
  const bool typed2 = solved ? a == "n2" : (a == "equal" || a == "sep");
  return {typed1 ? gvi(n.n1) : dvi(n.n1), typed2 ? gvi(n.n2) : dvi(n.n2)};
}
/// "N1 = 60, N2 = 80", marked as group_sizes marks them.
inline std::string group_sizes_text(const Scenario& s, const blocks::TwoN& n) {
  const auto [a, b] = group_sizes(s, n);
  return "N1 = " + a + ", N2 = " + b;
}

/// How the group sizes follow from the allocation where they were neither solved for nor both typed: N2 = ⌈r × N1⌉,
/// or a percentage's split of N.
inline void group_sizes_step(const Scenario& s, const blocks::TwoN& n, Steps& out) {
  if (s.solve == "n") return;
  const std::string& a = s.s("alloc");
  if (a == "ratio")
    out.push_back({"Group sizes", "N2 = ⌈r × N1⌉", "⌈" + gv(s("r")) + " × " + gvi(n.n1) + "⌉ = " + dvi(n.n2)});
  else if (a == "pct")
    out.push_back({"Group sizes", "N1 = N × pct/100, rounded, N2 = N − N1",
                   "N1 = " + gvi(s("ntotal")) + " × " + gv(s("pct1")) + "/100 = " + dvi(n.n1) + ", N2 = " +
                       gvi(s("ntotal")) + " − " + dvi(n.n1) + " = " + dvi(n.n2)});
}

/// The unknown of a two-group size solve, as find_size takes it: two_group_unknown, set off by a comma where it
/// carries a clause ("N1, with N2 = ⌈r × N1⌉ (r = 2),").
inline std::string sizes_unknown(const Scenario& s) {
  return two_group_unknown(s) + (s.is("alloc", "equal") ? "" : ",");
}

/// The rule of a two-group search under the scenario's allocation (as derive::two_group_search words it).
inline std::string two_group_rule(const Scenario& s) {
  const std::string& a = s.s("alloc");
  return a == "equal"   ? "N1 = N2, the smallest whole number whose power reaches 1 − β"
         : a == "ratio" ? "N1 = the smallest whole number whose power reaches 1 − β, with N2 = ⌈r × N1⌉"
         : a == "n1"    ? "N2 = the smallest whole number whose power reaches 1 − β, with N1 fixed"
         : a == "n2"    ? "N1 = the smallest whole number whose power reaches 1 − β, with N2 fixed"
                        : "N = the smallest total whose power reaches 1 − β, a fixed percentage of it (rounded) in "
                          "group 1";
}

/// derive::two_group_search with the fixed group of an allocation marked given, for a power that is costly (`at`, the
/// power at the answer, is the row's; only the smaller design's is computed) or saw-toothed (`sawtooth`).
inline Step two_group_step(const Scenario& s, const std::function<double(double, double)>& power,
                           const blocks::TwoN& n, double at, double target, double min_n, bool sawtooth) {
  const blocks::TwoN less = blocks::smaller_two_group(s, n);
  const std::string rule = sawtooth ? sawtooth_rule(two_group_rule(s)) : two_group_rule(s);
  if (less.n1 < min_n || less.n2 < min_n)
    return size_step("Smallest sizes", rule, group_sizes_text(s, n), at, "", NAN, target);
  return size_step("Smallest sizes", rule, group_sizes_text(s, n), at, group_sizes_text(s, less),
                   power(less.n1, less.n2), target);
}

/// The actual size of an exact test, where the row has it, to follow its power in an answer: " (actual α* = 0.0328)".
inline std::string alpha_star(const Row& r) {
  return r.has("actual_alpha") && std::isfinite(r.num("actual_alpha"))
             ? " (actual α* = " + dvp(r.num("actual_alpha")) + ")"
             : "";
}

/// The usual closed form for the size of a comparison of two proportions (Fleiss; Schlesselman's for cases and
/// controls), for comparison with the answer: group a's size, with group b r times as large (`r_sym` names r; empty
/// for equal groups, which the formula then leaves out), for the pooled test or the unpooled one, and with Fleiss's
/// continuity correction (`cc`: the corrected test, and Fisher's exact test, which it approximates). `pa_text` and
/// `pb_text` are the proportions as marked; `unknown` names the size ("N1 = N2").
inline Step two_prop_formula(Alt a, double alpha, double power, double pa, double pb, double r,
                             const std::string& r_sym, bool pooled, bool cc, const std::string& ga,
                             const std::string& gb, const std::string& pa_text, const std::string& pb_text,
                             const std::string& unknown) {
  const double za = stats::qnorm(tail_alpha(alpha, a), false), zb = stats::qnorm(power), d = pa - pb;
  const bool equal = r_sym.empty();
  const std::string Pa = "P" + ga, Pb = "P" + gb, R = equal ? "" : gv(r);
  const std::string lead = std::string("z(") + level_symbol(a) + ")", beta = "z(1 − β)";
  const std::string va = Pa + "(1 − " + Pa + ")", vb = Pb + "(1 − " + Pb + ")" + (equal ? "" : "/" + r_sym);
  const std::string na = pa_text + " × " + dv(1 - pa), nb = pb_text + " × " + dv(1 - pb) + (equal ? "" : "/" + R);
  const double v1 = pa * (1 - pa) + pb * (1 - pb) / r;
  std::string formula, numbers;
  double n;
  if (pooled) {
    const double pbar = (pa + r * pb) / (1 + r), v0 = (1 + 1 / r) * pbar * (1 - pbar);
    n = std::pow(za * std::sqrt(v0) + zb * std::sqrt(v1), 2) / (d * d);
    formula = unknown + " ≈ (" + lead + " √(" + (equal ? "2" : "(1 + 1/" + r_sym + ")") + " p̄(1 − p̄)) + " + beta +
              " √(" + va + " + " + vb + "))² / δ², p̄ = " +
              (equal ? "(" + Pa + " + " + Pb + ")/2" : "(" + Pa + " + " + r_sym + " " + Pb + ")/(1 + " + r_sym + ")");
    numbers = "p̄ = " +
              (equal ? "(" + pa_text + " + " + pb_text + ")/2"
                     : "(" + pa_text + " + " + R + " × " + pb_text + ")/(1 + " + R + ")") +
              " = " + dv(pbar) + ": (" + dv(za) + " × √(" + (equal ? "2" : "(1 + 1/" + R + ")") + " × " + dv(pbar) +
              " × " + dv(1 - pbar) + ") + " + dvs(zb) + " × √(" + na + " + " + nb + "))² / " + dvs(d) + "² = " + dv(n);
  } else {
    n = (za + zb) * (za + zb) * v1 / (d * d);
    formula = unknown + " ≈ (" + lead + " + " + beta + ")² (" + va + " + " + vb + ") / δ²";
    numbers = "(" + dv(za) + " + " + dvs(zb) + ")² × (" + na + " + " + nb + ") / " + dvs(d) + "² = " + dv(n);
  }
  formula += ", δ = " + Pa + " − " + Pb;
  if (cc) {
    // Fleiss's correction of the size for the continuity correction (1/N1 + 1/N2)/2
    const double k = equal ? 4 : 2 * (r + 1) / r, ncc = n / 4 * std::pow(1 + std::sqrt(1 + k / (n * std::fabs(d))), 2);
    formula += "; with Fleiss's continuity correction, n/4 × (1 + √(1 + " +
               (equal ? std::string("4/(n |δ|)") : "2(" + r_sym + " + 1)/(" + r_sym + " n |δ|)") +
               "))², n the size before it";
    numbers += "; " + dv(n) + "/4 × (1 + √(1 + " +
               (equal ? "4/(" + dv(n) + " × " + dv(std::fabs(d)) + ")"
                      : "2 × (" + R + " + 1)/(" + R + " × " + dv(n) + " × " + dv(std::fabs(d)) + ")") +
               "))² = " + dv(ncc);
  }
  formula += ", for comparison; the steps below give the answer";
  return {"Normal approximation", formula, numbers};
}

}  // namespace ssc::procs::derive
