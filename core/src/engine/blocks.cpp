#include "engine/blocks.hpp"

#include <cmath>

#include "util/numtext.hpp"

namespace ssc::blocks {

ParamDef alternative(const std::string& two, const std::string& upper, const std::string& lower,
                     const std::string& def) {
  std::vector<Option> opts;
  if (!two.empty()) opts.push_back({"two", two, {}});
  if (!upper.empty()) opts.push_back({"upper", upper, {}});
  if (!lower.empty()) opts.push_back({"lower", lower, {}});
  return select("alternative", "Alternative Hypothesis", opts, def)
      .grp(G_TEST)
      .hint("Direction of the alternative hypothesis. One-sided tests put all of alpha in one tail.");
}

Alt alt(const Scenario& s) {
  auto it = s.sel.find("alternative");
  if (it == s.sel.end() || it->second == "two") return Alt::Two;
  return it->second == "upper" ? Alt::Upper : Alt::Lower;
}

std::string sidedness(Alt a) { return a == Alt::Two ? "two-sided" : "one-sided"; }

ParamDef power_param(const std::string& def) {
  return num("power", "Power", def)
      .grp(G_POWER)
      .sym("Target Power")
      .range(0, 1, true, true)
      .sugg({"0.80", "0.90", "0.80 0.90", "0.70 to 0.95 by 0.05"})
      .hint("Probability of rejecting the null hypothesis when the alternative is true (1 − β).")
      .show(when(not_solving({"power"})));
}

ParamDef alpha_param(const std::string& def, const std::string& label) {
  return num("alpha", label, def)
      .grp(G_POWER)
      .sym("Alpha")
      .range(0, 1, true, true)
      .sugg({"0.05", "0.01", "0.025", "0.01 0.05", "0.10"})
      .hint("Probability of rejecting a true null hypothesis (type I error rate).");
}

ParamDef dropout_param() {
  return num("dropout", "Dropout Rate", "0")
      .grp(G_N)
      .sym("Dropout")
      .range(0, 1, false, true)
      .sugg({"0", "0.10", "0.15", "0.20"})
      .hint("Expected proportion of subjects without usable data. Sample sizes are inflated by 1/(1 − rate).");
}

ParamDef n_param(const std::string& id, const std::string& label, const std::string& def, double min_n) {
  return num(id, label, def)
      .grp(G_N)
      .at_least(min_n)
      .whole()
      .sugg({def, "10 to 100 by 10", "20 40 60 80 100"})
      .show(when(not_solving({"n"})));
}

std::vector<ParamDef> two_group_params(double min_n, const std::string& n1_default) {
  std::vector<ParamDef> v;
  v.push_back(select("alloc", "Group Allocation",
                     {{"equal", "Equal (N1 = N2)", {}},
                      {"ratio", "Ratio R = N2/N1", {}},
                      {"n1", "Enter N1, solve for N2", when(solving({"n"}))},
                      {"n2", "Enter N2, solve for N1", when(solving({"n"}))},
                      {"pct", "Percentage of total in group 1", {}},
                      {"sep", "Enter N1 and N2 individually", when(not_solving({"n"}))}},
                     "equal")
                  .grp(G_N)
                  .hint("How the subjects are divided between the two groups."));
  v.push_back(num("n1", "N1 (Group 1 Sample Size)", n1_default)
                  .grp(G_N)
                  .sym("N1")
                  .at_least(min_n)
                  .whole()
                  .sugg({n1_default, "10 to 100 by 10", "20 40 60 80 100"})
                  .hint("Number of subjects in group 1 (with equal allocation, also in group 2).")
                  .show(when_any({{not_solving({"n"}), is("alloc", {"equal", "ratio", "sep"})},
                                  {solving({"n"}), is("alloc", {"n1"})}})));
  v.push_back(num("n2", "N2 (Group 2 Sample Size)", n1_default)
                  .grp(G_N)
                  .sym("N2")
                  .at_least(min_n)
                  .whole()
                  .sugg({n1_default, "10 to 100 by 10"})
                  .hint("Number of subjects in group 2, used when the two group sizes are entered separately or "
                        "when group 2's size is fixed and N1 is found.")
                  .show(when_any({{not_solving({"n"}), is("alloc", {"sep"})}, {solving({"n"}), is("alloc", {"n2"})}})));
  v.push_back(num("r", "R (Allocation Ratio N2/N1)", "1")
                  .grp(G_N)
                  .sym("R")
                  .above(0)
                  .sugg({"1", "2", "0.5", "1 2 3"})
                  .hint("Size of group 2 relative to group 1: N2 = R × N1, rounded up to a whole number.")
                  .show(when(is("alloc", {"ratio"}))));
  v.push_back(num("ntotal", "N (Total Sample Size)", "40")
                  .grp(G_N)
                  .sym("N")
                  .at_least(2 * min_n)
                  .whole()
                  .sugg({"40", "50 to 200 by 50"})
                  .hint("Total number of subjects in both groups; the percentage in group 1 decides how they are "
                        "split.")
                  .show(when({not_solving({"n"}), is("alloc", {"pct"})})));
  v.push_back(num("pct1", "Percent in Group 1", "50")
                  .grp(G_N)
                  .sym("% in G1")
                  .range(0, 100, true, true)
                  .sugg({"50", "33.33", "66.67", "25 50 75"})
                  .hint("The share of all subjects placed in group 1: N1 = N × percent / 100 rounded to the nearest "
                        "whole number; N2 = N − N1.")
                  .show(when(is("alloc", {"pct"}))));
  return v;
}

namespace {
double ratio_n2(double n1, double r) { return solve::round_up(r * n1); }
TwoN split_pct(double N, double pct) {
  double n1 = std::floor(N * pct / 100 + 0.5);
  return {n1, N - n1};
}
}  // namespace

TwoN two_group_n(const Scenario& s) {
  const std::string& a = s.s("alloc");
  if (a == "equal") return {s("n1"), s("n1")};
  if (a == "ratio") return {s("n1"), ratio_n2(s("n1"), s("r"))};
  if (a == "sep") return {s("n1"), s("n2")};
  if (a == "pct") return split_pct(s("ntotal"), s("pct1"));
  throw CalcError("unsupported allocation for this calculation: " + a);
}

TwoN solve_two_group_n(const Scenario& s, const std::function<double(double, double)>& power, double target,
                       double min_n, double guess_n1, bool monotone, const RunContext& ctx, double max_n) {
  const std::string& a = s.s("alloc");
  solve::NSearch o;
  o.monotone = monotone;
  o.hi = max_n;
  if (a == "equal") {
    o.lo = min_n;
    o.guess = guess_n1;
    const double n = solve::min_n([&](double n1) { return power(n1, n1); }, target, o, &ctx);
    return {n, n};
  }
  if (a == "ratio") {
    const double r = s("r");
    o.lo = solve::first_whole([&](double x) { return ratio_n2(x, r) >= min_n; }, min_n, max_n,
                              "the allocation ratio is too small: group 2 would have too few subjects");
    o.guess = guess_n1;
    const double n1 = solve::min_n([&](double x) { return power(x, ratio_n2(x, r)); }, target, o, &ctx);
    return {n1, ratio_n2(n1, r)};
  }
  // with one group fixed, the power is capped by that group however large the other: say so
  auto capped = [&](const char* fixed, double size, const char* other) {
    return CalcError("the target power cannot be reached with " + std::string(fixed) + " = " + fnum(size) +
                     ", however large " + other + " is: the fixed group limits the power");
  };
  if (a == "n1") {
    const double n1 = s("n1");
    o.lo = min_n;
    try {
      return {n1, solve::min_n([&](double x) { return power(n1, x); }, target, o, &ctx)};
    } catch (const CalcError& e) {
      if (std::string(e.what()).find("cannot be reached") == std::string::npos) throw;
      throw capped("N1", n1, "N2");
    }
  }
  if (a == "n2") {
    const double n2 = s("n2");
    o.lo = min_n;
    try {
      return {solve::min_n([&](double x) { return power(x, n2); }, target, o, &ctx), n2};
    } catch (const CalcError& e) {
      if (std::string(e.what()).find("cannot be reached") == std::string::npos) throw;
      throw capped("N2", n2, "N1");
    }
  }
  if (a == "pct") {
    const double pct = s("pct1");
    // from the total at which the smaller share reaches min_n (counting up from 2 min_n took 10^7 steps, a second,
    // at a percentage of 1e-300); both groups grow with the total, so the first total that holds is the same
    const double share = std::fmin(pct, 100 - pct) / 100;
    double lo = std::fmax(2 * min_n, std::floor((min_n - 0.5) / share) - 1);
    if (!(lo <= 1e7)) throw CalcError("the percentage leaves one group empty");
    while (true) {
      const TwoN t = split_pct(lo, pct);
      if (t.n1 >= min_n && t.n2 >= min_n) break;
      lo += 1;
      if (lo > 1e7) throw CalcError("the percentage leaves one group empty");
    }
    o.lo = lo;
    if (std::isfinite(guess_n1)) o.guess = guess_n1 * 100 / pct;
    const double N = solve::min_n(
        [&](double x) {
          const TwoN t = split_pct(x, pct);
          return power(t.n1, t.n2);
        },
        target, o, &ctx);
    return split_pct(N, pct);
  }
  throw CalcError("unsupported allocation: " + a);
}

void put_two_group(Row& r, const TwoN& n) {
  r.set("n1", n.n1);
  r.set("n2", n.n2);
  r.set("n", n.total());
}

std::vector<Column> two_group_columns() {
  return {col("n1", "N1", "int", "Sample size of group 1"), col("n2", "N2", "int", "Sample size of group 2"),
          col("n", "N", "int", "Total sample size N1 + N2")};
}

Column col(const std::string& key, const std::string& label, const std::string& fmt, const std::string& tip) {
  return Column{key, label, fmt, tip};
}

std::vector<Column> power_columns(const Scenario& s) {
  if (s.solve == "n")
    return {col("target_power", "Target Power", "prob", "Desired power"),
            col("power", "Actual Power", "prob", "Power achieved with the reported (whole-number) sample size")};
  return {col("power", "Power", "prob", "Probability of rejecting H0 when H1 is true")};
}

std::string fnum(double x, int max_decimals) {
  if (!std::isfinite(x)) return std::isnan(x) ? "—" : (x > 0 ? "∞" : "−∞");
  if (x == std::floor(x) && std::fabs(x) < 1e15) return text::format_fixed(x, 0);
  const double ax = std::fabs(x);
  if (ax >= 1e6 || ax < 1e-4) return text::format_g(x, 4);
  std::string s = text::format_fixed(x, max_decimals);
  if (s.find('.') != std::string::npos) {
    while (!s.empty() && s.back() == '0') s.pop_back();
    if (!s.empty() && s.back() == '.') s.pop_back();
  }
  if (s == "-0" || s.empty()) s = "0";
  return s;
}

std::string fpct(double p, int decimals) { return text::format_fixed(100 * p, decimals) + "%"; }

std::string fint(double x) { return text::format_fixed(x, 0); }

std::string power_text(double p) {
  // never print 100.0% for a power that is not exactly one: up to four decimals, then "more than 99.9999%"
  int dec = 1;
  while (dec < 4 && p < 1 && text::format_fixed(100 * p, dec) == "100." + std::string(dec, '0')) ++dec;
  if (p < 1 && text::format_fixed(100 * p, dec) == "100." + std::string(dec, '0')) return ">99.9999%";
  return fpct(p, dec);
}

}  // namespace ssc::blocks
