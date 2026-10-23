// More proportion designs: non-inferiority / superiority / equivalence for two
// proportions (Farrington–Manning and Miettinen–Nurminen score tests), paired proportions (McNemar),
// chi-square tests (effect size w) and cluster-randomized proportions.
#include <limits>
#include <unordered_map>

#include "procedures/binomial_tools.hpp"
#include "procedures/common.hpp"
#include "procedures/mcnemar.hpp"
#include "util/parallel.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// =====================================================================================
// Farrington–Manning restricted maximum likelihood
// =====================================================================================
// Difference: maximize the two-sample binomial likelihood under p1 − p2 = s.
// Closed-form root of the score cubic (Farrington & Manning 1990), with a
// golden-section fallback when the trigonometric solution is degenerate.
void rmle_diff(double ph1, double ph2, double theta, double s, double& pt1, double& pt2) {
  const double a = 1 + theta;
  const double b = -(1 + theta + ph1 + theta * ph2 + s * (theta + 2));
  const double c = s * s + s * (2 * ph1 + theta + 1) + ph1 + theta * ph2;
  const double d = -ph1 * s * (1 + s);
  const double v = b * b * b / (27 * a * a * a) - b * c / (6 * a * a) + d / (2 * a);
  const double uu = b * b / (9 * a * a) - c / (3 * a);
  bool ok = uu > 0;
  if (ok) {
    const double u = std::copysign(std::sqrt(uu), v);
    double arg = v / (u * u * u);
    arg = std::max(-1.0, std::min(1.0, arg));
    const double w = (stats::PI + std::acos(arg)) / 3;
    pt1 = 2 * u * std::cos(w) - b / (3 * a);
    pt2 = pt1 - s;
    const double lo = std::max(0.0, s), hi = std::min(1.0, 1 + s);
    ok = std::isfinite(pt1) && pt1 >= lo - 1e-12 && pt1 <= hi + 1e-12;
    if (ok) {
      pt1 = std::min(std::max(pt1, lo), hi);
      pt2 = pt1 - s;
    }
  }
  if (!ok) {  // golden-section maximization (likelihood is concave in p1)
    const double lo = std::max(0.0, s) + 1e-12, hi = std::min(1.0, 1 + s) - 1e-12;
    auto L = [&](double p1) {
      const double p2 = p1 - s;
      auto lg = [](double x) { return x > 0 ? std::log(x) : -1e300; };
      return ph1 * lg(p1) + (1 - ph1) * lg(1 - p1) + theta * (ph2 * lg(p2) + (1 - ph2) * lg(1 - p2));
    };
    double x0 = lo, x3 = hi;
    const double g = 0.6180339887498949;
    double x1 = x3 - g * (x3 - x0), x2 = x0 + g * (x3 - x0);
    double f1 = L(x1), f2 = L(x2);
    for (int it = 0; it < 200 && x3 - x0 > 1e-14; ++it) {
      if (f1 < f2) {
        x0 = x1; x1 = x2; f1 = f2; x2 = x0 + g * (x3 - x0); f2 = L(x2);
      } else {
        x3 = x2; x2 = x1; f2 = f1; x1 = x3 - g * (x3 - x0); f1 = L(x1);
      }
    }
    pt1 = 0.5 * (x0 + x3);
    pt2 = pt1 - s;
  }
}

// Ratio: maximize under p1 = phi · p2 (quadratic in p2).
void rmle_ratio(double ph1, double ph2, double theta, double phi, double& pt1, double& pt2) {
  const double A = phi * (1 + theta);
  const double B = -(phi + ph1 + theta + phi * theta * ph2);
  const double C = ph1 + theta * ph2;
  const double disc = std::max(0.0, B * B - 4 * A * C);
  pt2 = (-B - std::sqrt(disc)) / (2 * A);
  pt2 = std::min(std::max(pt2, 0.0), std::min(1.0, 1 / phi));
  pt1 = phi * pt2;
}

// Odds ratio: maximize under p1/(1 − p1) = psi · p2/(1 − p2).  The score equations keep the expected number of
// responses n1·p1 + n2·p2 at the observed one, a quadratic in p2 (Miettinen & Nurminen 1985), whose root in [0, 1]
// is taken without cancellation (the quadratic is linear when psi = 1).
void rmle_odds(double ph1, double ph2, double theta, double psi, double& pt1, double& pt2) {
  const double m = ph1 + theta * ph2;  // responses per subject of group 1
  const double a = theta * (psi - 1), b = psi + theta - m * (psi - 1), c = -m;
  const double q = -0.5 * (b + std::copysign(std::sqrt(std::max(0.0, b * b - 4 * a * c)), b));
  double root = q != 0 ? c / q : 0;
  if (!(root >= -1e-12 && root <= 1 + 1e-12) && a != 0) root = q / a;
  pt2 = std::min(std::max(root, 0.0), 1.0);
  pt1 = psi * pt2 / (1 + (psi - 1) * pt2);
}

enum class Scale { Diff, Ratio, Odds };

struct NiSpec {
  Scale scale;
  double bound;  // D0, R0 or OR0
  bool upper;    // H1: the difference, ratio or odds ratio exceeds the bound
  double alpha;
};

// A score statistic for observed proportions: its numerator, the numerator's variance at the boundary of H0 (at the
// restricted maximum likelihood estimates), and how the numerator moves with the two proportions (for its variance
// under H1).  For the odds ratio the numerator is group 1's score ph1 − pt1 (Miettinen & Nurminen 1985), which moves
// with the proportions through pt1 as well.
struct Score {
  double num, var0, g1, g2;
};

Score fm_score(const NiSpec& S, double ph1, double ph2, double n1, double n2) {
  const double theta = n2 / n1;
  double pt1, pt2;
  if (S.scale == Scale::Ratio) {
    rmle_ratio(ph1, ph2, theta, S.bound, pt1, pt2);
    return {ph1 - S.bound * ph2, pt1 * (1 - pt1) / n1 + S.bound * S.bound * pt2 * (1 - pt2) / n2, 1, -S.bound};
  }
  if (S.scale == Scale::Odds) {
    rmle_odds(ph1, ph2, theta, S.bound, pt1, pt2);
    const double a = n1 * pt1 * (1 - pt1), b = n2 * pt2 * (1 - pt2);  // each group's information
    if (!(a + b > 0)) return {ph1 - pt1, 0, 1, 0};                     // every response alike: no information
    return {ph1 - pt1, a * b / (a + b) / (n1 * n1), b / (a + b), -a / (a + b) * theta};
  }
  rmle_diff(ph1, ph2, theta, S.bound, pt1, pt2);
  return {ph1 - ph2 - S.bound, pt1 * (1 - pt1) / n1 + pt2 * (1 - pt2) / n2, 1, -1};
}

// Score statistic for observed proportions.
double fm_z(const NiSpec& S, double ph1, double ph2, double n1, double n2) {
  const Score c = fm_score(S, ph1, ph2, n1, n2);
  if (!(c.var0 > 0)) return c.num == 0 ? 0.0 : std::copysign(INFINITY, c.num);
  return c.num / std::sqrt(c.var0);
}

// Normal-approximation power for one bound: the numerator has the mean and variance it has under the true
// proportions, and the critical value its variance at the boundary.
double fm_power_normal(const NiSpec& S, double n1, double n2, double p1, double p2) {
  using namespace stats;
  const Score c = fm_score(S, p1, p2, n1, n2);
  const double v1 = c.g1 * c.g1 * p1 * (1 - p1) / n1 + c.g2 * c.g2 * p2 * (1 - p2) / n2;
  const double z = qnorm(S.alpha, false);
  return pnorm(((S.upper ? c.num : -c.num) - z * std::sqrt(c.var0)) / std::sqrt(v1));
}

// Normal-approximation power of two one-sided score tests that must both reject (equivalence): L that the measure
// exceeds the lower limit, U that it falls below the upper one.  On the difference scale both numerators are the
// same estimate; on the ratio and odds-ratio scales they are jointly normal.
double tost_power_normal(const NiSpec& L, const NiSpec& U, double n1, double n2, double p1, double p2) {
  using namespace stats;
  const Score l = fm_score(L, p1, p2, n1, n2), u = fm_score(U, p1, p2, n1, n2);
  const double w1 = p1 * (1 - p1) / n1, w2 = p2 * (1 - p2) / n2;
  const double sl = std::sqrt(l.g1 * l.g1 * w1 + l.g2 * l.g2 * w2), su = std::sqrt(u.g1 * u.g1 * w1 + u.g2 * u.g2 * w2);
  const double z = qnorm(L.alpha, false);
  if (L.scale == Scale::Diff)
    return std::max(0.0, pnorm((-u.num - z * std::sqrt(u.var0)) / su) - pnorm((z * std::sqrt(l.var0) - l.num) / sl));
  const double r = (l.g1 * u.g1 * w1 + l.g2 * u.g2 * w2) / (sl * su);
  // P(num_L > z·sd0_L and num_U < −z·sd0_U), with the second written as −num_U > z·sd0_U
  return clamp01(bvn_upper((z * std::sqrt(l.var0) - l.num) / sl, (z * std::sqrt(u.var0) + u.num) / su, -r));
}

constexpr double MAX_CELLS = 4e7;

// A table's score tests cost about seven times the cell of a z-test, the unit of the exact budgets (as measured: 3.5
// against 0.5 ns a cell; the two tests of an equivalence share the loop over tables and cost about as much as one).
constexpr double FM_CELL_COST = 7;

// Exact power of a set of score tests that must all reject (one for NI, two for TOST); NaN when the
// enumeration is too large (the callers then use the normal approximation).
double fm_power_exact(const std::vector<NiSpec>& tests, double n1d, double n2d, double p1, double p2,
                      const RunContext& ctx) {
  // the counts are ints (binom_window computes n + 1); their sums are taken as doubles
  constexpr double most = std::numeric_limits<int>::max() - 1.0;
  if (!(n1d >= 0 && n2d >= 0 && n1d <= most && n2d <= most)) return NAN;
  const int n1 = static_cast<int>(n1d), n2 = static_cast<int>(n2d);
  const BinomWindow w1 = binom_window(n1, p1), w2 = binom_window(n2, p2);
  const double cells = static_cast<double>(w1.hi - w1.lo + 1) * (w2.hi - w2.lo + 1);
  if (cells > MAX_CELLS) return NAN;
  std::vector<double> zc;
  for (const auto& t : tests) zc.push_back(stats::qnorm(t.alpha, false));
  const std::size_t nr = static_cast<std::size_t>(w1.hi - w1.lo + 1);
  std::vector<double> part(nr, 0.0);
  auto body = [&](std::size_t k) {
    const int x1 = w1.lo + static_cast<int>(k);
    double s = 0;
    for (int x2 = w2.lo; x2 <= w2.hi; ++x2) {
      bool rej = true;
      for (std::size_t i = 0; i < tests.size() && rej; ++i) {
        const double z = fm_z(tests[i], static_cast<double>(x1) / n1, static_cast<double>(x2) / n2, n1, n2);
        rej = tests[i].upper ? z > zc[i] : z < -zc[i];
      }
      if (rej) s += w2.pmf[static_cast<std::size_t>(x2 - w2.lo)];
    }
    part[k] = w1.pmf[k] * s;
  };
  if (cells > 1e5 && ctx.threads > 1) {
    parallel_for(nr, ctx.threads, body, ctx.token);
  } else {
    for (std::size_t k = 0; k < nr; ++k) {
      ctx.check();  // so that Stop ends it (a step is a whole row: thousands of cells)
      body(k);
    }
  }
  double s = 0;
  for (double v : part) s += v;
  return stats::clamp01(s);
}

double p1_from(const Scenario& s, double p2) {
  double p1;
  if (s.is("es", "diff"))
    p1 = p2 + s("d1");
  else if (s.is("es", "ratio"))
    p1 = p2 * s("r1");
  else if (s.is("es", "or")) {
    const double o = s("or1") * p2 / (1 - p2);
    p1 = std::isinf(o) ? 1 : o / (1 + o);  // odds past the largest double: P1 is 1
  } else
    p1 = s("p1");
  if (!(p1 > 0 && p1 < 1)) throw CalcError("the implied P1 = " + fnum(p1) + " is not strictly between 0 and 1");
  return p1;
}

Scale scale_of(const Scenario& s) {
  return s.is("scale", "ratio") ? Scale::Ratio : (s.is("scale", "or") ? Scale::Odds : Scale::Diff);
}

// P1 on a boundary of the hypotheses: P2 plus the difference, P2 times the ratio, or the P1 whose odds are the odds
// ratio times those of P2.
double p1_at(Scale sc, double bound, double p2) {
  if (sc == Scale::Ratio) return bound * p2;
  if (sc == Scale::Odds) return bound * p2 / (1 - p2 + bound * p2);
  return p2 + bound;
}

// The comparison of P1 with P2 on a scale.
double measure(Scale sc, double p1, double p2) {
  if (sc == Scale::Ratio) return p1 / p2;
  if (sc == Scale::Odds) return p1 / (1 - p1) / (p2 / (1 - p2));
  return p1 - p2;
}

// Result keys and words for each scale: the boundary, the actual value, the quantity.
const char* actual_key(Scale sc) { return sc == Scale::Ratio ? "r1" : (sc == Scale::Odds ? "or1" : "d1"); }
std::string quantity(Scale sc) { return sc == Scale::Ratio ? "P1/P2" : (sc == Scale::Odds ? "the odds ratio" : "P1 − P2"); }

// A starting sample size per group for the searches: the normal formula on the scale's natural metric (the log scale
// for ratios), with `gap` the distance of the actual value from the boundary on that metric.
double group_guess(Scale sc, double zf, double p1, double p2, double gap) {
  double v;
  if (sc == Scale::Ratio)
    v = (1 - p1) / p1 + (1 - p2) / p2;
  else if (sc == Scale::Odds)
    v = 1 / (p1 * (1 - p1)) + 1 / (p2 * (1 - p2));
  else
    v = p1 * (1 - p1) + p2 * (1 - p2);
  return std::ceil(zf * v / (gap * gap));
}
double metric(Scale sc, double x) { return sc == Scale::Diff ? x : std::log(x); }

std::vector<ParamDef> actual_effect_params() {
  return {select("es", "Actual Effect Input",
                 {{"props", "P1", {}},
                  {"diff", "Difference D1 = P1 − P2", {}},
                  {"ratio", "Ratio R1 = P1/P2", {}},
                  {"or", "Odds ratio OR1", {}}},
                 "props")
              .grp(G_EFFECT)
              .hint("How you give the true proportion of the treatment group: directly as P1, or through its "
                    "relation to the reference proportion P2."),
          num("p1", "P1 (Actual Treatment Proportion)", "0.6")
              .grp(G_EFFECT)
              .sym("P1")
              .range(0, 1, true, true)
              .sugg({"0.6", "0.55 0.6"})
              .hint("The proportion of the treatment group (group 1) with the outcome that you expect to be true, "
                    "between 0 and 1.")
              .show(when(is("es", {"props"}))),
          num("d1", "D1 (Actual Difference P1 − P2)", "0")
              .grp(G_EFFECT)
              .sym("D1")
              .range(-1, 1, true, true)
              .sugg({"0", "0.05"})
              .hint("The true difference P1 − P2 you expect (0 when the treatments are equally effective); "
                    "P1 = P2 + D1.")
              .show(when(is("es", {"diff"}))),
          num("r1", "R1 (Actual Ratio P1/P2)", "1")
              .grp(G_EFFECT)
              .sym("R1")
              .above(0)
              .sugg({"1", "1.1"})
              .hint("The true ratio P1/P2 you expect (1 when the treatments are equally effective); P1 = R1 × P2.")
              .show(when(is("es", {"ratio"}))),
          num("or1", "OR1 (Actual Odds Ratio)", "1")
              .grp(G_EFFECT)
              .sym("OR1")
              .above(0)
              .sugg({"1", "1.2"})
              .hint("The true odds ratio you expect (1 when the treatments are equally effective): the odds "
                    "P1/(1 − P1) divided by the odds P2/(1 − P2).")
              .show(when(is("es", {"or"}))),
          num("p2", "P2 (Reference Proportion)", "0.6")
              .grp(G_EFFECT)
              .sym("P2")
              .range(0, 1, true, true)
              .sugg({"0.6", "0.5 0.6 0.7"})
              .hint("The proportion of the reference group (group 2, e.g. the standard treatment) expected to have "
                    "the outcome, between 0 and 1.")};
}

const char* FM_REF =
    "Farrington, C. P. & Manning, G. (1990). Test statistics and sample size formulae for comparative binomial trials "
    "with null hypothesis of non-zero risk difference or non-unity relative risk. Statistics in Medicine 9, "
    "1447–1454.";
const char* MN_REF =
    "Miettinen, O. & Nurminen, M. (1985). Comparative analysis of two rates. Statistics in Medicine 4, 213–226.";

// =====================================================================================
// Non-inferiority / superiority by a margin
// =====================================================================================
Procedure two_prop_ni() {
  Procedure p;
  p.id = "prop.two.ni";
  p.name = "Non-Inferiority Tests for Two Proportions";
  p.category = "Proportions";
  p.subcategory = "Two Independent Proportions";
  p.order = 2;
  p.summary =
      "One-sided score test that the difference, ratio or odds ratio of two proportions exceeds — or falls below — a "
      "margin (Farrington–Manning for the difference and ratio, Miettinen–Nurminen for the odds ratio), covering "
      "non-inferiority and superiority-by-a-margin designs. Power by exact enumeration or the normal approximation.";
  p.keywords = {"non-inferiority", "noninferiority", "farrington-manning", "score test", "margin", "risk difference",
                "relative risk", "risk ratio", "odds ratio", "miettinen-nurminen", "active control", "two proportions"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      select("scale", "Comparison Scale",
             {{"diff", "Difference P1 − P2", {}}, {"ratio", "Ratio P1/P2", {}}, {"or", "Odds ratio", {}}}, "diff")
          .grp(G_TEST)
          .hint("Compare the proportions by their difference P1 − P2, their ratio P1/P2 (the relative risk) or their "
                "odds ratio; the margin is given on the same scale."),
      select("h1", "Alternative Hypothesis",
             {{"upper", "H1: P1 − P2 > D0  (or P1/P2 > R0, OR > OR0)", {}},
              {"lower", "H1: P1 − P2 < D0  (or P1/P2 < R0, OR < OR0)", {}}},
             "upper")
          .grp(G_TEST)
          .hint("Upper: higher proportions are better (non-inferiority uses D0 < 0, or R0 or OR0 below 1). Lower: "
                "higher proportions are worse."),
      select("method", "Power Calculation",
             {{"normal", "Normal approximation", {}}, {"exact", "Exact (binomial enumeration)", {}}}, "normal")
          .grp(G_TEST)
          .hint("The normal approximation is quick and usually close; enumeration gives the exact power of the "
                "score test but is slow for large groups."),
      power_param(),
      alpha_param("0.025", "Alpha (One-Sided)").range(0, 0.5, true, true)};
  for (auto& q : two_group_params(2, "100")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  p.params.push_back(num("d0", "D0 (Difference Boundary under H0)", "-0.1")
                         .grp(G_EFFECT)
                         .sym("D0")
                         .range(-1, 1, true, true)
                         .sugg({"-0.1", "-0.05 -0.1", "0.05"})
                         .hint("The difference P1 − P2 at the boundary of the null hypothesis, e.g. −0.10 for a "
                               "non-inferiority margin of 10 percentage points when higher is better.")
                         .show(when(is("scale", {"diff"}))));
  p.params.push_back(num("r0", "R0 (Ratio Boundary under H0)", "0.8")
                         .grp(G_EFFECT)
                         .sym("R0")
                         .above(0)
                         .sugg({"0.8", "0.9", "1.25"})
                         .hint("The ratio P1/P2 at the boundary of the null hypothesis, e.g. 0.80 when the treatment "
                               "may be at most 20% worse (relative) and higher is better.")
                         .show(when(is("scale", {"ratio"}))));
  p.params.push_back(num("or0", "OR0 (Odds Ratio Boundary under H0)", "0.67")
                         .grp(G_EFFECT)
                         .sym("OR0")
                         .above(0)
                         .sugg({"0.67", "0.5 0.67", "1.5"})
                         .hint("The odds ratio at the boundary of the null hypothesis, e.g. 0.67 when the treatment's "
                               "odds of the outcome may be at most a third lower and higher is better.")
                         .show(when(is("scale", {"or"}))));
  for (auto& q : actual_effect_params()) p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Scale sc = scale_of(s);
    const char* bkey = sc == Scale::Ratio ? "r0" : (sc == Scale::Odds ? "or0" : "d0");
    const NiSpec S{sc, s(bkey), s.is("h1", "upper"), s("alpha")};
    const double p2 = s("p2"), p1 = p1_from(s, p2);
    const double p10 = p1_at(sc, S.bound, p2);  // P1 on the boundary of H0
    if (!(p10 > 0 && p10 < 1))
      throw CalcError("the boundary of the null hypothesis, P1 = " + fnum(p10) +
                      (sc == Scale::Ratio  ? " (R0 × P2)"
                       : sc == Scale::Odds ? " (odds OR0 times those of P2)"
                                           : " (P2 + D0)") +
                      ", must lie strictly between 0 and 1");
    bool exact = s.is("method", "exact");
    bool fell_back = false;
    ExactSearch search;
    auto P = [&](double n1, double n2) {
      // a point of the power curve is enumerated only while all the curves stay within the exact search's budget,
      // judged at its curve's largest point (in either group), at the cost of its score test
      const TwoN last = ctx.curve_last ? two_group_n(*ctx.curve_last) : TwoN{n1, n2};
      if (exact && !on_curve_too_long(ctx, FM_CELL_COST * exact_power_work(last.n1, last.n2, p1, p2))) {
        search.count(FM_CELL_COST * exact_power_work(n1, n2, p1, p2));
        const double v = fm_power_exact({S}, n1, n2, p1, p2, ctx);
        if (!std::isnan(v)) return v;
        fell_back = true;
      }
      return fm_power_normal(S, n1, n2, p1, p2);
    };
    TwoN n;
    if (s.solve == "n") {
      const double act = measure(sc, p1, p2);
      if ((S.upper && act <= S.bound) || (!S.upper && act >= S.bound))
        throw CalcError("the actual effect does not lie inside the alternative hypothesis");
      const double zf = z_factor(S.alpha, s("power"), Alt::Upper);
      double g;
      if (sc == Scale::Odds) {
        g = group_guess(sc, zf, p1, p2, std::log(act) - std::log(S.bound));
      } else {
        const double eff = sc == Scale::Ratio ? (p1 - S.bound * p2) : (p1 - p2 - S.bound);
        const double v = p1 * (1 - p1) + (sc == Scale::Ratio ? S.bound * S.bound : 1.0) * p2 * (1 - p2);
        g = zf * v / (eff * eff);
      }
      // the groups the search tries: with a ratio allocation, group 2 is R times group 1
      const double ratio = s.is("alloc", "ratio") ? s("r") : 1.0, g1 = g * (1 + 1 / ratio) / 2;
      if (exact && FM_CELL_COST * exact_search_work(g1, ratio * g1, p1, p2) > EXACT_SEARCH_BUDGET)
        exact = false, fell_back = true;
      n = search.run(exact, fell_back,
                     [&] { return solve_two_group_n(s, P, s("power"), 2, std::ceil(g), !exact, ctx, 1e7); });
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    if (exact) fell_back = false;  // the answer is labeled by its own evaluation, not by the search's steps
    put_power(s, r, P(n.n1, n.n2));
    r.set("p1", p1);
    r.set("p2", p2);
    r.set(bkey, S.bound);
    r.set(actual_key(sc), measure(sc, p1, p2));
    r.set("alpha", S.alpha);
    if (exact && !fell_back && !ctx.curve_points) {  // a curve uses only the power
      const double size = fm_power_exact({S}, n.n1, n.n2, std::min(std::max(p10, 1e-9), 1 - 1e-9), p2, ctx);
      if (!std::isnan(size)) r.set("actual_alpha", size);
    }
    if (fell_back) r.notes.push_back("Samples too large for exact enumeration; the normal approximation was used.");
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    const Scale sc = scale_of(s);
    c.push_back(col("p1", "P1", "num"));
    c.push_back(col("p2", "P2", "num"));
    if (sc == Scale::Ratio) {
      c.push_back(col("r0", "R0", "num", "Ratio under H0"));
      c.push_back(col("r1", "R1", "num", "Actual ratio"));
    } else if (sc == Scale::Odds) {
      c.push_back(col("or0", "OR0", "num", "Odds ratio under H0"));
      c.push_back(col("or1", "OR1", "num", "Actual odds ratio"));
    } else {
      c.push_back(col("d0", "D0", "num", "Difference under H0"));
      c.push_back(col("d1", "D1", "num", "Actual difference"));
    }
    c.push_back(col("alpha", "Alpha", "alpha"));
    if (s.is("method", "exact")) c.push_back(col("actual_alpha", "Actual Alpha", "prob", "Exact size at the H0 boundary"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const Scale sc = scale_of(s);
    const std::string q = sc == Scale::Ratio ? "P1/P2" : (sc == Scale::Odds ? "OR" : "P1 − P2");
    const std::string b = fnum(r.num(sc == Scale::Ratio ? "r0" : (sc == Scale::Odds ? "or0" : "d0")));
    const bool up = s.is("h1", "upper");
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
           " for the one-sided " + (sc == Scale::Odds ? "Miettinen–Nurminen" : "Farrington–Manning") +
           " score test of H0: " + q + (up ? " ≤ " : " ≥ ") + b + " versus H1: " + q + (up ? " > " : " < ") + b +
           (sc == Scale::Odds ? " (OR: the odds ratio of group 1 to group 2)" : "") + " at α = " + fnum(r.num("alpha")) +
           ", when the true proportions are P1 = " + fnum(r.num("p1")) + " and P2 = " + fnum(r.num("p2")) + ".";
  };
  p.definitions = {{"D0 / R0 / OR0", "Boundary of the null hypothesis (non-inferiority or superiority margin)."},
                   {"P1, P2", "Actual proportions in the treatment and reference groups."},
                   {"OR", "Odds ratio: the odds P1/(1 − P1) divided by the odds P2/(1 − P2)."},
                   {"Actual Alpha", "Exact rejection probability at the boundary of H0 (P2 fixed)."}};
  p.references = {FM_REF, MN_REF, REF_CHOW, REF_JULIOUS};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"d0", "-0.1"}, {"p1", "0.6"}, {"p2", "0.6"}};
  return p;
}
SSC_REGISTER(two_prop_ni);

// =====================================================================================
// Equivalence (TOST with two score tests)
// =====================================================================================
Procedure two_prop_equiv() {
  Procedure p;
  p.id = "prop.two.equiv";
  p.name = "Equivalence Tests for Two Proportions";
  p.category = "Proportions";
  p.subcategory = "Two Independent Proportions";
  p.order = 3;
  p.summary =
      "Two one-sided score tests showing that the difference, ratio or odds ratio of two proportions lies within "
      "equivalence limits (Farrington–Manning for the difference and ratio, Miettinen–Nurminen for the odds ratio). "
      "Power by the normal approximation or exact enumeration.";
  p.keywords = {"equivalence", "tost", "two proportions", "farrington-manning", "biosimilar", "therapeutic equivalence",
                "relative risk", "risk ratio", "odds ratio", "miettinen-nurminen"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {select("scale", "Comparison Scale",
                     {{"diff", "Difference P1 − P2", {}}, {"ratio", "Ratio P1/P2", {}}, {"or", "Odds ratio", {}}},
                     "diff")
                  .grp(G_TEST)
                  .hint("Compare the proportions by their difference P1 − P2, their ratio P1/P2 (the relative risk) "
                        "or their odds ratio; the equivalence limits are given on the same scale."),
              select("method", "Power Calculation",
                     {{"normal", "Normal approximation", {}}, {"exact", "Exact (binomial enumeration)", {}}}, "normal")
                  .grp(G_TEST)
                  .hint("The normal approximation is quick and usually close; enumeration gives the exact power of "
                        "the two score tests but is slow for large groups."),
              power_param(), alpha_param("0.05", "Alpha (Each One-Sided Test)").range(0, 0.5, true, true)};
  for (auto& q : two_group_params(2, "100")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  p.params.push_back(num("dl", "DL (Lower Equivalence Limit for P1 − P2)", "-0.1")
                         .grp(G_EFFECT)
                         .sym("DL")
                         .range(-1, 0, true, true)
                         .sugg({"-0.1", "-0.15"})
                         .hint("The most negative difference P1 − P2 that still counts as equivalent, e.g. −0.10 for "
                               "10 percentage points.")
                         .show(when(is("scale", {"diff"}))));
  p.params.push_back(num("du", "DU (Upper Equivalence Limit for P1 − P2)", "0.1")
                         .grp(G_EFFECT)
                         .sym("DU")
                         .range(0, 1, true, true)
                         .sugg({"0.1", "0.15"})
                         .hint("The largest positive difference P1 − P2 that still counts as equivalent, e.g. 0.10.")
                         .show(when(is("scale", {"diff"}))));
  p.params.push_back(num("rl", "RL (Lower Equivalence Limit for P1/P2)", "0.8")
                         .grp(G_EFFECT)
                         .sym("RL")
                         .range(0, 1, true, true)
                         .sugg({"0.8", "0.9"})
                         .hint("The smallest ratio P1/P2 that still counts as equivalent, e.g. 0.80.")
                         .show(when(is("scale", {"ratio"}))));
  p.params.push_back(num("ru", "RU (Upper Equivalence Limit for P1/P2)", "1.25")
                         .grp(G_EFFECT)
                         .sym("RU")
                         .above(1)
                         .sugg({"1.25", "1.1"})
                         .hint("The largest ratio P1/P2 that still counts as equivalent, e.g. 1.25 (= 1/0.80).")
                         .show(when(is("scale", {"ratio"}))));
  p.params.push_back(num("orl", "ORL (Lower Equivalence Limit for the Odds Ratio)", "0.67")
                         .grp(G_EFFECT)
                         .sym("ORL")
                         .range(0, 1, true, true)
                         .sugg({"0.67", "0.8"})
                         .hint("The smallest odds ratio that still counts as equivalent, e.g. 0.67.")
                         .show(when(is("scale", {"or"}))));
  p.params.push_back(num("oru", "ORU (Upper Equivalence Limit for the Odds Ratio)", "1.5")
                         .grp(G_EFFECT)
                         .sym("ORU")
                         .above(1)
                         .sugg({"1.5", "1.25"})
                         .hint("The largest odds ratio that still counts as equivalent, e.g. 1.5 (= 1/0.67).")
                         .show(when(is("scale", {"or"}))));
  for (auto& q : actual_effect_params()) p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Scale sc = scale_of(s);
    const char* lkey = sc == Scale::Ratio ? "rl" : (sc == Scale::Odds ? "orl" : "dl");
    const char* ukey = sc == Scale::Ratio ? "ru" : (sc == Scale::Odds ? "oru" : "du");
    const double p2 = s("p2"), p1 = p1_from(s, p2), alpha = s("alpha");
    const double lo = s(lkey), hi = s(ukey), act = measure(sc, p1, p2);
    if (!(act > lo && act < hi))
      throw CalcError("the actual " + std::string(sc == Scale::Ratio ? "ratio P1/P2" : (sc == Scale::Odds ? "odds ratio"
                                                                                                        : "difference")) +
                      " must lie inside the equivalence limits");
    if (sc == Scale::Ratio && !(hi * p2 < 1))
      throw CalcError("the upper equivalence limit, P1 = RU × P2 = " + fnum(hi * p2) + ", must lie below 1");
    const NiSpec L{sc, lo, true, alpha}, U{sc, hi, false, alpha};
    bool exact = s.is("method", "exact");
    bool fell_back = false;
    ExactSearch search;
    auto P = [&](double n1, double n2) {
      // a point of the power curve is enumerated only while all the curves stay within the exact search's budget,
      // judged at its curve's largest point (in either group), at the cost of its score tests
      const TwoN last = ctx.curve_last ? two_group_n(*ctx.curve_last) : TwoN{n1, n2};
      if (exact && !on_curve_too_long(ctx, FM_CELL_COST * exact_power_work(last.n1, last.n2, p1, p2))) {
        search.count(FM_CELL_COST * exact_power_work(n1, n2, p1, p2));
        const double v = fm_power_exact({L, U}, n1, n2, p1, p2, ctx);
        if (!std::isnan(v)) return v;
        fell_back = true;
      }
      return tost_power_normal(L, U, n1, n2, p1, p2);
    };
    TwoN n;
    if (s.solve == "n") {
      const double m = metric(sc, act), ml = metric(sc, lo), mu = metric(sc, hi);
      const double gap = std::min(m - ml, mu - m);
      const double zf = stats::qnorm(alpha, false) +
                        stats::qnorm(std::fabs(m - (ml + mu) / 2) < 1e-12 ? 1 - (1 - s("power")) / 2 : s("power"));
      const double g = group_guess(sc, zf * zf, p1, p2, gap);
      // the groups the search tries: with a ratio allocation, group 2 is R times group 1
      const double ratio = s.is("alloc", "ratio") ? s("r") : 1.0, g1 = g * (1 + 1 / ratio) / 2;
      if (exact && FM_CELL_COST * exact_search_work(g1, ratio * g1, p1, p2) > EXACT_SEARCH_BUDGET)
        exact = false, fell_back = true;
      n = search.run(exact, fell_back,
                     [&] { return solve_two_group_n(s, P, s("power"), 2, g, !exact, ctx, 1e7); });
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    if (exact) fell_back = false;  // the answer is labeled by its own evaluation, not by the search's steps
    put_power(s, r, P(n.n1, n.n2));
    r.set("p1", p1);
    r.set("p2", p2);
    r.set(actual_key(sc), act);
    r.set(lkey, lo);
    r.set(ukey, hi);
    r.set("alpha", alpha);
    if (fell_back) r.notes.push_back("Samples too large for exact enumeration; the normal approximation was used.");
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    const Scale sc = scale_of(s);
    c.push_back(col("p1", "P1", "num"));
    c.push_back(col("p2", "P2", "num"));
    if (sc == Scale::Ratio) {
      for (auto x : {col("r1", "R1", "num", "Actual ratio"), col("rl", "RL", "num"), col("ru", "RU", "num")})
        c.push_back(x);
    } else if (sc == Scale::Odds) {
      for (auto x : {col("or1", "OR1", "num", "Actual odds ratio"), col("orl", "ORL", "num"), col("oru", "ORU", "num")})
        c.push_back(x);
    } else {
      for (auto x : {col("d1", "D1", "num", "Actual difference"), col("dl", "DL", "num"), col("du", "DU", "num")})
        c.push_back(x);
    }
    c.push_back(col("alpha", "Alpha", "alpha"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const Scale sc = scale_of(s);
    const char* lkey = sc == Scale::Ratio ? "rl" : (sc == Scale::Odds ? "orl" : "dl");
    const char* ukey = sc == Scale::Ratio ? "ru" : (sc == Scale::Odds ? "oru" : "du");
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
           " to conclude that " + quantity(sc) + " lies between " + fnum(r.num(lkey)) + " and " + fnum(r.num(ukey)) +
           " using two one-sided score tests at α = " + fnum(r.num("alpha")) + " each, when P1 = " + fnum(r.num("p1")) +
           " and P2 = " + fnum(r.num("p2")) + ".";
  };
  p.definitions = {{"DL, DU / RL, RU / ORL, ORU", "Equivalence limits for P1 − P2, P1/P2 or the odds ratio."},
                   {"D1 / R1 / OR1", "Actual difference, ratio or odds ratio."}};
  p.references = {FM_REF, MN_REF, REF_CHOW};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"p1", "0.6"}, {"p2", "0.6"}};
  return p;
}
SSC_REGISTER(two_prop_equiv);

// =====================================================================================
// McNemar test for paired proportions
// =====================================================================================
}  // namespace

double mcnemar_power(const Scenario& s, int n, double p10, double p01, McNemarCache& cache, const RunContext& ctx) {
  using namespace stats;
  using namespace blocks;
  const double psi = p10 + p01, delta = p10 - p01, alpha = s("alpha");
  const Alt a = alt(s);
  if (s.is("method", "normal")) {  // Connor (1987)
    const double sd1 = std::sqrt(psi - delta * delta);
    const double rn = std::sqrt(static_cast<double>(n));
    if (a == Alt::Two) {
      const double z = qnorm(half_alpha(alpha), false);
      return std::min(1.0, pnorm((rn * delta - z * std::sqrt(psi)) / sd1) + pnorm((-rn * delta - z * std::sqrt(psi)) / sd1));
    }
    const double z = qnorm(alpha, false);
    return pnorm(((a == Alt::Upper ? 1 : -1) * rn * delta - z * std::sqrt(psi)) / sd1);
  }
  // exact conditional binomial test, unconditional power: sum over the number m of discordant pairs
  const BinomWindow wm = binom_window(n, psi);
  const double pi = p10 / psi;
  const double side = tail_alpha(alpha, a);
  double total = 0;
  int cu = -1, cl = -2;  // critical values of the previous m (−1, −2: not known)
  for (int m = wm.lo; m <= wm.hi; ++m) {
    ctx.check();
    if (m == 0) continue;
    auto hit = cache.find(m);
    if (hit == cache.end()) {
      // critical values move by at most one per step in m
      if (a != Alt::Lower) {
        if (cu < 0) cu = upper_critical(m, 0.5, side);
        while (cu <= m && pbinom(cu - 1, m, 0.5, false) > side) ++cu;
        while (cu > 0 && pbinom(cu - 2, m, 0.5, false) <= side) --cu;
      }
      if (a != Alt::Upper) {
        if (cl < -1) cl = lower_critical(m, 0.5, side);
        while (cl + 1 <= m && pbinom(cl + 1, m, 0.5) <= side) ++cl;
        while (cl >= 0 && pbinom(cl, m, 0.5) > side) --cl;
      }
      hit = cache.emplace(m, two_tail_prob(m, pi, a == Alt::Upper ? -1 : cl, a == Alt::Lower ? m + 1 : cu)).first;
    } else {
      cu = -1;
      cl = -2;
    }
    total += wm.at(m) * hit->second;
  }
  return clamp01(total);
}

namespace {

Procedure mcnemar() {
  Procedure p;
  p.id = "prop.paired.mcnemar";
  p.name = "Tests for Two Correlated Proportions (McNemar Test)";
  p.category = "Proportions";
  p.subcategory = "Paired Proportions";
  p.order = 1;
  p.summary =
      "Sample size and power for McNemar's test of paired binary outcomes (matched pairs or pre/post designs), from "
      "the two discordant-cell proportions. Exact power sums the conditional exact binomial test over the number of "
      "discordant pairs; Connor's normal approximation is also available.";
  p.keywords = {"mcnemar", "paired binary", "matched pairs", "discordant", "pre-post", "correlated proportions", "discordant pairs", "matched case-control", "paired proportions"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Pairs)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: P10 ≠ P01)", "One-sided (H1: P10 > P01)", "One-sided (H1: P10 < P01)"),
      select("method", "Test / Power Calculation",
             {{"exact", "Exact conditional test, exact power", {}}, {"normal", "Normal approximation (Connor)", {}}},
             "exact")
          .grp(G_TEST)
          .hint("The exact conditional (binomial) McNemar test with its exact power, or Connor's normal "
                "approximation, which is smoother but approximate."),
      power_param(),
      alpha_param(),
      n_param("n", "N (Number of Pairs)", "100", 1)
          .sym("N")
          .hint("Number of pairs, such as subjects with a yes-or-no outcome measured twice, or matched pairs of "
                "subjects."),
      dropout_param(),
      select("es", "Effect Size Input",
             {{"props", "Discordant proportions P10 and P01", {}}, {"or", "Discordance ψ and odds ratio P10/P01", {}}},
             "props")
          .grp(G_EFFECT)
          .hint("Describe the effect by the shares of the two kinds of discordant pairs (pairs whose outcomes "
                "differ), or by their total share ψ and their ratio."),
      num("p10", "P10 (Proportion: Yes then No)", "0.25")
          .grp(G_EFFECT)
          .sym("P10")
          .range(0, 1, true, true)
          .sugg({"0.25", "0.2 0.25"})
          .hint("Expected share of pairs whose first outcome is yes and second is no, between 0 and 1.")
          .show(when(is("es", {"props"}))),
      num("p01", "P01 (Proportion: No then Yes)", "0.1")
          .grp(G_EFFECT)
          .sym("P01")
          .range(0, 1, true, true)
          .sugg({"0.1", "0.05 0.1"})
          .hint("Expected share of pairs whose first outcome is no and second is yes; the test compares it with "
                "P10.")
          .show(when(is("es", {"props"}))),
      num("psi", "ψ (Discordant Proportion P10 + P01)", "0.35")
          .grp(G_EFFECT)
          .sym("ψ")
          .range(0, 1, true, false)
          .sugg({"0.35", "0.2 0.3 0.4"})
          .hint("Expected share of discordant pairs, whose two outcomes differ; only these pairs carry "
                "information for the test.")
          .show(when(is("es", {"or"}))),
      num("orr", "OR = P10/P01", "2.5")
          .grp(G_EFFECT)
          .sym("OR")
          .above(0)
          .sugg({"2", "2.5", "1.5 2 3"})
          .hint("Expected ratio of the two kinds of discordant pairs, P10/P01; 1 means no effect.")
          .show(when(is("es", {"or"}))),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    double p10, p01;
    if (s.is("es", "or")) {
      const double psi = s("psi"), o = s("orr");
      p10 = psi * o / (1 + o);
      p01 = psi / (1 + o);
    } else {
      p10 = s("p10");
      p01 = s("p01");
    }
    if (!(p10 + p01 <= 1)) throw CalcError("P10 + P01 cannot exceed 1");
    if (p10 == p01) throw CalcError("P10 equals P01; there is no effect");
    const Alt a = alt(s);
    McNemarCache cache;
    // the exact power sums over ~√n discordant counts, each a few binomial tails (0.25 s at a million pairs)
    const double n_max = s.is("method", "normal") ? 1e7 : 1e6;
    int n;
    if (s.solve == "n") {
      if ((a == Alt::Upper && p10 < p01) || (a == Alt::Lower && p10 > p01))
        throw CalcError("the effect has the sign opposite to the one-sided alternative");
      const double psi = p10 + p01, d = p10 - p01;
      const double za = stats::qnorm(tail_alpha(s("alpha"), a), false), zb = stats::qnorm(s("power"));
      solve::NSearch o;
      o.lo = 1;
      o.guess = std::ceil(std::pow(za * std::sqrt(psi) + zb * std::sqrt(psi - d * d), 2) / (d * d));
      o.monotone = s.is("method", "normal");
      o.hi = n_max;
      n = static_cast<int>(solve::min_n(
          [&](double x) { return mcnemar_power(s, static_cast<int>(x), p10, p01, cache, ctx); }, s("power"), o, &ctx));
    } else {
      if (!s.is("method", "normal") && s("n") > n_max)
        throw CalcError("the exact test is computed for up to 1,000,000 pairs; use the normal approximation for "
                        "more");
      n = static_cast<int>(s("n"));
    }
    r.set("n", n);
    put_power(s, r, mcnemar_power(s, n, p10, p01, cache, ctx));
    r.set("p10", p10);
    r.set("p01", p01);
    r.set("psi", p10 + p01);
    r.set("delta", p10 - p01);
    r.set("orr", p10 / p01);
    r.set("alpha", s("alpha"));
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int", "Number of pairs"), col("p10", "P10", "num"), col("p01", "P01", "num"),
                   col("psi", "ψ", "num", "Discordant proportion"), col("delta", "P10 − P01", "num"),
                   col("orr", "OR", "num", "P10/P01"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return fint(r.num("n")) + " pairs " + achieves(s, r) + " to detect a difference between the discordant "
           "proportions P10 = " + fnum(r.num("p10")) + " and P01 = " + fnum(r.num("p01")) + " (odds ratio " +
           fnum(r.num("orr"), 3) + ") using a " + alt_words(alt(s)) + " McNemar test" +
           (s.is("method", "exact") ? " (exact conditional binomial)" : " (normal approximation)") + " at α = " +
           fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"N", "Number of pairs."},
                   {"P10", "Proportion of pairs with response in the first but not the second condition."},
                   {"P01", "Proportion of pairs with response in the second but not the first condition."},
                   {"ψ", "Proportion of discordant pairs P10 + P01."}};
  p.references = {"Connor, R. J. (1987). Sample size for testing differences in proportions for the paired-sample "
                  "design. Biometrics 43, 207–211.",
                  "Lachin, J. M. (1992). Power and sample size evaluation for the McNemar test with application to "
                  "matched case-control studies. Statistics in Medicine 11, 1239–1251.",
                  REF_FLEISS};
  p.curve_param = "n";
  p.curve_label = "N (Pairs)";
  p.curve_min = 1;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"p10", "0.25"}, {"p01", "0.1"}};
  return p;
}
SSC_REGISTER(mcnemar);

// =====================================================================================
// Chi-square tests (effect size w)
// =====================================================================================
double chisq_w(const Scenario& s, double& df) {
  if (s.is("es", "w")) {
    df = s("df");
    return s("w");
  }
  auto p0 = s.v("p0s"), p1 = s.v("p1s");
  if (p0.size() != p1.size() || p0.size() < 2)
    throw CalcError("the H0 and H1 proportion lists must have the same length (at least 2)");
  double s0 = 0, s1 = 0;
  for (double v : p0) s0 += v;
  for (double v : p1) s1 += v;
  double w2 = 0;
  for (std::size_t i = 0; i < p0.size(); ++i) {
    const double a = p0[i] / s0, b = p1[i] / s1;
    if (!(a > 0)) throw CalcError("H0 proportions must be positive");
    w2 += (b - a) * (b - a) / a;
  }
  df = static_cast<double>(p0.size() - 1);
  return std::sqrt(w2);
}

double chisq_power(double n, double w, double df, double alpha) {
  const double crit = stats::qchisq(alpha, df, false);
  return stats::pnchisq(crit, df, n * w * w, false);
}

Procedure chisq() {
  Procedure p;
  p.id = "prop.chisq";
  p.name = "Chi-Square Tests (Goodness-of-Fit and Contingency Tables)";
  p.category = "Proportions";
  p.subcategory = "Contingency Tables";
  p.order = 1;
  p.summary =
      "Sample size and power for chi-square tests, from Cohen's effect size w and the degrees of freedom (any r × c "
      "table: df = (r − 1)(c − 1)) or from null and alternative cell proportions (goodness of fit). Power uses the "
      "noncentral chi-square distribution with λ = N w².";
  p.keywords = {"chi-square", "contingency table", "goodness of fit", "cohen w", "independence", "multinomial", "r x c table", "cross-tabulation", "categorical", "test of independence"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}, {"w", "Effect Size (w)", "w", "Detectable w"}};
  p.params = {
      power_param(),
      alpha_param(),
      n_param("n", "N (Total Sample Size)", "100", 1)
          .sym("N")
          .hint("Total number of subjects classified into the categories or table cells."),
      dropout_param(),
      select("es", "Effect Size Input",
             {{"w", "Cohen's w and degrees of freedom", {}}, {"gof", "Cell proportions under H0 and H1", {}}}, "w")
          .grp(G_EFFECT)
          .hint("Give Cohen's effect size w and the degrees of freedom, or the category proportions under H0 and "
                "H1 (goodness of fit), from which w is computed.")
          .show(when(not_solving({"w"}))),
      num("w", "w (Cohen's Effect Size)", "0.3").grp(G_EFFECT).sym("w").above(0).sugg({"0.1", "0.3", "0.5", "0.1 0.3 0.5"})
          .hint("Cohen's w: how far the category proportions under H1 depart from those under H0. By Cohen's "
                "conventions 0.1 is small, 0.3 medium and 0.5 large.")
          .show(when({not_solving({"w"}), is("es", {"w"})})),
      num("df", "Degrees of Freedom", "2").grp(G_EFFECT).sym("DF").at_least(1).whole().sugg({"1", "2", "4"})
          .hint("Degrees of freedom of the test: (rows − 1) × (columns − 1) for a table, or the number of "
                "categories − 1 for goodness of fit.")
          .show(when_any({{solving({"w"})}, {is("es", {"w"})}})),
      vec("p0s", "Proportions under H0", "0.25 0.25 0.25 0.25").grp(G_EFFECT).sym("P0").above(0)
          .hint("The share of each category under the null hypothesis, one value per category; they are rescaled "
                "to sum to 1.")
          .show(when({not_solving({"w"}), is("es", {"gof"})})),
      vec("p1s", "Proportions under H1", "0.2 0.2 0.3 0.3").grp(G_EFFECT).sym("P1").at_least(0)
          .hint("The share of each category that you expect to be true, in the same order as under H0; they are "
                "rescaled to sum to 1.")
          .show(when({not_solving({"w"}), is("es", {"gof"})})),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double alpha = s("alpha");
    // a test on df degrees of freedom compares at least df + 1 categories, each of which needs a subject (a huge w
    // otherwise gave a total sample size of 1)
    auto check_n = [&](double n, double df) {
      if (n < df + 1)
        throw CalcError("N must be at least " + fint(df + 1) + ": a chi-square test on " + fint(df) +
                        (df == 1 ? " degree" : " degrees") + " of freedom compares at least " + fint(df + 1) +
                        " categories");
    };
    double df = 0, w = 0, n;
    if (s.solve == "w") {
      df = s("df");
      n = s("n");
      check_n(n, df);
      w = solve::solve_increasing([&](double x) { return chisq_power(n, x, df, alpha); }, s("power"), 1e-9, 0.3, 1e6);
    } else {
      w = chisq_w(s, df);
      if (!(w > 0)) throw CalcError("the effect size w is zero");
      if (!std::isfinite(w * w)) throw CalcError("the effect size w is too large to compute with");
      if (s.solve == "n") {
        solve::NSearch o;
        o.lo = df + 1;
        o.guess = (z_factor(alpha, s("power"), Alt::Two) + 2 * (df - 1)) / (w * w);
        n = solve::min_n([&](double x) { return chisq_power(x, w, df, alpha); }, s("power"), o, &ctx);
      } else {
        n = s("n");
        check_n(n, df);
      }
    }
    r.set("n", n);
    put_power(s, r, chisq_power(n, w, df, alpha));
    r.set("w", w);
    r.set("df", df);
    r.set("lambda", n * w * w);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("w", "w", "num", "Cohen's effect size"), col("df", "DF", "int"),
                   col("lambda", "λ = N w²", "num", "Noncentrality"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A total sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) +
           " to detect an effect size of w = " + fnum(r.num("w")) + " with a chi-square test on " + fint(r.num("df")) +
           " degrees of freedom at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"w", "Cohen's effect size √Σ(P1ᵢ − P0ᵢ)²/P0ᵢ."},
                   {"DF", "Degrees of freedom: k − 1 for goodness of fit, (r − 1)(c − 1) for an r × c table."},
                   {"λ", "Noncentrality parameter N w²."}};
  p.references = {REF_COHEN, REF_BK,
                  "Guenther, W. C. (1977). Power and sample size for approximate chi-square tests. The American "
                  "Statistician 31, 83–85."};
  p.curve_param = "n";
  p.curve_label = "N (Total Sample Size)";
  p.curve_min = 1;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"w", "0.1 0.3 0.5"}, {"df", "2"}};
  return p;
}
SSC_REGISTER(chisq);

// =====================================================================================
// Cluster-randomized two proportions
// =====================================================================================
Procedure cluster_props() {
  Procedure p;
  p.id = "prop.two.cluster";
  p.name = "Two Proportions in a Cluster-Randomized Design";
  p.category = "Cluster-Randomized";
  p.subcategory = "Proportions";
  p.order = 2;
  p.summary =
      "Number of clusters per arm (or power) for comparing two proportions when clusters are randomized, using the "
      "pooled z-test with variances inflated by the design effect 1 + (m − 1)ρ.";
  p.keywords = {"cluster", "group randomized", "icc", "design effect", "binary", "community", "cluster trial", "binary outcome"};
  p.solve_for = {{"k", "Number of Clusters", "k1", "Clusters in Group 1"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: P1 ≠ P2)", "One-sided (H1: P1 > P2)", "One-sided (H1: P1 < P2)"),
      power_param(),
      alpha_param(),
      num("k1", "k1 (Clusters in Group 1)", "10")
          .grp(G_N)
          .sym("k1")
          .at_least(1)
          .whole()
          .sugg({"10", "5 to 30 by 5"})
          .hint("Number of clusters (such as clinics, schools or villages) randomized to group 1.")
          .show(when(not_solving({"k"}))),
      num("kr", "Cluster Ratio k2/k1", "1")
          .grp(G_N)
          .sym("k2/k1")
          .above(0)
          .sugg({"1"})
          .hint("Clusters in group 2 per cluster in group 1; k2 = ratio × k1, rounded up. Use 1 for equal "
                "numbers of clusters."),
      num("m", "m (Subjects per Cluster)", "50")
          .grp(G_N)
          .sym("m")
          .at_least(1)
          .whole()
          .sugg({"50", "20 50 100"})
          .hint("Number of subjects observed in each cluster, assumed the same for every cluster."),
      num("icc", "ρ (Intracluster Correlation)", "0.02")
          .grp(G_EFFECT)
          .sym("ICC")
          .range(0, 1, false, true)
          .sugg({"0.01", "0.02", "0.05"})
          .hint("The share of the outcome's variance that lies between clusters, from 0 to 1; values from 0.01 to "
                "0.05 are common. Take it from earlier studies of similar clusters."),
      num("p1", "P1 (Group 1 Proportion)", "0.3")
          .grp(G_EFFECT)
          .sym("P1")
          .range(0, 1, true, true)
          .sugg({"0.3", "0.25 0.3"})
          .hint("The share of subjects in group 1 (e.g. the intervention) expected to have the outcome, between "
                "0 and 1."),
      num("p2", "P2 (Group 2 Proportion)", "0.2")
          .grp(G_EFFECT)
          .sym("P2")
          .range(0, 1, true, true)
          .sugg({"0.2"})
          .hint("The share of subjects in group 2 (e.g. the control) expected to have the outcome, between 0 "
                "and 1."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    using namespace stats;
    const double alpha = s("alpha"), m = s("m"), icc = s("icc"), kr = s("kr"), p1 = s("p1"), p2 = s("p2");
    const double de = 1 + (m - 1) * icc;
    const Alt a = alt(s);
    auto k2_of = [&](double k1) { return std::ceil(kr * k1 - 1e-9); };
    auto P = [&](double k1) {
      const double k2 = k2_of(k1);
      const double n1 = k1 * m, n2 = k2 * m;
      const double pb = (n1 * p1 + n2 * p2) / (n1 + n2);
      const double s0 = std::sqrt(de * pb * (1 - pb) * (1 / n1 + 1 / n2));
      const double s1 = std::sqrt(de * (p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2));
      const double d = p1 - p2;
      if (a == Alt::Two) {
        const double z = qnorm(half_alpha(alpha), false);
        return std::min(1.0, pnorm((d - z * s0) / s1) + pnorm((-d - z * s0) / s1));
      }
      const double z = qnorm(alpha, false);
      return pnorm(((a == Alt::Upper ? d : -d) - z * s0) / s1);
    };
    double k1;
    if (s.solve == "k") {
      if (p1 == p2) throw CalcError("P1 equals P2");
      if ((a == Alt::Upper && (p1 - p2) < 0) || (a == Alt::Lower && (p1 - p2) > 0))
        throw CalcError("P1 − P2 has the sign opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = solve::first_whole([&](double x) { return k2_of(x) >= 1; }, 1, o.hi,
                                "the cluster ratio is too small: group 2 would have no clusters");
      o.guess = z_factor(alpha, s("power"), a) * de * (p1 * (1 - p1) + p2 * (1 - p2) / kr) / (m * (p1 - p2) * (p1 - p2)) + 1;
      k1 = solve::min_n(P, s("power"), o, &ctx);
    } else {
      k1 = s("k1");
    }
    const double k2 = k2_of(k1);
    r.set("k1", k1);
    r.set("k2", k2);
    r.set("m", m);
    r.set("n", (k1 + k2) * m);
    if (s.solve == "k") r.set("target_power", s("power"));
    r.set("power", P(k1));
    r.set("p1", p1);
    r.set("p2", p2);
    r.set("icc", icc);
    r.set("de", de);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c;
    if (s.solve == "k") c.push_back(col("target_power", "Target Power", "prob"));
    c.push_back(col("power", s.solve == "k" ? "Actual Power" : "Power", "prob"));
    for (auto x : {col("k1", "k1", "int"), col("k2", "k2", "int"), col("m", "m", "num"), col("n", "N", "int"),
                   col("p1", "P1", "num"), col("p2", "P2", "num"), col("icc", "ICC", "num"),
                   col("de", "Design Effect", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return fint(r.num("k1")) + " and " + fint(r.num("k2")) + " clusters of " + fnum(r.num("m")) +
           (r.num("m") == 1 ? " subject (" : " subjects (") +
           fint(r.num("n")) + " in total) achieve " + power_text(r.num("power")) + " power" +
           (std::isfinite(r.num("target_power")) ? " (target " + power_text(r.num("target_power")) + ")" : "") +
           " to detect P1 = " +
           fnum(r.num("p1")) + " versus P2 = " + fnum(r.num("p2")) + " with ICC " + fnum(r.num("icc")) +
           " (design effect " + fnum(r.num("de")) + ") using a " + alt_words(alt(s)) + " z-test at α = " +
           fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"k1, k2", "Clusters per group."}, {"m", "Subjects per cluster."},
                   {"Design Effect", "1 + (m − 1)ρ."}};
  p.references = {"Donner, A. & Klar, N. (2000). Design and Analysis of Cluster Randomization Trials in Health "
                  "Research. Arnold.",
                  "Hayes, R. J. & Moulton, L. H. (2017). Cluster Randomised Trials, 2nd ed. CRC Press."};
  p.curve_param = "k1";
  p.curve_label = "k1 (Clusters in Group 1)";
  p.curve_min = 1;
  p.example_solve = "k";
  p.example = {{"power", "0.8 0.9"}, {"m", "20 50 100"}, {"icc", "0.02"}, {"p1", "0.3"}, {"p2", "0.2"}};
  return p;
}
SSC_REGISTER(cluster_props);

}  // namespace
}  // namespace ssc::procs
