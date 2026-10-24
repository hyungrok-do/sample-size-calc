// One proportion: tests (exact binomial and z-tests, exact or approximate
// power) and confidence-interval precision.
#include "procedures/apart.hpp"
#include "procedures/binomial_tools.hpp"
#include "procedures/common.hpp"
#include "procedures/proportion_ci.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// =====================================================================================
// Tests for one proportion
// =====================================================================================
// Every supported statistic is monotone in the count x, so each rejection
// region is a lower tail {x ≤ cL} and/or an upper tail {x ≥ cU}.
struct Region {
  int cl = -1;  // reject if x ≤ cl
  int cu;       // reject if x ≥ cu
};

double z_stat(const std::string& test, int x, int n, double p0) {
  const double ph = static_cast<double>(x) / n;
  double num = ph - p0;
  if (test == "zp0cc" || test == "zphatcc") {
    const double a = std::fabs(num) - 0.5 / n;
    num = a <= 0 ? 0.0 : std::copysign(a, num);
  }
  double se;
  if (test == "zp0" || test == "zp0cc")
    se = std::sqrt(p0 * (1 - p0) / n);
  else
    se = std::sqrt(ph * (1 - ph) / n);
  if (se == 0) return num == 0 ? 0.0 : std::copysign(INFINITY, num);
  return num / se;
}

Region one_prop_region(const std::string& test, int n, double p0, double alpha, Alt a) {
  Region r;
  r.cu = n + 1;
  const double side_alpha = tail_alpha(alpha, a);
  if (test == "exact") {
    if (a != Alt::Lower) r.cu = upper_critical(n, p0, side_alpha);
    if (a != Alt::Upper) r.cl = lower_critical(n, p0, side_alpha);
    return r;
  }
  const double zc = stats::qnorm(side_alpha, false);
  if (a != Alt::Lower) {  // smallest x with z > zc (z is non-decreasing in x)
    int lo = -1, hi = n + 1;
    while (hi - lo > 1) {
      const int mid = lo + (hi - lo) / 2;
      if (z_stat(test, mid, n, p0) > zc)
        hi = mid;
      else
        lo = mid;
    }
    r.cu = hi;
  }
  if (a != Alt::Upper) {  // largest x with z < −zc
    int lo = -1, hi = n + 1;
    while (hi - lo > 1) {
      const int mid = lo + (hi - lo) / 2;
      if (z_stat(test, mid, n, p0) < -zc)
        lo = mid;
      else
        hi = mid;
    }
    r.cl = lo;
  }
  return r;
}

double one_prop_power_normal(const std::string& test, double n, double p0, double p1, double alpha, Alt a) {
  using namespace stats;
  const double q0 = 1 - p0, q1 = 1 - p1;
  const double cc = (test == "zp0cc" || test == "zphatcc") ? 0.5 / n : 0;
  const bool score = (test == "zp0" || test == "zp0cc" || test == "exact");
  const double s0 = std::sqrt((score ? p0 * q0 : p1 * q1) / n);
  const double s1 = std::sqrt(p1 * q1 / n);
  const double d = p1 - p0;
  auto upper = [&](double z) { return pnorm((d - cc - z * s0) / s1); };
  auto lower = [&](double z) { return pnorm((-d - cc - z * s0) / s1); };
  if (a == Alt::Two) {
    const double z = qnorm(half_alpha(alpha), false);
    return std::min(1.0, upper(z) + lower(z));
  }
  const double z = qnorm(alpha, false);
  return a == Alt::Upper ? upper(z) : lower(z);
}

struct OnePropResult {
  double power, actual_alpha;
  Region region;
};

OnePropResult one_prop_eval(const Scenario& s, int n, double p1) {
  const std::string& test = s.s("test");
  const double p0 = s("p0"), alpha = s("alpha");
  const Alt a = alt(s);
  OnePropResult res;
  if (test != "exact" && s.is("method", "normal")) {
    res.power = one_prop_power_normal(test, n, p0, p1, alpha, a);
    res.actual_alpha = NAN;
    res.region = Region{};
    return res;
  }
  res.region = one_prop_region(test, n, p0, alpha, a);
  res.power = two_tail_prob(n, p1, res.region.cl, res.region.cu);
  res.actual_alpha = two_tail_prob(n, p0, res.region.cl, res.region.cu);
  return res;
}

// non-breaking spaces keep each bound on one line in reports
std::string region_text(const Region& r, int n) {
  std::string t;
  if (r.cl >= 0) t += "X\u00a0≤\u00a0" + std::to_string(r.cl);
  if (r.cu <= n) t += std::string(t.empty() ? "" : " or ") + "X\u00a0≥\u00a0" + std::to_string(r.cu);
  return t.empty() ? "never rejects" : t;
}

Procedure one_prop() {
  Procedure p;
  p.id = "prop.one";
  p.name = "Tests for One Proportion";
  p.category = "Proportions";
  p.subcategory = "One Proportion";
  p.order = 1;
  p.summary =
      "Sample size and power for testing a single binomial proportion against a hypothesized value with the exact "
      "binomial test or z-tests (score or Wald, with or without continuity correction). Power is computed exactly by "
      "enumerating the binomial distribution, which also gives the actual significance level.";
  p.keywords = {"binomial", "single proportion", "exact test", "response rate", "z-test", "phase ii", "one-sample proportion", "success rate", "historical control"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}, {"p1", "Effect Size (P1)", "p1", "Detectable P1"}};
  p.params = {
      alternative("Two-sided (H1: P ≠ P0)", "One-sided (H1: P > P0)", "One-sided (H1: P < P0)"),
      select("test", "Test Type",
             {{"exact", "Exact binomial test", {}},
              {"zp0", "Z-test using P0 (score)", {}},
              {"zphat", "Z-test using P̂ (Wald)", {}},
              {"zp0cc", "Z-test using P0 with continuity correction", {}},
              {"zphatcc", "Z-test using P̂ with continuity correction", {}}},
             "exact")
          .grp(G_TEST)
          .hint("The test of the analysis. The exact binomial test never exceeds alpha; two-sided, it puts at "
                "most α/2 in each tail (equal-tailed, unlike R's binom.test). The z-tests are large-sample "
                "approximations whose standard error uses P0 (score) or the observed proportion (Wald)."),
      select("method", "Power Calculation",
             {{"exact", "Exact (binomial enumeration)", {}}, {"normal", "Normal approximation", {}}}, "exact")
          .grp(G_TEST)
          .hint("Binomial enumeration gives the exact power and the actual alpha of the chosen test.")
          .show(when(is_not("test", {"exact"}))),
      power_param(),
      alpha_param(),
      n_param("n", "N (Sample Size)", "50", 1).sym("N").hint("Number of subjects whose yes-or-no outcome is observed."),
      dropout_param(),
      num("p0", "P0 (Proportion under H0)", "0.5")
          .grp(G_EFFECT)
          .sym("P0")
          .range(0, 1, true, true)
          .sugg({"0.5", "0.2", "0.05"})
          .hint("The proportion assumed by the null hypothesis, between 0 and 1, such as a historical or standard "
                "response rate."),
      num("p1", "P1 (Proportion under H1)", "0.7")
          .grp(G_EFFECT)
          .sym("P1")
          .range(0, 1, true, true)
          .sugg({"0.7", "0.6 0.7 0.8"})
          .hint("The true proportion you expect, between 0 and 1; the study is sized to tell it apart from P0.")
          .show(when(not_solving({"p1"}))),
      select("p1dir", "Search Direction for P1", {{"above", "P1 > P0", {}}, {"below", "P1 < P0", {}}}, "above")
          .grp(G_EFFECT)
          .hint("A two-sided test can detect a P1 above or below P0; choose which of the two to find.")
          .show(when({solving({"p1"}), is("alternative", {"two"})})),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double p0 = s("p0");
    const Alt a = alt(s);
    const bool exact_power = s.is("test", "exact") || s.is("method", "exact");
    int n;
    double p1;
    if (s.solve == "n") {
      p1 = s("p1");
      if (p1 == p0) throw CalcError("P1 equals P0");
      if ((a == Alt::Upper && p1 < p0) || (a == Alt::Lower && p1 > p0))
        throw CalcError("P1 lies on the side opposite to the one-sided alternative");
      const double zf = std::sqrt(z_factor(s("alpha"), s("power"), a));
      const double za = stats::qnorm(tail_alpha(s("alpha"), a), false), zb = zf - za;
      const double g = std::pow((za * std::sqrt(p0 * (1 - p0)) + zb * std::sqrt(p1 * (1 - p1))) / (p1 - p0), 2);
      solve::NSearch o;
      o.lo = 1;
      o.guess = std::ceil(g);
      o.monotone = !exact_power;
      if (exact_power) o.period = 1 / std::fmin(p0, 1 - p0);
      o.hi = 5e6;
      n = static_cast<int>(solve::min_n([&](double x) { return one_prop_eval(s, static_cast<int>(x), p1).power; },
                                        s("power"), o, &ctx));
    } else {
      n = static_cast<int>(s("n"));
      if (s.solve == "power") {
        p1 = s("p1");
      } else {
        const bool up = a == Alt::Upper || (a == Alt::Two && !s.is("p1dir", "below"));
        auto f = [&](double x) { return one_prop_eval(s, n, x).power; };
        const double target = s("power");
        if (f(p0) >= target)
          throw CalcError("the target power must exceed the power at P1 = P0 (the actual significance level)");
        if (up) {
          if (f(1 - 1e-12) < target) throw CalcError("the target power cannot be reached with this N for any P1 > P0");
          p1 = numerics::brent([&](double x) { return f(x) - target; }, p0, 1 - 1e-12, 1e-13);
        } else {
          if (f(1e-12) < target) throw CalcError("the target power cannot be reached with this N for any P1 < P0");
          p1 = numerics::brent([&](double x) { return f(x) - target; }, 1e-12, p0, 1e-13);
        }
      }
    }
    const auto res = one_prop_eval(s, n, p1);
    r.set("n", n);
    put_power(s, r, res.power);
    r.set("p0", p0);
    r.set("p1", p1);
    r.set("diff", p1 - p0);
    r.set("alpha", s("alpha"));
    if (exact_power) {
      r.set("actual_alpha", res.actual_alpha);
      r.set("region", region_text(res.region, n));
    }
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("p0", "P0", "num"), col("p1", "P1", "num"),
                   col("diff", "P1 − P0", "num"), col("alpha", "Target Alpha", "alpha")})
      c.push_back(x);
    if (s.is("test", "exact") || s.is("method", "exact")) {
      c.push_back(col("actual_alpha", "Actual Alpha", "prob", "Exact size of the test under P0"));
      c.push_back(col("region", "Reject H0 if", "text", "Rejection region for the number of successes X"));
    }
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    static const std::map<std::string, std::string> names = {
        {"exact", "exact binomial test"}, {"zp0", "z-test using P0"}, {"zphat", "z-test using P̂"},
        {"zp0cc", "continuity-corrected z-test using P0"}, {"zphatcc", "continuity-corrected z-test using P̂"}};
    const auto [p0, p1] = fnum_apart(r.num("p0"), r.num("p1"));
    std::string t = "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) +
                    " to detect a difference of " + fnum(r.num("diff")) + " between the null proportion P0 = " + p0 +
                    " and the actual proportion P1 = " + p1 + " using a " + alt_words(alt(s)) + " " +
                    names.at(s.s("test")) + " with a target significance level of " + fnum(r.num("alpha"));
    if (r.has("actual_alpha")) t += " (actual significance level " + fprob(r.num("actual_alpha")) + ")";
    return t + ".";
  };
  p.definitions = {{"N", "Number of subjects (Bernoulli trials)."},
                   {"P0", "Proportion under the null hypothesis."},
                   {"P1", "Actual proportion at which power is evaluated."},
                   {"Target Alpha", "Nominal significance level."},
                   {"Actual Alpha", "Exact probability of rejecting H0 when P = P0 (discreteness makes it ≤ target)."},
                   {"Reject H0 if", "Rejection region for the observed number of successes X."}};
  p.notes = {"Because the binomial distribution is discrete, exact power does not increase monotonically with N; "
             "the reported N is the smallest sample size whose power reaches the target."};
  p.references = {REF_CHOW, REF_FLEISS, REF_MACHIN,
                  "Chernick, M. R. & Liu, C. Y. (2002). The saw-toothed behavior of power versus sample size and "
                  "software solutions. The American Statistician 56, 149–155."};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.curve_min = 1;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"p0", "0.5"}, {"p1", "0.65 0.7 0.75"}};
  return p;
}
SSC_REGISTER(one_prop);

// =====================================================================================
// Confidence interval for one proportion
// =====================================================================================
double ci_width(const Scenario& s, double n) {
  const double p = s("p");
  const bool two = s.is("sides", "two");
  const Interval ci = prop_ci(s.s("method"), n, n * p, s("conf"), two);
  // limits that meet (at a confidence near 0) can cross by a rounding error
  if (two) return ci.hi < ci.lo ? 0.0 : ci.hi - ci.lo;
  return s.is("sides", "upper") ? ci.hi - p : p - ci.lo;
}

// A one-sided bound's distance from the estimate p within the rounding of the bound, as at a confidence a rounding
// error above 50%: its sign then means nothing.
bool near_estimate(double w, double p) { return std::fabs(w) <= 1e-13 * p; }

// Why a Jeffreys bound can lie on the wrong side of the estimate (a distance below 0), as at a confidence not far
// above 50% for a proportion away from 1/2; the other intervals contain the estimate.
const char* const JEFFREYS_WRONG_SIDE =
    "the Jeffreys bound is a quantile of the posterior distribution, whose middle lies between the estimate and 1/2, "
    "further from the estimate than the bound's margin reaches at this confidence";

Procedure one_prop_ci() {
  Procedure p;
  p.id = "prop.one.ci";
  p.name = "Confidence Interval for One Proportion";
  p.category = "Confidence Intervals";
  p.subcategory = "Proportions";
  p.order = 1;
  p.summary =
      "Sample size needed so that a confidence interval for a proportion has a given width (or the width achieved by "
      "a given N), for the Wald, Wilson, Wilson with continuity correction, Agresti–Coull, Jeffreys and exact "
      "Clopper–Pearson intervals, evaluated at the anticipated proportion.";
  p.keywords = {"precision", "margin of error", "clopper-pearson", "wilson", "prevalence", "survey", "estimate a proportion", "agresti-coull", "jeffreys", "single proportion"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"width", "Interval Width", "width", "Width"}};
  p.params = {
      select("method", "Interval Method",
             {{"exact", "Exact (Clopper–Pearson)", {}},
              {"wilson", "Wilson score", {}},
              {"wilsoncc", "Wilson score with continuity correction", {}},
              {"ac", "Agresti–Coull", {}},
              {"jeffreys", "Jeffreys", {}},
              {"wald", "Simple asymptotic (Wald)", {}}},
             "exact")
          .grp(G_TEST)
          .hint("The formula of the interval. Clopper–Pearson never falls short of the stated confidence; Wilson "
                "and Agresti–Coull are shorter and accurate; Wald is the textbook formula, poor near 0 or 1."),
      select("sides", "Interval Type",
             {{"two", "Two-sided", {}}, {"lower", "One-sided lower bound", {}}, {"upper", "One-sided upper bound", {}}},
             "two")
          .grp(G_TEST)
          .hint("A two-sided interval has a lower and an upper limit; a one-sided interval gives only a lower or "
                "only an upper bound."),
      num("conf", "Confidence Level (1 − α)", "0.95")
          .grp(G_TEST)
          .sym("Confidence")
          .range(0, 1, true, true)
          .sugg({"0.95", "0.90 0.95 0.99"})
          .hint("How often intervals computed this way contain the true value; 0.95 gives a 95% interval."),
      n_param("n", "N (Sample Size)", "100", 1)
          .sym("N")
          .hint("Number of subjects whose yes-or-no outcome is observed.")
          .show(when(not_solving({"n"}))),
      num("width", "Width (two-sided: upper − lower; one-sided: distance to bound)", "0.1")
          .grp(G_EFFECT)
          .sym("Width")
          .range(0, 1, true, true)
          .sugg({"0.1", "0.05 0.1 0.15", "0.2"})
          .hint("The precision you want, on the proportion scale (0.1 = 10 percentage points): the full width of "
                "a two-sided interval, or the distance from the estimate to a one-sided bound.")
          .show(when(solving({"n"}))),
      num("p", "P (Anticipated Proportion)", "0.3")
          .grp(G_EFFECT)
          .sym("P")
          .range(0, 1, true, true)
          .sugg({"0.5", "0.3", "0.1 0.3 0.5"})
          .hint("The proportion you expect to observe. The interval is widest at 0.5, so 0.5 is the safe choice "
                "when you are unsure."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    if (!s.is("sides", "two")) check_one_sided_conf(s("conf"));
    double n;
    if (s.solve == "n") {
      const double target = s("width");
      const bool two = s.is("sides", "two");
      solve::NSearch o;
      o.lo = 1;
      const double z = stats::qnorm((1 - s("conf")) / (two ? 2 : 1), false);
      const double pp = s("p");
      o.guess = std::ceil((two ? 4 : 1) * z * z * pp * (1 - pp) / (target * target));
      o.hi = 1e8;
      o.goal = "the target width";
      // a one-sided bound at the estimate, on its wrong side or within a rounding error of it is a degenerate size;
      // a two-sided interval at a confidence near 0 closes on the estimate itself (a width of 0, or a rounding error
      // either side of it) at every size
      auto ok = [&](double x) {
        const double w = ci_width(s, x);
        return (two || (w > 0 && !near_estimate(w, pp))) && w <= target;
      };
      try {
        n = solve::min_n_staying(ok, o, &ctx);
      } catch (const CalcError&) {
        // a one-sided bound that is still on the wrong side of the estimate, or within a rounding error of it, at
        // the largest size: the reason, not only that the width was not reached
        const double w = ci_width(s, o.hi);
        if (two || (w > 0 && !near_estimate(w, pp))) throw;
        if (near_estimate(w, pp))
          throw CalcError("the confidence is so close to 50% that the bound lies within a rounding error of the "
                          "estimate, where its distance cannot be computed");
        throw CalcError(std::string("the bound lies on the wrong side of the estimate for every sample size up to ") +
                        text::format_g(o.hi, 10) + ": " +
                        (s.is("method", "jeffreys") ? JEFFREYS_WRONG_SIDE : "at this confidence"));
      }
      if (!two && n > o.lo && s.is("method", "jeffreys") &&
          (ci_width(s, n - 1) <= 0 || near_estimate(ci_width(s, n - 1), pp)))
        r.notes.push_back("Smaller samples put the bound on the wrong side of the estimate (" +
                          std::string(JEFFREYS_WRONG_SIDE) + "); this is the smallest that puts it beyond.");
      r.set("target_width", target);
    } else {
      n = s("n");
    }
    const double pp = s("p");
    const Interval ci = prop_ci(s.s("method"), n, n * pp, s("conf"), s.is("sides", "two"));
    r.set("n", n);
    r.set("width", ci_width(s, n));
    r.set("p", pp);
    r.set("lower", s.is("sides", "upper") ? NAN : ci.lo);
    r.set("upper", s.is("sides", "lower") ? NAN : ci.hi);
    r.set("conf", s("conf"));
  };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c;
    c.push_back(col("conf", "Confidence", "num"));
    c.push_back(col("n", "N", "int"));
    if (s.solve == "n") c.push_back(col("target_width", "Target Width", "num"));
    c.push_back(col("width", s.solve == "n" ? "Actual Width" : "Width", "num"));
    c.push_back(col("p", "P", "num", "Anticipated proportion"));
    c.push_back(col("lower", "Lower Limit", "num"));
    c.push_back(col("upper", "Upper Limit", "num"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    static const std::map<std::string, std::string> names = {
        {"exact", "exact (Clopper–Pearson)"}, {"wilson", "Wilson score"}, {"wilsoncc", "continuity-corrected Wilson"},
        {"ac", "Agresti–Coull"}, {"jeffreys", "Jeffreys"}, {"wald", "Wald"}};
    const bool two = s.is("sides", "two");
    const double w = r.num("width");
    std::string where = (two ? "with a width of " : "whose bound lies at a distance of ") + fnum_fine(w);
    if (!two && near_estimate(w, r.num("p")))
      where = "whose bound lies within a rounding error of the estimate";
    else if (!two && w < 0)
      where = "whose bound lies " + fnum_fine(-w) + (s.is("sides", "upper") ? " below" : " above") +
              " the estimate, on the wrong side of it (" +
              (s.is("method", "jeffreys") ? JEFFREYS_WRONG_SIDE : "at this confidence") + ")";
    return "A sample of " + fint(r.num("n")) + " gives a " + (two ? "two-sided " : "one-sided ") +
           level_text(r.num("conf")) + " " + names.at(s.s("method")) + " confidence interval " + where +
           " when the observed proportion equals " + fnum_fine(r.num("p")) + ".";
  };
  p.definitions = {{"N", "Sample size."},
                   {"Width", "Upper minus lower limit (two-sided) or distance from P to the bound (one-sided)."},
                   {"P", "Anticipated sample proportion at which the interval is evaluated."}};
  p.references = {"Newcombe, R. G. (1998). Two-sided confidence intervals for the single proportion: comparison of "
                  "seven methods. Statistics in Medicine 17, 857–872.",
                  "Clopper, C. J. & Pearson, E. S. (1934). The use of confidence or fiducial limits illustrated in the "
                  "case of the binomial. Biometrika 26, 404–413.",
                  "Agresti, A. & Coull, B. A. (1998). Approximate is better than 'exact' for interval estimation of "
                  "binomial proportions. The American Statistician 52, 119–126.",
                  REF_FLEISS};
  p.example_solve = "n";
  p.example = {{"width", "0.05 0.1 0.15"}, {"p", "0.1 0.3 0.5"}};
  return p;
}
SSC_REGISTER(one_prop_ci);

}  // namespace
}  // namespace ssc::procs
