// Two independent proportions: z-tests (pooled / unpooled / continuity
// corrected) and Fisher's exact test, with exact power by enumeration of the
// two binomial distributions (restricted to the region carrying all but
// 1e-16 of the probability) or with normal approximations.
#include <atomic>
#include <mutex>

#include "procedures/binomial_tools.hpp"
#include "procedures/common.hpp"
#include "procedures/proportion_ci.hpp"
#include "util/parallel.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

struct TwoPropTest {
  std::string test;  // zpool, zunpool, zpoolcc, fisher
  Alt alt;
  double alpha;
};

double z_two(const std::string& test, int x1, int n1, int x2, int n2) {
  const double p1 = static_cast<double>(x1) / n1, p2 = static_cast<double>(x2) / n2;
  double d = p1 - p2;
  if (test == "zpoolcc") {
    const double a = std::fabs(d) - 0.5 * (1.0 / n1 + 1.0 / n2);
    d = a <= 0 ? 0.0 : std::copysign(a, d);
  }
  double se;
  if (test == "zunpool") {
    se = std::sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2);
    if (se == 0) return d == 0 ? 0.0 : std::copysign(INFINITY, d);
  } else {
    const double pb = static_cast<double>(x1 + x2) / (n1 + n2);
    se = std::sqrt(pb * (1 - pb) * (1.0 / n1 + 1.0 / n2));
    if (se == 0) return 0.0;
  }
  return d / se;
}

// Fisher's exact test: rejection flags for x1 = lo..hi given total t.
// Two-sided p-value = sum of conditional probabilities ≤ that of the observed table.
void fisher_region(int n1, int n2, int t, double alpha, Alt alt, int& lo, int& hi, std::vector<char>& rej,
                   std::vector<double>& f) {
  const int slo = std::max(0, t - n2), shi = std::min(n1, t);
  // window around the mode where the null pmf is non-negligible
  const double N = n1 + n2;
  const double mean = static_cast<double>(t) * n1 / N;
  const double var = N > 1 ? mean * (N - t) / N * (N - n1) / (N - 1) : 0;
  const double half = 40 * std::sqrt(var) + 20;
  lo = std::max(slo, static_cast<int>(std::floor(mean - half)));
  hi = std::min(shi, static_cast<int>(std::ceil(mean + half)));
  const int m = hi - lo + 1;
  f.assign(static_cast<std::size_t>(m), 0.0);
  rej.assign(static_cast<std::size_t>(m), 0);
  int mode = static_cast<int>(std::floor((t + 1.0) * (n1 + 1.0) / (N + 2)));
  mode = std::min(std::max(mode, lo), hi);
  f[static_cast<std::size_t>(mode - lo)] = stats::dhyper(mode, n1, n2, t);
  for (int x = mode; x < hi; ++x)
    f[static_cast<std::size_t>(x + 1 - lo)] =
        f[static_cast<std::size_t>(x - lo)] * (static_cast<double>(n1 - x) * (t - x)) /
        (static_cast<double>(x + 1) * (n2 - t + x + 1));
  for (int x = mode; x > lo; --x)
    f[static_cast<std::size_t>(x - 1 - lo)] = f[static_cast<std::size_t>(x - lo)] *
                                              (static_cast<double>(x) * (n2 - t + x)) /
                                              (static_cast<double>(n1 - x + 1) * (t - x + 1));
  const double fuzz = 1 + 1e-7;
  if (alt == Alt::Upper) {  // reject large x1: P(X ≥ x) ≤ α
    double cum = 0;
    for (int i = m - 1; i >= 0; --i) {
      cum += f[static_cast<std::size_t>(i)];
      if (cum <= alpha * fuzz) rej[static_cast<std::size_t>(i)] = 1; else break;
    }
    return;
  }
  if (alt == Alt::Lower) {
    double cum = 0;
    for (int i = 0; i < m; ++i) {
      cum += f[static_cast<std::size_t>(i)];
      if (cum <= alpha * fuzz) rej[static_cast<std::size_t>(i)] = 1; else break;
    }
    return;
  }
  // two-sided: walk inwards from both ends in increasing order of probability
  int i = 0, j = m - 1;
  double cum = 0;
  while (i <= j) {
    const double v = std::min(f[static_cast<std::size_t>(i)], f[static_cast<std::size_t>(j)]);
    const double thr = v * fuzz;
    int gi = i, gj = j;
    double add = 0;
    while (gi <= gj && f[static_cast<std::size_t>(gi)] <= thr) add += f[static_cast<std::size_t>(gi++)];
    while (gj >= gi && f[static_cast<std::size_t>(gj)] <= thr) add += f[static_cast<std::size_t>(gj--)];
    cum += add;
    if (cum > alpha * fuzz) break;
    for (int k = i; k < gi; ++k) rej[static_cast<std::size_t>(k)] = 1;
    for (int k = gj + 1; k <= j; ++k) rej[static_cast<std::size_t>(k)] = 1;
    i = gi;
    j = gj;
  }
}

double two_prop_normal(const TwoPropTest& T, double n1, double n2, double p1, double p2) {
  using namespace stats;
  const double q1 = 1 - p1, q2 = 1 - p2;
  const double pb = (n1 * p1 + n2 * p2) / (n1 + n2), qb = 1 - pb;
  const double s1 = std::sqrt(p1 * q1 / n1 + p2 * q2 / n2);
  const double s0 = T.test == "zunpool" ? s1 : std::sqrt(pb * qb * (1 / n1 + 1 / n2));
  const double cc = (T.test == "zpoolcc" || T.test == "fisher") ? 0.5 * (1 / n1 + 1 / n2) : 0;
  const double d = p1 - p2;
  auto up = [&](double z) { return pnorm((d - cc - z * s0) / s1); };
  auto lo = [&](double z) { return pnorm((-d - cc - z * s0) / s1); };
  if (T.alt == Alt::Two) {
    const double z = qnorm(T.alpha / 2, false);
    return std::min(1.0, up(z) + lo(z));
  }
  const double z = qnorm(T.alpha, false);
  return T.alt == Alt::Upper ? up(z) : lo(z);
}

constexpr double MAX_CELLS = 6e7;

// Exact power by enumeration over the joint binomial window.
double two_prop_exact(const TwoPropTest& T, int n1, int n2, double p1, double p2, const RunContext& ctx) {
  const BinomWindow w1 = binom_window(n1, p1), w2 = binom_window(n2, p2);
  const double cells = static_cast<double>(w1.hi - w1.lo + 1) * (w2.hi - w2.lo + 1);
  if (cells > MAX_CELLS) return NAN;
  const bool par = cells > 2e5 && ctx.threads > 1;  // nested in the loop over scenarios: the idle workers help
  if (T.test == "fisher") {
    const int tlo = w1.lo + w2.lo, thi = w1.hi + w2.hi;
    const std::size_t nt = static_cast<std::size_t>(thi - tlo + 1);
    std::vector<double> part(nt, 0.0);
    auto body = [&](std::size_t k) {
      const int t = tlo + static_cast<int>(k);
      int lo, hi;
      std::vector<char> rej;
      std::vector<double> f;
      fisher_region(n1, n2, t, T.alpha, T.alt, lo, hi, rej, f);
      double s = 0;
      const int a = std::max({lo, w1.lo, t - w2.hi}), b = std::min({hi, w1.hi, t - w2.lo});
      for (int x1 = a; x1 <= b; ++x1)
        if (rej[static_cast<std::size_t>(x1 - lo)]) s += w1.at(x1) * w2.at(t - x1);
      // tables beyond the window (more than 40 SD from the null mean) are more extreme than any inside it: when the
      // window's edge on that side is rejected, so are they (with a large effect they carry nearly all the power)
      const int A = std::max(w1.lo, t - w2.hi), B = std::min(w1.hi, t - w2.lo);
      if (!rej.empty() && rej.back())
        for (int x1 = std::max(A, hi + 1); x1 <= B; ++x1) s += w1.at(x1) * w2.at(t - x1);
      if (!rej.empty() && rej.front())
        for (int x1 = A; x1 <= std::min(B, lo - 1); ++x1) s += w1.at(x1) * w2.at(t - x1);
      part[k] = s;
    };
    if (par)
      parallel_for(nt, ctx.threads, body, ctx.token);
    else
      for (std::size_t k = 0; k < nt; ++k) body(k);
    double s = 0;
    for (double v : part) s += v;
    return stats::clamp01(s);
  }
  const double zc = stats::qnorm(T.alt == Alt::Two ? T.alpha / 2 : T.alpha, false);
  const std::size_t nr = static_cast<std::size_t>(w1.hi - w1.lo + 1);
  std::vector<double> part(nr, 0.0);
  auto body = [&](std::size_t k) {
    const int x1 = w1.lo + static_cast<int>(k);
    const double b1 = w1.pmf[k];
    double s = 0;
    for (int x2 = w2.lo; x2 <= w2.hi; ++x2) {
      const double z = z_two(T.test, x1, n1, x2, n2);
      const bool rj = T.alt == Alt::Two ? std::fabs(z) > zc : (T.alt == Alt::Upper ? z > zc : z < -zc);
      if (rj) s += w2.pmf[static_cast<std::size_t>(x2 - w2.lo)];
    }
    part[k] = b1 * s;
  };
  if (par)
    parallel_for(nr, ctx.threads, body, ctx.token);
  else
    for (std::size_t k = 0; k < nr; ++k) body(k);
  double s = 0;
  for (double v : part) s += v;
  return stats::clamp01(s);
}

double p1_from_effect(const Scenario& s, double p2) {
  const std::string& es = s.s("es");
  double p1;
  if (es == "diff")
    p1 = p2 + s("d1");
  else if (es == "ratio")
    p1 = p2 * s("r1");
  else if (es == "or") {
    const double o = s("or1") * p2 / (1 - p2);
    p1 = o / (1 + o);
  } else
    p1 = s("p1");
  if (!(p1 > 0 && p1 < 1)) throw CalcError("the implied P1 = " + fnum(p1) + " is not between 0 and 1");
  return p1;
}

Procedure two_props() {
  Procedure p;
  p.id = "prop.two";
  p.name = "Tests for Two Proportions";
  p.category = "Proportions";
  p.subcategory = "Two Independent Proportions";
  p.order = 1;
  p.summary =
      "Sample size and power for comparing two independent proportions with z-tests (pooled, unpooled, continuity "
      "corrected; the pooled z-test equals the chi-square test) or Fisher's exact test. The effect can be given as "
      "proportions, a difference, a ratio or an odds ratio. Exact power enumerates both binomial distributions.";
  p.keywords = {"two proportions", "chi-square", "fisher exact", "risk difference", "relative risk", "odds ratio",
                "binary outcome", "response rate", "risk ratio", "binary endpoint", "percentages", "difference in proportions"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"},
                 {"power", "Power", "power", "Power"},
                 {"p1", "Effect Size (P1)", "p1", "Detectable P1"}};
  p.params = {
      alternative("Two-sided (H1: P1 ≠ P2)", "One-sided (H1: P1 > P2)", "One-sided (H1: P1 < P2)"),
      select("test", "Test Type",
             {{"zpool", "Z-test, pooled variance (= chi-square)", {}},
              {"zunpool", "Z-test, unpooled variance", {}},
              {"zpoolcc", "Z-test with continuity correction (Yates)", {}},
              {"fisher", "Fisher's exact test", {}}},
             "zpool")
          .grp(G_TEST)
          .hint("The test of the analysis. The pooled z-test is the usual chi-square test; the unpooled one uses "
                "separate variances; Yates' correction and Fisher's exact test are more conservative."),
      select("method", "Power Calculation",
             {{"exact", "Exact (binomial enumeration)", {}}, {"normal", "Normal approximation", {}}}, "exact")
          .grp(G_TEST)
          .hint("Enumeration gives the exact power of the chosen test; the normal approximation is the classical "
                "closed-form formula."),
      power_param(),
      alpha_param()};
  for (auto& q : two_group_params(2, "50")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  for (auto q :
       {select("es", "Effect Size Input",
               {{"props", "Proportions (P1, P2)", {}},
                {"diff", "Difference D1 = P1 − P2", {}},
                {"ratio", "Ratio R1 = P1/P2", {}},
                {"or", "Odds ratio OR1", {}}},
               "props")
            .grp(G_EFFECT)
            .hint("How you describe the effect: give P1 directly, or its difference from, ratio to or odds ratio "
                  "against P2.")
            .show(when(not_solving({"p1"}))),
        num("p1", "P1 (Group 1 Proportion under H1)", "0.6")
            .grp(G_EFFECT)
            .sym("P1")
            .range(0, 1, true, true)
            .sugg({"0.6", "0.5 0.6 0.7"})
            .hint("The share of subjects in group 1 (e.g. treatment) expected to have the outcome under the "
                  "alternative, between 0 and 1.")
            .show(when({not_solving({"p1"}), is("es", {"props"})})),
        num("d1", "D1 (Difference P1 − P2)", "0.2")
            .grp(G_EFFECT)
            .sym("D1")
            .range(-1, 1, true, true)
            .sugg({"0.1", "0.2", "0.1 0.15 0.2"})
            .hint("Expected difference between the group proportions (0.1 = 10 percentage points); P1 = P2 + D1.")
            .show(when({not_solving({"p1"}), is("es", {"diff"})})),
        num("r1", "R1 (Ratio P1/P2)", "1.5")
            .grp(G_EFFECT)
            .sym("R1")
            .above(0)
            .sugg({"1.5", "1.25 1.5 2"})
            .hint("Expected ratio of the group proportions (the relative risk); P1 = R1 × P2.")
            .show(when({not_solving({"p1"}), is("es", {"ratio"})})),
        num("or1", "OR1 (Odds Ratio)", "2")
            .grp(G_EFFECT)
            .sym("OR1")
            .above(0)
            .sugg({"2", "1.5 2 3"})
            .hint("Expected odds ratio of group 1 to group 2: the odds P1/(1 − P1) divided by the odds "
                  "P2/(1 − P2).")
            .show(when({not_solving({"p1"}), is("es", {"or"})})),
        num("p2", "P2 (Group 2 / Reference Proportion)", "0.4")
            .grp(G_EFFECT)
            .sym("P2")
            .range(0, 1, true, true)
            .sugg({"0.4", "0.1 0.2 0.3"})
            .hint("The share of subjects in group 2 (the control or reference group) expected to have the "
                  "outcome, between 0 and 1."),
        select("p1dir", "Search Direction for P1", {{"above", "P1 > P2", {}}, {"below", "P1 < P2", {}}}, "above")
            .grp(G_EFFECT)
            .hint("A two-sided test can detect a P1 above or below P2; choose which of the two to find.")
            .show(when({solving({"p1"}), is("alternative", {"two"})}))})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double p2 = s("p2");
    const Alt a = alt(s);
    const TwoPropTest T{s.s("test"), a, s("alpha")};
    bool exact = s.is("method", "exact");
    bool fell_back = false;
    auto power = [&](double n1, double n2, double x) {
      if (exact) {
        const double pw = two_prop_exact(T, static_cast<int>(n1), static_cast<int>(n2), x, p2, ctx);
        if (!std::isnan(pw)) return pw;
        fell_back = true;
      }
      return two_prop_normal(T, n1, n2, x, p2);
    };
    TwoN n;
    double p1;
    if (s.solve == "n") {
      p1 = p1_from_effect(s, p2);
      if (p1 == p2) throw CalcError("P1 equals P2");
      if ((a == Alt::Upper && p1 < p2) || (a == Alt::Lower && p1 > p2))
        throw CalcError("P1 − P2 has the sign opposite to the one-sided alternative");
      const double za = stats::qnorm(a == Alt::Two ? s("alpha") / 2 : s("alpha"), false), zb = stats::qnorm(s("power"));
      const double pb = (p1 + p2) / 2;
      const double g = std::pow(za * std::sqrt(2 * pb * (1 - pb)) + zb * std::sqrt(p1 * (1 - p1) + p2 * (1 - p2)), 2) /
                       ((p1 - p2) * (p1 - p2));
      if (exact && exact_search_work(g, g, p1, p2) > EXACT_SEARCH_BUDGET) exact = false, fell_back = true;
      n = solve_two_group_n(
          s, [&](double a1, double a2) { return power(a1, a2, p1); }, s("power"), 2, std::ceil(g), !exact, ctx, 1e7);
    } else {
      n = two_group_n(s);
      if (s.solve == "power") {
        p1 = p1_from_effect(s, p2);
      } else {
        const bool up = a == Alt::Upper || (a == Alt::Two && !s.is("p1dir", "below"));
        auto f = [&](double x) { return power(n.n1, n.n2, x) - s("power"); };
        const double lo = up ? p2 : 1e-9, hi = up ? 1 - 1e-9 : p2;
        if (f(p2) >= 0) throw CalcError("the target power must exceed the power at P1 = P2 (the significance level)");
        if (f(up ? hi : lo) < 0) throw CalcError("the target power cannot be reached with these sample sizes");
        p1 = numerics::brent(f, lo, hi, 1e-11);
      }
    }
    put_two_group(r, n);
    put_power(s, r, power(n.n1, n.n2, p1));
    r.set("p1", p1);
    r.set("p2", p2);
    r.set("diff", p1 - p2);
    r.set("ratio", p1 / p2);
    r.set("or", p1 / (1 - p1) / (p2 / (1 - p2)));
    r.set("alpha", s("alpha"));
    if (exact && !fell_back)
      r.set("actual_alpha", two_prop_exact(T, static_cast<int>(n.n1), static_cast<int>(n.n2), p2, p2, ctx));
    if (fell_back) r.notes.push_back("The samples are too large for exact enumeration; the normal approximation was used.");
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("p1", "P1", "num"), col("p2", "P2", "num"), col("diff", "D1 = P1 − P2", "num"),
                   col("ratio", "R1 = P1/P2", "num"), col("or", "OR1", "num"), col("alpha", "Target Alpha", "alpha")})
      c.push_back(x);
    if (s.is("method", "exact")) c.push_back(col("actual_alpha", "Actual Alpha", "prob", "Exact size when P1 = P2"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    static const std::map<std::string, std::string> names = {{"zpool", "pooled z-test (chi-square test)"},
                                                             {"zunpool", "unpooled z-test"},
                                                             {"zpoolcc", "continuity-corrected z-test"},
                                                             {"fisher", "Fisher's exact test"}};
    std::string t = "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
                    " to detect a difference between the group proportions P1 = " + fnum(r.num("p1")) + " and P2 = " +
                    fnum(r.num("p2")) + " (difference " + fnum(r.num("diff")) + ", ratio " + fnum(r.num("ratio"), 3) +
                    ", odds ratio " + fnum(r.num("or"), 3) + ") using a " + alt_words(alt(s)) + " " +
                    names.at(s.s("test")) + " at a target significance level of " + fnum(r.num("alpha"));
    if (r.has("actual_alpha")) t += " (actual " + text::format_fixed(r.num("actual_alpha"), 4) + ")";
    // the exact size is reported only when power was enumerated (large samples use the normal approximation)
    const bool enumerated = s.is("method", "exact") && r.has("actual_alpha");
    return t + (enumerated ? "; power was computed by exact enumeration." : "; power uses the normal approximation.");
  };
  p.definitions = {{"N1, N2", "Sample sizes of groups 1 and 2."},
                   {"P1, P2", "Response proportions in groups 1 (treatment) and 2 (control/reference)."},
                   {"D1, R1, OR1", "Difference, ratio and odds ratio of P1 relative to P2."},
                   {"Actual Alpha", "Exact probability of rejecting H0 when both groups have proportion P2."}};
  p.notes = {"With exact enumeration, power is a saw-toothed function of N; the reported N is the first one that "
             "reaches the target power."};
  p.references = {REF_FLEISS, REF_CHOW, REF_MACHIN,
                  "Casagrande, J. T., Pike, M. C. & Smith, P. G. (1978). An improved approximate formula for "
                  "calculating sample sizes for comparing two binomial distributions. Biometrics 34, 483–486.",
                  "Fisher, R. A. (1935). The logic of inductive inference. Journal of the Royal Statistical Society 98, "
                  "39–54."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"p1", "0.5 0.6"}, {"p2", "0.3"}};
  return p;
}
SSC_REGISTER(two_props);

// =====================================================================================
// Confidence interval for the difference of two proportions
// =====================================================================================
// Limits for P1 − P2 evaluated at the anticipated counts x_i = n_i·p_i:
// Newcombe's hybrid score intervals combine the Wilson limits of each
// proportion (Newcombe 1998, methods 10 and 11); Agresti–Caffo adds one
// success and one failure to each group; Wald is the simple asymptotic interval.
Interval diff_ci(const std::string& method, double n1, double x1, double n2, double x2, double conf, bool two) {
  const double p1 = x1 / n1, p2 = x2 / n2, d = p1 - p2;
  if (method == "newcombe" || method == "newcombecc") {
    const std::string one = method == "newcombe" ? "wilson" : "wilsoncc";
    const Interval c1 = prop_ci(one, n1, x1, conf, two), c2 = prop_ci(one, n2, x2, conf, two);
    return {d - std::hypot(p1 - c1.lo, c2.hi - p2), d + std::hypot(c1.hi - p1, p2 - c2.lo)};
  }
  const double z = stats::qnorm(two ? (1 - conf) / 2 : 1 - conf, false);
  if (method == "ac") {
    const double t1 = (x1 + 1) / (n1 + 2), t2 = (x2 + 1) / (n2 + 2);
    const double h = z * std::sqrt(t1 * (1 - t1) / (n1 + 2) + t2 * (1 - t2) / (n2 + 2));
    return {std::max(-1.0, t1 - t2 - h), std::min(1.0, t1 - t2 + h)};
  }
  double h = z * std::sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2);
  if (method == "waldcc") h += 0.5 * (1 / n1 + 1 / n2);
  return {std::max(-1.0, d - h), std::min(1.0, d + h)};
}

struct DiffCi {
  double n2;
  Interval ci;
  double width;
};

DiffCi diff_ci_at(const Scenario& s, double n1) {
  const double p1 = s("p1"), p2 = s("p2");
  DiffCi c;
  c.n2 = std::ceil(s("r") * n1 - 1e-9);
  c.ci = diff_ci(s.s("method"), n1, n1 * p1, c.n2, c.n2 * p2, s("conf"), s.is("sides", "two"));
  if (s.is("sides", "two"))
    c.width = c.ci.hi - c.ci.lo;
  else
    c.width = s.is("sides", "upper") ? c.ci.hi - (p1 - p2) : (p1 - p2) - c.ci.lo;
  return c;
}

Procedure two_props_ci() {
  Procedure p;
  p.id = "prop.two.ci";
  p.name = "Confidence Interval for the Difference of Two Proportions";
  p.category = "Confidence Intervals";
  p.subcategory = "Proportions";
  p.order = 2;
  p.summary =
      "Group sample sizes needed so that a confidence interval for the difference between two independent "
      "proportions has a given width (or the width achieved by given sample sizes), for Newcombe's hybrid score "
      "intervals with or without continuity correction, the Agresti–Caffo interval and the Wald interval, evaluated "
      "at the anticipated proportions.";
  p.keywords = {"precision", "risk difference", "margin of error", "newcombe", "agresti-caffo", "two groups", "difference in proportions", "two proportions"};
  p.solve_for = {{"n", "Sample Size", "n1", "N1"}, {"width", "Interval Width", "width", "Width"}};
  p.params = {
      select("method", "Interval Method",
             {{"newcombe", "Newcombe hybrid score (Wilson)", {}},
              {"newcombecc", "Newcombe hybrid score with continuity correction", {}},
              {"ac", "Agresti–Caffo", {}},
              {"wald", "Simple asymptotic (Wald)", {}},
              {"waldcc", "Wald with continuity correction", {}}},
             "newcombe")
          .grp(G_TEST)
          .hint("The formula of the interval for P1 − P2. Newcombe's score intervals and Agresti–Caffo are "
                "accurate; Wald is the textbook formula, poor with small groups or proportions near 0 or 1."),
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
      num("n1", "N1 (Group 1 Sample Size)", "100")
          .grp(G_N)
          .sym("N1")
          .at_least(1)
          .whole()
          .sugg({"100", "50 to 200 by 50"})
          .hint("Number of subjects in group 1; group 2 gets R × N1, rounded up.")
          .show(when(not_solving({"n"}))),
      num("r", "R (Allocation Ratio N2/N1)", "1")
          .grp(G_N)
          .sym("R")
          .above(0)
          .sugg({"1", "2"})
          .hint("Size of group 2 relative to group 1 (N2 = R × N1, rounded up); 1 gives equal groups."),
      num("width", "Width (two-sided: upper − lower; one-sided: distance to bound)", "0.2")
          .grp(G_EFFECT)
          .sym("Width")
          .range(0, 2, true, true)
          .sugg({"0.2", "0.1 0.15 0.2"})
          .hint("The precision you want for P1 − P2 (0.1 = 10 percentage points): the full width of a two-sided "
                "interval, or the distance from the estimate to a one-sided bound.")
          .show(when(solving({"n"}))),
      num("p1", "P1 (Anticipated Proportion, Group 1)", "0.6")
          .grp(G_EFFECT)
          .sym("P1")
          .range(0, 1, true, true)
          .sugg({"0.6", "0.5 0.6 0.7"})
          .hint("The proportion you expect to observe in group 1, between 0 and 1; the width depends on it."),
      num("p2", "P2 (Anticipated Proportion, Group 2)", "0.4")
          .grp(G_EFFECT)
          .sym("P2")
          .range(0, 1, true, true)
          .sugg({"0.4", "0.3 0.4"})
          .hint("The proportion you expect to observe in group 2, between 0 and 1."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    double n1;
    if (s.solve == "n") {
      const double target = s("width");
      const bool two = s.is("sides", "two");
      const double z = stats::qnorm((1 - s("conf")) / (two ? 2 : 1), false);
      const double p1 = s("p1"), p2 = s("p2");
      const double e = two ? target / 2 : target;
      solve::NSearch o;
      o.lo = 1;
      o.guess = std::ceil(z * z * (p1 * (1 - p1) + p2 * (1 - p2) / s("r")) / (e * e));
      o.hi = 1e8;
      o.goal = "the target width";
      n1 = solve::min_n_staying(
          [&](double x) { const double w = diff_ci_at(s, x).width; return w > 0 && w <= target; }, o, &ctx);
      r.set("target_width", target);
    } else {
      n1 = s("n1");
    }
    const DiffCi c = diff_ci_at(s, n1);
    r.set("n1", n1);
    r.set("n2", c.n2);
    r.set("n", n1 + c.n2);
    r.set("width", c.width);
    r.set("p1", s("p1"));
    r.set("p2", s("p2"));
    r.set("diff", s("p1") - s("p2"));
    r.set("lower", s.is("sides", "upper") ? NAN : c.ci.lo);
    r.set("upper", s.is("sides", "lower") ? NAN : c.ci.hi);
    r.set("conf", s("conf"));
  };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c = {col("conf", "Confidence", "num"), col("n1", "N1", "int"), col("n2", "N2", "int"),
                             col("n", "N", "int")};
    if (s.solve == "n") c.push_back(col("target_width", "Target Width", "num"));
    c.push_back(col("width", s.solve == "n" ? "Actual Width" : "Width", "num"));
    c.push_back(col("p1", "P1", "num", "Anticipated proportion in group 1"));
    c.push_back(col("p2", "P2", "num", "Anticipated proportion in group 2"));
    c.push_back(col("diff", "P1 − P2", "num"));
    c.push_back(col("lower", "Lower Limit", "num"));
    c.push_back(col("upper", "Upper Limit", "num"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    static const std::map<std::string, std::string> names = {
        {"newcombe", "Newcombe hybrid score"}, {"newcombecc", "continuity-corrected Newcombe hybrid score"},
        {"ac", "Agresti–Caffo"}, {"wald", "Wald"}, {"waldcc", "continuity-corrected Wald"}};
    const bool two = s.is("sides", "two");
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " give a " +
           (two ? "two-sided " : "one-sided ") + fpct(r.num("conf"), 0) + " " + names.at(s.s("method")) +
           " confidence interval for the difference of two proportions " +
           (two ? "with a width of " : "whose bound lies at a distance of ") + fnum(r.num("width")) +
           " when the observed proportions equal " + fnum(r.num("p1")) + " and " + fnum(r.num("p2")) + ".";
  };
  p.definitions = {{"N1, N2", "Group sample sizes; N2 = R·N1 rounded up."},
                   {"Width", "Upper minus lower limit (two-sided) or distance from P1 − P2 to the bound (one-sided)."},
                   {"P1, P2", "Anticipated proportions at which the interval is evaluated."}};
  p.references = {"Newcombe, R. G. (1998). Interval estimation for the difference between independent proportions: "
                  "comparison of eleven methods. Statistics in Medicine 17, 873–890.",
                  "Agresti, A. & Caffo, B. (2000). Simple and effective confidence intervals for proportions and "
                  "differences of proportions result from adding two successes and two failures. The American "
                  "Statistician 54, 280–288.",
                  REF_FLEISS};
  p.example_solve = "n";
  p.example = {{"width", "0.1 0.2"}, {"p1", "0.6"}, {"p2", "0.4 0.5"}};
  return p;
}
SSC_REGISTER(two_props_ci);

}  // namespace
}  // namespace ssc::procs
