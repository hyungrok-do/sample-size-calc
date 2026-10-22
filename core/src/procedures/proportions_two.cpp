// Two independent proportions: z-tests (pooled / unpooled / continuity
// corrected) and Fisher's exact test, with exact power by enumeration of the
// two binomial distributions (restricted to the region carrying all but
// 1e-16 of the probability) or with normal approximations.  Also the
// confidence interval for their difference, and case-control studies:
// unmatched (the exposure of cases against controls) and matched (Dupont).
#include <atomic>
#include <limits>
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
const char* const FELL_BACK = "The samples are too large for exact enumeration; the normal approximation was used.";
// Fisher's test also builds the null distribution of every table total: about 7 to 9 times the cells when
// neither proportion is near 0 or 1, and far more when one is (few totals, each with a wide null distribution)
constexpr double MAX_FISHER_WORK = 12 * MAX_CELLS;

// A loop over k = 0..n-1, on the idle workers when it is large (it is nested in the loop over scenarios), and
// one that Stop can end either way.
template <class F>
void enumerate(std::size_t n, double work, const RunContext& ctx, F&& body) {
  if (work > 2e5 && ctx.threads > 1) {
    parallel_for(n, ctx.threads, body, ctx.token);
    return;
  }
  for (std::size_t k = 0; k < n; ++k) {
    ctx.check();  // a step is a whole row or table total: thousands of cells
    body(k);
  }
}

// Exact power by enumeration over the joint binomial window; NaN when that is too large (the callers then use
// the normal approximation and say so).
double two_prop_exact(const TwoPropTest& T, double n1d, double n2d, double p1, double p2, const RunContext& ctx) {
  // the counts, their totals and the sizes' sum are ints
  if (!(n1d >= 0 && n2d >= 0 && n1d + n2d <= std::numeric_limits<int>::max())) return NAN;
  const int n1 = static_cast<int>(n1d), n2 = static_cast<int>(n2d);
  const BinomWindow w1 = binom_window(n1, p1), w2 = binom_window(n2, p2);
  const double cells = static_cast<double>(w1.hi - w1.lo + 1) * (w2.hi - w2.lo + 1);
  if (cells > MAX_CELLS) return NAN;
  if (T.test == "fisher") {
    const int tlo = w1.lo + w2.lo, thi = w1.hi + w2.hi;
    const std::size_t nt = static_cast<std::size_t>(thi - tlo + 1);
    // fisher_region's window is 80 null SDs wide, the widest at the total nearest N/2
    const double N = n1d + n2d, tm = std::clamp(N / 2, static_cast<double>(tlo), static_cast<double>(thi));
    const double var = N > 1 ? tm * n1d / N * (N - tm) / N * (N - n1d) / (N - 1) : 0;
    const double work = static_cast<double>(nt) * (80 * std::sqrt(var) + 41);
    if (work > MAX_FISHER_WORK) return NAN;
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
    enumerate(nt, work, ctx, body);
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
  enumerate(nr, cells, ctx, body);
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
      // a point of the power curve is enumerated only while all the curves stay within the exact search's budget
      if (exact && !on_curve_too_long(ctx, exact_power_work(n1 * ctx.curve_scale, n2 * ctx.curve_scale, x, p2,
                                                             T.test == "fisher"))) {
        const double pw = two_prop_exact(T, n1, n2, x, p2, ctx);
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
      // the groups the search tries: with a ratio allocation, group 2 is R times group 1
      const double ratio = s.is("alloc", "ratio") ? s("r") : 1.0, g1 = g * (1 + 1 / ratio) / 2;
      if (exact && exact_search_work(g1, ratio * g1, p1, p2) > EXACT_SEARCH_BUDGET) exact = false, fell_back = true;
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
    if (exact) fell_back = false;  // the answer is labeled by its own evaluation, not by the search's steps
    put_power(s, r, power(n.n1, n.n2, p1));
    r.set("p1", p1);
    r.set("p2", p2);
    r.set("diff", p1 - p2);
    r.set("ratio", p1 / p2);
    r.set("or", p1 / (1 - p1) / (p2 / (1 - p2)));
    r.set("alpha", s("alpha"));
    if (exact && !fell_back && !ctx.curve_points) {  // a curve uses only the power
      // the size can need a larger enumeration than the power (when P1 is the more extreme proportion)
      const double size = two_prop_exact(T, n.n1, n.n2, p2, p2, ctx);
      if (!std::isnan(size)) r.set("actual_alpha", size);
    }
    if (fell_back) r.notes.push_back(FELL_BACK);
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
    // large samples use the normal approximation, and say so in a note
    const bool enumerated = s.is("method", "exact") && std::find(r.notes.begin(), r.notes.end(), FELL_BACK) == r.notes.end();
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

// =====================================================================================
// Case-control studies
// =====================================================================================
// Exposure among the cases whose odds are the odds ratio times those among the controls.
double case_exposure(double odds_ratio, double p0) { return odds_ratio * p0 / (1 - p0 + odds_ratio * p0); }
double odds_ratio_of(double p1, double p0) { return p1 / (1 - p1) / (p0 / (1 - p0)); }

ParamDef exposure_param() {
  return num("p0", "P0 (Exposure Prevalence among Controls)", "0.3")
      .grp(G_EFFECT)
      .sym("P0")
      .range(0, 1, true, true)
      .sugg({"0.3", "0.1 0.2 0.3"})
      .hint("The share of controls (people without the disease) who were exposed, between 0 and 1; it stands for "
            "the exposure in the population the cases arose from.");
}

ParamDef odds_ratio_param(const char* hint) {
  return num("or", "OR (Odds Ratio of Exposure)", "2").grp(G_EFFECT).sym("OR").above(0).sugg({"2", "1.5 2 3"}).hint(hint);
}

ParamDef or_direction_param() {
  return select("ordir", "Search Direction for the Odds Ratio", {{"above", "OR > 1", {}}, {"below", "OR < 1", {}}},
                "above")
      .grp(G_EFFECT)
      .hint("A two-sided test can detect an odds ratio above or below 1; choose which of the two to find.")
      .show(when({solving({"or"}), is("alternative", {"two"})}));
}

const char* REF_SCHLESSELMAN = "Schlesselman, J. J. (1982). Case-Control Studies: Design, Conduct, Analysis. Oxford "
                               "University Press.";

Procedure case_control() {
  Procedure p;
  p.id = "prop.cc";
  p.name = "Unmatched Case-Control Study (Odds Ratio)";
  p.category = "Proportions";
  p.subcategory = "Case-Control Studies";
  p.order = 1;
  p.summary =
      "Numbers of cases and controls needed to detect an odds ratio of exposure in an unmatched case-control study, "
      "from the exposure prevalence among the controls and the number of controls per case, comparing the exposure "
      "of cases and controls with a chi-square test (with or without continuity correction) or Fisher's exact test. "
      "Exact power enumerates both binomial distributions.";
  p.keywords = {"case-control study", "case-control", "case control", "odds ratio", "exposure", "controls per case", "unmatched",
                "epidemiology", "risk factor", "retrospective study"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"},
                 {"power", "Power", "power", "Power"},
                 {"or", "Effect Size (OR)", "or", "Detectable Odds Ratio"}};
  p.params = {
      alternative("Two-sided (H1: OR ≠ 1)", "One-sided (H1: OR > 1)", "One-sided (H1: OR < 1)"),
      select("test", "Test Type",
             {{"zpool", "Chi-square test (pooled z-test)", {}},
              {"zpoolcc", "Chi-square test with continuity correction (Yates)", {}},
              {"zunpool", "Z-test, unpooled variance", {}},
              {"fisher", "Fisher's exact test", {}}},
             "zpool")
          .grp(G_TEST)
          .hint("The test comparing the exposure of cases and controls. The chi-square test is the usual choice; "
                "Yates' correction and Fisher's exact test are more conservative."),
      select("method", "Power Calculation",
             {{"exact", "Exact (binomial enumeration)", {}}, {"normal", "Normal approximation", {}}}, "exact")
          .grp(G_TEST)
          .hint("Enumeration gives the exact power of the chosen test; the normal approximation is the classical "
                "closed-form formula (Schlesselman; Fleiss for the continuity correction)."),
      power_param(),
      alpha_param(),
      n_param("ncase", "Number of Cases", "100", 2).sym("Cases").hint("Number of cases (people with the disease)."),
      num("m", "M (Controls per Case)", "1")
          .grp(G_N)
          .sym("M")
          .above(0)
          .sugg({"1", "2", "4", "1 2 3 4"})
          .hint("Number of controls per case: controls = M × cases, rounded up. More than four controls per case "
                "adds little power."),
      select("es", "Effect Size Input",
             {{"or", "Odds ratio OR", {}}, {"p1", "Exposure among cases P1", {}}}, "or")
          .grp(G_EFFECT)
          .hint("Describe the effect by the odds ratio of exposure, or by the share of cases expected to have been "
                "exposed.")
          .show(when(not_solving({"or"}))),
      odds_ratio_param("The odds ratio of exposure to detect: the odds of exposure among cases divided by that among "
                       "controls (about the relative risk when the disease is rare).")
          .show(when({not_solving({"or"}), is("es", {"or"})})),
      num("p1", "P1 (Exposure Prevalence among Cases)", "0.46")
          .grp(G_EFFECT)
          .sym("P1")
          .range(0, 1, true, true)
          .sugg({"0.46", "0.4 0.5"})
          .hint("The share of cases expected to have been exposed, between 0 and 1.")
          .show(when({not_solving({"or"}), is("es", {"p1"})})),
      exposure_param(),
      or_direction_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double p0 = s("p0"), m = s("m");
    const Alt a = alt(s);
    const TwoPropTest T{s.s("test"), a, s("alpha")};
    bool exact = s.is("method", "exact");
    bool fell_back = false;
    auto power = [&](double n1, double n2, double x) {
      // a point of the power curve is enumerated only while all the curves stay within the exact search's budget
      if (exact && !on_curve_too_long(ctx, exact_power_work(n1 * ctx.curve_scale, n2 * ctx.curve_scale, x, p0,
                                                             T.test == "fisher"))) {
        const double pw = two_prop_exact(T, n1, n2, x, p0, ctx);
        if (!std::isnan(pw)) return pw;
        fell_back = true;
      }
      return two_prop_normal(T, n1, n2, x, p0);
    };
    auto controls = [&](double cases) { return std::ceil(m * cases - 1e-9); };
    TwoN n;
    double p1;
    if (s.solve == "n") {
      p1 = s.is("es", "p1") ? s("p1") : case_exposure(s("or"), p0);
      if (p1 == p0) throw CalcError("the odds ratio is 1: cases and controls have the same exposure");
      if ((a == Alt::Upper && p1 < p0) || (a == Alt::Lower && p1 > p0))
        throw CalcError("the odds ratio is on the side opposite to the one-sided alternative");
      // Schlesselman's formula for the number of cases
      const double za = stats::qnorm(a == Alt::Two ? s("alpha") / 2 : s("alpha"), false), zb = stats::qnorm(s("power"));
      const double pb = (p1 + m * p0) / (1 + m);
      const double g = std::pow(za * std::sqrt((1 + 1 / m) * pb * (1 - pb)) + zb * std::sqrt(p1 * (1 - p1) + p0 * (1 - p0) / m), 2) /
                       ((p1 - p0) * (p1 - p0));
      if (exact && exact_search_work(g, m * g, p1, p0) > EXACT_SEARCH_BUDGET) exact = false, fell_back = true;
      Scenario t = s;  // the two-group search with controls = M × cases
      t.sel["alloc"] = "ratio";
      t.num["r"] = m;
      n = solve_two_group_n(
          t, [&](double c, double k) { return power(c, k, p1); }, s("power"), 2, std::ceil(g), !exact, ctx, 1e7);
    } else {
      n = {s("ncase"), controls(s("ncase"))};
      if (s.solve == "power") {
        p1 = s.is("es", "p1") ? s("p1") : case_exposure(s("or"), p0);
      } else {
        const bool up = a == Alt::Upper || (a == Alt::Two && !s.is("ordir", "below"));
        auto f = [&](double x) { return power(n.n1, n.n2, x) - s("power"); };
        const double lo = up ? p0 : 1e-9, hi = up ? 1 - 1e-9 : p0;
        if (f(p0) >= 0) throw CalcError("the target power must exceed the power at OR = 1 (the significance level)");
        if (f(up ? hi : lo) < 0) throw CalcError("the target power cannot be reached with these numbers of cases and controls");
        p1 = numerics::brent(f, lo, hi, 1e-11);
      }
    }
    r.set("ncase", n.n1);
    r.set("nctrl", n.n2);
    r.set("n", n.total());
    if (exact) fell_back = false;  // the answer is labeled by its own evaluation, not by the search's steps
    put_power(s, r, power(n.n1, n.n2, p1));
    r.set("m", m);
    r.set("p1", p1);
    r.set("p0", p0);
    r.set("or", odds_ratio_of(p1, p0));
    r.set("alpha", s("alpha"));
    if (exact && !fell_back && !ctx.curve_points) {  // a curve uses only the power
      // the size can need a larger enumeration than the power (when P1 is the more extreme proportion)
      const double size = two_prop_exact(T, n.n1, n.n2, p0, p0, ctx);
      if (!std::isnan(size)) r.set("actual_alpha", size);
    }
    if (fell_back) r.notes.push_back(FELL_BACK);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("ncase", "Cases", "int"), col("nctrl", "Controls", "int"), col("n", "N", "int", "Cases plus controls"),
                   col("m", "M", "num", "Controls per case"), col("or", "OR", "num", "Odds ratio of exposure"),
                   col("p1", "P1", "num", "Exposure among cases"), col("p0", "P0", "num", "Exposure among controls"),
                   col("alpha", "Target Alpha", "alpha")})
      c.push_back(x);
    if (s.is("method", "exact")) c.push_back(col("actual_alpha", "Actual Alpha", "prob", "Exact size when OR = 1"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    static const std::map<std::string, std::string> names = {{"zpool", "chi-square test"},
                                                             {"zunpool", "unpooled z-test"},
                                                             {"zpoolcc", "continuity-corrected chi-square test"},
                                                             {"fisher", "Fisher's exact test"}};
    std::string t = fint(r.num("ncase")) + " cases and " + fint(r.num("nctrl")) + " controls " + achieves(s, r) +
                    " to detect an odds ratio of exposure of " + fnum(r.num("or"), 3) + " (exposure " +
                    fpct(r.num("p1")) + " among cases and " + fpct(r.num("p0")) + " among controls) using a " +
                    alt_words(alt(s)) + " " + names.at(s.s("test")) + " at a target significance level of " +
                    fnum(r.num("alpha"));
    if (r.has("actual_alpha")) t += " (actual " + text::format_fixed(r.num("actual_alpha"), 4) + ")";
    const bool enumerated = s.is("method", "exact") && std::find(r.notes.begin(), r.notes.end(), FELL_BACK) == r.notes.end();
    return t + (enumerated ? "; power was computed by exact enumeration." : "; power uses the normal approximation.");
  };
  p.definitions = {{"Cases, Controls", "People with and without the disease; controls = M × cases, rounded up."},
                   {"P1, P0", "Shares of cases and of controls who were exposed."},
                   {"OR", "Odds ratio of exposure, P1/(1 − P1) divided by P0/(1 − P0)."},
                   {"Actual Alpha", "Exact probability of rejecting H0 when cases and controls share the exposure P0."}};
  p.notes = {"With exact enumeration, power is a saw-toothed function of the number of cases; the reported number "
             "is the first one that reaches the target power.",
             "Cases and controls are sampled independently; for individually matched controls use the matched "
             "case-control procedure."};
  p.references = {REF_SCHLESSELMAN, REF_FLEISS,
                  "Kelsey, J. L., Whittemore, A. S., Evans, A. S. & Thompson, W. D. (1996). Methods in Observational "
                  "Epidemiology, 2nd ed. Oxford University Press."};
  p.curve_param = "ncase";
  p.curve_label = "Number of Cases";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"or", "2"}, {"p0", "0.2 0.3"}, {"m", "1 2"}};
  return p;
}
SSC_REGISTER(case_control);

// Matched case-control studies (Dupont 1988): sets of one case and M controls whose exposures are correlated (phi)
// through the matching.  Given t exposed members in a set, the case is the exposed one with probability
// t·OR/(t·OR + M − t + 1) (t/(M + 1) under H0); sets with no or all members exposed carry no information.
struct Matched {
  double delta = 0, v0 = 0, v1 = 0;  // per set: mean of (case exposed − expected under H0), null and true variance
};

Matched matched_moments(double odds_ratio, double p0, double phi, int M) {
  const double p1 = case_exposure(odds_ratio, p0), q1 = 1 - p1, q0 = 1 - p0;
  // exposure of a control given that its case is exposed (plus) or not (minus)
  const double pp = p0 + phi * std::sqrt(q1 * p0 * q0 / p1), pm = p0 - phi * std::sqrt(p1 * p0 * q0 / q1);
  if (!(pp >= 0 && pp <= 1 && pm >= 0 && pm <= 1))
    throw CalcError("the correlation φ is too large in size for these exposures: the exposure of a control given its "
                    "case's would fall outside 0 to 1");
  Matched mo;
  for (int t = 1; t <= M; ++t) {
    // the chance that a set has t exposed members: the case and t − 1 controls, or t controls
    const double pt = p1 * std::exp(stats::lchoose(M, t - 1)) * std::pow(pp, t - 1) * std::pow(1 - pp, M - t + 1) +
                      q1 * std::exp(stats::lchoose(M, t)) * std::pow(pm, t) * std::pow(1 - pm, M - t);
    const double e = t * odds_ratio / (t * odds_ratio + M - t + 1), e0 = t / (M + 1.0);
    mo.delta += pt * (e - e0);
    mo.v0 += pt * e0 * (1 - e0);
    mo.v1 += pt * e * (1 - e);
  }
  return mo;
}

double matched_power(double n, double odds_ratio, double p0, double phi, int M, double alpha, Alt a) {
  using namespace stats;
  const Matched mo = matched_moments(odds_ratio, p0, phi, M);
  const double rn = std::sqrt(n), s0 = std::sqrt(mo.v0), s1 = std::sqrt(mo.v1);
  if (!(s1 > 0)) return NAN;
  auto tail = [&](double d, double z) { return pnorm((rn * d - z * s0) / s1); };
  if (a == Alt::Two) {
    const double z = qnorm(alpha / 2, false);
    return std::min(1.0, tail(mo.delta, z) + tail(-mo.delta, z));
  }
  return tail(a == Alt::Upper ? mo.delta : -mo.delta, qnorm(alpha, false));
}

Procedure matched_case_control() {
  Procedure p;
  p.id = "prop.cc.matched";
  p.name = "Matched Case-Control Study (1:M Matching)";
  p.category = "Proportions";
  p.subcategory = "Case-Control Studies";
  p.order = 2;
  p.summary =
      "Number of matched sets — each one case and M individually matched controls — needed to detect an odds ratio "
      "of exposure with the conditional (Mantel–Haenszel) test, allowing for the correlation of exposure within "
      "matched sets (Dupont 1988). With one control per case and no correlation it is the classical matched-pairs "
      "formula.";
  p.keywords = {"matched case-control", "matching", "matched pairs", "1:m matching", "odds ratio",
                "conditional logistic regression", "mantel-haenszel", "dupont", "epidemiology", "exposure"};
  p.solve_for = {{"n", "Sample Size", "n", "Matched Sets"},
                 {"power", "Power", "power", "Power"},
                 {"or", "Effect Size (OR)", "or", "Detectable Odds Ratio"}};
  p.params = {
      alternative("Two-sided (H1: OR ≠ 1)", "One-sided (H1: OR > 1)", "One-sided (H1: OR < 1)"),
      power_param(),
      alpha_param(),
      n_param("n", "Number of Matched Sets (Cases)", "100", 1)
          .sym("Sets")
          .hint("Number of matched sets, each one case with its M controls."),
      num("mctl", "M (Controls per Case)", "1")
          .grp(G_N)
          .sym("M")
          .at_least(1)
          .whole()
          .range(1, 20)
          .sugg({"1", "2", "4", "1 2 3 4"})
          .hint("Number of controls matched to each case, a whole number; beyond four the gain is small."),
      odds_ratio_param("The odds ratio of exposure to detect, as estimated by conditional logistic regression.")
          .show(when(not_solving({"or"}))),
      exposure_param(),
      num("phi", "φ (Correlation of Exposure within Matched Sets)", "0.2")
          .grp(G_EFFECT)
          .sym("φ")
          .range(-1, 1, true, true)
          .sugg({"0.2", "0 0.1 0.2"})
          .hint("Correlation between the exposures of a case and its matched control, from matching on factors "
                "related to the exposure; 0 when matching is unrelated to it. Dupont suggests 0.2 when unknown."),
      or_direction_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double p0 = s("p0"), phi = s("phi"), alpha = s("alpha");
    const int M = static_cast<int>(s("mctl"));
    const Alt a = alt(s);
    double n, odds_ratio;
    if (s.solve == "n") {
      odds_ratio = s("or");
      if (odds_ratio == 1) throw CalcError("the odds ratio is 1: there is no effect to detect");
      if ((a == Alt::Upper && odds_ratio < 1) || (a == Alt::Lower && odds_ratio > 1))
        throw CalcError("the odds ratio is on the side opposite to the one-sided alternative");
      const Matched mo = matched_moments(odds_ratio, p0, phi, M);
      const double za = stats::qnorm(a == Alt::Two ? alpha / 2 : alpha, false), zb = stats::qnorm(s("power"));
      solve::NSearch o;
      o.lo = 1;
      o.hi = 1e8;
      o.guess = std::ceil(std::pow(za * std::sqrt(mo.v0) + zb * std::sqrt(mo.v1), 2) / (mo.delta * mo.delta));
      n = solve::min_n([&](double x) { return matched_power(x, odds_ratio, p0, phi, M, alpha, a); }, s("power"), o, &ctx);
    } else {
      n = s("n");
      if (s.solve == "power") {
        odds_ratio = s("or");
      } else {
        const bool up = a == Alt::Upper || (a == Alt::Two && !s.is("ordir", "below"));
        // the search runs on the log odds ratio, up to where φ still gives valid exposures
        auto f = [&](double lor) { return matched_power(n, std::exp(lor), p0, phi, M, alpha, a) - s("power"); };
        if (f(0) >= 0) throw CalcError("the target power must exceed the power at OR = 1 (the significance level)");
        // widen the bracket until the power reaches the target; where φ gives no valid exposures (too strong an
        // effect), close in on the last valid point instead
        double good = 0, far = up ? 1 : -1, fail = NAN;
        while (true) {
          bool valid = true;
          double v = NAN;
          try {
            v = f(far);
          } catch (const CalcError&) {
            valid = false;
          }
          if (valid && v >= 0) break;
          if (valid)
            good = far;
          else
            fail = far;
          if (std::isnan(fail)) {
            if (std::fabs(far) > 30) throw CalcError("the target power cannot be reached with this number of matched sets");
            far *= 2;
          } else {
            if (std::fabs(fail - good) < 1e-10)
              throw CalcError("the target power cannot be reached with this number of matched sets: the correlation φ "
                              "allows no stronger effect");
            far = (good + fail) / 2;
          }
        }
        odds_ratio = std::exp(numerics::brent(f, std::fmin(good, far), std::fmax(good, far), 1e-11));
      }
    }
    r.set("n", n);
    r.set("nctrl", n * M);
    r.set("ntotal", n * (M + 1));
    put_power(s, r, matched_power(n, odds_ratio, p0, phi, M, alpha, a));
    r.set("mctl", M);
    r.set("or", odds_ratio);
    r.set("p1", case_exposure(odds_ratio, p0));
    r.set("p0", p0);
    r.set("phi", phi);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "Sets", "int", "Matched sets (cases)"), col("nctrl", "Controls", "int"),
                   col("ntotal", "N", "int", "Cases plus controls"), col("mctl", "M", "int", "Controls per case"),
                   col("or", "OR", "num", "Odds ratio of exposure"), col("p1", "P1", "num", "Exposure among cases"),
                   col("p0", "P0", "num", "Exposure among controls"),
                   col("phi", "φ", "num", "Correlation of exposure within matched sets"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const int M = static_cast<int>(r.num("mctl"));
    return fint(r.num("n")) + " matched sets of one case and " + (M == 1 ? "one control" : fint(M) + " controls") +
           " (" + fint(r.num("ntotal")) + " subjects) " + achieves(s, r) + " to detect an odds ratio of exposure of " +
           fnum(r.num("or"), 3) + " (exposure " + fpct(r.num("p1")) + " among cases and " + fpct(r.num("p0")) +
           " among controls, correlation of exposure within sets " + fnum(r.num("phi")) + ") using a " +
           alt_words(alt(s)) + " conditional test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"Sets", "Matched sets: one case with M matched controls."},
                   {"P1, P0", "Shares of cases and of controls who were exposed."},
                   {"OR", "Odds ratio of exposure, P1/(1 − P1) divided by P0/(1 − P0)."},
                   {"φ", "Correlation between the exposures of a case and a matched control."}};
  p.notes = {"Power uses Dupont's normal approximation to the conditional test on the matched sets in which some but "
             "not all members are exposed."};
  p.references = {"Dupont, W. D. (1988). Power calculations for matched case-control studies. Biometrics 44, "
                  "1157–1168.",
                  REF_SCHLESSELMAN,
                  "Breslow, N. E. & Day, N. E. (1980). Statistical Methods in Cancer Research, Vol. I: The Analysis of "
                  "Case-Control Studies. IARC."};
  p.curve_param = "n";
  p.curve_label = "Matched Sets";
  p.curve_min = 1;
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"or", "2"}, {"p0", "0.3"}, {"mctl", "1 2 4"}};
  return p;
}
SSC_REGISTER(matched_case_control);

}  // namespace
}  // namespace ssc::procs
