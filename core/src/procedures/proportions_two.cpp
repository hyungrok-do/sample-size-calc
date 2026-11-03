// Two independent proportions: z-tests (pooled / unpooled / continuity
// corrected) and Fisher's exact test, with exact power by enumeration of the
// two binomial distributions (restricted to the region carrying all but
// 1e-16 of the probability) or with normal approximations.  Also the
// confidence interval for their difference, and case-control studies:
// unmatched (the exposure of cases against controls) and matched (Dupont).
#include <atomic>
#include <limits>
#include <mutex>

#include "procedures/apart.hpp"
#include "procedures/binomial_tools.hpp"
#include "procedures/common.hpp"
#include "procedures/derive_props.hpp"
#include "procedures/proportion_ci.hpp"
#include "procedures/proportion_effect.hpp"
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

// The normal approximation's power and its parts (which the derivation shows): the pooled proportion, the standard
// errors under H0 and H1, the continuity correction, the critical value and each tail.
struct TwoPropNormal {
  double d, pb, s0, s1, cc, z, up, lo, power;
};

TwoPropNormal two_prop_normal_parts(const TwoPropTest& T, double n1, double n2, double p1, double p2) {
  using namespace stats;
  TwoPropNormal r;
  const double q1 = 1 - p1, q2 = 1 - p2;
  r.pb = (n1 * p1 + n2 * p2) / (n1 + n2);
  const double qb = 1 - r.pb;
  r.s1 = std::sqrt(p1 * q1 / n1 + p2 * q2 / n2);
  r.s0 = T.test == "zunpool" ? r.s1 : std::sqrt(r.pb * qb * (1 / n1 + 1 / n2));
  r.cc = (T.test == "zpoolcc" || T.test == "fisher") ? 0.5 * (1 / n1 + 1 / n2) : 0;
  r.d = p1 - p2;
  // variances P(1 − P)/n that underflow to 0 (proportions near 5e-324) would make a difference between them
  // infinitely significant, a power of 1 (an infinite count, which makes them no number, is said by the runner)
  if (r.d != 0 && (r.s0 == 0 || r.s1 == 0)) throw CalcError("the proportions are too small to compute with");
  r.z = qnorm(tail_alpha(T.alpha, T.alt), false);
  r.up = pnorm((r.d - r.cc - r.z * r.s0) / r.s1);
  r.lo = pnorm((-r.d - r.cc - r.z * r.s0) / r.s1);
  r.power = T.alt == Alt::Two ? std::min(1.0, r.up + r.lo) : (T.alt == Alt::Upper ? r.up : r.lo);
  return r;
}

double two_prop_normal(const TwoPropTest& T, double n1, double n2, double p1, double p2) {
  return two_prop_normal_parts(T, n1, n2, p1, p2).power;
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
  const double zc = stats::qnorm(tail_alpha(T.alpha, T.alt), false);
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

// ------------------------------------------------------------------ derivation
// The groups' symbols, 1 and 2 (treatment and control) or 1 and 0 (cases and controls), and their sizes and
// proportions as written, each marked given or derived.
struct Groups {
  std::string a, b;
  std::string n1, n2, p1, p2;
};

// The power as the runs compute it: by enumeration where asked for and not too large, else the normal approximation.
double two_prop_power(const TwoPropTest& T, double n1, double n2, double p1, double p2, bool enumerated,
                      const RunContext& ctx) {
  if (enumerated) {
    const double pw = two_prop_exact(T, n1, n2, p1, p2, ctx);
    if (!std::isnan(pw)) return pw;
  }
  return two_prop_normal(T, n1, n2, p1, p2);
}

// The z statistic of a table (as z_two), in symbols and with the table's numbers.
std::string z_two_formula(const std::string& test, const Groups& g) {
  const std::string h1 = "p̂" + g.a, h2 = "p̂" + g.b, n1 = "N" + g.a, n2 = "N" + g.b;
  const std::string props = ", " + h1 + " = x" + g.a + "/" + n1 + ", " + h2 + " = x" + g.b + "/" + n2;
  const std::string pooled = "√(p̄(1 − p̄)(1/" + n1 + " + 1/" + n2 + ")), p̄ = (x" + g.a + " + x" + g.b + ")/(" + n1 +
                             " + " + n2 + ")";
  if (test == "zunpool")
    return "Z = (" + h1 + " − " + h2 + ")/√(" + h1 + "(1 − " + h1 + ")/" + n1 + " + " + h2 + "(1 − " + h2 + ")/" + n2 +
           ")" + props;
  if (test == "zpoolcc")
    return "Z = (|" + h1 + " − " + h2 + "| − (1/" + n1 + " + 1/" + n2 + ")/2)/" + pooled + ", with the sign of " + h1 +
           " − " + h2 + props;
  return "Z = (" + h1 + " − " + h2 + ")/" + pooled + props;
}

std::string z_two_numbers(const std::string& test, int x1, int n1, int x2, int n2, const Groups& g) {
  using namespace derive;
  const double h1 = static_cast<double>(x1) / n1, h2 = static_cast<double>(x2) / n2;
  const std::string f = dvi(x1) + "/" + g.n1 + " − " + dvi(x2) + "/" + g.n2;
  const double pb = static_cast<double>(x1 + x2) / (n1 + n2);
  const std::string se =
      test == "zunpool"
          ? "√(" + dv(h1) + " × " + dv(1 - h1) + "/" + g.n1 + " + " + dv(h2) + " × " + dv(1 - h2) + "/" + g.n2 + ")"
          : "√(" + dv(pb) + " × " + dv(1 - pb) + " × (1/" + g.n1 + " + 1/" + g.n2 + "))";
  const std::string top = test == "zpoolcc" ? "(|" + f + "| − " + dv(0.5 * (1.0 / n1 + 1.0 / n2)) + ")" : "(" + f + ")";
  return top + "/" + se + " = " + dv(z_two(test, x1, n1, x2, n2));
}

// Fisher's conditional test at the table total t: the values of x1 it rejects (as fisher_region finds them) and
// their conditional size.
std::string fisher_text(int n1, int n2, int t, double alpha, Alt a, const std::string& x) {
  using namespace derive;
  int lo, hi;
  std::vector<char> rej;
  std::vector<double> f;
  fisher_region(n1, n2, t, alpha, a, lo, hi, rej, f);
  const int m = hi - lo + 1;
  int i = 0, j = m - 1;  // rejected: [0, i) at the bottom and (j, m) at the top
  while (i < m && rej[static_cast<std::size_t>(i)]) ++i;
  while (j >= i && rej[static_cast<std::size_t>(j)]) --j;
  double size = 0;
  for (int k = 0; k < m; ++k)
    if (rej[static_cast<std::size_t>(k)]) size += f[static_cast<std::size_t>(k)];
  std::string region;
  if (i > 0) region = x + " ≤ " + dvi(lo + i - 1);
  if (j < m - 1) region += (region.empty() ? "" : " or ") + x + " ≥ " + dvi(lo + j + 1);
  if (region.empty()) return "no value of " + x + " is rejected, a conditional size of " + dvi(0);
  return "reject when " + region + ", a conditional size of " + dvp(size);
}

// How the power of a test of two independent proportions comes about at the sizes given or found, by enumeration
// (`enumerated`) or by the normal approximation (`fell_back`: in place of an enumeration too large), and the
// enumeration's size where the row has it. The difference δ, which the normal approximation uses, comes before.
void two_prop_power_steps(const TwoPropTest& T, double n1, double n2, double p1, double p2, bool enumerated,
                          bool fell_back, const Row& r, const Groups& g, derive::Steps& out) {
  using namespace derive;
  const std::string N1 = "N" + g.a, N2 = "N" + g.b, P1 = "P" + g.a, P2 = "P" + g.b, X1 = "x" + g.a, X2 = "x" + g.b;
  if (!enumerated) {
    const TwoPropNormal q = two_prop_normal_parts(T, n1, n2, p1, p2);
    const std::string s1 = "√(" + P1 + "(1 − " + P1 + ")/" + N1 + " + " + P2 + "(1 − " + P2 + ")/" + N2 + ")";
    const std::string v1 = "√(" + g.p1 + " × " + dv(1 - p1) + "/" + g.n1 + " + " + g.p2 + " × " + dv(1 - p2) + "/" +
                           g.n2 + ") = " + dsig(q.s1);
    if (T.test == "zunpool") {
      out.push_back({"Standard error", "σ0 = σ1 = " + s1, v1});
    } else {
      out.push_back({"Standard error under H0",
                     "σ0 = √(p̄(1 − p̄)(1/" + N1 + " + 1/" + N2 + ")), p̄ = (" + N1 + " " + P1 + " + " + N2 + " " + P2 +
                         ")/(" + N1 + " + " + N2 + ")",
                     "p̄ = (" + g.n1 + " × " + g.p1 + " + " + g.n2 + " × " + g.p2 + ")/(" + g.n1 + " + " + g.n2 +
                         ") = " + dv(q.pb) + ", σ0 = √(" + dv(q.pb) + " × " + dv(1 - q.pb) + " × (1/" + g.n1 + " + 1/" +
                         g.n2 + ")) = " + dsig(q.s0)});
      out.push_back({"Standard error under H1", "σ1 = " + s1, v1});
    }
    if (q.cc > 0)
      out.push_back({"Continuity correction",
                     "c = (1/" + N1 + " + 1/" + N2 + ")/2" +
                         (T.test == "fisher" ? ", Fisher's exact test approximated by the corrected z-test" : ""),
                     "(1/" + g.n1 + " + 1/" + g.n2 + ")/2 = " + dsig(q.cc)});
    out.push_back(z_critical(T.alpha, T.alt));
    Step pw = tails_power(T.alt, "δ", q.cc > 0 ? "c" : "", "σ0", "σ1", q.d, q.cc, q.z, q.s0, q.s1, q.up, q.lo,
                          r.num("power"));
    if (fell_back) pw.formula += ", the normal approximation (the samples are too large for exact enumeration)";
    out.push_back(pw);
    return;
  }
  // the exact power: the probability, under P1 and P2, of the tables that the test rejects
  const int i1 = static_cast<int>(n1), i2 = static_cast<int>(n2);
  std::string over;
  if (T.test == "fisher") {
    const int t = static_cast<int>(std::lround(n1 * p1 + n2 * p2));
    const std::string given = "given the table's total t = " + X1 + " + " + X2 + ", " + X1 +
                              " is hypergeometric under H0, and the test rejects ";
    out.push_back(
        {"Conditional test",
         given + (T.alt == Alt::Two     ? "its least probable values, whose probabilities add up to at most α"
                  : T.alt == Alt::Upper ? "when P(X ≥ " + X1 + " | t) ≤ α"
                                        : "when P(X ≤ " + X1 + " | t) ≤ α"),
         "at the expected total t = " + dvi(t) + ": " + fisher_text(i1, i2, t, T.alpha, T.alt, X1)});
    over = "over the tables the conditional test rejects";
  } else {
    out.push_back(z_critical(T.alpha, T.alt));
    const int x1 = static_cast<int>(std::clamp(std::lround(n1 * p1), 0L, static_cast<long>(i1))),
              x2 = static_cast<int>(std::clamp(std::lround(n2 * p2), 0L, static_cast<long>(i2)));
    out.push_back({"Test statistic", z_two_formula(T.test, g),
                   "at the expected counts " + X1 + " = " + dvi(x1) + ", " + X2 + " = " + dvi(x2) + ": " +
                       z_two_numbers(T.test, x1, i1, x2, i2, g)});
    over = T.alt == Alt::Two ? "over the tables with |Z| > z" : (T.alt == Alt::Upper ? "over the tables with Z > z"
                                                                                     : "over the tables with Z < −z");
  }
  auto windows = [&](double a, double b) {
    return "Σ over " + X1 + " = " + window_text(i1, a) + ", " + X2 + " = " + window_text(i2, b);
  };
  if (r.has("actual_alpha"))
    out.push_back({"Actual alpha", "α* = the same sum with " + P1 + " = " + P2,
                   windows(p2, p2) + " = " + dvp(r.num("actual_alpha"))});
  out.push_back({"Power", table_sum(g.a, g.b, over), windows(p1, p2) + " = " + dvp(r.num("power"))});
}

// The difference the test is to detect, P1 − P2 (or P1 − P0), as the groups' proportions read.
Step diff_step(const Groups& g, double p1, double p2) {
  using namespace derive;
  return {"Difference to detect", "δ = P" + g.a + " − P" + g.b, g.p1 + " − " + g.p2 + " = " + dv(p1 - p2)};
}

// How a test of two proportions' answer comes about, following what was solved for: solving for the sizes, Fleiss's
// closed form for comparison, the power at the sizes found and at the next smaller design; solving for P1, the
// power at the P1 found, where it equals the target.
derive::Derivation two_prop_derivation(const Scenario& s, const Row& r) {
  using namespace derive;
  const TwoPropTest T{s.s("test"), alt(s), s("alpha")};
  const double n1 = r.num("n1"), n2 = r.num("n2"), p1 = r.num("p1"), p2 = r.num("p2"), target = s.get("power", NAN);
  const double power = r.num("power");
  const bool fell_back = has_note(r, FELL_BACK), enumerated = s.is("method", "exact") && !fell_back;
  const bool solved = s.solve == "n", effect = s.solve == "p1";
  const blocks::TwoN n{n1, n2};
  const auto [N1, N2] = group_sizes(s, n);
  const Groups g{"1", "2", N1, N2, mv(p1, p1_typed(s)), gv(p2)};
  Derivation out;
  if (solved) {
    out.find = find_size(sizes_unknown(s), target);
  } else if (effect) {
    const bool up = T.alt == Alt::Upper || (T.alt == Alt::Two && !s.is("p1dir", "below"));
    out.find = find_effect(std::string("P1 ") + (up ? "above" : "below") + " P2", target,
                           "with " + group_sizes_text(s, n));
  } else {
    out.find = find_power("with " + group_sizes_text(s, n));
  }
  group_sizes_step(s, n, out.steps);
  if (!effect) p1_step(s, p1, p2, out.steps);
  if (solved || !enumerated) out.steps.push_back(diff_step(g, p1, p2));
  const std::string& alloc = s.s("alloc");
  if (solved && (alloc == "equal" || alloc == "ratio"))
    out.steps.push_back(two_prop_formula(T.alt, T.alpha, target, p1, p2, alloc == "ratio" ? s("r") : 1.0,
                                         alloc == "ratio" ? "r" : "", T.test != "zunpool",
                                         T.test == "zpoolcc" || T.test == "fisher", "1", "2", g.p1, g.p2,
                                         alloc == "ratio" ? "N1" : "N1 = N2"));
  two_prop_power_steps(T, n1, n2, p1, p2, enumerated, fell_back, r, g, out.steps);
  if (solved) {
    RunContext ctx;
    out.steps.push_back(two_group_step(
        s, [&](double a1, double a2) { return two_prop_power(T, a1, a2, p1, p2, enumerated, ctx); }, n, power,
        target, 2, enumerated));
    out.answer = group_sizes_text(s, n) + " (N = " + dvi(n.total()) + "), with power " + dvp(power, target) +
                 alpha_star(r);
  } else if (effect) {
    out.steps.push_back(effect_search("P1", p1, power, target));
    out.answer = "P1 = " + dv(p1) + " (P1 − P2 = " + dv(p1 - p2) + ")";
  } else {
    out.answer = "power = " + dvp(power) + alpha_star(r);
  }
  return out;
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
    ExactSearch search;
    auto power = [&](double n1, double n2, double x) {
      // a point of the power curve is enumerated only while all the curves stay within the exact search's budget,
      // judged at its curve's largest point (in either group)
      const TwoN last = ctx.curve_last ? two_group_n(*ctx.curve_last) : TwoN{n1, n2};
      if (exact && !on_curve_too_long(ctx, exact_power_work(last.n1, last.n2, x, p2, T.test == "fisher"))) {
        search.count(exact_power_work(n1, n2, x, p2, T.test == "fisher"));
        const double pw = two_prop_exact(T, n1, n2, x, p2, ctx);
        if (!std::isnan(pw)) return pw;
        fell_back = true;
      }
      return two_prop_normal(T, n1, n2, x, p2);
    };
    TwoN n;
    double p1;
    if (s.solve == "n") {
      p1 = p1_from(s, p2);
      if (p1 == p2) throw CalcError("P1 equals P2");
      if ((a == Alt::Upper && p1 < p2) || (a == Alt::Lower && p1 > p2))
        throw CalcError("P1 − P2 has the sign opposite to the one-sided alternative");
      const double za = stats::qnorm(tail_alpha(s("alpha"), a), false), zb = stats::qnorm(s("power"));
      const double pb = (p1 + p2) / 2;
      const double g = std::pow(za * std::sqrt(2 * pb * (1 - pb)) + zb * std::sqrt(p1 * (1 - p1) + p2 * (1 - p2)), 2) /
                       ((p1 - p2) * (p1 - p2));
      // the groups the search tries: with a ratio allocation, group 2 is R times group 1
      const double ratio = s.is("alloc", "ratio") ? s("r") : 1.0, g1 = g * (1 + 1 / ratio) / 2;
      if (exact && exact_search_work(g1, ratio * g1, p1, p2, T.test == "fisher") > EXACT_SEARCH_BUDGET)
        exact = false, fell_back = true;
      n = search.run(exact, fell_back, [&] {
        return solve_two_group_n(
            s, [&](double a1, double a2) { return power(a1, a2, p1); }, s("power"), 2, std::ceil(g), !exact, ctx, 1e7);
      });
    } else {
      n = two_group_n(s);
      if (s.solve == "power") {
        p1 = p1_from(s, p2);
      } else {
        const bool up = a == Alt::Upper || (a == Alt::Two && !s.is("p1dir", "below"));
        const double lo = up ? p2 : 1e-9, hi = up ? 1 - 1e-9 : p2;
        // the search evaluates the power some 40 times: past the exact search's budget (judged at the approximate
        // answer), the normal approximation
        if (exact) {
          auto g = [&](double x) { return two_prop_normal(T, n.n1, n.n2, x, p2) - s("power"); };
          if (g(p2) < 0 && g(up ? hi : lo) >= 0 &&
              40 * exact_power_work(n.n1, n.n2, numerics::brent(g, lo, hi, 1e-9), p2, T.test == "fisher") >
                  EXACT_SEARCH_BUDGET)
            exact = false, fell_back = true;
        }
        auto f = [&](double x) { return power(n.n1, n.n2, x) - s("power"); };
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
    // a given odds ratio as given: recomputed from a P1 near 1, it would lose its digits with 1 − P1
    r.set("or", s.solve != "p1" && s.is("es", "or") ? s("or1") : p1 / (1 - p1) / (p2 / (1 - p2)));
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
    // the proportions apart from each other, the difference from 0 and the ratios from 1, their values of no effect;
    // an odds ratio whose P1 rounds to P2 (1.4 at P2 = 0.9999999999999999) as the 1 that the power is computed for,
    // beside a difference of 0 and a ratio of 1
    const auto [p1, p2] = fnum_apart(r.num("p1"), r.num("p2"));
    const double odds = r.num("p1") == r.num("p2") ? 1.0 : r.num("or");
    std::string t = "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
                    " to detect a difference between the group proportions P1 = " + p1 + " and P2 = " + p2 +
                    " (difference " + fnum_apart(r.num("diff"), 0).first + ", ratio " +
                    fnum_from(r.num("ratio"), 1, 3) + ", odds ratio " + fnum_from(odds, 1, 3) + ") using a " +
                    alt_words(alt(s)) + " " + names.at(s.s("test")) + " at a target significance level of " +
                    fnum(r.num("alpha"));
    if (r.has("actual_alpha")) t += " (actual " + fprob(r.num("actual_alpha")) + ")";
    // large samples use the normal approximation, and say so in a note
    const bool enumerated = s.is("method", "exact") && std::find(r.notes.begin(), r.notes.end(), FELL_BACK) == r.notes.end();
    return t + (enumerated ? "; power was computed by exact enumeration." : "; power uses the normal approximation.");
  };
  p.derivation = two_prop_derivation;
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
  // R × N1 rounded up, which gives group 2 a subject however small R is (rounded up after taking off 1e-9, a product
  // of 1e-9 or less left it none, and its interval was 0/0, which the limits' clamps pressed into −1 to 1)
  c.n2 = solve::round_up(s("r") * n1);
  c.ci = diff_ci(s.s("method"), n1, n1 * p1, c.n2, c.n2 * p2, s("conf"), s.is("sides", "two"));
  if (s.is("sides", "two"))
    c.width = c.ci.hi - c.ci.lo;
  else
    c.width = s.is("sides", "upper") ? c.ci.hi - (p1 - p2) : (p1 - p2) - c.ci.lo;
  return c;
}

// A one-sided bound's distance from the estimate within the rounding of the bound, as at a confidence a rounding
// error above 50%: its sign then means nothing.
bool near_estimate(double w, double p1, double p2) { return std::fabs(w) <= 1e-13 * std::fmax(p1, p2); }

// Why an Agresti–Caffo bound can lie on the wrong side of the estimate (a distance below 0), as at a confidence not
// far above 50%: each group's proportion is taken as (x + 1)/(n + 2), (1 − 2p)/(n + 2) nearer 1/2, and the center
// moves by the difference of the two, toward no difference with equal groups but either way with unequal ones
// (P1 = P2 = 0.9 with 10 and 100 subjects moves it from 0 to −0.059); the other intervals contain the estimate.
const char* const AC_WRONG_SIDE =
    "the success and failure that the Agresti–Caffo interval adds to each group pull each group's proportion toward "
    "1/2, which moves its center from the estimate by more than its margin at this confidence";

// How the limits for P1 − P2 come about at the anticipated counts, in the method's own formula, and the width;
// solving for N1, the Wald interval's closed form for comparison, the width at the sizes found and one size below.
derive::Derivation diff_ci_derivation(const Scenario& s, const Row& r) {
  using namespace derive;
  const std::string& m = s.s("method");
  const double n1 = r.num("n1"), n2 = r.num("n2"), p1 = r.num("p1"), p2 = r.num("p2"), conf = r.num("conf");
  const bool two = s.is("sides", "two"), up = s.is("sides", "upper"), solved = s.solve == "n";
  const double x1 = n1 * p1, x2 = n2 * p2, h1 = x1 / n1, h2 = x2 / n2, d = h1 - h2;
  const double a = two ? (1 - conf) / 2 : 1 - conf, z = stats::qnorm(a, false);
  const Interval ci = diff_ci(m, n1, x1, n2, x2, conf, two);
  const std::string N1 = size(n1, solved), N2 = dvi(n2), P1 = gv(p1), P2 = gv(p2), R = gv(s("r"));
  const std::string level = two ? "α/2" : "α";
  const double w = r.num("width"), target = s.get("width", NAN);
  const std::string what = two ? "the interval's width W" : "the bound's distance W from the estimate";
  // both limits, or the one a one-sided interval has
  auto limits = [&](const std::string& l, const std::string& u) { return two ? l + "; " + u : (up ? u : l); };
  const std::string lims = limits("L = " + dv(ci.lo), "U = " + dv(ci.hi));
  Derivation out;
  if (solved) {
    out.find = "the smallest N1, with N2 = ⌈R × N1⌉, from which " + what + " stays at or below " + gv(target);
    // the size of the Wald interval, beside the answer below (the search starts from it)
    out.steps.push_back(
        {"Normal approximation",
         std::string("N1 ≈ ") + (two ? "4 " : "") + "z² (P1(1 − P1) + P2(1 − P2)/R)/W², z = Φ⁻¹(1 − " + level +
             "), the Wald interval's size, for comparison; the steps below give the answer",
         (two ? "4 × " : "") + dv(z) + "² × (" + P1 + " × " + dv(1 - p1) + " + " + P2 + " × " + dv(1 - p2) + "/" + R +
             ") / " + gv(target) + "² = " +
             dv((two ? 4 : 1) * z * z * (p1 * (1 - p1) + p2 * (1 - p2) / s("r")) / (target * target))});
  } else {
    out.find = what + " with N1 = " + N1 + " and N2 = ⌈R × N1⌉";
  }
  out.steps.push_back({"Group 2 size", "N2 = ⌈R × N1⌉", "⌈" + R + " × " + N1 + "⌉ = " + N2});
  out.steps.push_back(
      {"Anticipated counts", "x1 = N1 P1, x2 = N2 P2",
       "x1 = " + N1 + " × " + P1 + " = " + dv(x1) + ", x2 = " + N2 + " × " + P2 + " = " + dv(x2)});
  out.steps.push_back({"Critical value", "z = Φ⁻¹(1 − " + level + "), α = 1 − confidence",
                       "Φ⁻¹(1 − " + dv(a) + ") = " + dv(z)});
  if (m == "newcombe" || m == "newcombecc") {
    const std::string one = m == "newcombe" ? "wilson" : "wilsoncc";
    const Interval c1 = prop_ci(one, n1, x1, conf, two), c2 = prop_ci(one, n2, x2, conf, two);
    out.steps.push_back({"Limits of each proportion",
                         std::string("[l1, u1] and [l2, u2], the Wilson score intervals") +
                             (m == "newcombecc" ? " with continuity correction" : "") +
                             " for P1 and P2 at confidence 1 − α",
                         "l1 = " + dv(c1.lo) + ", u1 = " + dv(c1.hi) + ", l2 = " + dv(c2.lo) + ", u2 = " + dv(c2.hi)});
    out.steps.push_back(
        {"Limits",
         limits("L = D − √((P̂1 − l1)² + (u2 − P̂2)²)", "U = D + √((u1 − P̂1)² + (P̂2 − l2)²)") +
             ", D = P̂1 − P̂2, P̂i = xi/Ni",
         limits("L = " + dvs(d) + " − √((" + dv(h1) + " − " + dv(c1.lo) + ")² + (" + dv(c2.hi) + " − " + dv(h2) +
                    ")²) = " + dv(ci.lo),
                "U = " + dvs(d) + " + √((" + dv(c1.hi) + " − " + dv(h1) + ")² + (" + dv(h2) + " − " + dv(c2.lo) +
                    ")²) = " + dv(ci.hi))});
  } else if (m == "ac") {
    const double t1 = (x1 + 1) / (n1 + 2), t2 = (x2 + 1) / (n2 + 2);
    const double h = z * std::sqrt(t1 * (1 - t1) / (n1 + 2) + t2 * (1 - t2) / (n2 + 2));
    out.steps.push_back({"Adjusted proportions",
                         "P̃i = (xi + 1)/(Ni + 2), a success and a failure added to each group",
                         "P̃1 = (" + dv(x1) + " + 1)/(" + N1 + " + 2) = " + dv(t1) + ", P̃2 = (" + dv(x2) + " + 1)/(" +
                             N2 + " + 2) = " + dv(t2)});
    out.steps.push_back({"Limits",
                         std::string(two ? "P̃1 − P̃2 ± H" : up ? "U = P̃1 − P̃2 + H" : "L = P̃1 − P̃2 − H") +
                             ", H = z √(P̃1(1 − P̃1)/(N1 + 2) + P̃2(1 − P̃2)/(N2 + 2)), kept within −1 and 1",
                         "H = " + dv(z) + " × √(" + dv(t1) + " × " + dv(1 - t1) + "/(" + N1 + " + 2) + " + dv(t2) +
                             " × " + dv(1 - t2) + "/(" + N2 + " + 2)) = " + dv(h) + ": " + lims});
  } else {
    const bool cc = m == "waldcc";
    double h = z * std::sqrt(h1 * (1 - h1) / n1 + h2 * (1 - h2) / n2);
    if (cc) h += 0.5 * (1 / n1 + 1 / n2);
    out.steps.push_back({"Limits",
                         std::string(two ? "D ± H" : up ? "U = D + H" : "L = D − H") +
                             ", H = z √(P̂1(1 − P̂1)/N1 + P̂2(1 − P̂2)/N2)" + (cc ? " + (1/N1 + 1/N2)/2" : "") +
                             ", D = P̂1 − P̂2, kept within −1 and 1",
                         "D = " + dv(d) + ", H = " + dv(z) + " × √(" + dv(h1) + " × " + dv(1 - h1) + "/" + N1 + " + " +
                             dv(h2) + " × " + dv(1 - h2) + "/" + N2 + ")" +
                             (cc ? " + (1/" + N1 + " + 1/" + N2 + ")/2" : "") + " = " + dv(h) + ": " + lims});
  }
  const std::string shown = solved ? prob(w, target) : num(w), wt = derived(shown);
  if (two)
    out.steps.push_back({"Width", "W = U − L", difference(ci.hi, ci.lo, shown) + " = " + wt});
  else if (up)
    out.steps.push_back({"Width", "W = U − (P1 − P2)", difference(ci.hi, p1 - p2, shown) + " = " + wt});
  else
    out.steps.push_back({"Width", "W = (P1 − P2) − L", difference(p1 - p2, ci.lo, shown) + " = " + wt});
  if (solved) {
    const std::string rule =
        "N1 = the smallest size from which the width stays at or below the target, with N2 = ⌈R × N1⌉";
    const std::string at = "width at N1 = " + dvi(n1) + ", N2 = " + N2 + ": " + wt;
    if (n1 <= 1) {
      out.steps.push_back({"Smallest sizes", rule,
                           "at N1 = " + dvi(n1) + ", N2 = " + N2 + ", the smallest sizes, the width is already " + wt +
                               " ≤ " + gv(target)});
    } else {
      const DiffCi below = diff_ci_at(s, n1 - 1);
      const std::string less = "N1 = " + dvi(n1 - 1) + ", N2 = " + dvi(below.n2);
      if (below.width > target)
        out.steps.push_back({"Smallest sizes", rule,
                             "width at " + less + ": " + dvp(below.width, target) + " > " + gv(target) + " ≥ " + at});
      else  // a one-sided bound on the wrong side of the estimate, or within a rounding error of it
        out.steps.push_back({"Smallest sizes", rule + ", with the bound beyond the estimate",
                             "at " + less + " the bound lies " +
                                 (near_estimate(below.width, p1, p2)
                                      ? "within a rounding error of the estimate"
                                      : "on the wrong side of the estimate (W = " + dv(below.width) + ")") +
                                 "; " + at + " ≤ " + gv(target)});
    }
    out.answer = "N1 = " + dvi(n1) + ", N2 = " + N2 + " (N = " + dvi(n1 + n2) + "), with W = " + wt;
  } else {
    out.answer = "W = " + wt + " (" + lims + ")";
  }
  return out;
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
    if (!s.is("sides", "two")) check_one_sided_conf(s("conf"));
    double n1;
    if (s.solve == "n") {
      const double target = s("width");
      const bool two = s.is("sides", "two");
      const double z = stats::qnorm((1 - s("conf")) / (two ? 2 : 1), false);
      const double p1 = s("p1"), p2 = s("p2");
      const double e = two ? target / 2 : target;
      solve::NSearch o;
      o.hi = 1e8;
      o.lo = 1;
      o.guess = std::ceil(z * z * (p1 * (1 - p1) + p2 * (1 - p2) / s("r")) / (e * e));
      o.goal = "the target width";
      // a one-sided bound at the estimate, on its wrong side or within a rounding error of it is a degenerate size;
      // a two-sided interval at a confidence near 0 closes on the estimate itself (a width of 0, or a rounding error
      // either side of it) at every size
      auto ok = [&](double x) {
        const double w = diff_ci_at(s, x).width;
        return (two || (w > 0 && !near_estimate(w, p1, p2))) && w <= target;
      };
      try {
        n1 = solve::min_n_staying(ok, o, &ctx);
      } catch (const CalcError&) {
        // a one-sided bound that is still on the wrong side of the estimate, or within a rounding error of it, at
        // the largest size: the reason, not only that the width was not reached (said of that size alone: a smaller
        // one can put the bound on the right side, as N1 = 1 at R = 0.5 does by giving the groups equal sizes)
        const double w = diff_ci_at(s, o.hi).width;
        if (two || (w > 0 && !near_estimate(w, p1, p2))) throw;
        if (near_estimate(w, p1, p2))
          throw CalcError("the confidence is so close to 50% that the bound lies within a rounding error of the "
                          "estimate, where its distance cannot be computed");
        throw CalcError(std::string("the bound still lies on the wrong side of the estimate with ") +
                        text::format_g(o.hi, 10) + " subjects in group 1, the largest sample size searched: " +
                        (s.is("method", "ac") ? AC_WRONG_SIDE : "at this confidence"));
      }
      // the size below the answer, where the bound lies on the wrong side or within a rounding error of the estimate
      // (with P1 = 0.989, P2 = 0.176 and 50.1%, 5.7e-14 beyond it at 50 digits: its sign in doubles means nothing,
      // and the reason need not hold)
      if (!two && n1 > o.lo && s.is("method", "ac")) {
        const DiffCi below = diff_ci_at(s, n1 - 1);
        const std::string sizes = "With " + fint(n1 - 1) + " and " + fint(below.n2) + " subjects (one fewer in group 1)";
        if (near_estimate(below.width, p1, p2))
          r.notes.push_back(sizes + " the bound lies within a rounding error of the estimate.");
        else if (below.width < 0)
          r.notes.push_back(sizes + " the bound lies on the wrong side of the estimate (" + std::string(AC_WRONG_SIDE) +
                            ").");
      }
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
    const auto [p1, p2] = fnum_apart(r.num("p1"), r.num("p2"));
    const double w = r.num("width");
    // the width with the digits that keep it to within 5% (0.00016, not 0.0002)
    std::string where = (two ? "with a width of " : "whose bound lies at a distance of ") + fnum_apart(w, 0).first;
    if (!two && near_estimate(w, r.num("p1"), r.num("p2")))
      where = "whose bound lies within a rounding error of the estimate";
    else if (!two && w < 0)
      where = "whose bound lies " + fnum_apart(-w, 0).first + (s.is("sides", "upper") ? " below" : " above") +
              " the estimate, on the wrong side of it (" + (s.is("method", "ac") ? AC_WRONG_SIDE : "at this confidence") +
              ")";
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " give a " +
           (two ? "two-sided " : "one-sided ") + level_text(r.num("conf")) + " " + names.at(s.s("method")) +
           " confidence interval for the difference of two proportions " + where + " when the observed proportions "
           "equal " + p1 + " and " + p2 + ".";
  };
  p.derivation = diff_ci_derivation;
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

// The exposures among cases and controls as percentages with one decimal, or with the decimals it takes to keep
// their difference to within 5% and to show neither as 0% or 100%, which it is not (99.9995% among cases and
// 99.999% among controls, not 100.0% for both).
std::pair<std::string, std::string> pct_apart(double a, double b) {
  auto pct = [](double p, int d) {
    const double x = 100 * p;
    return x != 0 && x < 1e-3 ? text::format_g(x, d) : text::format_fixed(x, d);
  };
  auto off = [](double p, double shown) { return (shown == 0 && p != 0) || (shown == 100 && p != 1); };
  std::string ta, tb;
  for (int d = 1; d <= 17; ++d) {
    ta = pct(a, d);
    tb = pct(b, d);
    double va, vb;
    if (text::parse_double(ta, va) && text::parse_double(tb, vb) && !off(a, va) && !off(b, vb) &&
        std::fabs((va - vb) - 100 * (a - b)) <= 0.05 * 100 * std::fabs(a - b))
      break;
  }
  return {ta + "%", tb + "%"};
}

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

// The cases as group 1 and the controls as group 0 of the test of two proportions, following what was solved for:
// solving for the cases, Schlesselman's closed form for comparison, the power at the numbers found and one case
// fewer; solving for the odds ratio, the exposure among cases at which the power equals the target, and from it the
// odds ratio.
derive::Derivation case_control_derivation(const Scenario& s, const Row& r) {
  using namespace derive;
  const TwoPropTest T{s.s("test"), alt(s), s("alpha")};
  const double nc = r.num("ncase"), nk = r.num("nctrl"), p1 = r.num("p1"), p0 = r.num("p0"), m = r.num("m");
  const double odds_ratio = r.num("or"), target = s.get("power", NAN), power = r.num("power");
  const bool fell_back = has_note(r, FELL_BACK), enumerated = s.is("method", "exact") && !fell_back;
  const bool solved = s.solve == "n", effect = s.solve == "or";
  const Groups g{"1", "0", size(nc, solved), dvi(nk), mv(p1, !effect && s.is("es", "p1")), gv(p0)};
  const std::string with =
      "with N1 = " + gvi(nc) + " cases and M = " + gv(m) + (m == 1 ? " control" : " controls") + " per case";
  Derivation out;
  if (solved) {
    out.find = find_size("number of cases N1, with N0 = ⌈M × N1⌉ controls,", target);
  } else if (effect) {
    const bool up = T.alt == Alt::Upper || (T.alt == Alt::Two && !s.is("ordir", "below"));
    out.find = find_effect(std::string("odds ratio OR ") + (up ? "above" : "below") + " 1", target, with);
  } else {
    out.find = find_power(with);
  }
  if (!effect && s.is("es", "or"))
    out.steps.push_back({"Exposure among cases", "P1 = OR P0/(1 − P0 + OR P0)",
                         gv(odds_ratio) + " × " + gv(p0) + "/(" + dv(1 - p0) + " + " + gv(odds_ratio) + " × " + gv(p0) +
                             ") = " + dv(p1)});
  out.steps.push_back({"Controls", "N0 = ⌈M × N1⌉ for N1 cases and M controls per case",
                       "⌈" + gv(m) + " × " + g.n1 + "⌉ = " + g.n2});
  if (solved || !enumerated) out.steps.push_back(diff_step(g, p1, p0));
  if (solved)
    out.steps.push_back(two_prop_formula(T.alt, T.alpha, target, p1, p0, m, "M", T.test != "zunpool",
                                         T.test == "zpoolcc" || T.test == "fisher", "1", "0", g.p1, g.p2, "N1"));
  two_prop_power_steps(T, nc, nk, p1, p0, enumerated, fell_back, r, g, out.steps);
  if (solved) {
    RunContext ctx;
    const double less = nc - 1, less_ctrl = solve::round_up(m * less);
    auto sizes = [](double a, double b) {
      return "N1 = " + dvi(a) + ", N0 = " + dvi(b) + " (N = " + dvi(a + b) + ")";
    };
    const std::string rule = "N1 = the smallest number of cases whose power reaches 1 − β, with N0 = ⌈M × N1⌉";
    out.steps.push_back(size_step(
        "Smallest numbers", enumerated ? sawtooth_rule(rule) : rule, sizes(nc, nk), power,
        less < 2 || less_ctrl < 2 ? "" : sizes(less, less_ctrl),
        less < 2 || less_ctrl < 2 ? NAN : two_prop_power(T, less, less_ctrl, p1, p0, enumerated, ctx), target));
    out.answer = "N1 = " + dvi(nc) + " cases and N0 = " + dvi(nk) + " controls (N = " + dvi(nc + nk) +
                 "), with power " + dvp(power, target) + alpha_star(r);
  } else if (effect) {
    // the search finds the exposure among cases; the odds ratio follows from it
    out.steps.push_back(effect_search("P1", p1, power, target));
    out.steps.push_back({"Odds ratio", "OR = (P1/(1 − P1)) / (P0/(1 − P0))",
                         "(" + dv(p1) + "/" + dv(1 - p1) + ") / (" + gv(p0) + "/" + dv(1 - p0) + ") = " +
                             dv(odds_ratio)});
    out.answer = "OR = " + dv(odds_ratio) + " (P1 = " + dv(p1) + ")";
  } else {
    out.answer = "power = " + dvp(power) + alpha_star(r);
  }
  return out;
}

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
    ExactSearch search;
    // M × cases rounded up, at least 1 however small M is
    auto controls = [&](double cases) { return solve::round_up(m * cases); };
    auto power = [&](double n1, double n2, double x) {
      // a point of the power curve is enumerated only while all the curves stay within the exact search's budget,
      // judged at its curve's largest point
      const double last1 = ctx.curve_last ? (*ctx.curve_last)("ncase") : n1;
      const double last2 = ctx.curve_last ? controls(last1) : n2;
      if (exact && !on_curve_too_long(ctx, exact_power_work(last1, last2, x, p0, T.test == "fisher"))) {
        search.count(exact_power_work(n1, n2, x, p0, T.test == "fisher"));
        const double pw = two_prop_exact(T, n1, n2, x, p0, ctx);
        if (!std::isnan(pw)) return pw;
        fell_back = true;
      }
      return two_prop_normal(T, n1, n2, x, p0);
    };
    TwoN n;
    double p1;
    if (s.solve == "n") {
      p1 = s.is("es", "p1") ? s("p1") : case_exposure(s("or"), p0);
      if (p1 == p0) throw CalcError("the odds ratio is 1: cases and controls have the same exposure");
      if ((a == Alt::Upper && p1 < p0) || (a == Alt::Lower && p1 > p0))
        throw CalcError("the odds ratio is on the side opposite to the one-sided alternative");
      // Schlesselman's formula for the number of cases
      const double za = stats::qnorm(tail_alpha(s("alpha"), a), false), zb = stats::qnorm(s("power"));
      const double pb = (p1 + m * p0) / (1 + m);
      const double g = std::pow(za * std::sqrt((1 + 1 / m) * pb * (1 - pb)) + zb * std::sqrt(p1 * (1 - p1) + p0 * (1 - p0) / m), 2) /
                       ((p1 - p0) * (p1 - p0));
      if (exact && exact_search_work(g, m * g, p1, p0, T.test == "fisher") > EXACT_SEARCH_BUDGET)
        exact = false, fell_back = true;
      // the search starts where both groups have 2 (in the words of this procedure, not of a two-group allocation)
      if (controls(1e7) < 2) throw CalcError("M is too small: even 10000000 cases would have fewer than 2 controls");
      Scenario t = s;  // the two-group search with controls = M × cases
      t.sel["alloc"] = "ratio";
      t.num["r"] = m;
      n = search.run(exact, fell_back, [&] {
        return solve_two_group_n(
            t, [&](double c, double k) { return power(c, k, p1); }, s("power"), 2, std::ceil(g), !exact, ctx, 1e7);
      });
    } else {
      n = {s("ncase"), controls(s("ncase"))};
      if (s.solve == "power") {
        p1 = s.is("es", "p1") ? s("p1") : case_exposure(s("or"), p0);
      } else {
        const bool up = a == Alt::Upper || (a == Alt::Two && !s.is("ordir", "below"));
        const double lo = up ? p0 : 1e-9, hi = up ? 1 - 1e-9 : p0;
        // the search evaluates the power some 40 times: past the exact search's budget (judged at the approximate
        // answer), the normal approximation
        if (exact) {
          auto g = [&](double x) { return two_prop_normal(T, n.n1, n.n2, x, p0) - s("power"); };
          if (g(p0) < 0 && g(up ? hi : lo) >= 0 &&
              40 * exact_power_work(n.n1, n.n2, numerics::brent(g, lo, hi, 1e-9), p0, T.test == "fisher") >
                  EXACT_SEARCH_BUDGET)
            exact = false, fell_back = true;
        }
        auto f = [&](double x) { return power(n.n1, n.n2, x) - s("power"); };
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
    // a given odds ratio as given: recomputed from a P1 near 1, it would lose its digits with 1 − P1 (and be
    // infinite once P1 rounds to 1)
    r.set("or", s.solve != "or" && s.is("es", "or") ? s("or") : odds_ratio_of(p1, p0));
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
    const auto [e1, e0] = pct_apart(r.num("p1"), r.num("p0"));
    std::string t = fint(r.num("ncase")) + " cases and " + fint(r.num("nctrl")) + " controls " + achieves(s, r) +
                    " to detect an odds ratio of exposure of " + fnum_from(r.num("or"), 1, 3) + " (exposure " + e1 +
                    " among cases and " + e0 + " among controls) using a " + alt_words(alt(s)) + " " +
                    names.at(s.s("test")) + " at a target significance level of " + fnum(r.num("alpha"));
    if (r.has("actual_alpha")) t += " (actual " + fprob(r.num("actual_alpha")) + ")";
    const bool enumerated = s.is("method", "exact") && std::find(r.notes.begin(), r.notes.end(), FELL_BACK) == r.notes.end();
    return t + (enumerated ? "; power was computed by exact enumeration." : "; power uses the normal approximation.");
  };
  p.derivation = case_control_derivation;
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
  double pp = 0, pm = 0;             // a control's exposure when its case is exposed, and when not
};

Matched matched_moments(double odds_ratio, double p0, double phi, int M) {
  // 1 − P1 and 1 − e below in the forms that keep their digits when the odds ratio is large (from about 1e16 they
  // rounded to 0: the variance under H1 was 0, and the conditional exposure pm no number)
  const double q0 = 1 - p0, p1 = case_exposure(odds_ratio, p0), q1 = q0 / (q0 + odds_ratio * p0);
  // exposure of a control given that its case is exposed (plus) or not (minus)
  const double pp = p0 + phi * std::sqrt(q1 * p0 * q0 / p1), pm = p0 - phi * std::sqrt(p1 * p0 * q0 / q1);
  if (!(pp >= 0 && pp <= 1 && pm >= 0 && pm <= 1))
    throw CalcError("the correlation φ is too large in size for these exposures: the exposure of a control given its "
                    "case's would fall outside 0 to 1");
  Matched mo;
  mo.pp = pp;
  mo.pm = pm;
  for (int t = 1; t <= M; ++t) {
    // the chance that a set has t exposed members: the case and t − 1 controls, or t controls
    const double pt = p1 * std::exp(stats::lchoose(M, t - 1)) * std::pow(pp, t - 1) * std::pow(1 - pp, M - t + 1) +
                      q1 * std::exp(stats::lchoose(M, t)) * std::pow(pm, t) * std::pow(1 - pm, M - t);
    const double e = t * odds_ratio / (t * odds_ratio + M - t + 1), e0 = t / (M + 1.0);
    mo.delta += pt * (e - e0);
    mo.v0 += pt * e0 * (1 - e0);
    mo.v1 += pt * e * (M - t + 1) / (t * odds_ratio + M - t + 1);
  }
  return mo;
}

// The power and its parts (which the derivation shows): the per-set moments, √n, their SDs, the critical value and
// each tail.
struct MatchedPower {
  Matched mo;
  double rn, s0, s1, z, up = NAN, lo = NAN, power = NAN;
};

MatchedPower matched_power_parts(double n, double odds_ratio, double p0, double phi, int M, double alpha, Alt a) {
  using namespace stats;
  MatchedPower r;
  r.mo = matched_moments(odds_ratio, p0, phi, M);
  r.rn = std::sqrt(n);
  r.s0 = std::sqrt(r.mo.v0);
  r.s1 = std::sqrt(r.mo.v1);
  r.z = qnorm(tail_alpha(alpha, a), false);
  if (!(r.s1 > 0)) return r;
  r.up = pnorm((r.rn * r.mo.delta - r.z * r.s0) / r.s1);
  r.lo = pnorm((r.rn * -r.mo.delta - r.z * r.s0) / r.s1);
  r.power = a == Alt::Two ? std::min(1.0, r.up + r.lo) : (a == Alt::Upper ? r.up : r.lo);
  return r;
}

double matched_power(double n, double odds_ratio, double p0, double phi, int M, double alpha, Alt a) {
  return matched_power_parts(n, odds_ratio, p0, phi, M, alpha, a).power;
}

// Dupont's moments of a matched set and the normal power of the conditional test, following what was solved for:
// solving for the sets, Dupont's closed form for comparison, the power at the number found and one fewer; solving for
// the odds ratio, the power at the odds ratio found, where it equals the target.
derive::Derivation matched_derivation(const Scenario& s, const Row& r) {
  using namespace derive;
  const Alt a = alt(s);
  const double n = r.num("n"), odds_ratio = r.num("or"), p0 = r.num("p0"), p1 = r.num("p1"), phi = r.num("phi");
  const double alpha = s("alpha"), target = s.get("power", NAN), power = r.num("power");
  const int M = static_cast<int>(r.num("mctl"));
  const bool solved = s.solve == "n", effect = s.solve == "or";
  const MatchedPower q = matched_power_parts(n, odds_ratio, p0, phi, M, alpha, a);
  const std::string OR = mv(odds_ratio, !effect), P0 = gv(p0), P1 = dv(p1), F = gvs(phi);
  const std::string with =
      "with n = " + gvi(n) + " matched sets of one case and M = " + gvi(M) + (M == 1 ? " control" : " controls");
  Derivation out;
  if (solved) {
    out.find = find_size("number of matched sets n", target);
  } else if (effect) {
    const bool up = a == Alt::Upper || (a == Alt::Two && !s.is("ordir", "below"));
    out.find = find_effect(std::string("odds ratio OR ") + (up ? "above" : "below") + " 1", target, with);
  } else {
    out.find = find_power(with);
  }
  out.steps.push_back({"Exposure among cases", "P1 = OR P0/(1 − P0 + OR P0)",
                       OR + " × " + P0 + "/(" + dv(1 - p0) + " + " + OR + " × " + P0 + ") = " + P1});
  out.steps.push_back(
      {"Exposure of a control given its case's",
       "p₊ = P0 + φ √((1 − P1) P0 (1 − P0)/P1) when the case is exposed, p₋ = P0 − φ √(P1 P0 (1 − P0)/(1 − P1)) "
       "when not",
       "p₊ = " + P0 + " + " + F + " × √(" + dv(1 - p1) + " × " + P0 + " × " + dv(1 - p0) + "/" + P1 + ") = " +
           dv(q.mo.pp) + ", p₋ = " + P0 + " − " + F + " × √(" + P1 + " × " + P0 + " × " + dv(1 - p0) + "/" +
           dv(1 - p1) + ") = " + dv(q.mo.pm)});
  out.steps.push_back(
      {"Per-set moments",
       "δ = Σₜ πₜ (eₜ − t/(M + 1)), v0 = Σₜ πₜ t/(M + 1) × (1 − t/(M + 1)), v1 = Σₜ πₜ eₜ (M − t + 1)/(t OR + M − t + "
       "1), for sets with t = 1 to M exposed members (πₜ their chance from P1, p₊ and p₋), eₜ = t OR/(t OR + M − t + 1) "
       "the chance that the case is one of them",
       "M = " + gvi(M) + ": δ = " + dsig(q.mo.delta) + ", v0 = " + dsig(q.mo.v0) + ", v1 = " + dsig(q.mo.v1)});
  if (solved) {
    // Dupont's closed form, beside the answer below (the search starts from it)
    const double za = stats::qnorm(tail_alpha(alpha, a), false), zb = stats::qnorm(target);
    out.steps.push_back(
        {"Normal approximation",
         std::string("n ≈ (z(") + level_symbol(a) +
             ") √v0 + z(1 − β) √v1)² / δ², for comparison; the steps below give the answer",
         "(" + dv(za) + " × √" + dsig(q.mo.v0) + " + " + dvs(zb) + " × √" + dsig(q.mo.v1) + ")² / " +
             dsig(q.mo.delta) + "² = " +
             dv(std::pow(za * std::sqrt(q.mo.v0) + zb * std::sqrt(q.mo.v1), 2) / (q.mo.delta * q.mo.delta))});
  }
  out.steps.push_back(z_critical(alpha, a));
  out.steps.push_back(tails_power(a, "√n δ", "", "√v0", "√v1", q.rn * q.mo.delta, 0, q.z, q.s0, q.s1, q.up, q.lo,
                                  power));
  if (solved) {
    out.steps.push_back(
        n_search([&](double x) { return matched_power(x, odds_ratio, p0, phi, M, alpha, a); }, n, target, 1, "n"));
    out.answer =
        "n = " + dvi(n) + " matched sets (" + dvi(n * (M + 1)) + " subjects), with power " + dvp(power, target);
  } else if (effect) {
    out.steps.push_back(effect_search("OR", odds_ratio, power, target));
    out.answer = "OR = " + dv(odds_ratio) + " (P1 = " + P1 + ")";
  } else {
    out.answer = "power = " + dvp(power);
  }
  return out;
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
      const double za = stats::qnorm(tail_alpha(alpha, a), false), zb = stats::qnorm(s("power"));
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
    const auto [e1, e0] = pct_apart(r.num("p1"), r.num("p0"));
    return fint(r.num("n")) + " matched sets of one case and " + (M == 1 ? "one control" : fint(M) + " controls") +
           " (" + fint(r.num("ntotal")) + " subjects) " + achieves(s, r) + " to detect an odds ratio of exposure of " +
           fnum_from(r.num("or"), 1, 3) + " (exposure " + e1 + " among cases and " + e0 +
           " among controls, correlation of exposure within sets " + fnum_fine(r.num("phi")) + ") using a " +
           alt_words(alt(s)) + " conditional test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.derivation = matched_derivation;
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
