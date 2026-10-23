// Simulation-based power for comparing an ordinal outcome between two groups
// (sim.ordinal): ordered categories such as the modified Rankin Scale, a
// pain grade or a Likert item.  The control group's category probabilities
// are given; the treatment group's follow from a common odds ratio
// (proportional odds) or are given too, of any shape.
//
// Each simulated study is analyzed with the proportional odds model
// (cumulative logit, fitted by Newton's method; Wald or likelihood-ratio
// test of the group effect, as MASS::polr) or the Mann–Whitney test with
// ties.  The quick approximation is Whitehead's (1993) formula for the
// proportional odds model and Noether's formula with ties for given
// probabilities.
#include "procedures/sim_ordinal.hpp"

#include <algorithm>
#include <array>
#include <cstdint>

#include "procedures/sim_common.hpp"

namespace ssc::procs {
namespace {

double expit(double t) { return 1 / (1 + std::exp(-t)); }
double logit(double p) { return std::log(p / (1 - p)); }

}  // namespace

namespace ordinal {

PoFit po_fit(const double* group1, const double* group2, int k) {
  PoFit f;
  std::array<std::array<double, CMAX>, 2> n{};
  int m = 0;  // categories with observations
  for (int j = 0; j < k; ++j)
    if (group1[j] + group2[j] > 0) {
      n[0][m] = group1[j];
      n[1][m] = group2[j];
      ++m;
    }
  double n1 = 0, n2 = 0, total = 0;
  for (int j = 0; j < m; ++j) {
    n1 += n[0][j];
    n2 += n[1][j];
  }
  total = n1 + n2;
  // the null model: one multinomial for both groups
  for (int j = 0; j < m; ++j) {
    const double c = n[0][j] + n[1][j];
    f.loglik0 += c * std::log(c / total);
  }
  if (m < 2 || n1 == 0 || n2 == 0) {
    f.single = true;
    return f;
  }
  // separation: one group entirely at or above the other (overlapping in one category at most); the estimate
  // is infinite and the likelihood's supremum is the two groups' own multinomials
  int first[2] = {m, m}, last[2] = {-1, -1};
  for (int g = 0; g < 2; ++g)
    for (int j = 0; j < m; ++j)
      if (n[g][j] > 0) {
        first[g] = std::min(first[g], j);
        last[g] = j;
      }
  if (last[0] <= first[1] || last[1] <= first[0]) {
    f.separated = true;
    f.ok = true;
    f.beta = last[0] <= first[1] ? -INFINITY : INFINITY;
    for (int g = 0; g < 2; ++g)
      for (int j = 0; j < m; ++j)
        if (n[g][j] > 0) f.loglik += n[g][j] * std::log(n[g][j] / (g == 0 ? n1 : n2));
    return f;
  }
  const int p = m;  // θ_1 … θ_{m−1} and β
  std::array<double, CMAX + 1> par{};
  double cum = 0;
  for (int j = 0; j + 1 < m; ++j) {
    cum += n[0][j] + n[1][j];
    par[j] = logit(cum / total);
  }
  par[m - 1] = 0;
  auto loglik = [&](const std::array<double, CMAX + 1>& q) -> double {
    double ll = 0;
    for (int g = 0; g < 2; ++g) {
      const double eta = g == 0 ? q[m - 1] : 0;
      double below = 0;
      for (int j = 0; j < m; ++j) {
        const double up = j + 1 < m ? expit(q[j] - eta) : 1;
        const double prob = up - below;
        if (n[g][j] > 0) {
          if (!(prob > 0)) return -INFINITY;
          ll += n[g][j] * std::log(prob);
        }
        below = up;
      }
    }
    return ll;
  };
  double ll = loglik(par);
  std::array<double, (CMAX + 1) * (CMAX + 1)> h{};
  std::array<double, CMAX + 1> grad{};
  bool converged = false;
  for (int iter = 0; iter < 100; ++iter) {
    h.fill(0);
    grad.fill(0);
    for (int g = 0; g < 2; ++g) {
      const double x = g == 0 ? 1 : 0, eta = x * par[m - 1];
      for (int j = 0; j < m; ++j) {
        const double c = n[g][j];
        if (c == 0) continue;
        // P = F(a) − F(b), a = θ_j − η (∞ for the top category), b = θ_{j−1} − η (−∞ for the bottom one)
        const bool top = j + 1 == m, bottom = j == 0;
        const double Fa = top ? 1 : expit(par[j] - eta), Fb = bottom ? 0 : expit(par[j - 1] - eta);
        const double fa = top ? 0 : Fa * (1 - Fa), fb = bottom ? 0 : Fb * (1 - Fb);
        const double da = top ? 0 : fa * (1 - 2 * Fa), db = bottom ? 0 : fb * (1 - 2 * Fb);
        const double P = Fa - Fb, P2 = P * P;
        const int ia = j, ib = j - 1, ibeta = m - 1;
        if (!top) grad[ia] += c * fa / P;
        if (!bottom) grad[ib] -= c * fb / P;
        grad[ibeta] -= c * x * (fa - fb) / P;
        auto add = [&](int r, int s, double v) {
          h[r * p + s] += c * v;
          if (r != s) h[s * p + r] += c * v;
        };
        if (!top) add(ia, ia, da / P - fa * fa / P2);
        if (!bottom) add(ib, ib, -db / P - fb * fb / P2);
        if (!top && !bottom) add(ia, ib, fa * fb / P2);
        if (x != 0) {
          if (!top) add(ia, ibeta, -(da / P - fa * (fa - fb) / P2));
          if (!bottom) add(ib, ibeta, -(-db / P + fb * (fa - fb) / P2));
          add(ibeta, ibeta, (da - db) / P - (fa - fb) * (fa - fb) / P2);
        }
      }
    }
    // Newton step: solve (−H) d = grad by Cholesky
    std::array<double, (CMAX + 1) * (CMAX + 1)> l{};
    for (int i = 0; i < p * p; ++i) l[i] = -h[i];
    bool pd = true;
    for (int j = 0; j < p && pd; ++j) {
      double s = l[j * p + j];
      for (int q = 0; q < j; ++q) s -= l[j * p + q] * l[j * p + q];
      if (!(s > 0)) {
        pd = false;
        break;
      }
      l[j * p + j] = std::sqrt(s);
      for (int i = j + 1; i < p; ++i) {
        double v = l[i * p + j];
        for (int q = 0; q < j; ++q) v -= l[i * p + q] * l[j * p + q];
        l[i * p + j] = v / l[j * p + j];
      }
    }
    if (!pd) return f;
    if (converged) {
      // covariance of the estimates at the maximum: the inverse of −H, from its Cholesky factor (the β entry)
      std::array<double, CMAX + 1> e{};
      e[m - 1] = 1;
      for (int i = 0; i < p; ++i) {
        for (int q = 0; q < i; ++q) e[i] -= l[i * p + q] * e[q];
        e[i] /= l[i * p + i];
      }
      for (int i = p - 1; i >= 0; --i) {
        for (int q = i + 1; q < p; ++q) e[i] -= l[q * p + i] * e[q];
        e[i] /= l[i * p + i];
      }
      f.ok = std::fabs(par[m - 1]) < 50 && e[m - 1] > 0;
      f.beta = par[m - 1];
      f.se = std::sqrt(e[m - 1]);
      f.loglik = ll;
      return f;
    }
    std::array<double, CMAX + 1> d = grad;
    for (int i = 0; i < p; ++i) {
      for (int q = 0; q < i; ++q) d[i] -= l[i * p + q] * d[q];
      d[i] /= l[i * p + i];
    }
    for (int i = p - 1; i >= 0; --i) {
      for (int q = i + 1; q < p; ++q) d[i] -= l[q * p + i] * d[q];
      d[i] /= l[i * p + i];
    }
    double step = 1, gain = 0;
    for (int i = 0; i < p; ++i) gain += grad[i] * d[i];
    std::array<double, CMAX + 1> next{};
    double ll_next = -INFINITY;
    for (int half = 0; half < 40; ++half, step /= 2) {
      for (int i = 0; i < p; ++i) next[i] = par[i] + step * d[i];
      bool ordered = true;
      for (int j = 1; j + 1 < m; ++j) ordered = ordered && next[j] > next[j - 1];
      ll_next = ordered ? loglik(next) : -INFINITY;
      if (ll_next >= ll - 1e-12) break;
    }
    if (!(ll_next >= ll - 1e-12)) return f;
    par = next;
    const double improvement = ll_next - ll;
    ll = ll_next;
    // converged: one more pass computes the Hessian at the maximum for the standard error
    converged = gain < 1e-10 && improvement < 1e-10;
  }
  return f;
}

}  // namespace ordinal

namespace {

using namespace blocks;
using sim::Rng;
using ordinal::CMAX;
using Table2 = std::array<std::array<double, CMAX>, 2>;  // counts by group (0: group 1) and category

struct Design {
  int k = 0;
  std::array<double, CMAX> p1{}, p2{};  // category probabilities (group 1, group 2), lowest category first
  bool po = true;                       // the treatment probabilities come from a common odds ratio
  double odds_ratio = 1;
};

std::array<double, CMAX> probs_of(const std::vector<double>& v, const char* what, int k) {
  if (static_cast<int>(v.size()) != k)
    throw CalcError(std::string(what) + " needs one probability per category (" + std::to_string(k) + ")");
  double total = 0;
  for (double x : v) {
    if (!(x >= 0)) throw CalcError(std::string(what) + ": probabilities cannot be negative");
    total += x;
  }
  if (std::fabs(total - 1) > 0.005)
    throw CalcError(std::string(what) + " must add up to 1 (they add up to " + text::format_g(total, 4) + ")");
  std::array<double, CMAX> p{};
  for (int j = 0; j < k; ++j) p[j] = v[j] / total;
  return p;
}

Design design_of(const Scenario& s) {
  Design d;
  const std::vector<double>& p2 = s.v("p2");
  if (p2.size() < 2) throw CalcError("enter at least two category probabilities for group 2");
  if (p2.size() > CMAX) throw CalcError("at most 20 categories");
  d.k = static_cast<int>(p2.size());
  d.p2 = probs_of(p2, "the probabilities of group 2", d.k);
  d.po = !s.is("effect", "probs");
  if (d.po) {
    d.odds_ratio = s("or");
    double below = 0, cum2 = 0;
    for (int j = 0; j < d.k; ++j) {
      cum2 += d.p2[j];
      const double up = j + 1 == d.k ? 1 : (cum2 >= 1 ? 1 : expit(logit(cum2) - std::log(d.odds_ratio)));
      d.p1[j] = std::fmax(0.0, up - below);
      below = up;
    }
  } else {
    d.p1 = probs_of(s.v("p1"), "the probabilities of group 1", d.k);
  }
  return d;
}

/// P(Y1 > Y2) + P(Y1 = Y2)/2, the probability that group 1 does better (higher categories) counting ties as half.
double win_probability(const Design& d) {
  double below = 0, w = 0;
  for (int j = 0; j < d.k; ++j) {
    w += d.p1[j] * (below + d.p2[j] / 2);
    below += d.p2[j];
  }
  return w;
}

enum Test { WALD = 0, LR = 1, RANKS = 2 };

// Bits of sim::Outcome::also: each test's rejection (bit WALD, LR or RANKS), then why a test could not be computed
// (SEPARATED is also the share of separated studies).
constexpr unsigned SINGLE = 3, SEPARATED = 4, NO_FIT = 5;

const std::vector<simc::Cause>& causes_of(Test test) {
  using simc::Failures;
  static const simc::Cause single{SINGLE, Failures::Design, "every outcome fell in one category"};
  static const simc::Cause no_fit{NO_FIT, Failures::Numerical, "the proportional odds fit did not converge"};
  static const std::vector<simc::Cause> wald = {
      single,
      {SEPARATED, Failures::Design,
       "one group's outcomes all lay at or above the other's, where the estimate is infinite and the Wald test "
       "cannot reject, though the likelihood-ratio and Mann–Whitney tests can"},
      no_fit};
  static const std::vector<simc::Cause> lr = {single, no_fit}, ranks = {single};
  return test == WALD ? wald : (test == LR ? lr : ranks);
}

sim::Outcome ordinal_study(const Design& d, long n1, long n2, long rep, std::uint64_t seed, Test test, double alpha,
                           Alt a) {
  Table2 t{};
  std::array<std::array<double, CMAX>, 2> cum{};
  for (int g = 0; g < 2; ++g) {
    double c = 0;
    for (int j = 0; j < d.k; ++j) cum[g][j] = c += (g == 0 ? d.p1[j] : d.p2[j]);
  }
  const long n[2] = {n1, n2};
  for (int g = 0; g < 2; ++g) {
    Rng rng(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g + 1));
    for (long i = 0; i < n[g]; ++i) {
      const double u = rng.uniform();
      int j = 0;
      while (j + 1 < d.k && u > cum[g][j]) ++j;
      t[g][j] += 1;
    }
  }
  sim::Outcome o;
  const ordinal::PoFit f = ordinal::po_fit(t[0].data(), t[1].data(), d.k);
  const unsigned fit_failed = f.ok ? 0u : (f.single ? SINGLE : NO_FIT);
  bool wald = false, lr = false;
  if (f.ok) {
    if (!f.separated) {
      wald = simc::p_value(f.beta / f.se, INFINITY, a) < alpha;
      o.estimate = f.beta;
      o.se = f.se;
    }
    const double stat = std::fmax(0.0, 2 * (f.loglik - f.loglik0));
    const double z = (f.beta > 0 ? 1 : -1) * std::sqrt(stat);
    lr = a == Alt::Two ? stats::pchisq(stat, 1, false) < alpha : simc::p_value(z, INFINITY, a) < alpha;
  }
  // Mann–Whitney on the table: midranks of each category
  const double a1 = static_cast<double>(n1), a2 = static_cast<double>(n2);
  double before = 0, r1 = 0, ties = 0;
  for (int j = 0; j < d.k; ++j) {
    const double c = t[0][j] + t[1][j];
    r1 += t[0][j] * (before + (c + 1) / 2);
    ties += c * c * c - c;
    before += c;
  }
  const bool ranks = simc::rank_sum_p(r1 - a1 * (a1 + 1) / 2, a1, a2, ties, a) < alpha;
  // every outcome in one category leaves the ranks no variance, as it leaves the model nothing to fit; separation
  // leaves the Wald test no finite estimate
  unsigned failed = fit_failed;
  if (test == RANKS) failed = f.single ? SINGLE : 0u;
  if (test == WALD && f.separated) failed = SEPARATED;
  o.also = (wald ? 1u : 0u) | (lr ? 2u : 0u) | (ranks ? 4u : 0u) | (fit_failed ? 1u << fit_failed : 0u) |
           (f.separated ? 1u << SEPARATED : 0u);
  if (failed != 0) {
    o.failed = true;
  } else {
    o.reject = test == WALD ? wald : (test == LR ? lr : ranks);
  }
  return o;
}

Procedure ordinal() {
  Procedure p;
  p.id = "sim.ordinal";
  p.name = "Two Groups with an Ordinal Outcome (Proportional Odds Model or Mann–Whitney), by Simulation";
  p.category = "Simulation";
  p.subcategory = "Ordinal Outcomes";
  p.order = 1;
  p.simulation = true;
  p.summary =
      "Sample size (or power) for comparing two groups on ordered categories (the modified Rankin Scale, a pain "
      "grade, a Likert item), estimated by simulating studies and analyzing each with the proportional odds "
      "model (Wald or likelihood-ratio test) or the Mann–Whitney test. The treatment effect is a common odds "
      "ratio, or the treatment group's own category probabilities.";
  p.keywords = {"simulation", "monte carlo", "ordinal", "ordered categories", "proportional odds",
                "cumulative logit", "ordinal logistic regression", "polr", "common odds ratio", "shift analysis",
                "modified rankin scale", "mrs", "likert", "grades", "pain grade", "severity score",
                "mann-whitney", "wilcoxon", "whitehead", "categories"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total Sample Size)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: OR ≠ 1)", "One-sided (H1: group 1 higher, OR > 1)",
                  "One-sided (H1: group 1 lower, OR < 1)"),
      select("test", "Test",
             {{"wald", "Proportional odds model, Wald test", {}},
              {"lr", "Proportional odds model, likelihood-ratio test", {}},
              {"mw", "Mann–Whitney (Wilcoxon rank-sum) test with ties", {}}},
             "wald")
          .grp(G_TEST)
          .hint("The test that decides the power and the sample size; the other two are run on the same simulated "
                "studies and reported alongside."),
      power_param(),
      alpha_param(),
      num("n1", "N1 (Subjects in Group 1)", "100")
          .grp(G_N)
          .sym("N1")
          .whole()
          .at_least(2)
          .sugg({"100", "50 to 200 by 50"})
          .hint("Number of subjects in group 1; group 2 gets N1 × the allocation ratio, rounded up.")
          .show(when(not_solving({"n"}))),
      num("r", "Allocation Ratio N2/N1", "1")
          .grp(G_N)
          .sym("N2/N1")
          .above(0)
          .sugg({"1", "2"})
          .hint("Size of group 2 relative to group 1 (N2 = ratio × N1, rounded up); 1 gives equal groups."),
      vec("p2", "Category Probabilities of Group 2 (Lowest First)", "0.2 0.25 0.25 0.2 0.1")
          .grp(G_EFFECT)
          .sym("P2")
          .range(0, 1)
          .hint("The share of group 2 (e.g. the control) in each category, from the lowest to the highest, adding "
                "up to 1; from 2 to 20 categories."),
      select("effect", "Treatment Effect",
             {{"or", "A common odds ratio (proportional odds)", {}},
              {"probs", "The category probabilities of group 1", {}}},
             "or")
          .grp(G_EFFECT)
          .hint("A common odds ratio shifts every cumulative split of the scale by the same amount; giving group "
                "1's probabilities allows any shape, and shows how the tests fare when the odds are not "
                "proportional."),
      num("or", "OR (Odds of a Higher Category, Group 1 vs 2)", "1.8")
          .grp(G_EFFECT)
          .sym("OR")
          .above(0)
          .sugg({"1.8", "1.5 2"})
          .hint("The odds of being above any split of the scale in group 1, divided by those in group 2; above 1 "
                "group 1 tends to higher categories. When the higher categories are worse, an effective "
                "treatment has an odds ratio below 1.")
          .show(when(is("effect", {"or"}))),
      vec("p1", "Category Probabilities of Group 1 (Lowest First)", "0.1 0.2 0.3 0.25 0.15")
          .grp(G_EFFECT)
          .sym("P1")
          .range(0, 1)
          .hint("The share of group 1 (e.g. the treatment) in each category, from the lowest to the highest, adding "
                "up to 1.")
          .show(when(is("effect", {"probs"}))),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double alpha = s("alpha"), ratio = s("r");
    const Design d = design_of(s);
    const Test test = s.is("test", "lr") ? LR : (s.is("test", "mw") ? RANKS : WALD);
    const double win = win_probability(d);
    auto n2_of = [&](double n1) { return std::ceil(ratio * n1 - 1e-9); };
    // Whitehead's variance of log OR̂ for proportional odds; Noether's with ties for given probabilities
    auto approx = [&](double n1) {
      const double n2 = n2_of(n1), total = n1 + n2, pi = n1 / total;
      double cube = 0;
      for (int j = 0; j < d.k; ++j) {
        const double avg = pi * d.p1[j] + (1 - pi) * d.p2[j];
        cube += avg * avg * avg;
      }
      if (d.po) {
        const double se = std::sqrt(3 / (total * pi * (1 - pi) * (1 - cube)));
        return z_power(std::log(d.odds_ratio) / se, alpha, a);
      }
      const double se = std::sqrt((total + 1) * (1 - cube) / (12 * n1 * n2));
      return z_power((win - 0.5) / se, alpha, a);
    };
    int lo1 = d.k, hi1 = -1, lo2 = d.k, hi2 = -1;
    for (int j = 0; j < d.k; ++j) {
      if (d.p1[j] > 0) lo1 = std::min(lo1, j), hi1 = j;
      if (d.p2[j] > 0) lo2 = std::min(lo2, j), hi2 = j;
    }
    if (lo1 == hi1 && lo2 == hi2 && lo1 == lo2)
      throw CalcError("every subject of both groups falls in the same category, so there is nothing to compare");
    if (test == WALD) {  // one group certain to lie at or above the other: the Wald test can never reject
      if (hi2 <= lo1 || hi1 <= lo2)
        throw CalcError("one group lies entirely at or above the other, so the proportional odds model's estimate is "
                        "infinite and its Wald test cannot reject; use the likelihood-ratio or Mann–Whitney test");
    }
    if (s.solve == "n") {
      const double effect = win - 0.5;
      if (std::fabs(effect) < 1e-12) throw CalcError("the two groups have the same distribution over the categories");
      if ((a == Alt::Upper && effect < 0) || (a == Alt::Lower && effect > 0))
        throw CalcError("the effect has the direction opposite to the one-sided alternative");
    }
    const double lo = s.solve != "n" ? 2 : solve::first_whole([&](double x) { return n2_of(x) >= 2; }, 2, 1e5,
                                                           "the allocation ratio is too small: group 2 would have "
                                                           "fewer than 2 subjects");
    double n1;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = lo;
      o.hi = 1e5;
      n1 = solve::min_n(approx, s("power"), o, &ctx);
    } else {
      n1 = s("n1");
      if (n2_of(n1) < 2) throw CalcError("group 2 needs at least 2 subjects");
    }
    if (ctx.approximate) {
      simc::put_approximation(s, r, approx(n1));
    } else {
      const long reps = static_cast<long>(s("reps"));
      const auto seed = static_cast<std::uint64_t>(s("seed"));
      auto power_at = [&](double m1, long R, double p0, double p1) {
        simc::check_study_size(m1 + n2_of(m1), "subjects");
        const long k1 = static_cast<long>(m1), k2 = static_cast<long>(n2_of(m1));
        return simc::simulate_tests(
            R, [&](long rep) { return ordinal_study(d, k1, k2, rep, seed, test, alpha, a); }, ctx, p0, p1);
      };
      sim::Summary sm;
      if (s.solve == "n") {
        sim::SizeSearch o;
        o.lo = lo;
        o.hi = 1e5;
        o.guess = n1;
        o.reps = reps;
        o.explore_reps = simc::explore_reps(reps);
        try {
          const sim::SizeResult res = sim::search_size(power_at, s("power"), o, ctx);
          n1 = res.n;
          sm = res.s;
        } catch (const CalcError& e) {
          // a large effect that the Wald test still misses: separation (one group at or above the other)
          if (test == WALD && approx(o.hi) >= s("power"))
            throw CalcError(std::string(e.what()) + "; in most simulated studies one group lies entirely at or "
                            "above the other, where the Wald test cannot reject: use the likelihood-ratio or "
                            "Mann–Whitney test");
          throw;
        }
      } else {
        sm = power_at(n1, reps, 0, 1);
      }
      simc::put_summary(s, r, sm, causes_of(test));
      r.set("approx_power", approx(n1));
      r.set("power_wald", simc::other_power(sm, WALD, causes_of(WALD)));
      r.set("power_lr", simc::other_power(sm, LR, causes_of(LR)));
      r.set("power_mw", simc::other_power(sm, RANKS, causes_of(RANKS)));
      if (std::isfinite(sm.mean_estimate)) r.set("mean_logor", sm.mean_estimate);
      if (sm.also[SEPARATED] > 0) r.set("separated", sm.also_power(SEPARATED));
    }
    put_two_group(r, {n1, n2_of(n1)});
    std::vector<double> p1(d.p1.begin(), d.p1.begin() + d.k), p2(d.p2.begin(), d.p2.begin() + d.k);
    for (double& x : p1) x = std::round(x * 1e4) / 1e4;
    r.set("p1", list_text(p1));
    r.set("p2", list_text(p2));
    r.set("k", static_cast<double>(d.k));
    if (d.po) r.set("or", d.odds_ratio);
    r.set("win", win);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("k", "Categories", "int", "Number of ordered categories"),
                   col("or", "OR", "num", "Common odds ratio of a higher category, group 1 vs 2"),
                   col("mean_logor", "Mean log OR̂", "num", "Average estimated log odds ratio (proportional odds model)"),
                   col("win", "P(Better)", "prob", "P(Y1 > Y2) + P(Y1 = Y2)/2: the chance that a subject of group 1 "
                                                    "is in a higher category, counting ties as half"),
                   col("p1", "P1", "text", "Category probabilities of group 1"),
                   col("p2", "P2", "text", "Category probabilities of group 2"),
                   col("power_wald", "Power (Wald)", "prob", "Simulated power of the proportional odds Wald test"),
                   col("power_lr", "Power (LR)", "prob", "Simulated power of the likelihood-ratio test"),
                   col("power_mw", "Power (M–W)", "prob", "Simulated power of the Mann–Whitney test"),
                   col("separated", "Separated", "prob", "Share of studies with complete separation (one group "
                                                         "entirely above the other)"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns()) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string test = s.is("test", "lr")   ? "likelihood-ratio test of the proportional odds model"
                             : s.is("test", "mw") ? "Mann–Whitney test with ties"
                                                  : "Wald test of the proportional odds model";
    std::string effect = r.has("or") && std::isfinite(r.num("or"))
                             ? "a common odds ratio of " + fnum(r.num("or")) + " (category probabilities " +
                                   text_of(r, "p2") + " in group 2 and " + text_of(r, "p1") + " in group 1)"
                             : "category probabilities of " + text_of(r, "p1") + " in group 1 against " +
                                   text_of(r, "p2") + " in group 2";
    std::string t = "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " (total " +
                    fint(r.num("n")) + ") " + simc::power_phrase(s, r) + " to detect " + effect +
                    " over " + fint(r.num("k")) + " ordered categories, using a " + alt_words(alt(s)) + " " + test +
                    " at α = " + fnum(r.num("alpha")) + ".";
    t += simc::others_sentence(r, "studies", {{"Wald test", "power_wald"}, {"likelihood-ratio test", "power_lr"},
                                             {"Mann–Whitney test", "power_mw"}});
    return t + simc::method_note(r);
  };
  p.definitions = {{"N1, N2", "Number of subjects in groups 1 and 2; N = N1 + N2."},
                   {"OR", "Common odds ratio: the odds of a higher category in group 1 divided by those in group 2, "
                          "the same at every split of the scale (proportional odds)."},
                   {"P(Better)", "P(Y1 > Y2) + P(Y1 = Y2)/2, the probability that a subject of group 1 is in a higher "
                                 "category than one of group 2, counting ties as half; 0.5 with no difference."},
                   {"Power", "Share of the simulated studies whose test was significant, with its 95% interval "
                             "(Wilson) and Monte Carlo standard error."}};
  p.notes = {"The proportional odds model logit P(Y ≤ j) = θj − β·group is fitted by maximum likelihood (Newton's "
             "method; as MASS::polr), with the Wald z-test of β or the likelihood-ratio test; categories nobody "
             "falls in are merged away.",
             "With separation (every subject of one group at or above every subject of the other) the estimate "
             "is infinite: the Wald test then cannot reject, while the likelihood-ratio and Mann–Whitney tests "
             "can. Such a study counts as a failed analysis of the Wald test (not rejecting), as a study with every "
             "outcome in one category does for every test.",
             "The quick approximation is Whitehead's (1993) formula for a common odds ratio, and Noether's formula "
             "with the tie correction when group 1's probabilities are given.",
             "Giving group 1's probabilities directly shows the power when the odds are not proportional, for "
             "example when a treatment helps only the most severe patients."};
  p.references = {"Whitehead, J. (1993). Sample size calculations for ordered categorical data. Statistics in "
                  "Medicine 12, 2257–2271.",
                  "McCullagh, P. (1980). Regression models for ordinal data. Journal of the Royal Statistical "
                  "Society B 42, 109–142.",
                  "Roozenbeek, B. et al. (2011). The added value of ordinal analysis in clinical trials: an example "
                  "in traumatic brain injury. Critical Care 15, R127.",
                  simc::REF_SIM};
  p.curve_param = "n1";
  p.curve_label = "N1 (Subjects in Group 1)";
  p.example_solve = "n";
  p.example = {{"power", "0.8"}, {"p2", "0.2 0.25 0.25 0.2 0.1"}, {"or", "1.8"}, {"reps", "1000"}};
  return p;
}
SSC_REGISTER(ordinal);

}  // namespace
}  // namespace ssc::procs
