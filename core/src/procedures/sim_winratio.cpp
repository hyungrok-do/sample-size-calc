// Simulation-based power for a hierarchical composite endpoint analyzed by
// pairwise comparisons (sim.winratio): every patient of the treatment group
// is compared with every patient of the control group, first on death, then
// (if neither died during the time both were followed) on a non-fatal event
// such as a hospitalization, and last on a measured outcome; each comparison
// is a win, a loss or a tie (Pocock et al. 2012).
//
// Deaths and non-fatal events follow exponential distributions with hazard
// ratios between the groups; patients enter uniformly during accrual and are
// followed until the analysis or loss to follow-up.  Each study is analyzed
// with the Finkelstein–Schoenfeld test (the net benefit, wins minus losses,
// with its permutation variance from all pairs) or the win ratio's Wald test
// (log of wins over losses with the two-sample U-statistic variance; Bebu &
// Lachin 2016).  The quick approximation is Yu & Ganju's (2022) formula with
// the probabilities of a win, a loss and a tie computed from 200,000 pairs.
#include "procedures/sim_winratio.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <vector>

#include "procedures/sim_common.hpp"
#include "util/memo.hpp"

namespace ssc::procs {

namespace winratio {
namespace {

/// Fenwick tree counting entries by rank.
class Fenwick {
 public:
  explicit Fenwick(std::size_t n) : t_(n + 1, 0) {}
  void add(std::size_t rank) {
    for (std::size_t i = rank + 1; i < t_.size(); i += i & (~i + 1)) t_[i] += 1;
  }
  /// Entries with a rank below `rank`.
  long below(std::size_t rank) const {
    long s = 0;
    for (std::size_t i = rank; i > 0; i -= i & (~i + 1)) s += t_[i];
    return s;
  }

 private:
  std::vector<long> t_;
};

struct Point {
  double x, y;
};
struct Query {
  double x, y;
  std::size_t k;  // the patient asking
};
enum class X { AtLeast, AtMost };  // a point counts when x ≥ X (or x ≤ X) …
enum class Y { Below, Above };     // … and y < Y (or y > Y)

/// Adds to out[q.k] the number of points with x (≥ or ≤) q.x and y (< or >) q.y, by sweeping points and
/// queries in order of x through a Fenwick tree over the points' y.
void dominance(std::vector<Point>& pts, std::vector<Query>& qs, X xr, Y yr, std::vector<double>& out) {
  if (pts.empty() || qs.empty()) return;
  std::vector<double> ys(pts.size());
  for (std::size_t i = 0; i < pts.size(); ++i) ys[i] = pts[i].y;
  std::sort(ys.begin(), ys.end());
  ys.erase(std::unique(ys.begin(), ys.end()), ys.end());
  const bool down = xr == X::AtLeast;
  std::sort(pts.begin(), pts.end(), [&](const Point& a, const Point& b) { return down ? a.x > b.x : a.x < b.x; });
  std::sort(qs.begin(), qs.end(), [&](const Query& a, const Query& b) { return down ? a.x > b.x : a.x < b.x; });
  Fenwick f(ys.size());
  std::size_t next = 0;
  long inserted = 0;
  for (const Query& q : qs) {
    while (next < pts.size() && (down ? pts[next].x >= q.x : pts[next].x <= q.x)) {
      f.add(static_cast<std::size_t>(std::lower_bound(ys.begin(), ys.end(), pts[next].y) - ys.begin()));
      ++next;
      ++inserted;
    }
    if (yr == Y::Below)
      out[q.k] += static_cast<double>(
          f.below(static_cast<std::size_t>(std::lower_bound(ys.begin(), ys.end(), q.y) - ys.begin())));
    else
      out[q.k] += static_cast<double>(
          inserted - f.below(static_cast<std::size_t>(std::upper_bound(ys.begin(), ys.end(), q.y) - ys.begin())));
  }
}

}  // namespace

void pair_counts(const std::vector<Patient>& p, long n1, double margin, Counts& out) {
  const std::size_t n = p.size(), m1 = static_cast<std::size_t>(n1);
  for (int g = 0; g < 2; ++g) {
    out.win[g].assign(n, 0);
    out.loss[g].assign(n, 0);
  }
  std::vector<double> sorted;
  std::vector<Point> pts;
  std::vector<Query> qs;
  std::vector<std::pair<double, std::size_t>> runs;
  // dead patients at exactly time t (or with the event at exactly t): the few pairs a sweep cannot separate
  auto equal_run = [&](double t) {
    auto lo = std::lower_bound(runs.begin(), runs.end(), std::make_pair(t, std::size_t{0}));
    auto hi = lo;
    while (hi != runs.end() && hi->first == t) ++hi;
    return std::make_pair(lo, hi);
  };
  for (int g = 0; g < 2; ++g) {
    const std::size_t from = g == 0 ? 0 : m1, to = g == 0 ? m1 : n;
    std::vector<double>& win = out.win[g];
    std::vector<double>& loss = out.loss[g];
    // --- death: k beats l when l died while k was known to be alive, l beats k when k died while l was
    sorted.clear();
    for (std::size_t l = from; l < to; ++l) sorted.push_back(p[l].death);
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t k = 0; k < n; ++k)
      win[k] += static_cast<double>(std::lower_bound(sorted.begin(), sorted.end(), p[k].alive) - sorted.begin());
    sorted.clear();
    for (std::size_t l = from; l < to; ++l) sorted.push_back(p[l].alive);
    std::sort(sorted.begin(), sorted.end());
    for (std::size_t k = 0; k < n; ++k)
      loss[k] += static_cast<double>(sorted.end() - std::upper_bound(sorted.begin(), sorted.end(), p[k].death));
    // --- the non-fatal event, for pairs that death leaves tied: with k alive, every l with l.death ≥ k.alive;
    // with k dead at t, every l alive with l.alive ≤ t, and every l dead at exactly t
    pts.clear();
    qs.clear();
    for (std::size_t l = from; l < to; ++l) pts.push_back({p[l].death, p[l].event});
    for (std::size_t k = 0; k < n; ++k)
      if (std::isinf(p[k].death)) qs.push_back({p[k].alive, p[k].free, k});
    dominance(pts, qs, X::AtLeast, Y::Below, win);  // l.event < k.free
    pts.clear();
    for (std::size_t l = from; l < to; ++l) pts.push_back({p[l].death, p[l].free});
    for (auto& q : qs) q.y = p[q.k].event;
    dominance(pts, qs, X::AtLeast, Y::Above, loss);  // l.free > k.event
    pts.clear();
    qs.clear();
    for (std::size_t l = from; l < to; ++l)
      if (std::isinf(p[l].death)) pts.push_back({p[l].alive, p[l].event});
    for (std::size_t k = 0; k < n; ++k)
      if (!std::isinf(p[k].death)) qs.push_back({p[k].death, p[k].free, k});
    dominance(pts, qs, X::AtMost, Y::Below, win);
    pts.clear();
    for (std::size_t l = from; l < to; ++l)
      if (std::isinf(p[l].death)) pts.push_back({p[l].alive, p[l].free});
    for (auto& q : qs) q.y = p[q.k].event;
    dominance(pts, qs, X::AtMost, Y::Above, loss);
    runs.clear();
    for (std::size_t l = from; l < to; ++l)
      if (!std::isinf(p[l].death)) runs.push_back({p[l].death, l});
    std::sort(runs.begin(), runs.end());
    for (std::size_t k = 0; k < n; ++k) {
      if (std::isinf(p[k].death)) continue;
      const auto r = equal_run(p[k].death);
      for (auto it = r.first; it != r.second; ++it) {
        win[k] += p[it->second].event < p[k].free;
        loss[k] += p[it->second].free > p[k].event;
      }
    }
    // --- the measured outcome, between patients alive at the end of follow-up (the others have no value) that
    // death and the event leave tied: with k free of the event, every such l with l.event ≥ k.free; with k's
    // event at t, every such l free of the event with l.free ≤ t, and every l with its event at exactly t
    pts.clear();
    qs.clear();
    for (std::size_t l = from; l < to; ++l)
      if (!std::isnan(p[l].value)) pts.push_back({p[l].event, p[l].value});
    if (pts.empty()) continue;
    for (std::size_t k = 0; k < n; ++k)
      if (!std::isnan(p[k].value) && std::isinf(p[k].event)) qs.push_back({p[k].free, p[k].value - margin, k});
    dominance(pts, qs, X::AtLeast, Y::Below, win);  // l.value < k.value − margin
    for (auto& q : qs) q.y = p[q.k].value + margin;
    std::vector<Point> again = pts;
    dominance(again, qs, X::AtLeast, Y::Above, loss);  // l.value > k.value + margin
    pts.clear();
    qs.clear();
    for (std::size_t l = from; l < to; ++l)
      if (!std::isnan(p[l].value) && std::isinf(p[l].event)) pts.push_back({p[l].free, p[l].value});
    for (std::size_t k = 0; k < n; ++k)
      if (!std::isnan(p[k].value) && !std::isinf(p[k].event)) qs.push_back({p[k].event, p[k].value - margin, k});
    again = pts;
    dominance(again, qs, X::AtMost, Y::Below, win);
    for (auto& q : qs) q.y = p[q.k].value + margin;
    dominance(pts, qs, X::AtMost, Y::Above, loss);
    runs.clear();
    for (std::size_t l = from; l < to; ++l)
      if (!std::isnan(p[l].value) && !std::isinf(p[l].event)) runs.push_back({p[l].event, l});
    std::sort(runs.begin(), runs.end());
    for (std::size_t k = 0; k < n; ++k) {
      if (std::isnan(p[k].value) || std::isinf(p[k].event)) continue;
      const auto r = equal_run(p[k].event);
      for (auto it = r.first; it != r.second; ++it) {
        win[k] += p[it->second].value < p[k].value - margin;
        loss[k] += p[it->second].value > p[k].value + margin;
      }
    }
  }
}

void pair_counts_direct(const std::vector<Patient>& p, long n1, double margin, Counts& out) {
  const std::size_t n = p.size(), m1 = static_cast<std::size_t>(n1);
  for (int g = 0; g < 2; ++g) {
    out.win[g].assign(n, 0);
    out.loss[g].assign(n, 0);
  }
  for (std::size_t k = 0; k < n; ++k)
    for (std::size_t l = 0; l < n; ++l) {
      const int c = compare(p[k], p[l], margin);
      const int g = l < m1 ? 0 : 1;
      out.win[g][k] += c > 0;
      out.loss[g][k] += c < 0;
    }
}

}  // namespace winratio

namespace {

using namespace blocks;
using sim::Rng;

constexpr double LN2 = 0.69314718055994530942;

struct Design {
  double lam_d = 0, hr_d = 1;   // control hazard of death, hazard ratio (treatment / control)
  bool event = true;            // a non-fatal event is the second priority
  double lam_h = 0, hr_h = 1;
  bool score = false;           // a measured outcome is the last priority
  double delta = 0, margin = 0; // its difference in SD units (treatment − control) and the smallest that counts
  double accrual = 24, follow = 12, eta = 0;
};

using winratio::Patient;

Patient draw(const Design& d, bool treated, Rng& rng) {
  const double entry = rng.uniform() * d.accrual;
  const double death = rng.exponential() / (d.lam_d * (treated ? d.hr_d : 1));
  const double event = d.event ? rng.exponential() / (d.lam_h * (treated ? d.hr_h : 1)) : INFINITY;
  const double lost = d.eta > 0 ? rng.exponential() / d.eta : INFINITY;
  const double value = d.score ? rng.normal() + (treated ? d.delta : 0) : 0;
  const double end = std::fmin(d.accrual + d.follow - entry, lost);
  Patient p;
  p.death = death <= end ? death : INFINITY;
  p.alive = std::fmin(death, end);
  p.event = event <= end && event < death ? event : INFINITY;
  p.free = d.event ? std::fmin(event, p.alive) : INFINITY;
  p.value = d.score && death > end ? value : NAN;
  return p;
}

struct Probabilities {
  double win = 0, loss = 0, tie = 0, deaths = 0;
};

/// Probabilities of a win, a loss and a tie for a random treatment–control pair (200,000 pairs from a fixed
/// stream, so the approximation does not depend on the seed), and the share of patients who died.
Probabilities pair_probabilities(const Design& d) {
  constexpr long PAIRS = 200000;
  // half the pairs draw the treatment patient from the first stream, half from the second, so that the Monte
  // Carlo error cancels between the groups: identical groups give exactly as many wins as losses
  long w = 0, l = 0, deaths = 0;
  for (int half = 0; half < 2; ++half) {
    Rng rt(0x5eed, half == 0 ? 1 : 2), rc(0x5eed, half == 0 ? 2 : 1);
    for (long i = 0; i < PAIRS / 2; ++i) {
      const Patient a = draw(d, true, rt), b = draw(d, false, rc);
      const int c = winratio::compare(a, b, d.margin);
      w += c > 0;
      l += c < 0;
      deaths += std::isfinite(a.death) + std::isfinite(b.death);
    }
  }
  Probabilities p;
  p.win = static_cast<double>(w) / PAIRS;
  p.loss = static_cast<double>(l) / PAIRS;
  p.tie = 1 - p.win - p.loss;
  p.deaths = static_cast<double>(deaths) / (2.0 * PAIRS);
  return p;
}

/// pair_probabilities(d), drawn once per design: the power curve and the size search ask for the same design
/// at every sample size, and the 200,000 pairs take tens of milliseconds.
Probabilities cached_pair_probabilities(const Design& d) {
  using Key = std::array<double, 11>;  // every field of the design
  static OnceCache<Key, Probabilities> cache;
  const Key key{d.lam_d, d.hr_d, d.event ? 1.0 : 0.0, d.lam_h, d.hr_h, d.score ? 1.0 : 0.0,
                d.delta, d.margin, d.accrual, d.follow, d.eta};
  return *cache.get(key, [&] { return pair_probabilities(d); });
}

enum Test { FS = 0, WR = 1 };

// Bits of sim::Outcome::also: each test's rejection (bit FS or WR), then why a test could not be computed.
constexpr unsigned TIED = 2, BALANCED = 3, ONE_SIDED = 4, FLAT = 5;

const std::vector<simc::Cause>& causes_of(Test test) {
  using simc::Failures;
  static const std::vector<simc::Cause> fs = {
      {TIED, Failures::Design, "every pair of patients tied"},
      {BALANCED, Failures::Design, "every patient won as often as lost, so the permutation variance was 0"}};
  static const std::vector<simc::Cause> wr = {
      {ONE_SIDED, Failures::Design, "each had no wins or no losses, so no finite win ratio"},
      {FLAT, Failures::Design, "the win ratio's standard error was 0"}};
  return test == FS ? fs : wr;
}

sim::Outcome winratio_study(const Design& d, long n1, long n2, long rep, std::uint64_t seed, Test test, double alpha,
                            Alt a) {
  // per study, freed with it (nothing outlives the simulation)
  std::vector<Patient> pts;
  winratio::Counts counts;
  std::vector<double> u, w1, l1, w2, l2;
  const long n = n1 + n2;
  pts.resize(static_cast<std::size_t>(n));
  for (int g = 0; g < 2; ++g) {
    Rng rng(seed, static_cast<std::uint64_t>(rep), static_cast<std::uint64_t>(g + 1));
    const long start = g == 0 ? 0 : n1, count = g == 0 ? n1 : n2;
    for (long i = 0; i < count; ++i) pts[static_cast<std::size_t>(start + i)] = draw(d, g == 0, rng);
  }
  winratio::pair_counts(pts, n1, d.margin, counts);
  // U_k over all patients (for the Finkelstein–Schoenfeld variance); each patient's wins and losses against
  // the other group (for the win ratio)
  u.resize(static_cast<std::size_t>(n));
  w1.resize(static_cast<std::size_t>(n1));
  l1.resize(static_cast<std::size_t>(n1));
  w2.resize(static_cast<std::size_t>(n2));
  l2.resize(static_cast<std::size_t>(n2));
  double wins = 0, losses = 0, decided = 0;
  for (long k = 0; k < n; ++k) {
    const auto i = static_cast<std::size_t>(k);
    u[i] = counts.win[0][i] + counts.win[1][i] - counts.loss[0][i] - counts.loss[1][i];
    decided += counts.win[0][i] + counts.win[1][i];
    if (k < n1) {
      w1[i] = counts.win[1][i];
      l1[i] = counts.loss[1][i];
      wins += w1[i];
      losses += l1[i];
    } else {
      w2[static_cast<std::size_t>(k - n1)] = counts.win[0][i];
      l2[static_cast<std::size_t>(k - n1)] = counts.loss[0][i];
    }
  }
  const double a1 = static_cast<double>(n1), a2 = static_cast<double>(n2), nn = a1 + a2;
  sim::Outcome o;
  // Finkelstein–Schoenfeld: T = wins − losses, permutation variance n1 n2 / (N (N − 1)) Σ U_k²
  double uu = 0;
  for (double x : u) uu += x * x;
  const double var_t = a1 * a2 / (nn * (nn - 1)) * uu;
  const unsigned fs_failed = var_t > 0 ? 0u : (decided == 0 ? TIED : BALANCED);
  const bool fs = var_t > 0 && simc::p_value((wins - losses) / std::sqrt(var_t), INFINITY, a) < alpha;
  // win ratio: log(θW/θL) with the U-statistic covariance of (θW, θL)
  bool wr = false;
  double logwr = NAN, se = NAN;
  unsigned wr_failed = ONE_SIDED;
  if (wins > 0 && losses > 0) {
    const double tw = wins / (a1 * a2), tl = losses / (a1 * a2);
    // The delta-method variance vw/θW² + vl/θL² − 2cwl/(θW θL), from the U-statistic covariance of (θW, θL), is
    // within each group the variance of a patient's share of the other group beaten over θW less the share lost to
    // over θL, divided by the group's size.  Times wins · losses / (n1 n2) those are the whole numbers
    // w · losses − l · wins (exact below 2⁵³, so up to the 200,000 patients allowed), whose mean is 0: their sum of
    // squares is 0 exactly where the variance is, where the difference of the three rounded terms could fall to or
    // below 0, or stay just above it and reject.
    double ss1 = 0, ss2 = 0;
    for (std::size_t i = 0; i < w1.size(); ++i) {
      const double e = w1[i] * losses - l1[i] * wins;
      ss1 += e * e;
    }
    for (std::size_t i = 0; i < w2.size(); ++i) {
      const double e = l2[i] * losses - w2[i] * wins;
      ss2 += e * e;
    }
    const double v = (a1 / (a1 - 1) * ss1 + a2 / (a2 - 1) * ss2) / (wins * losses) / (wins * losses);
    wr_failed = FLAT;
    if (v > 0) {
      wr_failed = 0;
      logwr = std::log(tw / tl);
      se = std::sqrt(v);
      wr = simc::p_value(logwr / se, INFINITY, a) < alpha;
    }
  }
  o.also = (fs ? 1u : 0u) | (wr ? 2u : 0u) | (fs_failed ? 1u << fs_failed : 0u) |
           (wr_failed ? 1u << wr_failed : 0u);
  if ((test == FS ? fs_failed : wr_failed) != 0) {
    o.failed = true;  // no variance, or no wins or no losses (the win ratio is 0 or infinite)
  } else {
    o.reject = test == FS ? fs : wr;
  }
  o.estimate = logwr;
  o.se = se;
  o.extra = (wins - losses) / (a1 * a2);
  return o;
}

Design design_of(const Scenario& s) {
  Design d;
  d.lam_d = LN2 / s("med_d");
  // the hazard of a median of about 5e-324 is infinite, and every simulated time no number
  if (!std::isfinite(d.lam_d)) throw CalcError("the median time to death is too short to compute with");
  d.hr_d = s("hr_d");
  d.event = !s.is("tier2", "none");
  if (d.event) {
    d.lam_h = LN2 / s("med_h");
    if (!std::isfinite(d.lam_h)) throw CalcError("the median time to the non-fatal event is too short to compute with");
    d.hr_h = s("hr_h");
  }
  d.score = s.is("tier3", "score");
  if (d.score) {
    d.delta = s("delta_s");
    d.margin = s("margin_s");
  }
  d.accrual = s("accrual");
  d.follow = s("follow");
  d.eta = s("eta");
  return d;
}

Procedure winratio() {
  Procedure p;
  p.id = "sim.winratio";
  p.name = "Win Ratio for a Hierarchical Composite Endpoint (Death, Then a Non-Fatal Event, Then a Score), by "
           "Simulation";
  p.category = "Simulation";
  p.subcategory = "Composite Endpoints";
  p.order = 1;
  p.simulation = true;
  p.summary =
      "Sample size (or power) for a hierarchical composite endpoint analyzed by comparing every treatment patient "
      "with every control patient, first on death, then on a non-fatal event such as a hospitalization, then on a "
      "measured outcome: the Finkelstein–Schoenfeld test or the win ratio, by simulating whole trials with "
      "staggered entry and censoring.";
  p.keywords = {"simulation", "monte carlo", "win ratio", "win odds", "net benefit", "hierarchical composite",
                "composite endpoint", "prioritized outcomes", "finkelstein-schoenfeld", "pocock", "buyse",
                "generalized pairwise comparisons", "death and hospitalization", "heart failure", "cardiology",
                "kccq", "time to event", "survival", "pairwise comparisons"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total Sample Size)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: win ratio ≠ 1)", "One-sided (H1: the treatment wins more, win ratio > 1)",
                  "One-sided (H1: the control wins more, win ratio < 1)"),
      select("test", "Test",
             {{"fs", "Finkelstein–Schoenfeld (net benefit, permutation variance)", {}},
              {"wr", "Win ratio (log win ratio, U-statistic variance)", {}}},
             "fs")
          .grp(G_TEST)
          .hint("The Finkelstein–Schoenfeld test compares wins with losses using their variance when the groups "
                "do not differ; the win ratio's Wald test uses the variance of the estimated win ratio. Both are "
                "run on each simulated trial and reported."),
      power_param(),
      alpha_param(),
      num("n1", "N1 (Patients in the Treatment Group)", "300")
          .grp(G_N)
          .sym("N1")
          .whole()
          .at_least(2)
          .sugg({"300", "100 to 500 by 100"})
          .hint("Number of patients in the treatment group; the control group gets N1 × the allocation ratio, "
                "rounded up.")
          .show(when(not_solving({"n"}))),
      num("r", "Allocation Ratio N2/N1", "1")
          .grp(G_N)
          .sym("N2/N1")
          .above(0)
          .sugg({"1", "2"})
          .hint("Size of the control group relative to the treatment group; 1 gives equal groups."),
      num("med_d", "Median Time to Death (Control)", "60")
          .grp("First Priority: Death")
          .sym("Median (Death)")
          .above(0)
          .sugg({"60", "36 60 120"})
          .hint("Median survival in the control group, in the time unit of accrual and follow-up (e.g. months)."),
      num("hr_d", "HR for Death (Treatment vs Control)", "0.8")
          .grp(G_EFFECT)
          .sym("HR (Death)")
          .above(0)
          .sugg({"0.8", "0.7 0.8 0.9 1"})
          .hint("Hazard ratio of death; below 1 the treatment saves lives."),
      select("tier2", "Second Priority",
             {{"event", "A non-fatal event (the time to the first one)", {}}, {"none", "None", {}}}, "event")
          .grp("Second Priority: A Non-Fatal Event")
          .hint("Pairs in which neither patient died during the time both were followed are compared on the "
                "time to a non-fatal event, such as a hospitalization for heart failure."),
      num("med_h", "Median Time to the Non-Fatal Event (Control)", "24")
          .grp("Second Priority: A Non-Fatal Event")
          .sym("Median (Event)")
          .above(0)
          .sugg({"24", "12 24"})
          .hint("Median time to the first non-fatal event in the control group (without death).")
          .show(when(is("tier2", {"event"}))),
      num("hr_h", "HR for the Non-Fatal Event", "0.75")
          .grp(G_EFFECT)
          .sym("HR (Event)")
          .above(0)
          .sugg({"0.75", "0.7 0.8"})
          .hint("Hazard ratio of the non-fatal event between the groups.")
          .show(when(is("tier2", {"event"}))),
      select("tier3", "Last Priority",
             {{"none", "None", {}}, {"score", "A measured outcome at the end of follow-up (higher is better)", {}}},
             "none")
          .grp("Last Priority: A Measured Outcome")
          .hint("Pairs still tied are compared on a measured outcome (such as a quality-of-life score) among "
                "patients alive at the end of their follow-up."),
      num("delta_s", "Difference in the Measured Outcome (in SDs)", "0.2")
          .grp(G_EFFECT)
          .sym("Δ (SDs)")
          .sugg({"0.2", "0.1 0.2 0.3"})
          .hint("Mean of the treatment group minus the control group, in standard deviations; positive favors "
                "the treatment.")
          .show(when(is("tier3", {"score"}))),
      num("margin_s", "Smallest Difference That Counts (in SDs)", "0")
          .grp("Last Priority: A Measured Outcome")
          .sym("Margin")
          .at_least(0)
          .sugg({"0", "0.2"})
          .hint("Differences between two patients smaller than this count as ties.")
          .show(when(is("tier3", {"score"}))),
      num("accrual", "Accrual Time (Uniform Entry)", "24")
          .grp("Accrual and Follow-Up")
          .sym("Accrual")
          .at_least(0)
          .sugg({"24", "12 24"})
          .hint("Length of the enrollment period; patients enter uniformly during it."),
      num("follow", "Additional Follow-Up after Accrual", "12")
          .grp("Accrual and Follow-Up")
          .sym("Follow-Up")
          .above(0)
          .sugg({"12", "12 24"})
          .hint("Time from the end of accrual to the analysis: the last patient is followed this long."),
      num("eta", "Loss to Follow-Up Hazard Rate", "0")
          .grp("Accrual and Follow-Up")
          .sym("η")
          .at_least(0)
          .sugg({"0", "0.005"})
          .hint("Rate at which patients are lost to follow-up (per time unit), in both groups."),
      simc::reps_param(),
      simc::seed_param()};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double alpha = s("alpha"), ratio = s("r");
    const Design d = design_of(s);
    const Test test = s.is("test", "wr") ? WR : FS;
    const Probabilities pr = cached_pair_probabilities(d);
    const double nb = pr.win - pr.loss;
    auto n2_of = [&](double n1) { return std::ceil(ratio * n1 - 1e-9); };
    // Yu & Ganju (2022): Var(net benefit) ≈ (1 − p_tie²)/(3 κ(1 − κ) N), Var(log WR) ≈ 4(1 + p_tie)/(3 κ(1 − κ) N (1 − p_tie))
    auto approx = [&](double n1) {
      const double n2 = n2_of(n1), total = n1 + n2, kappa = n1 / total;
      if (test == FS) return z_power(nb / std::sqrt((1 - pr.tie * pr.tie) / (3 * kappa * (1 - kappa) * total)), alpha, a);
      if (!(pr.win > 0 && pr.loss > 0)) return static_cast<double>(NAN);
      const double v = 4 * (1 + pr.tie) / (3 * kappa * (1 - kappa) * total * (1 - pr.tie));
      return z_power(std::log(pr.win / pr.loss) / std::sqrt(v), alpha, a);
    };
    if (s.solve == "n") {
      // the same distribution in both groups (every hazard ratio 1, no score difference): wins and losses are
      // equally likely, which the net benefit, estimated from random pairs, shows only approximately
      const bool same = d.hr_d == 1 && (!d.event || d.hr_h == 1) && (!d.score || d.delta == 0);
      if (same || std::fabs(nb) < 1e-12) throw CalcError("wins and losses are equally likely; power cannot exceed alpha");
      if ((a == Alt::Upper && nb < 0) || (a == Alt::Lower && nb > 0))
        throw CalcError("the effect has the direction opposite to the one-sided alternative");
    }
    // at most 50,000 patients per trial (far beyond any trial analyzed this way; each costs O(N log N))
    const double cap = std::fmax(2, std::floor(50000 / (1 + ratio)));
    const double lo = s.solve != "n" ? 2 : solve::first_whole([&](double x) { return n2_of(x) >= 2; }, 2, cap,
                                                           "the allocation ratio is too small: the control group "
                                                           "would have fewer than 2 patients");
    const double n1_max = std::fmax(lo, cap);
    double n1;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = lo;
      o.hi = n1_max;
      n1 = solve::min_n(approx, s("power"), o, &ctx);
    } else {
      n1 = s("n1");
      if (n2_of(n1) < 2) throw CalcError("the control group needs at least 2 patients");
    }
    if (ctx.approximate) {
      simc::put_approximation(s, r, approx(n1));
    } else {
      const long reps = static_cast<long>(s("reps"));
      const auto seed = static_cast<std::uint64_t>(s("seed"));
      auto power_at = [&](double m1, long R, double p0, double p1) {
        simc::check_study_size(m1 + n2_of(m1), "patients", 2e5);  // O(N log N), about 0.5 s there
        const long k1 = static_cast<long>(m1), k2 = static_cast<long>(n2_of(m1));
        return simc::simulate_tests(
            R, [&](long rep) { return winratio_study(d, k1, k2, rep, seed, test, alpha, a); }, ctx, p0, p1);
      };
      sim::Summary sm;
      if (s.solve == "n") {
        sim::SizeSearch o;
        o.lo = lo;
        o.hi = n1_max;
        o.guess = n1;
        o.reps = reps;
        o.explore_reps = simc::explore_reps(reps);
        const sim::SizeResult res = sim::search_size(power_at, s("power"), o, ctx);
        n1 = res.n;
        sm = res.s;
      } else {
        sm = power_at(n1, reps, 0, 1);
      }
      simc::put_summary(s, r, sm, causes_of(test));
      r.set("approx_power", approx(n1));
      simc::put_other(r, sm, "power_fs", FS, causes_of(FS), test == FS);
      simc::put_other(r, sm, "power_wr", WR, causes_of(WR), test == WR);
      if (std::isfinite(sm.mean_estimate)) r.set("mean_wr", std::exp(sm.mean_estimate));
      if (std::isfinite(sm.mean_extra)) r.set("mean_nb", sm.mean_extra);
    }
    put_two_group(r, {n1, n2_of(n1)});
    r.set("pwin", pr.win);
    r.set("ploss", pr.loss);
    r.set("ptie", pr.tie);
    r.set("wr", pr.loss > 0 ? pr.win / pr.loss : NAN);
    r.set("wo", (pr.win + pr.tie / 2) / (pr.loss + pr.tie / 2));
    r.set("nb", nb);
    r.set("deaths", pr.deaths);
    r.set("med_d", s("med_d"));
    r.set("hr_d", d.hr_d);
    if (d.event) {
      r.set("med_h", s("med_h"));
      r.set("hr_h", d.hr_h);
    }
    if (d.score) r.set("delta_s", d.delta);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = simc::power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("wr", "Win Ratio", "num", "Wins divided by losses over all treatment–control pairs (expected)"),
                   col("wo", "Win Odds", "num", "Win odds: ties counted as half a win and half a loss"),
                   col("nb", "Net Benefit", "num", "Share of pairs won minus share lost (expected)"),
                   col("pwin", "Wins", "prob", "Probability that a treatment patient wins against a control patient"),
                   col("ploss", "Losses", "prob", "Probability that a treatment patient loses"),
                   col("ptie", "Ties", "prob", "Probability of a tie on every priority"),
                   col("mean_wr", "Mean WR̂", "num", "Geometric mean of the estimated win ratios"),
                   col("mean_nb", "Mean NB̂", "num", "Average estimated net benefit"),
                   col("deaths", "Deaths", "prob", "Share of patients who die during follow-up"),
                   col("hr_d", "HR (Death)", "num", "Hazard ratio of death"),
                   col("hr_h", "HR (Event)", "num", "Hazard ratio of the non-fatal event"),
                   col("delta_s", "Δ (SDs)", "num", "Difference in the measured outcome"),
                   col("power_fs", "Power (F–S)", "prob", "Simulated power of the Finkelstein–Schoenfeld test"),
                   col("power_wr", "Power (WR)", "prob", "Simulated power of the win ratio's Wald test"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    for (const auto& x : simc::run_columns()) c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string test = s.is("test", "wr") ? "Wald test of the win ratio" : "Finkelstein–Schoenfeld test";
    std::string comp = "death (median " + fnum(r.num("med_d")) + " in the control group, hazard ratio " +
                       fnum(r.num("hr_d")) + ")";
    if (r.has("hr_h") && std::isfinite(r.num("hr_h")))
      comp += ", then a non-fatal event (median " + fnum(r.num("med_h")) + ", hazard ratio " + fnum(r.num("hr_h")) + ")";
    if (r.has("delta_s") && std::isfinite(r.num("delta_s")))
      comp += ", then a measured outcome (difference " + fnum(r.num("delta_s")) + " SD)";
    std::string t = "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " (total " +
                    fint(r.num("n")) + ") " + simc::power_phrase(s, r) + " to detect a win ratio of " +
                    fnum(r.num("wr"), 3) + " (net benefit " + fnum(r.num("nb"), 3) + "; " + fpct(r.num("pwin")) +
                    " of the pairs won, " + fpct(r.num("ploss")) + " lost and " + fpct(r.num("ptie")) +
                    " tied) for the hierarchical composite of " + comp + ", using a " + alt_words(alt(s)) + " " + test +
                    " at α = " + fnum(r.num("alpha")) + ".";
    t += simc::others_sentence(r, "trials",
                               {{"Finkelstein–Schoenfeld test", "power_fs"}, {"win ratio test", "power_wr"}});
    return t + simc::method_note(r);
  };
  p.definitions = {{"N1, N2", "Number of patients in the treatment and control groups; N = N1 + N2."},
                   {"Win, loss, tie", "Comparing a treatment patient with a control patient over the time both were "
                                      "followed: the one who dies first loses; if neither died, the one with the "
                                      "earlier non-fatal event loses; if still undecided, the one with the better "
                                      "measured outcome wins; otherwise it is a tie."},
                   {"Win Ratio", "Wins divided by losses; above 1 favors the treatment."},
                   {"Net Benefit", "Share of pairs won minus share lost (Buyse's net benefit)."},
                   {"Power", "Share of the simulated trials whose test was significant, with its 95% interval "
                             "(Wilson) and Monte Carlo standard error."}};
  p.notes = {"Death and the non-fatal event are independent exponential times within each group; the non-fatal "
             "event counts only if it comes before death. Patients enter uniformly during accrual and are followed "
             "until the analysis, or until they are lost to follow-up.",
             "The Finkelstein–Schoenfeld statistic is wins minus losses; its variance comes from comparing every "
             "pair of patients in the trial (both groups pooled). The comparisons are counted exactly by sorting "
             "(O(N log N)) rather than pair by pair, so large trials simulate quickly.",
             "The win ratio's standard error is the two-sample U-statistic (delta-method) estimate; a trial with no "
             "wins or no losses has no finite win ratio and counts as a failed analysis for that test, as a trial "
             "in which every pair ties does for the Finkelstein–Schoenfeld test.",
             "The expected win, loss and tie probabilities, and the quick approximation (Yu & Ganju 2022), come "
             "from 200,000 simulated pairs drawn from a fixed stream."};
  p.references = {"Pocock, S. J., Ariti, C. A., Collier, T. J. & Wang, D. (2012). The win ratio: a new approach to "
                  "the analysis of composite endpoints in clinical trials based on clinical priorities. European "
                  "Heart Journal 33, 176–182.",
                  "Finkelstein, D. M. & Schoenfeld, D. A. (1999). Combining mortality and longitudinal measures in "
                  "clinical trials. Statistics in Medicine 18, 1341–1354.",
                  "Bebu, I. & Lachin, J. M. (2016). Large sample inference for a win ratio analysis of a composite "
                  "outcome based on prioritized components. Biostatistics 17, 178–187.",
                  "Yu, R. X. & Ganju, J. (2022). Sample size formula for a win ratio endpoint. Statistics in "
                  "Medicine 41, 950–963.",
                  simc::REF_SIM};
  p.curve_param = "n1";
  p.curve_label = "N1 (Patients in the Treatment Group)";
  p.example_solve = "n";
  p.example = {{"power", "0.8"}, {"hr_d", "0.8"}, {"hr_h", "0.75"}, {"reps", "500"}};
  return p;
}
SSC_REGISTER(winratio);

}  // namespace
}  // namespace ssc::procs
