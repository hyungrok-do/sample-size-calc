// Single-arm phase II trials: Simon's two-stage designs (optimal, minimax and
// admissible) found by exhaustive search over all designs, and the exact
// operating characteristics of a given two-stage design.
//
// A design (r1, n1, r, n) treats n1 patients, stops for futility when at most
// r1 respond, otherwise treats n − n1 more and declares the treatment
// promising when more than r of all n patients respond.
#include <algorithm>

#include "procedures/common.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

const char* REF_SIMON =
    "Simon, R. (1989). Optimal two-stage designs for phase II clinical trials. Controlled Clinical Trials 10, 1–10.";
const char* REF_JUNG =
    "Jung, S.-H., Lee, T., Kim, K. & George, S. L. (2004). Admissible two-stage designs for phase II cancer clinical "
    "trials. Statistics in Medicine 23, 561–569.";

// Binomial probability tables for every sample size m ≤ nmax.
struct BinomTables {
  std::vector<std::vector<double>> pmf;    // pmf[m][x]
  std::vector<std::vector<double>> upper;  // upper[m][k] = P(X ≥ k), k = 0..m
  BinomTables(int nmax, double p) : pmf(nmax + 1), upper(nmax + 1) {
    for (int m = 0; m <= nmax; ++m) {
      auto& f = pmf[m];
      f.resize(static_cast<std::size_t>(m) + 1);
      for (int x = 0; x <= m; ++x) f[x] = stats::dbinom_raw(x, m, p, 1 - p);
      auto& u = upper[m];
      u.assign(static_cast<std::size_t>(m) + 1, 0.0);
      double s = 0;  // summed from the top so that small tails keep their accuracy
      for (int k = m; k >= 0; --k) u[k] = (s += f[k]);
    }
  }
  double tail(int m, int k) const { return k <= 0 ? 1.0 : (k > m ? 0.0 : upper[m][k]); }
};

struct Design {
  int r1 = -1, n1 = 0, r = -1, n = 0;
  double alpha = NAN, power = NAN, pet0 = NAN, pet1 = NAN, en0 = NAN, en1 = NAN;
  bool found() const { return n > 0; }
};

// Exact characteristics at response probability p: probability of early
// termination, of declaring the treatment promising, and expected N.
struct OC {
  double pet, reject, en;
};

OC characteristics(int r1, int n1, int r, int n, double p) {
  using namespace stats;
  const int n2 = n - n1;
  OC oc;
  oc.pet = pbinom(r1, n1, p, true);
  double s = 0;
  for (int x1 = r1 + 1; x1 <= n1; ++x1) {
    const int need = r - x1;  // reject when X2 > need
    s += dbinom(x1, n1, p) * (need < 0 ? 1.0 : pbinom(need, n2, p, false));
  }
  oc.reject = clamp01(s);
  oc.en = n1 + (1 - oc.pet) * n2;
  return oc;
}

Design complete(int r1, int n1, int r, int n, double p0, double p1) {
  Design d;
  d.r1 = r1;
  d.n1 = n1;
  d.r = r;
  d.n = n;
  const OC a = characteristics(r1, n1, r, n, p0), b = characteristics(r1, n1, r, n, p1);
  d.alpha = a.reject;
  d.pet0 = a.pet;
  d.en0 = a.en;
  d.power = b.reject;
  d.pet1 = b.pet;
  d.en1 = b.en;
  return d;
}

// Smallest n whose most powerful (randomized) single-stage test at level α
// reaches the power: no two-stage design with a smaller maximum exists.
int lower_bound_n(double p0, double p1, double alpha, double power) {
  using namespace stats;
  for (int n = 1; n < 100000; ++n) {
    const int c = static_cast<int>(qbinom(1 - alpha, n, p0));  // P(X > c) ≤ α < P(X ≥ c)
    const double gt = pbinom(c, n, p0, false), at = dbinom(c, n, p0);
    const double g = at > 0 ? std::min(1.0, std::max(0.0, (alpha - gt) / at)) : 0.0;
    if (pbinom(c, n, p1, false) + g * dbinom(c, n, p1) >= power) return n;
  }
  throw CalcError("the required sample size exceeds 100,000");
}

// Smallest n of a non-randomized single-stage exact design (for the default search range).
int single_stage_n(double p0, double p1, double alpha, double power, int from) {
  using namespace stats;
  for (int n = std::max(1, from); n < 100000; ++n) {
    const int c = static_cast<int>(qbinom(1 - alpha, n, p0));
    int k = c + 1;  // reject when X ≥ k, the smallest k with P(X ≥ k | p0) ≤ α
    while (k > 0 && pbinom(k - 2, n, p0, false) <= alpha) --k;
    while (pbinom(k - 1, n, p0, false) > alpha) ++k;
    if (pbinom(k - 1, n, p1, false) >= power) return n;
  }
  throw CalcError("the required sample size exceeds 100,000");
}

// The design with stage 1 size n1 and maximum size n that has the smallest
// EN0 below `bound`, if any.  EN0 falls as r1 grows, so the largest r1 that
// meets both error rates wins; the search runs down from the largest r1 for
// which stopping early under p1 alone costs at most β.  The smallest cutoff r
// with type I error ≤ α only grows as r1 falls.
Design best_for_stage1(int n, int n1, double alpha, double beta, double bound, const BinomTables& T0,
                       const BinomTables& T1) {
  Design best;
  const int n2 = n - n1;
  const auto& f0 = T0.pmf[n1];
  const auto& f1 = T1.pmf[n1];
  auto reject = [&](const std::vector<double>& f, const BinomTables& T, int r1, int r) {
    double s = 0;
    for (int x1 = r1 + 1; x1 <= n1; ++x1) s += f[x1] * T.tail(n2, r - x1 + 1);
    return s;
  };
  std::vector<double> cdf0;  // P(X1 ≤ r1 | p0) for r1 = 0..top
  double c0 = 0, c1 = 0;
  for (int r1 = 0; r1 < n1 && c1 + f1[r1] <= beta; ++r1) {
    c1 += f1[r1];
    cdf0.push_back(c0 += f0[r1]);
  }
  const int top = static_cast<int>(cdf0.size()) - 1;
  auto en0_of = [&](int r1) { return n1 + (1 - cdf0[static_cast<std::size_t>(r1)]) * n2; };
  if (top < 0 || !(en0_of(top) < bound)) return best;
  int lo = top - 1, hi = n;  // r = n never rejects
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    (reject(f0, T0, top, mid) <= alpha ? hi : lo) = mid;
  }
  int r = hi;
  for (int r1 = top; r1 >= 0; --r1) {
    const double en0 = en0_of(r1);
    if (!(en0 < bound) || (best.found() && en0 > best.en0)) break;
    while (r < n && reject(f0, T0, r1, r) > alpha) ++r;
    if (r >= n) break;  // no cutoff keeps the type I error, for this r1 or any smaller one
    if (reject(f1, T1, r1, r) >= 1 - beta) {  // later hits have the same EN0 and a smaller r1
      best.r1 = r1;
      best.n1 = n1;
      best.r = r;
      best.n = n;
      best.en0 = en0;
    }
  }
  return best;
}

// The design with maximum size n that has the smallest EN0 below `bound` (ties:
// the smallest n1, then r1), if any.  EN0 ≥ n1, so stage 1 sizes from `bound`
// on cannot qualify.
Design search_n(int n, double alpha, double beta, double bound, const BinomTables& T0, const BinomTables& T1,
                const RunContext& ctx) {
  const int last = std::isfinite(bound) ? std::min(n - 1, static_cast<int>(std::ceil(bound)) - 1) : n - 1;
  if (last < 1) return Design{};
  constexpr int CHUNK = 16;
  const std::size_t chunks = static_cast<std::size_t>((last + CHUNK - 1) / CHUNK);
  std::vector<Design> part(chunks);
  // the bound from smaller n prunes most stage-1 sizes, so a step is a fraction of a millisecond until n is in
  // the hundreds: shared between threads only beyond that
  parallel_for(
      chunks, chunks > 32 ? ctx.threads : 1,
      [&](std::size_t c) {
        Design best;
        const int from = 1 + static_cast<int>(c) * CHUNK, to = std::min(last, from + CHUNK - 1);
        for (int n1 = from; n1 <= to; ++n1) {
          const Design d = best_for_stage1(n, n1, alpha, beta, best.found() ? std::min(bound, best.en0) : bound,
                                           T0, T1);
          if (d.found()) best = d;
        }
        part[c] = best;
      },
      ctx.token);
  Design best;
  for (const auto& d : part)
    if (d.found() && (!best.found() || d.en0 < best.en0)) best = d;
  return best;
}

struct SearchResult {
  std::vector<Design> per_n;  // best design for each n (n = 0 when none)
  int n_lo = 0, n_hi = 0;
};

// Maximum sizes are searched in increasing order: a design whose EN0 is not
// below the smallest EN0 at a smaller maximum size is dominated (neither
// minimax, optimal nor admissible), which bounds the search at each size.
SearchResult search(double p0, double p1, double alpha, double power, int nmax, const RunContext& ctx) {
  SearchResult res;
  res.n_lo = lower_bound_n(p0, p1, alpha, power);
  const int single = single_stage_n(p0, p1, alpha, power, res.n_lo);
  res.n_hi = nmax > 0 ? nmax : std::max(single + 15, static_cast<int>(std::ceil(1.6 * single)));
  if (res.n_hi < res.n_lo)
    throw CalcError("no design exists with N ≤ " + std::to_string(res.n_hi) + " (at least " +
                    std::to_string(res.n_lo) + " is needed)");
  if (res.n_hi > 2000) throw CalcError("the search range is too large (N above 2,000); check p0 and p1");
  const BinomTables T0(res.n_hi, p0), T1(res.n_hi, p1);
  res.per_n.assign(static_cast<std::size_t>(res.n_hi - res.n_lo + 1), Design{});
  double frontier = INFINITY;
  for (int n = res.n_lo; n <= res.n_hi; ++n) {
    ctx.check();
    const Design d = search_n(n, alpha, 1 - power, frontier, T0, T1, ctx);
    res.per_n[static_cast<std::size_t>(n - res.n_lo)] = d;
    if (d.found()) frontier = d.en0;
  }
  return res;
}

// Minimax, optimal and the admissible designs between them: the vertices of
// the lower convex hull of (n, EN0) with weights q·n + (1 − q)·EN0, q ∈ [0, 1].
struct Choice {
  Design d;
  std::string label;
  double q_lo = 0, q_hi = 1;
};

std::vector<Choice> admissible(const SearchResult& sr) {
  std::vector<Design> front;  // Pareto front: n increasing, EN0 decreasing
  for (const auto& d : sr.per_n)
    if (d.found() && (front.empty() || d.en0 < front.back().en0 - 1e-12)) front.push_back(d);
  std::vector<Design> hull;
  auto cross = [](const Design& a, const Design& b, const Design& c) {
    return (b.n - a.n) * (c.en0 - a.en0) - (b.en0 - a.en0) * (c.n - a.n);
  };
  for (const auto& d : front) {
    while (hull.size() >= 2 && cross(hull[hull.size() - 2], hull.back(), d) <= 0) hull.pop_back();
    hull.push_back(d);
  }
  std::vector<Choice> out;
  for (std::size_t i = 0; i < hull.size(); ++i) {
    Choice c;
    c.d = hull[i];
    auto q_between = [&](const Design& a, const Design& b) {
      const double s = (a.en0 - b.en0) / (b.n - a.n);
      return s / (1 + s);
    };
    c.q_hi = i == 0 ? 1.0 : q_between(hull[i - 1], hull[i]);
    c.q_lo = i + 1 == hull.size() ? 0.0 : q_between(hull[i], hull[i + 1]);
    c.label = i == 0 ? (hull.size() == 1 ? "Minimax and optimal" : "Minimax")
                     : (i + 1 == hull.size() ? "Optimal" : "Admissible");
    out.push_back(c);
  }
  return out;
}

std::string ratio(int a, int b) { return std::to_string(a) + "/" + std::to_string(b); }

void oc_table(Row& r, const Design& d, double p0, double p1) {
  Table t;
  t.title = "Operating Characteristics";
  t.columns = {{"p", "True Response Rate"}, {"pet", "P(Stop After Stage 1)"}, {"rej", "P(Declare Promising)"},
               {"en", "Expected N"}};
  t.formats = {"num", "prob", "prob", "num"};
  std::vector<double> ps;
  for (int k = 0; k <= 4; ++k) ps.push_back(p0 + (p1 - p0) * k / 4.0);
  const double step = (p1 - p0) / 4.0;
  if (p0 - step > 0) ps.insert(ps.begin(), p0 - step);
  if (p1 + step < 1) ps.push_back(p1 + step);
  for (double p : ps) {
    const OC oc = characteristics(d.r1, d.n1, d.r, d.n, p);
    t.rows.push_back({Cell(p), Cell(oc.pet), Cell(oc.reject), Cell(oc.en)});
  }
  r.tables.push_back(std::move(t));
}

void put_design(Row& r, const Design& d) {
  r.set("r1", d.r1);
  r.set("n1", d.n1);
  r.set("r", d.r);
  r.set("n", d.n);
  r.set("actual_alpha", d.alpha);
  r.set("pet0", d.pet0);
  r.set("en0", d.en0);
  r.set("pet1", d.pet1);
  r.set("en1", d.en1);
}

Procedure simon() {
  Procedure p;
  p.id = "prop.one.simon";
  p.name = "Simon's Two-Stage Designs (Single-Arm Phase II)";
  p.category = "Proportions";
  p.subcategory = "Single-Arm Phase II";
  p.order = 1;
  p.summary =
      "Optimal, minimax and admissible two-stage designs for a single-arm phase II trial testing H0: p ≤ p0 against "
      "H1: p ≥ p1, found by exhaustive search with exact binomial probabilities; or the exact operating "
      "characteristics of a given design.";
  p.keywords = {"simon", "two-stage", "phase ii", "single arm", "response rate", "optimal", "minimax", "futility",
                "oncology", "early stopping for futility", "cancer", "single-stage"};
  p.solve_for = {{"n", "Sample Size (Design Search)", "n", "N (Maximum)"}, {"power", "Power (Given Design)", "power",
                                                                          "Power"}};
  p.params = {
      select("criterion", "Design Criterion",
             {{"optimal", "Optimal (smallest expected N under p0)", {}},
              {"minimax", "Minimax (smallest maximum N)", {}}},
             "optimal")
          .grp(G_TEST)
          .hint("Optimal designs need the fewest patients on average when the treatment is not worth pursuing; "
                "minimax designs have the smallest maximum. The report lists every admissible design in between.")
          .show(when(solving({"n"}))),
      power_param(),
      alpha_param("0.05", "Alpha (One-Sided)").range(0, 0.5, true, true),
      num("p0", "P0 (Response Rate Not Worth Pursuing)", "0.1")
          .grp(G_EFFECT)
          .sym("P0")
          .range(0, 1, true, true)
          .sugg({"0.1", "0.2", "0.05 0.1 0.2"})
          .hint("A response rate so low that the treatment is not worth further study (often the rate of the "
                "standard treatment), between 0 and 1."),
      num("p1", "P1 (Response Rate Worth Pursuing)", "0.3")
          .grp(G_EFFECT)
          .sym("P1")
          .range(0, 1, true, true)
          .sugg({"0.3", "0.4", "0.25 0.3"})
          .hint("A response rate high enough to make the treatment worth pursuing; the power is the chance of "
                "calling it promising at this rate."),
      num("nmax", "Largest N to Search (0 = Automatic)", "0")
          .grp(G_N)
          .sym("N Max")
          .at_least(0)
          .whole()
          .single()
          .hint("The largest total sample size the design search considers. With 0 it is 1.6 times the "
                "single-stage sample size, but at least 15 more than that size.")
          .show(when(solving({"n"}))),
      num("r1", "r1 (Stop if Responses in Stage 1 ≤ r1)", "1")
          .grp(G_N)
          .sym("r1")
          .at_least(0)
          .whole()
          .hint("Stopping rule of stage 1: the trial stops early for futility if at most r1 of the first n1 "
                "patients respond.")
          .show(when(solving({"power"}))),
      num("n1", "n1 (Stage 1 Sample Size)", "10")
          .grp(G_N)
          .sym("n1")
          .at_least(1)
          .whole()
          .hint("Number of patients treated in stage 1, before the interim decision.")
          .show(when(solving({"power"}))),
      num("r", "r (Promising if Total Responses > r)", "5")
          .grp(G_N)
          .sym("r")
          .at_least(0)
          .whole()
          .hint("Final decision rule: the treatment is called promising if more than r of all n patients "
                "respond.")
          .show(when(solving({"power"}))),
      num("n", "n (Total Sample Size)", "29")
          .grp(G_N)
          .sym("n")
          .at_least(2)
          .whole()
          .hint("Total number of patients over both stages when the trial does not stop after stage 1.")
          .show(when(solving({"power"}))),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double p0 = s("p0"), p1 = s("p1"), alpha = s("alpha");
    if (!(p1 > p0)) throw CalcError("P1 must be larger than P0");
    Design d;
    if (s.solve == "n") {
      const double power = s("power");
      if (!(power > alpha)) throw CalcError("the power must exceed alpha");
      const SearchResult sr = search(p0, p1, alpha, power, static_cast<int>(s("nmax")), ctx);
      const auto choices = admissible(sr);
      if (choices.empty())
        throw CalcError("no two-stage design with N ≤ " + std::to_string(sr.n_hi) +
                        " meets the error rates; increase the largest N searched");
      const Choice& pick = s.is("criterion", "minimax") ? choices.front() : choices.back();
      d = complete(pick.d.r1, pick.d.n1, pick.d.r, pick.d.n, p0, p1);
      Table t;
      t.title = "Minimax, Admissible and Optimal Designs";
      t.columns = {{"design", "Design"}, {"q", "Weight q on N"}, {"s1", "r1/n1"}, {"s2", "r/n"},
                   {"en0", "E(N | p0)"},  {"pet0", "PET(p0)"},    {"alpha", "Alpha"}, {"power", "Power"}};
      t.formats = {"text", "text", "text", "text", "num", "prob", "prob", "prob"};
      for (const auto& c : choices) {
        const Design e = complete(c.d.r1, c.d.n1, c.d.r, c.d.n, p0, p1);
        t.rows.push_back({Cell(c.label), Cell(fnum(c.q_lo, 3) + "–" + fnum(c.q_hi, 3)), Cell(ratio(e.r1, e.n1)),
                          Cell(ratio(e.r, e.n)), Cell(e.en0), Cell(e.pet0), Cell(e.alpha), Cell(e.power)});
      }
      r.tables.push_back(std::move(t));
      if (choices.back().d.n >= sr.n_hi - 1)
        r.notes.push_back("The optimal design lies at the edge of the search range (N ≤ " + std::to_string(sr.n_hi) +
                          "); a larger range may find a design with a smaller expected sample size.");
      r.set("design", pick.label);
    } else {
      const int r1 = static_cast<int>(s("r1")), n1 = static_cast<int>(s("n1")), rr = static_cast<int>(s("r")),
                n = static_cast<int>(s("n"));
      if (n > 100000) throw CalcError("designs are evaluated for up to 100,000 patients, as in the design search");
      if (!(n1 < n)) throw CalcError("n1 must be smaller than n");
      if (!(r1 < n1)) throw CalcError("r1 must be smaller than n1");
      if (!(rr >= r1 && rr < n)) throw CalcError("r must lie between r1 and n − 1");
      d = complete(r1, n1, rr, n, p0, p1);
      r.set("design", "Given");
    }
    put_power(s, r, d.power);
    put_design(r, d);
    oc_table(r, d, p0, p1);
    r.set("p0", p0);
    r.set("p1", p1);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("r1", "r1", "int", "Stop after stage 1 if responses ≤ r1"),
                   col("n1", "n1", "int", "Stage 1 sample size"),
                   col("r", "r", "int", "Declare the treatment promising if total responses > r"),
                   col("n", "n", "int", "Maximum (total) sample size"),
                   col("en0", "E(N | p0)", "num", "Expected sample size when p = p0"),
                   col("pet0", "PET(p0)", "prob", "Probability of early termination when p = p0"),
                   col("actual_alpha", "Actual Alpha", "prob", "Exact type I error of the design"),
                   col("p0", "P0", "num"), col("p1", "P1", "num"), col("alpha", "Alpha", "alpha"),
                   col("design", "Design", "text")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string kind = s.solve == "n" ? (s.is("criterion", "minimax") ? "minimax" : "optimal") : "given";
    const int r1 = static_cast<int>(r.num("r1")), n1 = static_cast<int>(r.num("n1")),
              rr = static_cast<int>(r.num("r")), n = static_cast<int>(r.num("n"));
    return "For P0 = " + fnum(r.num("p0")) + " versus P1 = " + fnum(r.num("p1")) + ", the " + kind +
           " two-stage design enters " + fint(n1) + " patients in stage 1 and stops for futility if " +
           (r1 == 0 ? std::string("no patient responds") : fint(r1) + " or fewer respond") + "; otherwise " +
           fint(n - n1) +
           " more patients are entered (" + fint(n) + " in total) and the treatment is considered promising if more "
           "than " + fint(rr) + " respond. Its exact one-sided type I error is " + fnum(r.num("actual_alpha"), 4) +
           " (target " + fnum(r.num("alpha")) + ") and its power is " + fnum(r.num("power"), 4) +
           ". When P = P0 the trial stops after stage 1 with probability " + fnum(r.num("pet0"), 3) +
           " and the expected sample size is " + fnum(r.num("en0"), 1) + ".";
  };
  p.definitions = {
      {"r1, n1", "Stage 1: treat n1 patients; stop for futility if at most r1 respond."},
      {"r, n", "Treat n patients in total; the treatment is promising if more than r respond."},
      {"E(N | p0)", "Expected number of patients when the true response rate is P0."},
      {"PET(p0)", "Probability of stopping after stage 1 when the true response rate is P0."},
      {"Optimal", "The design with the smallest expected sample size under P0."},
      {"Minimax", "The design with the smallest maximum sample size (ties: smallest E(N | p0))."},
      {"Admissible", "Designs minimizing q·n + (1 − q)·E(N | p0) for some weight q between 0 and 1."}};
  p.references = {REF_SIMON, REF_JUNG};
  p.notes = {"Error rates are exact (binomial); the type I error is one-sided."};
  p.example_solve = "n";
  p.example = {{"p0", "0.1 0.2"}, {"p1", "0.3"}, {"alpha", "0.05"}, {"power", "0.8"}};
  return p;
}
SSC_REGISTER(simon);

}  // namespace
}  // namespace ssc::procs
