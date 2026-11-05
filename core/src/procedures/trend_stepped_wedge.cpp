// Cochran–Armitage test for a linear trend in proportions, and stepped-wedge
// cluster-randomized designs (Hussey & Hughes).
#include <sstream>

#include "procedures/apart.hpp"
#include "procedures/common.hpp"
#include "procedures/derive_props.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// =====================================================================================
// Cochran–Armitage trend test
// =====================================================================================
// For group sizes n_i, scores d_i and proportions p_i, the statistic
// T = Σ n_i (d_i − d̄) p̂_i has mean Σ n_i (d_i − d̄) p_i, null variance
// p̄ q̄ Σ n_i (d_i − d̄)² and alternative variance Σ n_i (d_i − d̄)² p_i q_i
// (Nam, 1987).
struct Trend {
  std::vector<double> p, d, w;  // proportions, scores, relative group sizes
};

Trend trend_setup(const Scenario& s) {
  Trend t;
  t.p = s.v("props");
  const std::size_t k = t.p.size();
  if (k < 2) throw CalcError("enter at least two proportions");
  for (double x : t.p)
    if (!(x > 0 && x < 1)) throw CalcError("proportions must lie strictly between 0 and 1");
  if (s.is("scores", "custom")) {
    t.d = s.v("dvals");
    if (t.d.size() != k) throw CalcError("enter one score per group (" + std::to_string(k) + ")");
  } else {
    for (std::size_t i = 0; i < k; ++i) t.d.push_back(static_cast<double>(i));
  }
  if (s.is("alloc", "ratio")) {
    t.w = s.v("weights");
    if (t.w.size() != k) throw CalcError("enter one relative group size per group (" + std::to_string(k) + ")");
    for (double x : t.w)
      if (!(x > 0)) throw CalcError("relative group sizes must be positive");
  } else {
    t.w.assign(k, 1.0);
  }
  return t;
}

// n times the relative size, rounded up (at least 1 however small the size)
std::vector<double> trend_sizes(const Trend& t, double n) {
  std::vector<double> ni;
  for (double w : t.w) ni.push_back(solve::round_up(n * w));
  return ni;
}

// The direction of the trend: the sign of Σ wᵢ(dᵢ − d̄)pᵢ with the relative group sizes w (0: no linear trend).
// Here and in the power the scores are taken from the first one (dᵢ − d₁): a mean of the scores themselves rounds
// when they are large next to their spread, and leaves a trend that is not there.
double trend_direction(const Trend& t) {
  double W = 0, dbar = 0;
  for (std::size_t i = 0; i < t.w.size(); ++i) {
    W += t.w[i];
    dbar += t.w[i] * (t.d[i] - t.d[0]);
  }
  dbar /= W;
  double e = 0, scale = 0;
  for (std::size_t i = 0; i < t.w.size(); ++i) {
    e += t.w[i] * (t.d[i] - t.d[0] - dbar) * t.p[i];
    scale += t.w[i] * std::abs(t.d[i] - t.d[0] - dbar);
  }
  return std::abs(e) <= 1e-12 * scale ? 0.0 : e;
}

// The power and its parts (which the derivation shows): the group sizes and their total, the mean score (from the
// first) and proportion, the statistic's mean and its variances under H0 and H1, the critical value and each tail.
struct TrendPower {
  std::vector<double> ni;
  double N, dbar, pbar, e, v0, v1, z, up, lo, power;
};

TrendPower trend_parts(const Trend& t, double n, double alpha, Alt a) {
  using namespace stats;
  TrendPower r;
  r.ni = trend_sizes(t, n);
  const auto& ni = r.ni;
  double N = 0, dbar = 0, pbar = 0;
  for (std::size_t i = 0; i < ni.size(); ++i) {
    N += ni[i];
    dbar += ni[i] * (t.d[i] - t.d[0]);
    pbar += ni[i] * t.p[i];
  }
  dbar /= N;
  pbar /= N;
  double e = 0, v0 = 0, v1 = 0;
  for (std::size_t i = 0; i < ni.size(); ++i) {
    const double dd = t.d[i] - t.d[0] - dbar;
    e += ni[i] * dd * t.p[i];
    v0 += ni[i] * dd * dd;
    v1 += ni[i] * dd * dd * t.p[i] * (1 - t.p[i]);
  }
  v0 *= pbar * (1 - pbar);
  const double z = qnorm(tail_alpha(alpha, a), false);
  r.up = pnorm((e - z * std::sqrt(v0)) / std::sqrt(v1));
  r.lo = pnorm((-e - z * std::sqrt(v0)) / std::sqrt(v1));
  double pw = 0;
  if (a != Alt::Lower) pw += r.up;
  if (a != Alt::Upper) pw += r.lo;
  r.N = N, r.dbar = dbar, r.pbar = pbar, r.e = e, r.v0 = v0, r.v1 = v1, r.z = z;
  r.power = stats::clamp01(pw);
  return r;
}

double trend_power(const Trend& t, double n, double alpha, Alt a) { return trend_parts(t, n, alpha, a).power; }

// How the trend test's answer comes about, following what was solved for: the group sizes, the scores' and
// proportions' means, the statistic's mean and variances and the power at n (each term listed for up to 6 groups);
// solving for n, Nam's closed form for comparison and the power one size below.
derive::Derivation trend_derivation(const Scenario& s, const Row& r) {
  using namespace derive;
  const Trend t = trend_setup(s);
  const Alt a = alt(s);
  const double n = r.num("n"), alpha = r.num("alpha"), power = r.num("power"), target = s.get("power", NAN);
  const bool solved = s.solve == "n", ratio = s.is("alloc", "ratio");
  const TrendPower q = trend_parts(t, n, alpha, a);
  const std::size_t k = q.ni.size();
  const bool few = k <= 6;
  const std::string Nn = size(n, solved);
  std::string sizes, e_terms, v_terms;
  for (std::size_t i = 0; i < k; ++i) {
    const double dd = t.d[i] - t.d[0] - q.dbar;
    sizes += (i ? ", " : "") + dvi(q.ni[i]);
    e_terms += (i ? " + " : "") + dvi(q.ni[i]) + " × " + dvs(dd) + " × " + gv(t.p[i]);
    v_terms += (i ? " + " : "") + dvi(q.ni[i]) + " × " + dvs(dd) + "²";
  }
  const std::string unit = ratio ? "n, the size of a group of relative size 1," : "n per group";
  Derivation out;
  out.find = solved ? find_size(unit, target)
                    : find_power("with n = " + Nn + (ratio ? " in a group of relative size 1" : " per group"));
  if (solved) {
    // Nam's closed form: the statistic's mean and variances per unit of n, with the relative sizes unrounded, beside
    // the answer below
    double W = 0, d1 = 0, p1 = 0;
    for (std::size_t i = 0; i < k; ++i) {
      W += t.w[i];
      d1 += t.w[i] * (t.d[i] - t.d[0]);
      p1 += t.w[i] * t.p[i];
    }
    d1 /= W;
    p1 /= W;
    double e1 = 0, v01 = 0, v11 = 0;
    for (std::size_t i = 0; i < k; ++i) {
      const double dd = t.d[i] - t.d[0] - d1;
      e1 += t.w[i] * dd * t.p[i];
      v01 += t.w[i] * dd * dd;
      v11 += t.w[i] * dd * dd * t.p[i] * (1 - t.p[i]);
    }
    v01 *= p1 * (1 - p1);
    const double za = stats::qnorm(tail_alpha(alpha, a), false), zb = stats::qnorm(target);
    out.steps.push_back(
        {"Normal approximation",
         std::string("n ≈ (z(") + level_symbol(a) +
             ") √V0₁ + z(1 − β) √V1₁)² / E₁² (Nam), E₁, V0₁ and V1₁ the mean and variances below with nᵢ = " +
             (ratio ? "wᵢ, unrounded" : "1") + ", for comparison; the steps below give the answer",
         "E₁ = " + dv(e1) + ", V0₁ = " + dv(v01) + ", V1₁ = " + dv(v11) + ": (" + dv(za) + " × √" + dv(v01) + " + " +
             dvs(zb) + " × √" + dv(v11) + ")² / " + dvs(e1) + "² = " +
             dv(std::pow(za * std::sqrt(v01) + zb * std::sqrt(v11), 2) / (e1 * e1))});
  }
  out.steps.push_back({"Group sizes", ratio ? "nᵢ = ⌈n wᵢ⌉ for the relative sizes wᵢ, N = Σ nᵢ" : "nᵢ = n, N = k n",
                       "n = " + Nn + ": " + sizes + ", N = " + dvi(q.N)});
  out.steps.push_back({"Means", "d̄ = Σ nᵢ (dᵢ − d₁)/N with the scores dᵢ measured from the first, p̄ = Σ nᵢ pᵢ/N",
                       "d̄ = " + dv(q.dbar) + ", p̄ = " + dv(q.pbar)});
  out.steps.push_back({"Mean of the trend statistic", "E(T) = Σ nᵢ (dᵢ − d₁ − d̄) pᵢ, T = Σ nᵢ (dᵢ − d̄) p̂ᵢ",
                       (few ? e_terms : "Σ over the " + dvi(k) + " groups") + " = " + dv(q.e)});
  out.steps.push_back({"Variances under H0 and H1",
                       "V0 = p̄(1 − p̄) Σ nᵢ (dᵢ − d₁ − d̄)², V1 = Σ nᵢ (dᵢ − d₁ − d̄)² pᵢ(1 − pᵢ)",
                       "V0 = " + dv(q.pbar) + " × " + dv(1 - q.pbar) + " × (" +
                           (few ? v_terms : "Σ over the " + dvi(k) + " groups") + ") = " + dv(q.v0) + ", V1 = " +
                           dv(q.v1)});
  out.steps.push_back(z_critical(alpha, a));
  out.steps.push_back(tails_power(a, "E(T)", "", "√V0", "√V1", q.e, 0, q.z, std::sqrt(q.v0), std::sqrt(q.v1), q.up,
                                  q.lo, power));
  if (solved) {
    out.steps.push_back(n_search([&](double x) { return trend_power(t, x, alpha, a); }, n, target, 2, "n"));
    out.answer =
        "n = " + dvi(n) + " (group sizes " + sizes + ", N = " + dvi(q.N) + "), with power " + dvp(power, target);
  } else {
    out.answer = "power = " + dvp(power);
  }
  return out;
}

Procedure trend() {
  Procedure p;
  p.id = "prop.trend";
  p.name = "Cochran–Armitage Test for a Trend in Proportions";
  p.category = "Proportions";
  p.subcategory = "Several Proportions";
  p.order = 1;
  p.summary =
      "Sample size and power for the Cochran–Armitage test of a linear trend in response proportions across k "
      "ordered groups (for example increasing doses), with equal or unequal group sizes and arbitrary dose scores "
      "(Nam, 1987).";
  p.keywords = {"cochran-armitage", "trend", "dose response", "ordered groups", "proportions", "linear trend", "ordinal exposure", "dose-response trend"};
  p.solve_for = {{"n", "Sample Size", "n", "n (Per Group)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (any trend)", "One-sided (increasing)", "One-sided (decreasing)"),
      power_param(),
      alpha_param(),
      n_param("n", "n (Size of a Group with Relative Size 1)", "50", 2)
          .sym("n")
          .hint("Number of subjects in each group; with relative sizes, the size of a group with relative size 1 "
                "(the others get n × their size, rounded up)."),
      select("alloc", "Group Sizes", {{"equal", "Equal", {}}, {"ratio", "Relative sizes", {}}}, "equal")
          .grp(G_N)
          .hint("Give every group the same number of subjects, or size the groups in proportion to relative "
                "sizes."),
      vec("weights", "Relative Group Sizes", "1 1 1 2")
          .grp(G_N)
          .sym("Relative Sizes")
          .above(0)
          .hint("One relative size per group, in the order of the proportions; 1 1 2 makes the last group twice "
                "as large as the others.")
          .show(when(is("alloc", {"ratio"}))),
      vec("props", "Response Proportions by Group", "0.1 0.15 0.2 0.25")
          .grp(G_EFFECT)
          .sym("Proportions")
          .range(0, 1, true, true)
          .sugg({"0.1 0.15 0.2 0.25", "0.2 0.2 0.3", "0.05 0.1 0.2"})
          .hint("The share of subjects you expect to respond in each group, listed in the order of the groups "
                "(e.g. by increasing dose)."),
      select("scores", "Dose Scores", {{"equal", "Equally spaced (0, 1, 2, …)", {}}, {"custom", "Custom scores", {}}}, "equal")
          .grp(G_EFFECT)
          .hint("The values the trend test assigns to the ordered groups: equally spaced, or your own scores such "
                "as the actual doses."),
      vec("dvals", "Scores (e.g. Doses)", "0 1 2 4")
          .grp(G_EFFECT)
          .sym("Scores")
          .hint("One score per group, in the order of the proportions, such as the dose given to each group.")
          .show(when(is("scores", {"custom"}))),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Trend t = trend_setup(s);
    const double alpha = s("alpha");
    const Alt a = alt(s);
    double n;
    if (s.solve == "n") {
      const double dir = trend_direction(t);
      if (dir == 0) throw CalcError("the proportions have no linear trend across the scores");
      if ((a == Alt::Upper && dir < 0) || (a == Alt::Lower && dir > 0))
        throw CalcError("the trend in the proportions runs opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = 2;
      // groups that tiny relative sizes keep at 1 subject up to the largest n hold the power down however large n
      // grows
      const std::string held = solve::single_groups(t.w, o.hi, "relative size", "relative sizes");
      if (!held.empty() && trend_power(t, o.hi, alpha, a) < s("power"))
        throw CalcError("the target power cannot be reached with a sample size up to " + text::format_g(o.hi, 10) +
                        ", where " + held);
      n = solve::min_n([&](double x) { return trend_power(t, x, alpha, a); }, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    const auto ni = trend_sizes(t, n);
    double total = 0;
    for (double x : ni) total += x;
    r.set("n", n);
    r.set("N", total);
    put_power(s, r, trend_power(t, n, alpha, a));
    r.set("sizes", list_text(ni));
    // never a proportion of 0.99999 as 1, nor 0.10001 as 0.1 next to 0.1, and the scores alike
    std::string props, scores;
    for (const auto& x : fnum_apart_all(t.p)) props += (props.empty() ? "" : ", ") + x;
    for (const auto& x : fnum_apart_all(t.d)) scores += (scores.empty() ? "" : ", ") + x;
    r.set("props", props);
    r.set("scores", scores);
    r.set("k", static_cast<double>(t.p.size()));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "n", "int", "Size of a group with relative size 1"), col("N", "N (Total)", "int"),
                   col("k", "Groups", "int"), col("sizes", "Group Sizes", "text"), col("props", "Proportions", "text"),
                   col("scores", "Scores", "text"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sizes of " + text_of(r, "sizes") + " (" + fint(r.num("N")) + " in total) " + achieves(s, r) +
           " to detect a linear trend in the response proportions " + text_of(r, "props") + " (scores " +
           text_of(r, "scores") + ") with a " + alt_words(alt(s)) + " Cochran–Armitage test at α = " +
           fnum(r.num("alpha")) + ".";
  };
  p.derivation = trend_derivation;
  p.definitions = {{"Scores", "Values assigned to the ordered groups (for example doses)."},
                   {"Group Sizes", "Subjects per group: n times the relative size, rounded up."}};
  p.notes = {"Power uses the normal approximation without continuity correction."};
  p.curve_param = "n";
  p.curve_label = "n (Per Group)";
  p.curve_min = 2;
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"props", "0.1 0.15 0.2 0.25"}};
  return p;
}
SSC_REGISTER(trend);

// =====================================================================================
// Stepped-wedge cluster-randomized design
// =====================================================================================
// Cross-sectional design with S steps: I = S·g clusters, T = S + 1 periods; the
// g clusters of step k switch from control to intervention at period k + 1.  With
// m subjects per cluster-period, σ² = σ²_within/m and τ² the between-cluster
// variance, Hussey & Hughes (2007) give
//   Var(θ̂) = I σ² (σ² + T τ²) / [(IU − W) σ² + (U² + ITU − TW − IV) τ²],
// U = Σ X_ij, W = Σ_j (Σ_i X_ij)², V = Σ_i (Σ_j X_ij)².
// The variance and its parts (which the derivation shows).
struct SwDesign {
  double I, T, tau2, s2, U, W, V, var;
};

SwDesign sw_design(int steps, int per_step, double m, double sigma2, double icc) {
  // in floating point: S·g clusters overflow an int for large designs
  const double g = per_step, I = static_cast<double>(steps) * g, T = steps + 1.0;
  const double tau2 = icc * sigma2, s2 = (1 - icc) * sigma2 / m;
  double U = 0, W = 0, V = 0;
  for (int k = 1; k <= steps; ++k) {
    // clusters of step k are treated in periods k + 1 .. T
    const double row = T - k;
    U += g * row;
    V += g * row * row;
  }
  // in period j the clusters of steps 1 .. j − 1 are treated: (j − 1)·g of them
  for (int j = 1; j <= steps + 1; ++j) W += (g * (j - 1)) * (g * (j - 1));
  const double var = I * s2 * (s2 + T * tau2) / ((I * U - W) * s2 + (U * U + I * T * U - T * W - I * V) * tau2);
  return {I, T, tau2, s2, U, W, V, var};
}

double sw_variance(int steps, int per_step, double m, double sigma2, double icc) {
  return sw_design(steps, per_step, m, sigma2, icc).var;
}

struct SwEffect {
  double delta, sigma2;
  double unit = 1;  // the unit of σ² (σ itself when σ² would underflow)
};

SwEffect sw_effect(const Scenario& s) {
  // with no effect the power is alpha, which neither more clusters nor larger ones raise
  const bool power = s.solve == "power";
  if (s.is("outcome", "binary")) {
    const double p1 = s("p1"), p2 = s("p2");
    if (p1 == p2 && !power) throw CalcError("P1 equals P2");
    return {p1 - p2, (p1 * (1 - p1) + p2 * (1 - p2)) / 2};
  }
  const double d = s("delta"), sd = s("sigma");
  if (d == 0 && !power) throw CalcError("δ = 0: there is no effect to detect");
  // σ² underflows to 0 for a σ below about 1e-154 (the variance would be 0/0): it is then 1 in units of σ
  if (sd * sd < stats::DBL_MIN_) return {d, 1.0, sd};
  return {d, sd * sd};
}

double sw_se(int steps, int per_step, double m, const SwEffect& e, double icc) {
  return std::sqrt(sw_variance(steps, per_step, m, e.sigma2, icc)) * e.unit;
}

double sw_power(int steps, int per_step, double m, const SwEffect& e, double icc, double alpha) {
  using namespace stats;
  const double se = sw_se(steps, per_step, m, e, icc);
  const double z = qnorm(half_alpha(alpha), false), d = std::fabs(e.delta) / se;
  return pnorm(d - z) + pnorm(-d - z);
}

// How the stepped-wedge design's power comes about from Hussey & Hughes's variance, following what was solved for:
// solving for the clusters per step or the subjects per cluster-period, the power at the number found and one fewer.
derive::Derivation sw_derivation(const Scenario& s, const Row& r) {
  using namespace derive;
  const int steps = static_cast<int>(r.num("steps")), g = static_cast<int>(r.num("g"));
  const double m = r.num("m"), icc = r.num("icc"), alpha = r.num("alpha"), power = r.num("power");
  const SwEffect e = sw_effect(s);
  const SwDesign w = sw_design(steps, g, m, e.sigma2, icc);
  const double se = r.num("se"), z = stats::qnorm(half_alpha(alpha), false), lambda = std::fabs(e.delta) / se;
  const bool for_g = s.solve == "g", for_m = s.solve == "m", binary = s.is("outcome", "binary");
  const std::string G = size(g, for_g), Mm = size(m, for_m), S = gvi(steps);
  // in units of σ² where σ² would underflow (sw_effect), and SE = √Var(θ̂) × σ
  const std::string unit = e.unit == 1 ? "" : " in units of σ²";
  const std::string se_text = e.unit == 1 ? dsig(se) : "√" + dsig(w.var) + " × " + gv(e.unit) + " = " + dsig(se);
  Derivation out;
  if (for_g)
    out.find = find_size("number of clusters g switching at each step", s("power"));
  else if (for_m)
    out.find = find_size("number of subjects m per cluster and period", s("power"));
  else
    out.find = find_power("with S = " + S + " steps, g = " + G + " clusters switching at each step and m = " + Mm +
                          " subjects per cluster and period");
  if (binary)
    out.steps.push_back({"Effect and subject variance", "δ = P1 − P2, σ² = (P1(1 − P1) + P2(1 − P2))/2",
                         "δ = " + gv(s("p1")) + " − " + gv(s("p2")) + " = " + dv(e.delta) + ", σ² = (" + gv(s("p1")) +
                             " × " + dv(1 - s("p1")) + " + " + gv(s("p2")) + " × " + dv(1 - s("p2")) + ")/2 = " +
                             dv(e.sigma2)});
  else
    out.steps.push_back({"Subject variance", "σ², the square of σ",
                         e.unit == 1 ? gv(s("sigma")) + "² = " + dsig(e.sigma2)
                                     : "σ² = " + dvi(1) + " in units of σ² (" + gv(s("sigma")) + "² underflows)"});
  out.steps.push_back({"Design", "I = S g clusters, T = S + 1 periods",
                       "I = " + S + " × " + G + " = " + dvi(w.I) + ", T = " + S + " + 1 = " + dvi(w.T)});
  out.steps.push_back(
      {"Variance components", "τ² = ρ σ² between clusters, σₑ² = (1 − ρ) σ²/m for a cluster-period's mean",
       "τ² = " + gv(icc) + " × " + dsig(e.sigma2) + " = " + dsig(w.tau2) + ", σₑ² = " + dv(1 - icc) + " × " +
           dsig(e.sigma2) + "/" + Mm + " = " + dsig(w.s2) + unit});
  out.steps.push_back({"Treatment counts",
                       "U = Σᵢⱼ Xᵢⱼ, W = Σⱼ (Σᵢ Xᵢⱼ)², V = Σᵢ (Σⱼ Xᵢⱼ)², Xᵢⱼ = 1 when cluster i is treated in period j",
                       "U = " + dvi(w.U) + ", W = " + dvi(w.W) + ", V = " + dvi(w.V)});
  const std::string I = dvi(w.I), T = dvi(w.T), U = dvi(w.U), W = dvi(w.W), V = dvi(w.V);
  out.steps.push_back({"Standard error of the effect",
                       "SE = √Var(θ̂), Var(θ̂) = I σₑ² (σₑ² + T τ²)/((I U − W) σₑ² + (U² + I T U − T W − I V) τ²) "
                       "(Hussey & Hughes)",
                       "Var(θ̂) = " + I + " × " + dsig(w.s2) + " × (" + dsig(w.s2) + " + " + T + " × " + dsig(w.tau2) +
                           ")/((" + I + " × " + U + " − " + W + ") × " + dsig(w.s2) + " + (" + U + "² + " + I + " × " +
                           T + " × " + U + " − " + T + " × " + W + " − " + I + " × " + V + ") × " + dsig(w.tau2) +
                           ") = " + dsig(w.var) + unit + ", SE = " + se_text});
  out.steps.push_back(z_critical(alpha, Alt::Two));
  const double up = stats::pnorm(lambda - z), lo = stats::pnorm(-lambda - z);
  // |δ| is the input itself only for a continuous effect that is not negative
  const std::string D = !binary && e.delta >= 0 ? gv(e.delta) : dv(std::fabs(e.delta));
  out.steps.push_back({"Power", "Φ(λ − z) + Φ(−λ − z), λ = |δ|/SE",
                       "λ = " + D + "/" + dsig(se) + " = " + dv(lambda) + "; Φ(" + dv(lambda) + " − " + dv(z) +
                           ") + Φ(" + dv(-lambda) + " − " + dv(z) + ") = " + dvp(up) + " + " + dvp(lo) + " = " +
                           dvp(power)});
  const std::string N = dvi(r.num("N"));
  if (for_g) {
    out.steps.push_back(n_search([&](double x) { return sw_power(steps, static_cast<int>(x), m, e, icc, alpha); }, g,
                                 s("power"), 1, "g"));
    out.answer = "g = " + dvi(g) + " clusters per step (I = " + I + " clusters, N = " + N +
                 " observations), with power " + dvp(power, s("power"));
  } else if (for_m) {
    out.steps.push_back(
        n_search([&](double x) { return sw_power(steps, g, x, e, icc, alpha); }, m, s("power"), 1, "m"));
    out.answer = "m = " + dvi(m) + " subjects per cluster and period (N = " + N + " observations), with power " +
                 dvp(power, s("power"));
  } else {
    out.answer = "power = " + dvp(power);
  }
  return out;
}

Procedure stepped_wedge() {
  Procedure p;
  p.id = "cluster.sw";
  p.name = "Stepped-Wedge Cluster-Randomized Designs";
  p.category = "Cluster-Randomized";
  p.subcategory = "Stepped Wedge";
  p.order = 1;
  p.summary =
      "Power and the number of clusters (or of subjects per cluster and period) for a cross-sectional stepped-wedge "
      "trial in which groups of clusters switch from control to intervention at successive steps, using the "
      "variance of Hussey & Hughes (2007) with period effects and a random cluster effect.";
  p.keywords = {"stepped wedge", "cluster randomized", "icc", "crossover", "implementation", "hussey hughes", "staggered rollout", "sw-crt"};
  p.solve_for = {{"power", "Power", "power", "Power"},
                 {"g", "Clusters per Step", "g", "Clusters per Step"},
                 {"m", "Subjects per Cluster-Period", "m", "m"}};
  p.params = {
      power_param(),
      alpha_param("0.05", "Alpha (Two-Sided)"),
      num("steps", "S (Number of Steps)", "4")
          .grp("Design")
          .sym("Steps")
          .range(2, 1000)
          .whole()
          .sugg({"3", "4", "4 6"})
          .hint("Number of times a new group of clusters switches from control to intervention; there are S + 1 "
                "measurement periods. With one step the effect could not be told apart from the period."),
      num("g", "g (Clusters Switching at Each Step)", "2")
          .grp("Design")
          .sym("Clusters/Step")
          .at_least(1)
          .whole()
          .sugg({"1", "2", "1 2 3"})
          .hint("Number of clusters that switch to the intervention at each step; the trial has S × g clusters.")
          .show(when(not_solving({"g"}))),
      num("m", "m (Subjects per Cluster per Period)", "20")
          .grp("Design")
          .sym("m")
          .at_least(1)
          .whole()
          .sugg({"10", "20", "10 20 50"})
          .hint("Number of subjects measured in each cluster in each period (different subjects in every "
                "period).")
          .show(when(not_solving({"m"}))),
      num("icc", "ICC (Intracluster Correlation)", "0.05")
          .grp(G_EFFECT)
          .sym("ICC")
          .range(0, 1, false, true)
          .sugg({"0.01", "0.05", "0.1"})
          .hint("The share of the outcome's variance that lies between clusters, from 0 to 1; values from 0.01 to "
                "0.1 are common."),
      select("outcome", "Outcome", {{"continuous", "Continuous (difference in means)", {}}, {"binary", "Binary (two proportions)", {}}},
             "continuous")
          .grp(G_EFFECT)
          .hint("Whether the outcome is a measurement compared by its mean, or a yes-or-no event compared by its "
                "proportion."),
      num("delta", "δ (Treatment Effect)", "0.3")
          .grp(G_EFFECT)
          .sym("δ")
          .sugg({"0.2", "0.3"})
          .hint("The difference between the mean outcome under intervention and under control that the trial "
                "should detect.")
          .show(when(is("outcome", {"continuous"}))),
      num("sigma", "σ (Total Standard Deviation)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .hint("Standard deviation of one subject's outcome, including the variation between clusters.")
          .show(when(is("outcome", {"continuous"}))),
      num("p1", "P1 (Proportion under Intervention)", "0.3")
          .grp(G_EFFECT)
          .sym("P1")
          .range(0, 1, true, true)
          .hint("The share of subjects with the outcome while their cluster receives the intervention, between 0 "
                "and 1.")
          .show(when(is("outcome", {"binary"}))),
      num("p2", "P2 (Proportion under Control)", "0.2")
          .grp(G_EFFECT)
          .sym("P2")
          .range(0, 1, true, true)
          .hint("The share of subjects with the outcome while their cluster is still under control, between 0 "
                "and 1.")
          .show(when(is("outcome", {"binary"}))),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const int steps = static_cast<int>(s("steps"));
    const double icc = s("icc"), alpha = s("alpha");
    const SwEffect e = sw_effect(s);
    int g;
    double m;
    if (s.solve == "g") {
      m = s("m");
      solve::NSearch o;
      o.lo = 1;
      o.hi = 1e5;
      g = static_cast<int>(solve::min_n([&](double x) { return sw_power(steps, static_cast<int>(x), m, e, icc, alpha); },
                                        s("power"), o, &ctx));
    } else if (s.solve == "m") {
      g = static_cast<int>(s("g"));
      // power is bounded as m → ∞ (the between-cluster variance remains)
      if (sw_power(steps, g, 1e12, e, icc, alpha) < s("power"))
        throw CalcError("the target power cannot be reached with " + std::to_string(g) +
                        " clusters per step, however large the clusters");
      solve::NSearch o;
      o.lo = 1;
      o.hi = 1e9;
      m = solve::min_n([&](double x) { return sw_power(steps, g, x, e, icc, alpha); }, s("power"), o, &ctx);
    } else {
      g = static_cast<int>(s("g"));
      m = s("m");
    }
    const double I = static_cast<double>(steps) * g, T = steps + 1.0;
    if (s.solve != "power") r.set("target_power", s("power"));
    r.set("power", sw_power(steps, g, m, e, icc, alpha));
    r.set("steps", steps);
    r.set("g", g);
    r.set("clusters", I);
    r.set("periods", T);
    r.set("m", m);
    if (s.is("outcome", "binary"))
      r.notes.push_back("For a binary outcome the subject-level variance σ² is the average of P1(1 − P1) and "
                        "P2(1 − P2).");
    r.set("N", I * T * m);
    r.set("delta", e.delta);
    r.set("icc", icc);
    r.set("se", sw_se(steps, g, m, e, icc));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c;
    if (s.solve != "power") c.push_back(col("target_power", "Target Power", "prob"));
    for (auto x : {col("power", s.solve == "power" ? "Power" : "Actual Power", "prob"), col("steps", "Steps", "int"),
                   col("g", "Clusters/Step", "int"), col("clusters", "Clusters", "int"), col("periods", "Periods", "int"),
                   col("m", "m", "int", "Subjects per cluster per period"),
                   col("N", "N (Total)", "int", "Clusters × periods × m"), col("delta", "δ", "num"),
                   col("icc", "ICC", "num"), col("se", "SE(θ̂)", "num", "Standard error of the treatment effect"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A stepped-wedge design with " + fint(r.num("clusters")) + " clusters (" + fint(r.num("g")) +
           " switching to the intervention at each of " + fint(r.num("steps")) + " steps, " + fint(r.num("periods")) +
           " periods) and " + fint(r.num("m")) + " subjects per cluster and period (" + fint(r.num("N")) +
           " observations in total) has power " + fnum(r.num("power"), 4) + " to detect a treatment effect of " +
           fnum(r.num("delta")) + " with an intracluster correlation of " + fnum(r.num("icc")) +
           ", using a two-sided test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.derivation = sw_derivation;
  p.definitions = {{"Steps", "Number of times a group of clusters switches from control to intervention."},
                   {"Periods", "Measurement periods: one baseline period plus one per step."},
                   {"ICC", "Correlation between two subjects in the same cluster (and period)."}};
  p.notes = {"Cross-sectional design: different subjects are measured in each period. Power uses the normal "
             "approximation of Hussey & Hughes."};
  p.curve_param = "";
  p.example_solve = "power";
  p.example = {{"steps", "3 4 6"}, {"g", "2"}, {"m", "20"}, {"icc", "0.05"}, {"delta", "0.3"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(stepped_wedge);

}  // namespace
}  // namespace ssc::procs
