// Designs for means beyond two parallel groups: one-way ANOVA, 2×2 crossover
// (difference and bioequivalence ratio), cluster-randomized and
// repeated-measures comparisons.
#include "procedures/common.hpp"
#include "procedures/derive_means.hpp"
#include "procedures/integrals.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// =====================================================================================
// One-way ANOVA
// =====================================================================================
struct AnovaDesign {
  std::vector<double> sizes;  // group sizes
  double lambda;              // noncentrality
  double sd_means;            // weighted SD of the group means
  double f;                   // Cohen's f
  int k;
};

std::vector<double> anova_weights(const Scenario& s, int k) {
  if (s.is("alloc", "ratio")) {
    const auto& r = s.v("ratios");
    if (static_cast<int>(r.size()) != k)
      throw CalcError("the allocation pattern must have one entry per group (" + std::to_string(k) + ")");
    for (double x : r)
      if (!(x > 0)) throw CalcError("allocation ratios must be positive");
    return r;
  }
  return std::vector<double>(static_cast<std::size_t>(k), 1.0);
}

int anova_k(const Scenario& s) {
  if (s.is("estype", "f")) return static_cast<int>(s("k"));
  const int k = static_cast<int>(s.v("means").size());
  if (k < 2) throw CalcError("enter at least two group means");
  return k;
}

// the group sizes for n per weight-1 group: n × each weight, rounded up
std::vector<double> anova_sizes(const Scenario& s, int k, double n) {
  std::vector<double> sizes;
  for (double wi : anova_weights(s, k)) sizes.push_back(solve::round_up(n * wi));
  return sizes;
}

// effect-size multiplier `scale` stretches the means about their mean (for solving the effect)
AnovaDesign anova_design(const Scenario& s, double n, double scale_f = NAN) {
  const int k = anova_k(s);
  AnovaDesign d;
  d.k = k;
  d.sizes = anova_sizes(s, k, n);
  double N = 0;
  for (double x : d.sizes) N += x;
  if (s.is("estype", "f")) {
    d.f = std::isfinite(scale_f) ? scale_f : s("f");
    d.lambda = N * d.f * d.f;
    d.sd_means = NAN;
    return d;
  }
  const auto& mu = s.v("means");
  const double sigma = s("sigma");
  // from the differences from the first mean: a mean of the means themselves rounds when they are large next to
  // their spread (or to σ), and equal means of 0.1 with σ = 1e-30, or of 1e30 with σ = 1, showed an effect of
  // f = 1e13 that is not there
  double mbar = 0;
  for (int i = 0; i < k; ++i) mbar += d.sizes[i] * (mu[i] - mu[0]);
  mbar /= N;
  double ss = 0;
  for (int i = 0; i < k; ++i) ss += d.sizes[i] * (mu[i] - mu[0] - mbar) * (mu[i] - mu[0] - mbar);
  double sm = std::sqrt(ss / N);
  if (std::isfinite(scale_f)) {  // rescale to the requested f
    sm = scale_f * sigma;
  }
  d.sd_means = sm;
  d.f = sm / sigma;
  d.lambda = N * d.f * d.f;
  return d;
}

double anova_power(const AnovaDesign& d, double alpha) {
  double N = 0;
  for (double x : d.sizes) N += x;
  const double df1 = d.k - 1, df2 = N - d.k;
  if (!(df2 > 0)) return NAN;
  const double fc = stats::qf(alpha, df1, df2, false);
  return stats::pnf(fc, df1, df2, d.lambda, false);
}

Procedure anova() {
  Procedure p;
  p.id = "means.anova";
  p.name = "One-Way Analysis of Variance (ANOVA) F-Test";
  p.category = "Means";
  p.subcategory = "Many Means (ANOVA)";
  p.order = 1;
  p.summary =
      "Sample size and power for the one-way ANOVA F-test comparing k group means, given the hypothesized means "
      "(or Cohen's effect size f), with equal or unequal group sizes. Power uses the noncentral F distribution.";
  p.keywords = {"anova", "f-test", "one-way", "k groups", "several means", "cohen f", "several groups", "multiple groups", "between-group", "effect size f"};
  p.solve_for = {{"n", "Sample Size", "n", "n (Per-Group Size)"},
                 {"power", "Power", "power", "Power"},
                 {"f", "Effect Size (f)", "f", "Detectable f"}};
  p.params = {
      power_param(),
      alpha_param(),
      select("alloc", "Group Allocation",
             {{"equal", "Equal group sizes", {}}, {"ratio", "Allocation pattern (relative sizes)", {}}}, "equal")
          .grp(G_N)
          .hint("Give every group the same number of subjects, or size the groups in proportion to a pattern of "
                "relative weights."),
      n_param("n", "n (Size of Each Group / of a Weight-1 Group)", "20", 2)
          .sym("n")
          .hint("Number of subjects in each group. With an allocation pattern, the size of a group with weight 1 "
                "(a group with weight 2 gets twice as many)."),
      vec("ratios", "Allocation Pattern", "1 1 2")
          .grp(G_N)
          .sym("Pattern")
          .above(0)
          .hint("Relative group sizes; group i gets n × weight_i subjects (rounded up).")
          .show(when(is("alloc", {"ratio"}))),
      dropout_param(),
      select("estype", "Effect Size Input",
             {{"means", "Group means and standard deviation", {}}, {"f", "Cohen's f and number of groups", {}}}, "means")
          .grp(G_EFFECT)
          .hint("Describe the effect by the group means you expect and the within-group standard deviation, or "
                "directly by Cohen's effect size f and the number of groups."),
      vec("means", "Group Means μ1, μ2, …", "0 0.5 1")
          .grp(G_EFFECT)
          .sym("Means")
          .sugg({"0 0.5 1", "0 0 1", "10 12 14 16"})
          .hint("One mean per group. Separate alternative sets of means with ';'.")
          .show(when({is("estype", {"means"})})),
      num("sigma", "σ (Within-Group Standard Deviation)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .sugg({"1", "1 1.5 2"})
          .hint("The standard deviation of the outcome among subjects in the same group, assumed equal in all "
                "groups; take it from earlier studies or a pilot.")
          .show(when(is("estype", {"means"}))),
      num("k", "k (Number of Groups)", "3")
          .grp(G_EFFECT)
          .sym("k")
          .range(2, 1000)
          .whole()
          .sugg({"3", "2 3 4 5"})
          .hint("How many groups are compared (at least 2).")
          .show(when(is("estype", {"f"}))),
      num("f", "f (Cohen's Effect Size σm/σ)", "0.25")
          .grp(G_EFFECT)
          .sym("f")
          .at_least(0)
          .sugg({"0.1", "0.25", "0.4", "0.1 0.25 0.4"})
          .hint("Cohen's f: the SD of the group means divided by the within-group SD. By Cohen's conventions 0.10 "
                "is small, 0.25 medium and 0.40 large.")
          .show(when({is("estype", {"f"}), not_solving({"f"})})),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double alpha = s("alpha");
    // the F-test's error has N − k degrees of freedom: none when tiny allocation ratios give every group 1 subject
    const auto w = anova_weights(s, anova_k(s));
    auto spare = [&](double n) {
      for (double wi : w)
        if (solve::round_up(n * wi) >= 2) return true;
      return false;
    };
    const std::string single = "every group has 1 subject (n × its allocation ratio, rounded up)";
    const std::string no_df = ", which leaves the F-test no degrees of freedom for the error; at least one group "
                              "needs 2";
    AnovaDesign d;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = 2;
      if (!spare(o.hi)) throw CalcError(single + " for any n up to " + text::format_g(o.hi, 10) + no_df);
      d = anova_design(s, o.lo);
      if (!(d.f > 0))
        throw CalcError(s.is("estype", "f") ? "f = 0: there is no effect to detect"
                                            : "the group means are all equal; there is no effect to detect");
      o.guess = (z_factor(alpha, s("power"), Alt::Two) + 2.0 * (d.k - 1)) / (d.k * d.f * d.f) + 1;
      // groups that tiny ratios keep at 1 subject up to the largest n hold the power down however large n grows
      const std::string held = solve::single_groups(w, o.hi, "allocation ratio", "allocation ratios");
      if (!held.empty() && anova_power(anova_design(s, o.hi), alpha) < s("power"))
        throw CalcError("the target power cannot be reached with a sample size up to " + text::format_g(o.hi, 10) +
                        ", where " + held);
      const double n = solve::min_n([&](double x) { return anova_power(anova_design(s, x), alpha); }, s("power"), o, &ctx);
      d = anova_design(s, n);
      r.set("n", n);
    } else {
      const double n = s("n");
      if (!spare(n)) throw CalcError(single + no_df);
      if (s.solve == "power") {
        d = anova_design(s, n);
      } else {
        const double f = solve::solve_increasing([&](double ff) { return anova_power(anova_design(s, n, ff), alpha); },
                                                 s("power"), 1e-9, 0.25, 1e30, 1e-12, "the effect size f");
        d = anova_design(s, n, f);
      }
      r.set("n", n);
    }
    double N = 0;
    for (double x : d.sizes) N += x;
    put_power(s, r, anova_power(d, alpha));
    r.set("N", N);
    r.set("k", d.k);
    if (std::isfinite(d.sd_means)) r.set("sm", d.sd_means);
    if (s.is("estype", "means")) r.set("sigma", s("sigma"));
    r.set("f", d.f);
    r.set("lambda", d.lambda);
    r.set("df1", d.k - 1);
    r.set("df2", N - d.k);
    r.set("alpha", alpha);
    if (s.is("alloc", "ratio")) {
      std::string sz;
      for (double x : d.sizes) sz += (sz.empty() ? "" : ", ") + fint(x);
      r.set("sizes", sz);
    }
  };
  p.curve_inputs = [](Scenario& sc, const Row& r) {  // the solved f, not the means it was solved next to
    sc.sel["estype"] = "f";
    sc.num["f"] = r.num("f");
    sc.num["k"] = r.num("k");
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    c.push_back(col("n", "n", "int", "Size of each group (or of a group with weight 1)"));
    if (s.is("alloc", "ratio")) c.push_back(col("sizes", "Group Sizes", "text"));
    for (auto x : {col("N", "N", "int", "Total sample size"), col("k", "k", "int", "Number of groups")})
      c.push_back(x);
    if (s.is("estype", "means")) {
      c.push_back(col("sm", "σm", "num", "Standard deviation of the group means"));
      c.push_back(col("sigma", "σ", "num", "Within-group standard deviation"));
    }
    for (auto x : {col("f", "f = σm/σ", "num", "Cohen's effect size"), col("lambda", "λ", "num", "Noncentrality"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A total of " + fint(r.num("N")) + " subjects in " + fint(r.num("k")) + " groups (n = " + fint(r.num("n")) +
           (s.is("alloc", "ratio") ? " for a weight-1 group" : " per group") + ") " + achieves_singular(s, r) +
           " to detect differences among the group means with effect size f = " + fnum(r.num("f")) +
           " using the one-way ANOVA F-test at significance level " + fnum(r.num("alpha")) + ".";
  };
  // Solving for n, the design at the n found and one below; for power, at the n given; for f, the target power gives
  // the noncentrality λ = N f², and λ the effect.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const double alpha = s("alpha"), n = r.num("n"), N = r.num("N"), f = r.num("f"), lambda = r.num("lambda");
    const double df1 = r.num("df1"), df2 = r.num("df2"), target = s.get("power", NAN);
    const int k = anova_k(s);
    const bool ratio = s.is("alloc", "ratio"), by_f = s.is("estype", "f"), solved = s.solve == "n";
    const std::vector<double> sizes = anova_sizes(s, k, n);
    const std::string ns = size(n, solved), Ns = dvi(N);
    const std::string per = ratio ? " for a group of weight 1" : " per group";
    Derivation out;
    if (solved)
      out.find = find_size(ratio ? "n (the size of a group of weight 1)" : "n per group", target);
    else if (s.solve == "power")
      out.find = find_power("with n = " + gvi(n) + per);
    else
      out.find = find_effect("effect size f", target, "with n = " + gvi(n) + per);
    if (ratio) {
      const std::vector<double> w = anova_weights(s, k);
      std::vector<std::string> wt, sz;
      for (int i = 0; i < k; ++i) {
        wt.push_back(gv(w[i]));
        sz.push_back(dvi(sizes[i]));
      }
      out.steps.push_back({"Group sizes", "nᵢ = ⌈n × wᵢ⌉ for the allocation weights wᵢ, N = Σ nᵢ",
                           "⌈" + ns + " × (" + joined(wt, ", ") + ")⌉ = (" + joined(sz, ", ") + "), N = " + Ns});
    } else {
      out.steps.push_back({"Total sample size", "N = k × n", gvi(k) + " × " + ns + " = " + Ns});
    }
    if (s.solve == "f") {  // the noncentrality the target power needs, and the f that gives it
      out.steps.push_back(f_critical(alpha, df1, df2, "ν1 = k − 1 and ν2 = N − k"));
      out.steps.push_back(means_f_ncp_needed(df1, df2, lambda, alpha, target));
      out.steps.push_back({"Effect size detected", "f = √(λ / N)", "√(" + dv(lambda) + " / " + Ns + ") = " + dv(f)});
      out.answer = "f = " + dv(f);
      if (!by_f) {
        out.steps.push_back({"SD of the group means detected", "σm = f σ",
                             dv(f) + " × " + gv(s("sigma")) + " = " + dv(r.num("sm"))});
        out.answer += " (σm = " + dv(r.num("sm")) + ")";
      }
      return out;
    }
    if (!by_f) {
      const auto& mu = s.v("means");
      // the mean of the means as the run takes it, from their differences from the first
      double mbar = 0;
      for (int i = 0; i < k; ++i) mbar += sizes[i] * (mu[i] - mu[0]);
      mbar = mu[0] + mbar / N;
      std::vector<std::string> sum, dev;
      for (int i = 0; i < k; ++i) {
        const std::string times = ratio ? dvi(sizes[i]) + " × " : "";
        sum.push_back(times + gvs(mu[i]));
        dev.push_back(times + "(" + gv(mu[i]) + " − " + dvs(mbar) + ")²");
      }
      const std::string over = ratio ? Ns : gvi(k);
      out.steps.push_back({"Mean of the group means", ratio ? "μ̄ = Σ nᵢμᵢ / N" : "μ̄ = Σ μᵢ / k",
                           "(" + joined(sum, " + ") + ") / " + over + " = " + dv(mbar)});
      out.steps.push_back({"SD of the group means", ratio ? "σm = √(Σ nᵢ(μᵢ − μ̄)² / N)" : "σm = √(Σ (μᵢ − μ̄)² / k)",
                           "√((" + joined(dev, " + ") + ") / " + over + ") = " + dv(r.num("sm"))});
      out.steps.push_back(
          {"Effect size", "f = σm / σ", dv(r.num("sm")) + " / " + gv(s("sigma")) + " = " + dv(f)});
    }
    out.steps.push_back({"Noncentrality", "λ = N f²", Ns + " × " + (by_f ? gv(f) : dv(f)) + "² = " + dv(lambda)});
    out.steps.push_back(f_critical(alpha, df1, df2, "ν1 = k − 1 and ν2 = N − k"));
    out.steps.push_back(f_power_step(df1, df2, lambda, alpha));
    const double power = r.num("power");
    if (solved) {
      out.steps.push_back(
          n_search([&](double x) { return anova_power(anova_design(s, x), alpha); }, n, target, 2, "n"));
      out.answer = "n = " + dvi(n) + per + " (N = " + Ns + "), with power " + dvp(power);
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
  };
  p.definitions = {{"n, N", "Group size and total sample size."},
                   {"k", "Number of groups."},
                   {"σm", "Standard deviation of the group means: √(Σ nᵢ(μᵢ − μ̄)²/N)."},
                   {"f", "Cohen's effect size σm/σ."},
                   {"λ", "Noncentrality parameter N f² of the F distribution with (k − 1, N − k) df."}};
  p.curve_param = "n";
  p.curve_label = "n (Per-Group Size)";
  p.dropout_groups = {"N"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"means", "0 0.5 1; 0 0 1"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(anova);

// =====================================================================================
// 2×2 crossover (difference)
// =====================================================================================
double crossover_sigma_w(const Scenario& s) {
  const double v = s("sigma");
  return s.is("sdtype", "diff") ? v / std::sqrt(2.0) : v;
}

// The standard error of the treatment difference with n subjects, σw √(2/n), and the power of its test.
double crossover_se(double sw, double n) { return sw * std::sqrt(2 / n); }

double crossover_power(const Scenario& s, double n, double d) {
  const double ncp = d / crossover_se(crossover_sigma_w(s), n), alpha = s("alpha");
  return s.is("test", "z") ? z_power(ncp, alpha, alt(s)) : t_power(ncp, n - 2, alpha, alt(s));
}

Procedure crossover_diff() {
  Procedure p;
  p.id = "means.crossover";
  p.name = "Two-Treatment 2×2 Cross-Over Design (Difference of Means)";
  p.category = "Means";
  p.subcategory = "Cross-Over Designs";
  p.order = 1;
  p.summary =
      "Sample size and power for comparing two treatments in an AB/BA cross-over trial with a t-test on the period "
      "differences, given the within-subject standard deviation. Sequences are balanced.";
  p.keywords = {"crossover", "cross-over", "2x2", "ab/ba", "within-subject", "period", "two-period", "two-sequence", "period differences"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"},
                 {"power", "Power", "power", "Power"},
                 {"delta", "Effect Size (δ)", "delta", "Detectable δ"}};
  p.params = {
      alternative("Two-sided (H1: δ ≠ 0)", "One-sided (H1: δ > 0)", "One-sided (H1: δ < 0)"),
      select("test", "Test Statistic", {{"t", "T-test", {}}, {"z", "Z-test (σ known)", {}}}, "t")
          .grp(G_TEST)
          .hint("The t-test uses the sample standard deviation; the z-test assumes σ is known."),
      power_param(),
      alpha_param(),
      n_param("n", "N (Total Subjects, Both Sequences)", "24", 4)
          .sym("N")
          .hint("Total number of subjects; half receive the treatments in the order AB and half in the order BA, "
                "so N must be even."),
      dropout_param(),
      num("delta", "δ (Treatment Difference)", "0.5")
          .grp(G_EFFECT)
          .sym("δ")
          .sugg({"0.5", "0.25 0.5 0.75"})
          .hint("The true difference between the two treatment means that the study should detect, in the "
                "units of the outcome.")
          .show(when(not_solving({"delta"}))),
      select("sdtype", "Variability Input",
             {{"within", "σw: within-subject SD (residual)", {}},
              {"diff", "σd: SD of the period differences Y1 − Y2", {}}},
             "within")
          .grp(G_EFFECT)
          .hint("Which standard deviation you enter as the SD value: the within-subject SD σw, or the SD σd of a "
                "subject's difference between the two periods (σd = √2 × σw)."),
      num("sigma", "SD Value", "1")
          .grp(G_EFFECT)
          .sym("SD")
          .above(0)
          .sugg({"1", "0.5 1"})
          .hint("The standard deviation of the kind chosen above, in the units of the outcome; take it from an "
                "earlier cross-over study or a pilot."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double sw = crossover_sigma_w(s), alpha = s("alpha");
    const Alt a = alt(s);
    auto P = [&](double n, double d) { return crossover_power(s, n, d); };
    double n, d;
    if (s.solve == "n") {
      d = s("delta");
      if (d == 0) throw CalcError("δ = 0: there is no treatment difference to detect");
      if ((a == Alt::Upper && d < 0) || (a == Alt::Lower && d > 0))
        throw CalcError("δ has the sign opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = 4;
      o.step = 2;
      o.guess = z_factor(alpha, s("power"), a) * 2 * sw * sw / (d * d) + 2;
      n = solve::min_n([&](double x) { return P(x, d); }, s("power"), o, &ctx);
    } else {
      n = s("n");
      if (std::fmod(n, 2) != 0) throw CalcError("N must be even for balanced sequences");
      if (s.solve == "power") {
        d = s("delta");
      } else {
        const double sign = a == Alt::Lower ? -1 : 1;
        d = sign * solve::solve_increasing([&](double dd) { return P(n, sign * dd); }, s("power"), 0, sw,
                                           std::fmax(1e9 * sw, solve::BOUND_SPAN));
      }
    }
    r.set("n", n);
    r.set("n_seq", n / 2);
    put_power(s, r, P(n, d));
    r.set("delta", d);
    r.set("sw", sw);
    r.set("sd", sw * std::sqrt(2.0));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int", "Total subjects"), col("n_seq", "Per Sequence", "int"),
                   col("delta", "δ", "num", "Treatment difference"), col("sw", "σw", "num", "Within-subject SD"),
                   col("sd", "σd", "num", "SD of period differences"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A 2×2 cross-over design with " + fint(r.num("n")) + " subjects (" + fint(r.num("n_seq")) +
           " per sequence) " + achieves_singular(s, r) + " to detect a treatment difference of " +
           fnum(r.num("delta")) + " when the within-subject standard deviation is " + fnum(r.num("sw")) +
           ", using a " + alt_words(alt(s)) + " " + (s.is("test", "z") ? "z-test" : "t-test") + " at α = " +
           fnum(r.num("alpha")) + ".";
  };
  // The test of the period differences has noncentrality λ = δ / SE, SE = σw √(2/N): solving for N, at the N found
  // and two below; solving for δ, the target power gives λ, and λ the difference.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const Alt a = alt(s);
    const bool z = s.is("test", "z"), solved = s.solve == "n", from_diff = s.is("sdtype", "diff");
    const double n = r.num("n"), d = r.num("delta"), sw = r.num("sw"), alpha = s("alpha"), target = s.get("power", NAN);
    const double se = crossover_se(sw, n), ncp = d / se;
    const std::string N = size(n, solved), sws = from_diff ? dv(sw) : gv(sw);
    Derivation out;
    if (from_diff)
      out.steps.push_back({"Within-subject SD", "σw = σd / √2", gv(s("sigma")) + " / √2 = " + dv(sw)});
    const auto se_step = [&] {
      out.steps.push_back({"Standard error", "SE = σw √(2/N)", sws + " × √(2/" + N + ") = " + dv(se)});
    };
    const auto critical = [&] {
      out.steps.push_back(z ? z_critical(alpha, a) : t_critical(alpha, n - 2, a, "N − 2"));
    };
    if (s.solve == "delta") {  // the target power gives λ, and λ the difference
      out.find = find_effect("treatment difference δ", target, "with N = " + gvi(n));
      se_step();
      critical();
      out.steps.push_back(z ? z_ncp_needed(ncp, alpha, a, target) : t_ncp_needed(ncp, n - 2, alpha, a, target));
      out.steps.push_back({"Difference detected", "δ = λ × SE", dvs(ncp) + " × " + dv(se) + " = " + dv(d)});
      out.answer = "δ = " + dv(d);
      return out;
    }
    if (solved) {
      out.find = find_size("even N (N/2 in each sequence)", target);
      // the closed form of a z-test, beside the exact answer below
      const double za = stats::qnorm(tail_alpha(alpha, a), false), zb = stats::qnorm(target);
      out.steps.push_back({"Normal approximation",
                           std::string("N ≈ 2 (z(") + level_symbol(a) +
                               ") + z(1 − β))² σw² / δ², for comparison; the steps below give the answer",
                           "2 × (" + dv(za) + " + " + dv(zb) + ")² × " + sws + "² / " + gvs(d) + "² = " +
                               dv(2 * (za + zb) * (za + zb) * sw * sw / (d * d))});
    } else {
      out.find = find_power("with N = " + gvi(n));
    }
    se_step();
    out.steps.push_back({"Noncentrality", "λ = δ / SE", gvs(d) + " / " + dv(se) + " = " + dv(ncp)});
    critical();
    out.steps.push_back(z ? z_power_step(ncp, alpha, a) : t_power_step(ncp, n - 2, alpha, a));
    const double power = crossover_power(s, n, d);
    if (solved) {
      out.steps.push_back(
          even_n_search([&](double x) { return crossover_power(s, x, d); }, n, target, 4, "N/2 in each sequence"));
      out.answer = "N = " + dvi(n) + " (" + dvi(n / 2) + " per sequence), with power " + dvp(power);
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
  };
  p.definitions = {{"N", "Total number of subjects; half receive AB and half BA."},
                   {"σw", "Within-subject (residual) standard deviation."},
                   {"σd", "SD of a subject's period difference, √2 σw."},
                   {"δ", "Difference between the treatment means."}};
  p.curve_param = "n";
  p.curve_label = "N (Total Subjects)";
  p.curve_min = 4;
  p.dropout_groups = {"n"};
  p.dropout_multiple = 2;  // two sequences of equal size
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"delta", "0.5 1"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(crossover_diff);

// =====================================================================================
// Bioequivalence of a ratio of means (log scale, TOST)
// =====================================================================================
// The SD on the log scale, √ln(1 + CV²): √ln(1 + CV²) = CV(1 − CV²/4 + …) is CV itself to double precision below
// 1e-8 (computed, CV² underflows to 0 below about 1e-162, and so did the log-scale SD).
double be_log_sd(double cv) { return cv < 1e-8 ? cv : std::sqrt(std::log1p(cv * cv)); }

// The standard error of the log ratio with n subjects: σ √(2/n) in a 2×2 cross-over, σ √(4/n) in parallel groups.
double be_se(const Scenario& s, double sw, double n) {
  return sw * std::sqrt((s.is("design", "2x2") ? 2.0 : 4.0) / n);
}

double be_power(const Scenario& s, double sw, double n) {
  return tost_power(std::log(s("theta1")), std::log(s("theta2")), std::log(s("theta0")), be_se(s, sw, n), n - 2,
                    s("alpha"));
}

Procedure bioequivalence() {
  Procedure p;
  p.id = "means.be.ratio";
  p.name = "Bioequivalence: Ratio of Two Means (TOST, Log Scale)";
  p.category = "Means";
  p.subcategory = "Cross-Over Designs";
  p.order = 2;
  p.summary =
      "Average bioequivalence of a test and reference formulation with two one-sided t-tests on log-transformed "
      "data (e.g. limits 0.80–1.25), for a 2×2 cross-over or a parallel design. Exact power (Owen's Q) by numerical "
      "integration.";
  p.keywords = {"bioequivalence", "be", "tost", "ratio", "log-normal", "cv", "crossover", "pharmacokinetics", "auc", "cmax", "geometric mean ratio", "gmr", "average bioequivalence", "generic drug", "80-125"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      select("design", "Design",
             {{"2x2", "2×2 cross-over (CV = within-subject)", {}}, {"parallel", "Parallel groups (CV = total)", {}}},
             "2x2")
          .grp(G_TEST)
          .hint("In a 2×2 cross-over every subject receives both formulations and the CV is the within-subject "
                "CV; in a parallel design each subject receives one, and the CV is the total CV."),
      power_param(),
      alpha_param("0.05", "Alpha (Each One-Sided Test)").range(0, 0.5, true, true),
      n_param("n", "N (Total Subjects)", "24", 4)
          .sym("N")
          .hint("Total number of subjects, split evenly between the two sequences (or groups), so N must be "
                "even."),
      dropout_param(),
      num("theta0", "θ0 (True Ratio Test/Reference)", "0.95")
          .grp(G_EFFECT)
          .sym("θ0")
          .above(0)
          .sugg({"0.95", "0.9 0.95 1", "1"})
          .hint("The ratio of the test to the reference mean that you expect to be true. It must lie between "
                "the limits; values from 0.95 to 1.05 are typical."),
      num("theta1", "Lower Limit θL", "0.8")
          .grp(G_EFFECT)
          .sym("θL")
          .range(0, 1, true, true)
          .sugg({"0.8", "0.9"})
          .hint("Lower bioequivalence limit for the ratio of means; 0.80 is the usual regulatory value."),
      num("theta2", "Upper Limit θU", "1.25")
          .grp(G_EFFECT)
          .sym("θU")
          .above(1)
          .sugg({"1.25", "1.111"})
          .hint("Upper bioequivalence limit for the ratio of means, usually 1/θL (1.25 when θL = 0.80, 1.111 when "
                "θL = 0.90)."),
      num("cv", "CV (Coefficient of Variation)", "0.3")
          .grp(G_EFFECT)
          .sym("CV")
          .above(0)
          .sugg({"0.2", "0.3", "0.2 0.25 0.3 0.35"})
          .hint("The SD divided by the mean on the original scale (0.3 for 30%): within-subject for a cross-over, "
                "total for a parallel design. The log-scale SD is √ln(1 + CV²)."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double alpha = s("alpha"), cv = s("cv"), th0 = s("theta0"), th1 = s("theta1"), th2 = s("theta2");
    // on a limit itself the power is the size of the tests, which no sample size raises above alpha
    if (s.solve == "n" ? !(th0 > th1 && th0 < th2) : !(th0 >= th1 && th0 <= th2))
      throw CalcError(s.solve == "n" ? "the true ratio must lie strictly between the limits"
                                     : "the true ratio cannot lie outside the limits");
    const double sw = be_log_sd(cv);
    const bool cross = s.is("design", "2x2");
    auto P = [&](double n) { return be_power(s, sw, n); };
    double n;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = 4;
      o.step = 2;
      const double gap = std::min(std::log(th0 / th1), std::log(th2 / th0));
      const double zf = stats::qnorm(alpha, false) + stats::qnorm(s("power"));
      o.guess = (cross ? 2.0 : 4.0) * zf * zf * sw * sw / (gap * gap) + 2;
      n = solve::min_n(P, s("power"), o, &ctx);
    } else {
      n = s("n");
      if (std::fmod(n, 2) != 0) throw CalcError("N must be even (balanced sequences or groups)");
    }
    r.set("n", n);
    put_power(s, r, P(n));
    r.set("theta0", th0);
    r.set("theta1", th1);
    r.set("theta2", th2);
    r.set("cv", cv);
    r.set("sw", sw);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int", "Total subjects"), col("theta0", "θ0", "num", "True ratio"),
                   col("theta1", "θL", "num"), col("theta2", "θU", "num"), col("cv", "CV", "num"),
                   col("sw", "σ (log)", "num", "SD on the log scale"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return std::string(s.is("design", "2x2") ? "A 2×2 cross-over study" : "A parallel-group study") + " with " +
           fint(r.num("n")) + " subjects " + achieves_singular(s, r) +
           " to conclude average bioequivalence with two one-sided t-tests at α = " + fnum(r.num("alpha")) +
           " each, when the true ratio of means is " + fnum(r.num("theta0")) + ", the equivalence limits are " +
           fnum(r.num("theta1")) + " and " + fnum(r.num("theta2")) + ", and the CV is " + fnum(r.num("cv")) + ".";
  };
  // The two one-sided tests on the log scale: solving for N, at the N found and two below; solving for power, at
  // the N given.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const double n = r.num("n"), sw = r.num("sw"), se = be_se(s, sw, n), alpha = s("alpha");
    const double target = s.get("power", NAN);
    const double l1 = std::log(r.num("theta1")), l2 = std::log(r.num("theta2")), l0 = std::log(r.num("theta0"));
    const bool cross = s.is("design", "2x2"), solved = s.solve == "n";
    const std::string halves = cross ? "N/2 in each sequence" : "N/2 in each group", N = size(n, solved);
    const TostSymbols sym{"ln θL", "ln θU", "ln θ0"};
    const TostShown shown{dv(l1), dv(l2), dv(l0)};
    Derivation out;
    out.find = solved ? find_size("even N (" + halves + ")", target) : find_power("with N = " + gvi(n));
    out.steps.push_back(
        {"SD on the log scale", "σ = √ln(1 + CV²)", "√ln(1 + " + gv(r.num("cv")) + "²) = " + dv(sw)});
    out.steps.push_back({"Limits and true ratio on the log scale", "ln θL, ln θU, ln θ0",
                         "ln " + gv(r.num("theta1")) + " = " + dv(l1) + ", ln " + gv(r.num("theta2")) + " = " +
                             dv(l2) + ", ln " + gv(r.num("theta0")) + " = " + dv(l0)});
    if (solved)
      out.steps.push_back(means_tost_approx(l1, l2, l0, alpha, target, sym, shown, "N", cross ? "2σ²" : "4σ²",
                                            (cross ? "2 × " : "4 × ") + dv(sw) + "²", (cross ? 2 : 4) * sw * sw));
    out.steps.push_back({"Standard error", cross ? "SE = σ √(2/N)" : "SE = σ √(4/N)",
                         dv(sw) + " × √(" + (cross ? "2" : "4") + "/" + N + ") = " + dv(se)});
    for (auto& st : tost_steps(l1, l2, l0, se, n - 2, alpha, "N − 2", sym, shown)) out.steps.push_back(st);
    const double power = be_power(s, sw, n);
    if (solved) {
      out.steps.push_back(even_n_search([&](double x) { return be_power(s, sw, x); }, n, target, 4, halves));
      out.answer = "N = " + dvi(n) + " (" + dvi(n / 2) + (cross ? " per sequence" : " per group") +
                   "), with power " + dvp(power);
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
  };
  p.definitions = {{"θ0", "True ratio of the test to the reference mean (geometric means)."},
                   {"θL, θU", "Bioequivalence acceptance limits for the ratio."},
                   {"CV", "Coefficient of variation; within-subject for cross-over designs."},
                   {"σ (log)", "√ln(1 + CV²), the standard deviation of log-transformed data."}};
  p.curve_param = "n";
  p.curve_label = "N (Total Subjects)";
  p.curve_min = 4;
  p.dropout_groups = {"n"};
  p.dropout_multiple = 2;  // two sequences of equal size
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"theta0", "0.95"}, {"cv", "0.2 0.25 0.3"}};
  return p;
}
SSC_REGISTER(bioequivalence);

// =====================================================================================
// Cluster-randomized two means
// =====================================================================================
// k1 and k2 = ⌈r × k1⌉ clusters of m subjects, whose variance the design effect 1 + (m − 1)ρ inflates.
struct ClusterMeans {
  double alpha, m, icc, sigma, kr, de;
  Alt a;
  bool z;
  explicit ClusterMeans(const Scenario& s)
      : alpha(s("alpha")), m(s("m")), icc(s("icc")), sigma(s("sigma")), kr(s("kr")), de(1 + (m - 1) * icc),
        a(alt(s)), z(s.is("test", "z")) {}
  double k2(double k1) const { return solve::round_up(kr * k1); }
  double se(double k1, double k2) const { return sigma * std::sqrt(de / m * (1 / k1 + 1 / k2)); }
  double power(double k1, double diff) const {
    const double clusters2 = k2(k1), error = se(k1, clusters2);
    return z ? z_power(diff / error, alpha, a) : t_power(diff / error, k1 + clusters2 - 2, alpha, a);
  }
};

Procedure cluster_means() {
  Procedure p;
  p.id = "means.two.cluster";
  p.name = "Two Means in a Cluster-Randomized Design";
  p.category = "Cluster-Randomized";
  p.subcategory = "Means";
  p.order = 1;
  p.summary =
      "Number of clusters per arm (or power) for comparing two means when whole clusters are randomized. The "
      "variance is inflated by the design effect 1 + (m − 1)ρ; the test uses k1 + k2 − 2 degrees of freedom.";
  p.keywords = {"cluster", "group randomized", "icc", "intraclass", "design effect", "community trial", "cluster trial", "continuous outcome"};
  p.solve_for = {{"k", "Number of Clusters", "k1", "Clusters in Group 1"},
                 {"power", "Power", "power", "Power"},
                 {"mu1", "Effect Size (μ1)", "mu1", "Detectable μ1"}};
  p.params = {
      alternative("Two-sided (H1: μ1 ≠ μ2)", "One-sided (H1: μ1 > μ2)", "One-sided (H1: μ1 < μ2)"),
      select("test", "Test Statistic", {{"t", "T-test on cluster means (df = k1 + k2 − 2)", {}}, {"z", "Z-test", {}}},
             "t")
          .grp(G_TEST)
          .hint("The t-test allows for the few degrees of freedom that a small number of clusters gives; the "
                "z-test is a large-sample approximation that overstates the power when clusters are few."),
      power_param(),
      alpha_param(),
      num("k1", "k1 (Clusters in Group 1)", "10")
          .grp(G_N)
          .sym("k1")
          .at_least(2)
          .whole()
          .sugg({"10", "5 to 30 by 5"})
          .hint("Number of clusters (such as clinics, schools or villages) randomized to group 1.")
          .show(when(not_solving({"k"}))),
      num("kr", "Cluster Ratio k2/k1", "1")
          .grp(G_N)
          .sym("k2/k1")
          .above(0)
          .sugg({"1", "2"})
          .hint("Clusters in group 2 per cluster in group 1; k2 = ratio × k1, rounded up. Use 1 for equal "
                "numbers of clusters."),
      num("m", "m (Subjects per Cluster)", "20")
          .grp(G_N)
          .sym("m")
          .at_least(1)
          .whole()
          .sugg({"20", "10 20 50"})
          .hint("Number of subjects measured in each cluster, assumed the same for every cluster."),
      num("icc", "ρ (Intracluster Correlation)", "0.05")
          .grp(G_EFFECT)
          .sym("ICC")
          .range(0, 1, false, true)
          .sugg({"0.01", "0.05", "0.01 0.05 0.1"})
          .hint("The share of the outcome's variance that lies between clusters, from 0 to 1; values from 0.01 to "
                "0.1 are common. Take it from earlier studies of similar clusters."),
      num("mu1", "μ1 (Mean of Group 1)", "0.5")
          .grp(G_EFFECT)
          .sym("μ1")
          .sugg({"0.5", "0.25 0.5"})
          .hint("Mean outcome in group 1 (e.g. the intervention) that you expect under the alternative.")
          .show(when(not_solving({"mu1"}))),
      num("mu2", "μ2 (Mean of Group 2)", "0")
          .grp(G_EFFECT)
          .sym("μ2")
          .sugg({"0"})
          .hint("Mean outcome in group 2 (e.g. the control); the difference μ1 − μ2 is the effect to detect."),
      num("sigma", "σ (Subject-Level SD)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .sugg({"1"})
          .hint("Standard deviation of one subject's outcome within a group, including the variation between "
                "clusters (the total SD)."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const ClusterMeans c(s);
    const double alpha = c.alpha, m = c.m, icc = c.icc, sigma = c.sigma, kr = c.kr, de = c.de;
    const Alt a = c.a;
    auto k2_of = [&](double k1) { return c.k2(k1); };
    auto P = [&](double k1, double diff) { return c.power(k1, diff); };
    double k1, diff;
    if (s.solve == "k") {
      diff = s("mu1") - s("mu2");
      if (diff == 0) throw CalcError("μ1 equals μ2");
      if ((a == Alt::Upper && diff < 0) || (a == Alt::Lower && diff > 0))
        throw CalcError("the difference μ1 − μ2 has the sign opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = solve::first_whole([&](double x) { return k2_of(x) >= 2; }, 2, o.hi,
                                "the cluster ratio is too small: group 2 would have fewer than 2 clusters");
      o.guess = z_factor(alpha, s("power"), a) * sigma * sigma * de / m * (1 + 1 / kr) / (diff * diff) + 1;
      k1 = solve::min_n([&](double x) { return P(x, diff); }, s("power"), o, &ctx);
    } else {
      k1 = s("k1");
      if (k2_of(k1) < 2)
        throw CalcError("group 2 would have " + fint(k2_of(k1)) + (k2_of(k1) == 1 ? " cluster" : " clusters") +
                        " (k1 × the cluster ratio, rounded up); each group needs at least 2 clusters");
      if (s.solve == "power") {
        diff = s("mu1") - s("mu2");
      } else {
        const double sign = a == Alt::Lower ? -1 : 1;
        diff = sign * solve::solve_increasing([&](double d) { return P(k1, sign * d); }, s("power"), 0, sigma,
                                              std::fmax(1e9 * sigma, solve::BOUND_SPAN));
        solve::offset(s("mu2"), diff, [&](double d) { return P(k1, d); }, s("power"), "μ1", "μ2", true);
      }
    }
    const double k2 = k2_of(k1);
    r.set("k1", k1);
    r.set("k2", k2);
    r.set("m", m);
    r.set("n1", k1 * m);
    r.set("n2", k2 * m);
    r.set("n", (k1 + k2) * m);
    put_power(s, r, P(k1, diff));
    if (s.solve == "k") r.set("target_power", s("power"));
    r.set("mu1", s("mu2") + diff);
    r.set("mu2", s("mu2"));
    r.set("delta", diff);
    r.set("sigma", sigma);
    r.set("icc", icc);
    r.set("de", de);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c;
    if (s.solve == "k") c.push_back(col("target_power", "Target Power", "prob"));
    c.push_back(col("power", s.solve == "k" ? "Actual Power" : "Power", "prob"));
    for (auto x : {col("k1", "k1", "int", "Clusters in group 1"), col("k2", "k2", "int", "Clusters in group 2"),
                   col("m", "m", "num", "Subjects per cluster"), col("n1", "N1", "int"), col("n2", "N2", "int"),
                   col("n", "N", "int", "Total subjects"), col("mu1", "μ1", "num"), col("mu2", "μ2", "num"),
                   col("sigma", "σ", "num"), col("icc", "ICC", "num"), col("de", "Design Effect", "num"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    std::string pw = "achieve " + power_text(r.num("power")) + " power";
    if (s.solve == "k") pw += " (target " + power_text(s("power")) + ")";
    return fint(r.num("k1")) + " and " + fint(r.num("k2")) + " clusters of " + fnum(r.num("m")) +
           " subjects (" + fint(r.num("n")) + " in total) " + pw + " to detect a difference of " +
           fnum(r.num("delta")) + " between the means (σ = " + fnum(r.num("sigma")) + ", ICC = " +
           fnum_fine(r.num("icc")) + ", design effect " + fnum(r.num("de")) + ") with a " + alt_words(alt(s)) + " " +
           (s.is("test", "z") ? "z-test" : "t-test") + " at α = " + fnum(r.num("alpha")) + ".";
  };
  // The test of the cluster means has noncentrality λ = δ / SE: solving for the clusters, at the k1 found and one
  // below; solving for power, at the k1 given; solving for μ1, the target power gives λ, and λ the difference.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const ClusterMeans c(s);
    const double k1 = r.num("k1"), k2 = r.num("k2"), diff = r.num("delta"), se = c.se(k1, k2), ncp = diff / se;
    const double target = s.get("power", NAN), df = k1 + k2 - 2;
    const bool solved = s.solve == "k", same = c.kr == 1;
    // k2 is k1 itself with a ratio of 1, and otherwise ⌈r × k1⌉, derived
    const std::string K1 = size(k1, solved), K2 = solved || !same ? dvi(k2) : gvi(k2);
    const std::string at = same ? "k1 = k2 = " + gvi(k1) : "k1 = " + gvi(k1) + ", k2 = ⌈r × k1⌉ (r = " + gv(c.kr) + ")";
    Derivation out;
    if (solved)
      out.find = find_size(same ? "k1 = k2" : "k1, with k2 = ⌈r × k1⌉ (r = " + gv(c.kr) + ")", target);
    else if (s.solve == "power")
      out.find = find_power("with " + at);
    else
      out.find = find_effect("μ1", target, "with " + at);
    out.steps.push_back({"Design effect", "DE = 1 + (m − 1)ρ",
                         "1 + (" + gv(c.m) + " − 1) × " + gv(c.icc) + " = " + dv(c.de)});
    if (s.solve != "mu1") {
      out.steps.push_back({"Difference to detect", "δ = μ1 − μ2",
                           gv(s("mu1")) + " − " + gvs(s("mu2")) + " = " + dv(diff)});
      if (solved) {
        // the closed form of a z-test, beside the exact answer below
        const double za = stats::qnorm(tail_alpha(c.alpha, c.a), false), zb = stats::qnorm(target);
        const double approx = (za + zb) * (za + zb) * c.sigma * c.sigma * c.de * (1 + 1 / c.kr) / (c.m * diff * diff);
        out.steps.push_back(
            {"Normal approximation",
             std::string(same ? "k1 = k2 ≈ 2" : "k1 ≈ ") + "(z(" + level_symbol(c.a) + ") + z(1 − β))² σ² DE" +
                 (same ? "" : " (1 + 1/r)") + " / (m δ²), for comparison; the steps below give the answer",
             std::string(same ? "2 × " : "") + "(" + dv(za) + " + " + dv(zb) + ")² × " + gv(c.sigma) + "² × " +
                 dv(c.de) + (same ? "" : " × (1 + 1/" + gv(c.kr) + ")") + " / (" + gv(c.m) + " × " + dvs(diff) +
                 "²) = " + dv(approx)});
      }
    }
    if (!same)
      out.steps.push_back(
          {"Clusters in group 2", "k2 = ⌈r × k1⌉", "⌈" + gv(c.kr) + " × " + K1 + "⌉ = " + dvi(k2)});
    out.steps.push_back({"Standard error", "SE = σ √(DE/m × (1/k1 + 1/k2))",
                         gv(c.sigma) + " × √(" + dv(c.de) + "/" + gv(c.m) + " × (1/" + K1 + " + 1/" + K2 +
                             ")) = " + dv(se)});
    const auto critical = [&] {
      out.steps.push_back(c.z ? z_critical(c.alpha, c.a) : t_critical(c.alpha, df, c.a, "k1 + k2 − 2"));
    };
    if (s.solve == "mu1") {  // the target power gives λ, λ the difference, the difference μ1
      critical();
      out.steps.push_back(c.z ? z_ncp_needed(ncp, c.alpha, c.a, target)
                              : t_ncp_needed(ncp, df, c.alpha, c.a, target));
      out.steps.push_back({"Difference detected", "δ = λ × SE", dvs(ncp) + " × " + dv(se) + " = " + dv(diff)});
      out.steps.push_back(
          {"Mean detected", "μ1 = μ2 + δ", gv(r.num("mu2")) + " + " + dvs(diff) + " = " + dv(r.num("mu1"))});
      out.answer = "μ1 = " + dv(r.num("mu1")) + " (δ = " + dv(diff) + ")";
      return out;
    }
    out.steps.push_back({"Noncentrality", "λ = δ / SE", dvs(diff) + " / " + dv(se) + " = " + dv(ncp)});
    critical();
    out.steps.push_back(c.z ? z_power_step(ncp, c.alpha, c.a) : t_power_step(ncp, df, c.alpha, c.a));
    const double power = c.power(k1, diff);
    if (!solved) {
      out.answer = "power = " + dvp(power);
      return out;
    }
    // the search's numbers of clusters: k1 from 2 on, with at least 2 clusters in group 2 too
    const std::string rule = same ? "k1 = k2, the smallest whole number whose power reaches 1 − β"
                                  : "k1 = the smallest whole number whose power reaches 1 − β, with k2 = ⌈r × k1⌉";
    const auto sizes = [&](double x) { return "k1 = " + dvi(x) + ", k2 = " + dvi(c.k2(x)); };
    if (k1 - 1 < 2 || c.k2(k1 - 1) < 2)
      out.steps.push_back({"Smallest numbers of clusters", rule,
                           sizes(k1) + " are the smallest numbers allowed, and their power is already " +
                               dvp(power, target) + " ≥ " + gv(target)});
    else
      out.steps.push_back({"Smallest numbers of clusters", rule,
                           "power at " + sizes(k1 - 1) + ": " + dvp(c.power(k1 - 1, diff), target) + " < " +
                               gv(target) + " ≤ power at " + sizes(k1) + ": " + dvp(power, target)});
    out.answer = sizes(k1) + " (" + dvi(r.num("n")) + " subjects), with power " + dvp(power);
    return out;
  };
  p.definitions = {{"k1, k2", "Number of clusters randomized to each group."},
                   {"m", "Number of subjects per cluster."},
                   {"ICC", "Intracluster correlation coefficient ρ."},
                   {"Design Effect", "Variance inflation 1 + (m − 1)ρ relative to individual randomization."}};
  p.curve_param = "k1";
  p.curve_label = "k1 (Clusters in Group 1)";
  p.example_solve = "k";
  p.example = {{"power", "0.8 0.9"}, {"m", "10 20 40"}, {"icc", "0.05"}, {"mu1", "0.5"}};
  return p;
}
SSC_REGISTER(cluster_means);

// =====================================================================================
// Repeated measures: difference in time-averaged means
// =====================================================================================
double avg_var_factor(const Scenario& s) {
  const double m = s("m"), rho = s("rho");
  if (s.is("struct", "ar1")) {
    double sum = m;
    if (m <= 1000) {
      for (int k = 1; k < static_cast<int>(m); ++k) sum += 2 * (m - k) * std::pow(rho, k);
    } else {
      // the sum of (m − k)·rho^k in closed form, so that its cost does not grow with m (up to 1e9);
      // rho ≤ 0.999 keeps 1 − rho away from 0
      sum += 2 * rho * (m - 1 - m * rho + std::pow(rho, m)) / ((1 - rho) * (1 - rho));
    }
    return sum / (m * m);
  }
  return (1 + (m - 1) * rho) / m;
}

// The standard error of the difference in time-averaged means, σ √(v (1/n1 + 1/n2)) with v the variance factor of a
// subject's average, and the power of its test.
double repeated_se(double sigma, double vf, double n1, double n2) { return sigma * std::sqrt(vf * (1 / n1 + 1 / n2)); }

double repeated_power(const Scenario& s, double vf, double n1, double n2, double diff) {
  const double se = repeated_se(s("sigma"), vf, n1, n2), alpha = s("alpha");
  return s.is("test", "z") ? z_power(diff / se, alpha, alt(s)) : t_power(diff / se, n1 + n2 - 2, alpha, alt(s));
}

Procedure repeated() {
  Procedure p;
  p.id = "means.two.repeated";
  p.name = "Two Means in a Repeated-Measures Design (Time-Averaged Response)";
  p.category = "Means";
  p.subcategory = "Repeated Measures";
  p.order = 1;
  p.summary =
      "Compares two groups on the average of m repeated measurements per subject, accounting for the correlation "
      "between measurements (compound symmetry or AR(1)), following the approach of Diggle et al.";
  p.keywords = {"repeated measures", "longitudinal", "correlated", "compound symmetry", "ar(1)", "diggle", "gee", "time-averaged", "serial measurements", "multiple measurements per subject"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"},
                 {"power", "Power", "power", "Power"},
                 {"mu1", "Effect Size (μ1)", "mu1", "Detectable μ1"}};
  p.params = {
      alternative("Two-sided (H1: μ1 ≠ μ2)", "One-sided (H1: μ1 > μ2)", "One-sided (H1: μ1 < μ2)"),
      select("test", "Test Statistic", {{"z", "Z-test (large sample)", {}}, {"t", "T-test (df = N − 2)", {}}}, "z")
          .grp(G_TEST)
          .hint("The z-test is the usual large-sample approximation; the t-test with N − 2 degrees of freedom is "
                "accurate for small samples too. With fewer than about 20 subjects per group the z-test overstates "
                "the power by a few points (about 5 at 6–10 per group, 1 at 40)."),
      power_param(),
      alpha_param()};
  for (auto& q : two_group_params(2)) p.params.push_back(q);
  p.params.push_back(dropout_param());
  for (auto q : {num("m", "m (Measurements per Subject)", "4")
                     .grp(G_EFFECT)
                     .sym("m")
                     .at_least(1)
                     .whole()
                     .sugg({"3", "4", "2 4 6"})
                     .hint("How many times each subject is measured; the groups are compared on each subject's "
                           "average over these measurements."),
                 select("struct", "Correlation Structure",
                        {{"cs", "Compound symmetry (exchangeable)", {}}, {"ar1", "AR(1)", {}}}, "cs")
                     .grp(G_EFFECT)
                     .hint("Compound symmetry gives every pair of a subject's measurements the same correlation ρ; "
                           "AR(1) gives ρ to neighbors and ρ^k to measurements k steps apart."),
                 num("rho", "ρ (Correlation)", "0.5")
                     .grp(G_EFFECT)
                     .sym("ρ")
                     .range(-0.99, 0.999)
                     .sugg({"0.3", "0.5", "0.7"})
                     .hint("Correlation between two measurements of the same subject (between neighboring ones for "
                           "AR(1)); the higher it is, the less repeated measurements add."),
                 num("mu1", "μ1 (Average Mean, Group 1)", "0.5")
                     .grp(G_EFFECT)
                     .sym("μ1")
                     .sugg({"0.5", "0.3 0.5"})
                     .hint("Mean of group 1 averaged over the measurement times, expected under the alternative.")
                     .show(when(not_solving({"mu1"}))),
                 num("mu2", "μ2 (Average Mean, Group 2)", "0")
                     .grp(G_EFFECT)
                     .sym("μ2")
                     .hint("Mean of group 2 (e.g. the control) averaged over the measurement times."),
                 num("sigma", "σ (SD of a Single Measurement)", "1")
                     .grp(G_EFFECT)
                     .sym("σ")
                     .above(0)
                     .hint("Standard deviation of one measurement among subjects in the same group, assumed the "
                           "same at every time.")})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double alpha = s("alpha"), sigma = s("sigma");
    const double vf = avg_var_factor(s);
    if (!(vf > 0)) throw CalcError("the correlation structure is not positive definite for these values");
    const Alt a = alt(s);
    auto P = [&](double n1, double n2, double diff) { return repeated_power(s, vf, n1, n2, diff); };
    TwoN n;
    double diff;
    if (s.solve == "n") {
      diff = s("mu1") - s("mu2");
      if (diff == 0) throw CalcError("μ1 equals μ2");
      if ((a == Alt::Upper && diff < 0) || (a == Alt::Lower && diff > 0))
        throw CalcError("the difference μ1 − μ2 has the sign opposite to the one-sided alternative");
      n = solve_two_group_n(s, [&](double a1, double a2) { return P(a1, a2, diff); }, s("power"), 2,
                            z_factor(alpha, s("power"), a) * 2 * sigma * sigma * vf / (diff * diff) + 1, true, ctx);
    } else {
      n = two_group_n(s);
      if (s.solve == "power") {
        diff = s("mu1") - s("mu2");
      } else {
        const double sign = a == Alt::Lower ? -1 : 1;
        diff = sign * solve::solve_increasing([&](double d) { return P(n.n1, n.n2, sign * d); }, s("power"), 0, sigma,
                                              std::fmax(1e9 * sigma, solve::BOUND_SPAN));
        solve::offset(s("mu2"), diff, [&](double d) { return P(n.n1, n.n2, d); }, s("power"), "μ1", "μ2", true);
      }
    }
    put_two_group(r, n);
    put_power(s, r, P(n.n1, n.n2, diff));
    r.set("mu1", s("mu2") + diff);
    r.set("mu2", s("mu2"));
    r.set("m", s("m"));
    r.set("rho", s("rho"));
    r.set("sigma", sigma);
    r.set("vf", vf * s("m"));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("mu1", "μ1", "num"), col("mu2", "μ2", "num"), col("m", "m", "int", "Measurements per subject"),
                   col("rho", "ρ", "num"), col("sigma", "σ", "num"),
                   col("vf", "Variance Factor", "num", "Var(average) × m/σ² (1 for independent measurements)"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + ", each subject measured " +
           fint(r.num("m")) + " times, " + achieves(s, r) + " to detect a difference of " +
           fnum(r.num("mu1") - r.num("mu2")) + " in the time-averaged means when σ = " + fnum(r.num("sigma")) +
           " and the correlation between measurements is " + fnum_fine(r.num("rho")) + " (" +
           (s.is("struct", "ar1") ? "AR(1)" : "compound symmetry") + "), using a " + alt_words(alt(s)) + " test at α = " +
           fnum(r.num("alpha")) + ".";
  };
  // The test of the subjects' averages has noncentrality λ = δ / SE, SE = σ √(v (1/N1 + 1/N2)): solving for N, at
  // the sizes found and the next smaller; for power, at the sizes given; for μ1, the target power gives λ, and λ the
  // difference.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const Alt a = alt(s);
    const bool z = s.is("test", "z"), solved = s.solve == "n";
    const double alpha = s("alpha"), sigma = s("sigma"), m = s("m"), rho = s("rho"), vf = avg_var_factor(s);
    const double target = s.get("power", NAN);
    const blocks::TwoN n{r.num("n1"), r.num("n2")};
    // the difference the power was computed at: the one given, or the one found
    const double diff = s.solve == "mu1" ? r.num("mu1") - r.num("mu2") : s("mu1") - s("mu2");
    const double se = repeated_se(sigma, vf, n.n1, n.n2), ncp = diff / se, df = n.n1 + n.n2 - 2;
    const auto [N1, N2] = means_group_sizes(s, n, solved);
    Derivation out;
    if (solved)
      out.find = find_size(two_group_unknown(s), target);
    else if (s.solve == "power")
      out.find = find_power("with " + means_sizes_at(s, n));
    else
      out.find = find_effect("μ1", target, "with " + means_sizes_at(s, n));
    if (s.is("struct", "ar1")) {
      std::string sum;
      if (m <= 6) {  // the terms of the sum; beyond, the sum itself
        std::vector<std::string> t;
        for (int k = 1; k < static_cast<int>(m); ++k) t.push_back(dvi(m - k) + " × " + gvs(rho) + superscript(k));
        sum = t.empty() ? "" : " + 2 × (" + joined(t, " + ") + ")";
      } else {
        sum = " + 2 Σₖ (" + gvi(m) + " − k) × " + gvs(rho) + "ᵏ";
      }
      out.steps.push_back({"Variance of a subject's average",
                           "v = Var(Ȳ)/σ² = (m + 2 Σₖ (m − k) ρᵏ) / m², k = 1, …, m − 1",
                           "(" + gvi(m) + sum + ") / " + gvi(m) + "² = " + dv(vf)});
    } else {
      out.steps.push_back({"Variance of a subject's average", "v = Var(Ȳ)/σ² = (1 + (m − 1)ρ) / m",
                           "(1 + (" + gvi(m) + " − 1) × " + gvs(rho) + ") / " + gvi(m) + " = " + dv(vf)});
    }
    if (s.solve != "mu1")
      out.steps.push_back(
          {"Difference to detect", "δ = μ1 − μ2", gv(s("mu1")) + " − " + gvs(s("mu2")) + " = " + dv(diff)});
    if (solved)
      means_alloc_approx(out.steps, s, a, {"σ² v", gv(sigma) + "² × " + dv(vf), sigma * sigma * vf},
                         {"δ", dvs(diff), diff});
    else
      means_alloc_step(out.steps, s, n);
    out.steps.push_back({"Standard error", "SE = σ √(v (1/N1 + 1/N2))",
                         gv(sigma) + " × √(" + dv(vf) + " × (1/" + N1 + " + 1/" + N2 + ")) = " + dv(se)});
    const auto critical = [&] {
      out.steps.push_back(z ? z_critical(alpha, a) : t_critical(alpha, df, a, "N1 + N2 − 2"));
    };
    if (s.solve == "mu1") {  // the target power gives λ, λ the difference, the difference μ1
      critical();
      out.steps.push_back(z ? z_ncp_needed(ncp, alpha, a, target) : t_ncp_needed(ncp, df, alpha, a, target));
      out.steps.push_back({"Difference detected", "δ = λ × SE", dvs(ncp) + " × " + dv(se) + " = " + dv(diff)});
      out.steps.push_back(
          {"Mean detected", "μ1 = μ2 + δ", gv(r.num("mu2")) + " + " + dvs(diff) + " = " + dv(r.num("mu1"))});
      out.answer = "μ1 = " + dv(r.num("mu1")) + " (δ = " + dv(diff) + ")";
      return out;
    }
    out.steps.push_back({"Noncentrality", "λ = δ / SE", dvs(diff) + " / " + dv(se) + " = " + dv(ncp)});
    critical();
    out.steps.push_back(z ? z_power_step(ncp, alpha, a) : t_power_step(ncp, df, alpha, a));
    const auto power = [&](double x1, double x2) { return repeated_power(s, vf, x1, x2, diff); };
    if (solved) {
      out.steps.push_back(means_two_search(s, power, n, target, 2));
      out.answer = means_two_answer(s, n, power(n.n1, n.n2));
    } else {
      out.answer = "power = " + dvp(power(n.n1, n.n2));
    }
    return out;
  };
  p.definitions = {{"m", "Number of repeated measurements per subject."},
                   {"ρ", "Correlation between measurements (between adjacent ones for AR(1))."},
                   {"Variance Factor", "m·Var(subject average)/σ² — 1 + (m − 1)ρ under compound symmetry."}};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8"}, {"m", "2 4 6"}, {"rho", "0.3 0.6"}, {"mu1", "0.5"}};
  return p;
}
SSC_REGISTER(repeated);

}  // namespace
}  // namespace ssc::procs
