// Designs for means beyond two parallel groups: one-way ANOVA, 2×2 crossover
// (difference and bioequivalence ratio), cluster-randomized and
// repeated-measures comparisons.
#include "procedures/common.hpp"
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

// effect-size multiplier `scale` stretches the means about their mean (for solving the effect)
AnovaDesign anova_design(const Scenario& s, double n, double scale_f = NAN) {
  const int k = anova_k(s);
  const auto w = anova_weights(s, k);
  AnovaDesign d;
  d.k = k;
  for (double wi : w) d.sizes.push_back(std::ceil(n * wi - 1e-9));
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
  double mbar = 0;
  for (int i = 0; i < k; ++i) mbar += d.sizes[i] * mu[i];
  mbar /= N;
  double ss = 0;
  for (int i = 0; i < k; ++i) ss += d.sizes[i] * (mu[i] - mbar) * (mu[i] - mbar);
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
          .above(0)
          .sugg({"0.1", "0.25", "0.4", "0.1 0.25 0.4"})
          .hint("Cohen's f: the SD of the group means divided by the within-group SD. By Cohen's conventions 0.10 "
                "is small, 0.25 medium and 0.40 large.")
          .show(when({is("estype", {"f"}), not_solving({"f"})})),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double alpha = s("alpha");
    AnovaDesign d;
    if (s.solve == "n") {
      d = anova_design(s, 2);
      if (!(d.f > 0)) throw CalcError("the group means are all equal; there is no effect to detect");
      solve::NSearch o;
      o.lo = 2;
      o.guess = (z_factor(alpha, s("power"), Alt::Two) + 2.0 * (d.k - 1)) / (d.k * d.f * d.f) + 1;
      const double n = solve::min_n([&](double x) { return anova_power(anova_design(s, x), alpha); }, s("power"), o, &ctx);
      d = anova_design(s, n);
      r.set("n", n);
    } else if (s.solve == "power") {
      d = anova_design(s, s("n"));
      r.set("n", s("n"));
    } else {
      const double n = s("n");
      const double f = solve::solve_increasing([&](double ff) { return anova_power(anova_design(s, n, ff), alpha); },
                                               s("power"), 1e-9, 0.25, 1e6, 1e-12, "the effect size f");
      d = anova_design(s, n, f);
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
  p.definitions = {{"n, N", "Group size and total sample size."},
                   {"k", "Number of groups."},
                   {"σm", "Standard deviation of the group means: √(Σ nᵢ(μᵢ − μ̄)²/N)."},
                   {"f", "Cohen's effect size σm/σ."},
                   {"λ", "Noncentrality parameter N f² of the F distribution with (k − 1, N − k) df."}};
  p.references = {REF_COHEN,
                  "Desu, M. M. & Raghavarao, D. (1990). Sample Size Methodology. Academic Press.",
                  "Fleiss, J. L. (1986). The Design and Analysis of Clinical Experiments. Wiley.", REF_BK};
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
    const bool z = s.is("test", "z");
    auto P = [&](double n, double d) {
      const double ncp = d / (sw * std::sqrt(2 / n));
      return z ? z_power(ncp, alpha, a) : t_power(ncp, n - 2, alpha, a);
    };
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
        d = sign * solve::solve_increasing([&](double dd) { return P(n, sign * dd); }, s("power"), 0, sw, 1e9 * sw);
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
  p.definitions = {{"N", "Total number of subjects; half receive AB and half BA."},
                   {"σw", "Within-subject (residual) standard deviation."},
                   {"σd", "SD of a subject's period difference, √2 σw."},
                   {"δ", "Difference between the treatment means."}};
  p.references = {REF_CHOW, REF_JULIOUS,
                  "Senn, S. (2002). Cross-over Trials in Clinical Research, 2nd ed. Wiley."};
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
    if (!(th0 > th1 && th0 < th2)) throw CalcError("the true ratio must lie strictly between the limits");
    // √ln(1 + CV²) = CV(1 − CV²/4 + …) is CV itself to double precision below 1e-8 (computed, CV² underflows to 0
    // below about 1e-162, and so did the log-scale SD)
    const double sw = cv < 1e-8 ? cv : std::sqrt(std::log1p(cv * cv));
    const bool cross = s.is("design", "2x2");
    auto P = [&](double n) {
      const double se = sw * std::sqrt((cross ? 2.0 : 4.0) / n);
      return tost_power(std::log(th1), std::log(th2), std::log(th0), se, n - 2, alpha);
    };
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
  p.definitions = {{"θ0", "True ratio of the test to the reference mean (geometric means)."},
                   {"θL, θU", "Bioequivalence acceptance limits for the ratio."},
                   {"CV", "Coefficient of variation; within-subject for cross-over designs."},
                   {"σ (log)", "√ln(1 + CV²), the standard deviation of log-transformed data."}};
  p.references = {"Schuirmann, D. J. (1987). A comparison of the two one-sided tests procedure and the power approach "
                  "for assessing the equivalence of average bioavailability. Journal of Pharmacokinetics and "
                  "Biopharmaceutics 15, 657–680.",
                  "Owen, D. B. (1965). A special case of a bivariate non-central t-distribution. Biometrika 52, 437–446.",
                  "Diletti, E., Hauschke, D. & Steinijans, V. W. (1991). Sample size determination for bioequivalence "
                  "assessment by means of confidence intervals. Int. J. Clin. Pharmacol. Ther. Toxicol. 29, 1–8.",
                  REF_CHOW};
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
    const double alpha = s("alpha"), m = s("m"), icc = s("icc"), sigma = s("sigma"), kr = s("kr");
    const double de = 1 + (m - 1) * icc;
    const Alt a = alt(s);
    const bool z = s.is("test", "z");
    auto k2_of = [&](double k1) { return std::ceil(kr * k1 - 1e-9); };
    auto P = [&](double k1, double diff) {
      const double k2 = k2_of(k1);
      const double se = sigma * std::sqrt(de / m * (1 / k1 + 1 / k2));
      return z ? z_power(diff / se, alpha, a) : t_power(diff / se, k1 + k2 - 2, alpha, a);
    };
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
        throw CalcError("group 2 would have " + fint(k2_of(k1)) + " cluster (k1 × the cluster ratio, rounded up); "
                        "each group needs at least 2 clusters");
      if (s.solve == "power") {
        diff = s("mu1") - s("mu2");
      } else {
        const double sign = a == Alt::Lower ? -1 : 1;
        diff = sign * solve::solve_increasing([&](double d) { return P(k1, sign * d); }, s("power"), 0, sigma,
                                              1e9 * sigma);
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
           fnum(r.num("icc")) + ", design effect " + fnum(r.num("de")) + ") with a " + alt_words(alt(s)) + " " +
           (s.is("test", "z") ? "z-test" : "t-test") + " at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"k1, k2", "Number of clusters randomized to each group."},
                   {"m", "Number of subjects per cluster."},
                   {"ICC", "Intracluster correlation coefficient ρ."},
                   {"Design Effect", "Variance inflation 1 + (m − 1)ρ relative to individual randomization."}};
  p.references = {"Donner, A. & Klar, N. (2000). Design and Analysis of Cluster Randomization Trials in Health "
                  "Research. Arnold.",
                  "Campbell, M. J. & Walters, S. J. (2014). How to Design, Analyse and Report Cluster Randomised "
                  "Trials in Medicine and Health Related Research. Wiley.",
                  REF_MACHIN};
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
    const bool z = s.is("test", "z");
    auto P = [&](double n1, double n2, double diff) {
      const double se = sigma * std::sqrt(vf * (1 / n1 + 1 / n2));
      return z ? z_power(diff / se, alpha, a) : t_power(diff / se, n1 + n2 - 2, alpha, a);
    };
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
                                              1e9 * sigma);
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
           " and the correlation between measurements is " + fnum(r.num("rho")) + " (" +
           (s.is("struct", "ar1") ? "AR(1)" : "compound symmetry") + "), using a " + alt_words(alt(s)) + " test at α = " +
           fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"m", "Number of repeated measurements per subject."},
                   {"ρ", "Correlation between measurements (between adjacent ones for AR(1))."},
                   {"Variance Factor", "m·Var(subject average)/σ² — 1 + (m − 1)ρ under compound symmetry."}};
  p.references = {"Diggle, P. J., Heagerty, P., Liang, K.-Y. & Zeger, S. L. (2002). Analysis of Longitudinal Data, "
                  "2nd ed. Oxford University Press.",
                  "Liu, G. & Liang, K.-Y. (1997). Sample size calculations for studies with correlated observations. "
                  "Biometrics 53, 937–947."};
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
