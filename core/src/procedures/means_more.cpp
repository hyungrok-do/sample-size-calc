// More designs for means: two groups adjusted for covariates (ANCOVA), the
// two-way factorial ANOVA and the one-way repeated-measures ANOVA.
#include "procedures/common.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

const char* REF_COHEN = "Cohen, J. (1988). Statistical Power Analysis for the Behavioral Sciences, 2nd ed. Lawrence "
                        "Erlbaum.";

// =====================================================================================
// ANCOVA: two means adjusted for covariates
// =====================================================================================
// Covariates explaining a fraction R² of the response variance reduce the
// residual standard deviation to σ√(1 − R²); each covariate costs one degree
// of freedom (Borm, Fransen & Lemmens, 2007).
double ancova_power(double n1, double n2, double diff, double sigma, double r2, double q, double alpha, Alt a) {
  const double df = n1 + n2 - 2 - q;
  if (!(df > 0)) return NAN;
  return t_power(diff / (sigma * std::sqrt((1 - r2) * (1 / n1 + 1 / n2))), df, alpha, a);
}

Procedure ancova() {
  Procedure p;
  p.id = "means.two.ancova";
  p.name = "Two Means Adjusted for Covariates (ANCOVA)";
  p.category = "Means";
  p.subcategory = "Two Independent Means";
  p.order = 5;
  p.summary =
      "Sample size and power for comparing two group means after adjusting for one or more baseline covariates "
      "(analysis of covariance). Covariates that explain a fraction R² of the outcome variance shrink the residual "
      "standard deviation to σ√(1 − R²).";
  p.keywords = {"ancova", "covariate", "baseline adjustment", "analysis of covariance", "regression adjustment",
                "two-sample", "baseline covariate", "pre-test", "adjusted analysis"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"},
                 {"power", "Power", "power", "Power"},
                 {"mu1", "Effect Size (μ1)", "mu1", "Detectable μ1"}};
  p.params = {alternative("Two-sided (H1: μ1 ≠ μ2)", "One-sided (H1: μ1 > μ2)", "One-sided (H1: μ1 < μ2)"),
              power_param(), alpha_param()};
  for (auto& q : two_group_params(3, "50")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  for (auto q : {num("mu1", "μ1 (Adjusted Mean of Group 1)", "0.5")
                     .grp(G_EFFECT)
                     .sym("μ1")
                     .sugg({"0.5", "0.3 0.5"})
                     .hint("Mean outcome of group 1 (e.g. treatment) after adjusting for the covariates, expected "
                           "under the alternative.")
                     .show(when(not_solving({"mu1"}))),
                 num("mu2", "μ2 (Adjusted Mean of Group 2)", "0")
                     .grp(G_EFFECT)
                     .sym("μ2")
                     .hint("Mean outcome of group 2 (e.g. control) after adjusting for the covariates; μ1 − μ2 is "
                           "the treatment effect to detect."),
                 num("sigma", "σ (Standard Deviation, Unadjusted)", "1")
                     .grp(G_EFFECT)
                     .sym("σ")
                     .above(0)
                     .hint("Standard deviation of the outcome within a group, before adjusting for the covariates; "
                           "take it from earlier studies or a pilot."),
                 num("r2", "R² (Variance Explained by the Covariates)", "0.25")
                     .grp(G_EFFECT)
                     .sym("R²")
                     .range(0, 1, false, true)
                     .sugg({"0.25", "0 0.25 0.5", "0.36"})
                     .hint("For one covariate, the squared correlation between the covariate and the outcome."),
                 num("ncov", "Number of Covariates", "1")
                     .grp(G_EFFECT)
                     .sym("Covariates")
                     .at_least(0)
                     .whole()
                     .single()
                     .hint("How many covariates (such as the baseline value of the outcome) the analysis adjusts "
                           "for; each one costs a degree of freedom.")})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double mu2 = s("mu2"), sigma = s("sigma"), r2 = s("r2"), q = s("ncov"), alpha = s("alpha");
    const Alt a = alt(s);
    TwoN n;
    double mu1;
    if (s.solve == "mu1") {
      n = two_group_n(s);
      if (!(n.n1 + n.n2 - 2 - q > 0)) throw CalcError("the number of covariates must be less than N − 2");
      const double sign = a == Alt::Lower ? -1.0 : 1.0;
      const double se = sigma * std::sqrt((1 - r2) * (1 / n.n1 + 1 / n.n2));
      auto power = [&](double dd) { return ancova_power(n.n1, n.n2, dd, sigma, r2, q, alpha, a); };
      const double d = solve::solve_increasing([&](double dd) { return power(sign * dd); }, s("power"), 0, se,
                                               std::fmax(1e6 * se, solve::BOUND_SPAN), 1e-12, "μ1");
      mu1 = solve::offset(mu2, sign * d, power, s("power"), "μ1", "μ2", true);
    } else {
      mu1 = s("mu1");
      const double diff = mu1 - mu2;
      if (diff == 0) throw CalcError("μ1 equals μ2");
      if ((a == Alt::Upper && diff < 0) || (a == Alt::Lower && diff > 0))
        throw CalcError("the difference has the sign opposite to the one-sided alternative");
      if (s.solve == "n") {
        const double guess = 2 * z_factor(alpha, s("power"), a) * sigma * sigma * (1 - r2) / (diff * diff) + 2;
        // at least 3 per group, as the inputs; designs with too few degrees of freedom for the covariates fall short
        n = solve_two_group_n(
            s,
            [&](double n1, double n2) {
              const double pw = ancova_power(n1, n2, diff, sigma, r2, q, alpha, a);
              return std::isnan(pw) ? 0.0 : pw;
            },
            s("power"), 3, guess, true, ctx);
      } else {
        n = two_group_n(s);
      }
    }
    if (!(n.n1 + n.n2 - 2 - q > 0)) throw CalcError("the number of covariates must be less than N − 2");
    put_two_group(r, n);
    put_power(s, r, ancova_power(n.n1, n.n2, mu1 - mu2, sigma, r2, q, alpha, a));
    r.set("mu1", mu1);
    r.set("mu2", mu2);
    r.set("sigma", sigma);
    r.set("r2", r2);
    r.set("sadj", sigma * std::sqrt(1 - r2));
    r.set("ncov", q);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("mu1", "μ1", "num"), col("mu2", "μ2", "num"), col("sigma", "σ", "num"),
                   col("r2", "R²", "num", "Variance explained by the covariates"),
                   col("sadj", "Residual σ", "num", "σ√(1 − R²)"), col("ncov", "Covariates", "int"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
           " to detect a difference of " + fnum(r.num("mu1") - r.num("mu2")) +
           " between the covariate-adjusted means with an ANCOVA t-test (" + fint(r.num("ncov")) + " covariate" +
           (r.num("ncov") == 1 ? "" : "s") +
           (r.num("r2") > 0 ? " explaining R² = " + fnum(r.num("r2")) + " of the variance, so that the residual "
                              "standard deviation is " + fnum(r.num("sadj")) + " instead of " + fnum(r.num("sigma"))
                            : " explaining none of the variance, so that the standard deviation stays " +
                                  fnum(r.num("sigma"))) +
           "), " + alt_words(alt(s)) + " at α = " + fnum(r.num("alpha")) + ".";
  };
  p.notes = {"The formula (Borm et al.) leaves out the chance imbalance of the covariates between the groups, so it "
             "slightly overstates the power: by less than half a point with one covariate, but by several points "
             "with many covariates and small groups (5 covariates and 10 per group: 0.83 instead of 0.70)."};
  p.definitions = {{"R²", "Fraction of the outcome variance explained by the covariates (pooled within groups)."},
                   {"Residual σ", "Standard deviation after covariate adjustment, σ√(1 − R²)."}};
  p.references = {"Borm, G. F., Fransen, J. & Lemmens, W. A. J. G. (2007). A simple sample size formula for analysis "
                  "of covariance in randomized clinical trials. Journal of Clinical Epidemiology 60, 1234–1238.",
                  REF_CHOW};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"mu1", "0.5"}, {"mu2", "0"}, {"sigma", "1"}, {"r2", "0 0.25 0.5"}};
  return p;
}
SSC_REGISTER(ancova);

// =====================================================================================
// Two-way factorial ANOVA
// =====================================================================================
struct Factorial {
  int a = 0, b = 0;
  double ssa = 0, ssb = 0, ssab = 0;  // Σ effects² (main effects summed over their own levels)
  bool equal = true;                  // every cell mean the same number
};

// A sum kept to about twice the digits of a double, hi + lo, with the rounding error of each addition carried in lo
// (Knuth's two-sum): an effect is a difference of sums of the cell means, and next to large means it is smaller than
// the digits that a sum in doubles keeps (a row of 1e6 shifted the effect of B in the row of 0 0 1e-9 by 5%).
struct Sum2 {
  double hi = 0, lo = 0;
  void add(double x) {
    const double t = hi + x, v = t - hi;
    lo += (hi - (t - v)) + (x - v);
    hi = t;
  }
  void add(double x, double w) {  // x × w, exactly for a whole w (the product's rounding error from fma)
    const double p = x * w;
    add(p);
    add(std::fma(x, w, -p));
  }
  void add(const Sum2& x, double w) {
    add(x.hi, w);
    add(x.lo, w);
  }
  double value() const { return hi + lo; }
};

Factorial factorial(const Scenario& s) {
  Factorial f;
  f.a = static_cast<int>(s("a"));
  f.b = static_cast<int>(s("b"));
  const int A = f.a, B = f.b;
  const auto& m = s.v("means");
  if (static_cast<int>(m.size()) != A * B)
    throw CalcError("enter " + std::to_string(A * B) + " cell means (" + std::to_string(A) + " rows of " +
                    std::to_string(B) + ")");
  auto cell = [&](int i, int j) { return m[static_cast<std::size_t>(i * B + j)]; };
  // The effects times A·B, from the sums of the rows R, the columns C and all cells G: A R_i − G for the rows,
  // B C_j − G for the columns and A·B m_ij − A R_i − B C_j + G for the interaction.  Each is the sum of the cells
  // with whole weights, and next to it the most that the rounding of the means as typed can make of it: the sum of
  // the cells' half spacings h of the doubles there, with the same weights in size.
  std::vector<Sum2> R(A), C(B);
  Sum2 G;
  std::vector<double> hr(A, 0.0), hc(B, 0.0);
  double H = 0;
  auto half = [](double x) { return 0.5 * (std::nextafter(std::fabs(x), HUGE_VAL) - std::fabs(x)); };
  for (int i = 0; i < A; ++i)
    for (int j = 0; j < B; ++j) {
      const double v = cell(i, j), h = half(v);
      R[i].add(v);
      C[j].add(v);
      G.add(v);
      hr[i] += h;
      hc[j] += h;
      H += h;
      if (v != m[0]) f.equal = false;
    }
  // an effect smaller than what the rounding of the means as typed can make is none: 0.1 0.2 0.3 0.4 are additive,
  // but as doubles they leave an interaction of 7e-18, 7e12 times σ = 1e-30, where their rounding can make 1.9e-17;
  // one at least as large as that is kept however small, since no means that round to these are without it
  const double cells = static_cast<double>(A) * B;
  std::vector<double> ea(A), eb(B), eab(static_cast<std::size_t>(A * B));
  bool real_a = false, real_b = false, real_ab = false;
  for (int i = 0; i < A; ++i) {
    Sum2 e;
    e.add(R[i], A);
    e.add(G, -1);
    ea[i] = e.value();
    if (!(std::fabs(ea[i]) < (A - 1) * hr[i] + (H - hr[i]))) real_a = true;
  }
  for (int j = 0; j < B; ++j) {
    Sum2 e;
    e.add(C[j], B);
    e.add(G, -1);
    eb[j] = e.value();
    if (!(std::fabs(eb[j]) < (B - 1) * hc[j] + (H - hc[j]))) real_b = true;
  }
  for (int i = 0; i < A; ++i)
    for (int j = 0; j < B; ++j) {
      Sum2 e;
      e.add(cell(i, j), cells);
      e.add(R[i], -A);
      e.add(C[j], -B);
      e.add(G, 1);
      const double h = half(cell(i, j)), v = e.value();
      eab[static_cast<std::size_t>(i * B + j)] = v;
      // weights (A − 1)(B − 1) on the cell, 1 − A on the rest of its row, 1 − B on the rest of its column, 1 elsewhere
      const double bound = (A - 1.0) * (B - 1) * h + (A - 1.0) * (hr[i] - h) + (B - 1.0) * (hc[j] - h) +
                           (H - hr[i] - hc[j] + h);
      if (!(std::fabs(v) < bound)) real_ab = true;
    }
  for (double e : ea)
    if (real_a) f.ssa += (e / cells) * (e / cells);
  for (double e : eb)
    if (real_b) f.ssb += (e / cells) * (e / cells);
  for (double e : eab)
    if (real_ab) f.ssab += (e / cells) * (e / cells);
  return f;
}

struct FactorialPower {
  double a, b, ab;
};

FactorialPower factorial_power(const Factorial& f, double n, double sigma, double alpha) {
  using namespace stats;
  const double s2 = sigma * sigma, dfe = f.a * f.b * (n - 1);
  auto pw = [&](double df1, double lambda) -> double {
    if (!(dfe > 0)) return NAN;
    return pnf(qf(alpha, df1, dfe, false), df1, dfe, lambda, false);
  };
  return {pw(f.a - 1, n * f.b * f.ssa / s2), pw(f.b - 1, n * f.a * f.ssb / s2),
          pw((f.a - 1) * (f.b - 1), n * f.ssab / s2)};
}

Procedure anova2() {
  Procedure p;
  p.id = "means.anova2";
  p.name = "Two-Way Factorial ANOVA (Main Effects and Interaction)";
  p.category = "Means";
  p.subcategory = "Many Means (ANOVA)";
  p.order = 2;
  p.summary =
      "Power and per-cell sample size for the F-tests of a balanced two-way factorial design: the main effects of "
      "factors A and B and their interaction, from the hypothesized cell means and the within-cell standard "
      "deviation (noncentral F distribution).";
  p.keywords = {"factorial", "two-way anova", "interaction", "main effects", "f-test", "2x2", "two factors", "factorial design", "2 by 2"};
  p.solve_for = {{"n", "Sample Size", "n", "n (Per Cell)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      select("goal", "Size the Study For",
             {{"all", "Every effect that is present", {}}, {"A", "Main effect of A", {}}, {"B", "Main effect of B", {}},
              {"AB", "A × B interaction", {}}},
             "all")
          .grp(G_TEST)
          .hint("Which F-test must reach the target power: one chosen effect, or every effect the cell means "
                "imply (then the least powerful of them decides).")
          .show(when(solving({"n"}))),
      power_param(),
      alpha_param(),
      n_param("n", "n (Subjects per Cell)", "10", 2)
          .sym("n")
          .hint("Number of subjects in each cell (each combination of a level of A with a level of B); the "
                "total is n × A × B."),
      num("a", "Levels of Factor A (Rows)", "2")
          .grp(G_EFFECT)
          .sym("A Levels")
          .range(2, 100)  // the cell count A × B stays a small whole number
          .whole()
          .single()
          .hint("Number of levels (categories) of the first factor, A; they form the rows of the table of cell "
                "means."),
      num("b", "Levels of Factor B (Columns)", "2")
          .grp(G_EFFECT)
          .sym("B Levels")
          .range(2, 100)
          .whole()
          .single()
          .hint("Number of levels of the second factor, B; they form the columns of the table of cell means."),
      vec("means", "Cell Means (Row by Row)", "0 0 0 0.5")
          .grp(G_EFFECT)
          .sym("Cell Means")
          .sugg({"0 0 0 0.5", "0 0.5 0.5 1", "0 0 0 0 0.5 1"})
          .hint("A × B means listed row by row: first all levels of B for the first level of A, and so on."),
      num("sigma", "σ (Within-Cell Standard Deviation)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .hint("Standard deviation of the outcome among subjects in the same cell, assumed equal in all cells."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Factorial f = factorial(s);
    const double sigma = s("sigma"), alpha = s("alpha");
    // the power that decides the design: of the chosen effect, or the least of the effects present (also when
    // solving for power, so that the power curve and the card have one power to show)
    const std::string goal = s.sel.count("goal") ? s.s("goal") : "all";
    const double eps = 1e-12 * sigma * sigma;
    const bool hasA = f.ssa > eps, hasB = f.ssb > eps, hasAB = f.ssab > eps;
    auto P = [&](double x) {
      const auto p3 = factorial_power(f, x, sigma, alpha);
      if (goal == "A") return p3.a;
      if (goal == "B") return p3.b;
      if (goal == "AB") return p3.ab;
      double m = 1;
      if (hasA) m = std::min(m, p3.a);
      if (hasB) m = std::min(m, p3.b);
      if (hasAB) m = std::min(m, p3.ab);
      return hasA || hasB || hasAB ? m : alpha;
    };
    double n;
    if (s.solve == "n") {
      if (f.equal) throw CalcError("the cell means are all equal; there is no effect to detect");
      if (f.ssa == 0 && f.ssb == 0 && f.ssab == 0)
        throw CalcError("the cell means differ by no more than the rounding of numbers of their size; there is no "
                        "effect to detect");
      if ((goal == "A" && !hasA) || (goal == "B" && !hasB) || (goal == "AB" && !hasAB) ||
          (goal == "all" && !(hasA || hasB || hasAB)))
        throw CalcError("the cell means have no effect of the chosen kind");
      solve::NSearch o;
      o.lo = 2;
      n = solve::min_n(P, s("power"), o, &ctx);
      r.set("target_power", s("power"));
      r.set("goal", goal == "all" ? std::string("all effects") : goal);
    } else {
      n = s("n");
    }
    const auto p3 = factorial_power(f, n, sigma, alpha);
    r.set("n", n);
    r.set("N", n * f.a * f.b);
    r.set("power", P(n));
    r.set("pa", p3.a);
    r.set("pb", p3.b);
    r.set("pab", p3.ab);
    const double cells = f.a * f.b;
    r.set("fa", std::sqrt(f.ssa / f.a) / sigma);
    r.set("fb", std::sqrt(f.ssb / f.b) / sigma);
    r.set("fab", std::sqrt(f.ssab / cells) / sigma);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c;
    if (s.solve == "n") {
      c.push_back(col("target_power", "Target Power", "prob"));
      c.push_back(col("goal", "Sized For", "text"));
    }
    for (auto x : {col("n", "n per Cell", "int"), col("N", "N (Total)", "int"),
                   col("pa", "Power (A)", "prob", "Power of the F-test for the main effect of A"),
                   col("pb", "Power (B)", "prob", "Power of the F-test for the main effect of B"),
                   col("pab", "Power (A×B)", "prob", "Power of the F-test for the interaction"),
                   col("fa", "f (A)", "num", "Cohen's effect size of A"), col("fb", "f (B)", "num"),
                   col("fab", "f (A×B)", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "With " + fint(r.num("n")) + " subjects in each cell (" + fint(r.num("N")) +
           " in total), the F-tests at α = " + fnum(r.num("alpha")) + " have power " + fnum(r.num("pa"), 4) +
           " for the main effect of A (f = " + fnum(r.num("fa"), 3) + "), " + fnum(r.num("pb"), 4) +
           " for the main effect of B (f = " + fnum(r.num("fb"), 3) + ") and " + fnum(r.num("pab"), 4) +
           " for the interaction (f = " + fnum(r.num("fab"), 3) + ").";
  };
  p.definitions = {{"n per Cell", "Subjects in each of the A × B cells (balanced design)."},
                   {"f", "Cohen's effect size: SD of the effects divided by the within-cell SD."},
                   {"Power (A), (B), (A×B)", "Power of the F-tests for the two main effects and the interaction."}};
  p.references = {REF_COHEN, "Montgomery, D. C. (2017). Design and Analysis of Experiments, 9th ed. Wiley."};
  p.curve_param = "n";
  p.curve_label = "n (Per Cell)";
  p.curve_min = 2;
  p.example_solve = "n";
  p.example = {{"goal", "all"}, {"power", "0.8"}, {"a", "2"}, {"b", "2"}, {"means", "0 0.5 0.5 1.5"}, {"sigma", "1"}};
  return p;
}
SSC_REGISTER(anova2);

// =====================================================================================
// One-way repeated-measures ANOVA (one within-subject factor)
// =====================================================================================
// With m measurements per subject, compound symmetry (correlation ρ) and the
// Greenhouse–Geisser sphericity factor ε, the F-test has noncentrality
// λ = ε·n·m·f²/(1 − ρ) and (m − 1)ε, (n − 1)(m − 1)ε degrees of freedom.
double rm_power(double n, double m, double f, double rho, double eps, double alpha) {
  using namespace stats;
  const double df1 = (m - 1) * eps, df2 = (n - 1) * (m - 1) * eps;
  if (!(df2 > 0)) return NAN;
  return pnf(qf(alpha, df1, df2, false), df1, df2, eps * n * m * f * f / (1 - rho), false);
}

Procedure rm_anova() {
  Procedure p;
  p.id = "means.rm.anova";
  p.name = "Repeated-Measures ANOVA (One Within-Subject Factor)";
  p.category = "Means";
  p.subcategory = "Repeated Measures";
  p.order = 2;
  p.summary =
      "Number of subjects for the F-test that m repeated measurements (conditions or time points) have equal means, "
      "assuming a common correlation ρ between measurements and an optional Greenhouse–Geisser correction ε for "
      "departures from sphericity.";
  p.keywords = {"repeated measures", "within-subject", "rm anova", "longitudinal", "sphericity",
                "greenhouse-geisser", "within-subjects anova", "time points", "conditions"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Subjects)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      power_param(),
      alpha_param(),
      n_param("n", "N (Number of Subjects)", "20", 2)
          .sym("N")
          .hint("Number of subjects; every subject is measured under all m conditions or at all m time points."),
      dropout_param(),
      select("estype", "Effect Size Input",
             {{"means", "Means of the measurements and σ", {}}, {"f", "Cohen's f and number of measurements", {}}},
             "means")
          .grp(G_EFFECT)
          .hint("Describe the effect by the means you expect for the measurements and their standard deviation, "
                "or directly by Cohen's f and the number of measurements."),
      vec("means", "Means of the m Measurements", "0 0.2 0.4")
          .grp(G_EFFECT)
          .sym("Means")
          .sugg({"0 0.2 0.4", "0 0 0.5", "10 11 12 12"})
          .hint("The mean you expect at each condition or time point, one value per measurement; the number of "
                "values sets m.")
          .show(when(is("estype", {"means"}))),
      num("sigma", "σ (Standard Deviation of Each Measurement)", "1")
          .grp(G_EFFECT)
          .sym("σ")
          .above(0)
          .hint("Standard deviation of a single measurement among subjects, assumed the same for every "
                "measurement.")
          .show(when(is("estype", {"means"}))),
      num("m", "m (Number of Measurements)", "3")
          .grp(G_EFFECT)
          .sym("m")
          .at_least(2)
          .whole()
          .hint("How many times each subject is measured (conditions or time points).")
          .show(when(is("estype", {"f"}))),
      num("f", "f (Effect Size)", "0.25")
          .grp(G_EFFECT)
          .sym("f")
          .above(0)
          .sugg({"0.1", "0.25", "0.4"})
          .hint("Cohen's f: the SD of the m means divided by the SD of one measurement. By Cohen's conventions "
                "0.10 is small, 0.25 medium and 0.40 large.")
          .show(when(is("estype", {"f"}))),
      num("rho", "ρ (Correlation Between Measurements)", "0.5")
          .grp(G_EFFECT)
          .sym("ρ")
          .range(-1, 1, true, true)
          .sugg({"0.5", "0.3 0.5 0.7"})
          .hint("Correlation between two measurements of the same subject, assumed equal for every pair; a "
                "higher correlation gives more power."),
      num("eps", "ε (Sphericity Correction, 1 = Sphericity)", "1")
          .grp(G_EFFECT)
          .sym("ε")
          .range(0, 1, true, false)
          .sugg({"1", "0.75", "0.5"})
          .hint("Greenhouse–Geisser correction for unequal variances or correlations: 1 under sphericity, at least "
                "1/(m − 1). Below 1, the degrees of freedom and noncentrality are both scaled by ε, an approximation "
                "that can overstate the corrected test's power by a few points."),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    double m, f;
    if (s.is("estype", "f")) {
      m = s("m");
      f = s("f");
    } else {
      const auto& mu = s.v("means");
      m = static_cast<double>(mu.size());
      if (m < 2) throw CalcError("enter at least two means");
      // from the differences from the first mean, as in the two-way ANOVA
      double mean = 0, ss = 0;
      for (double v : mu) mean += (v - mu[0]) / m;
      for (double v : mu) ss += (v - mu[0] - mean) * (v - mu[0] - mean);
      f = std::sqrt(ss / m) / s("sigma");
      if (!(f > 0)) throw CalcError("the means are all equal; there is no effect to detect");
    }
    const double rho = s("rho"), eps = s("eps"), alpha = s("alpha");
    if (!(rho > -1 / (m - 1)))
      throw CalcError("with " + fint(m) + " measurements, equal correlations must exceed −1/(m − 1) = " +
                      fnum(-1 / (m - 1), 4) + " (otherwise they cannot all hold at once)");
    if (eps < 1 / (m - 1) - 1e-12) throw CalcError("ε cannot be smaller than 1/(m − 1) = " + fnum(1 / (m - 1), 4));
    double n;
    if (s.solve == "n") {
      solve::NSearch o;
      o.lo = 2;
      n = solve::min_n([&](double x) { return rm_power(x, m, f, rho, eps, alpha); }, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    r.set("n", n);
    put_power(s, r, rm_power(n, m, f, rho, eps, alpha));
    r.set("m", m);
    r.set("f", f);
    r.set("rho", rho);
    r.set("eps", eps);
    r.set("lambda", eps * n * m * f * f / (1 - rho));
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int", "Number of subjects"), col("m", "m", "int", "Measurements per subject"),
                   col("f", "f", "num", "Cohen's effect size (SD of the means / σ)"), col("rho", "ρ", "num"),
                   col("eps", "ε", "num"), col("lambda", "λ", "num", "Noncentrality parameter"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A sample of " + fint(r.num("n")) + " subjects, each measured " + fint(r.num("m")) + " times, " +
           achieves_singular(s, r) + " to detect differences among the means with effect size f = " +
           fnum(r.num("f"), 4) + ", assuming a correlation of " + fnum_fine(r.num("rho")) +
           " between measurements and a sphericity correction ε = " + fnum(r.num("eps")) +
           ", using the repeated-measures F-test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"f", "SD of the m means (divisor m) divided by the SD of a single measurement."},
                   {"ρ", "Correlation between two measurements of the same subject (compound symmetry)."},
                   {"ε", "Greenhouse–Geisser factor applied to the degrees of freedom and the noncentrality."}};
  p.references = {REF_COHEN,
                  "Muller, K. E. & Barton, C. N. (1989). Approximate power for repeated-measures ANOVA lacking "
                  "sphericity. Journal of the American Statistical Association 84, 549–555."};
  p.curve_param = "n";
  p.curve_label = "N (Subjects)";
  p.curve_min = 2;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"means", "0 0.2 0.4"}, {"sigma", "1"}, {"rho", "0.3 0.6"}};
  return p;
}
SSC_REGISTER(rm_anova);

// =====================================================================================
// Difference between two mean slopes in a longitudinal design
// =====================================================================================
// Every subject is measured at the same times t_1, ..., t_m.  With random
// subject slopes (SD σb) and independent residuals (SD σe), the least-squares
// slope of one subject has variance σb² + σe²/Sxx, Sxx = Σ(t_j − t̄)², whatever
// the variance of the intercepts.  The groups are compared with a two-sample
// t-test on these slopes, which for complete, balanced data is the mixed-model
// (generalized least squares) analysis of the slope difference.
double slope_sxx(const Scenario& s) {
  if (s.is("timing", "custom")) {
    const auto& t = s.v("times");
    if (t.size() < 2) throw CalcError("enter at least two measurement times");
    double mean = 0;
    for (double x : t) mean += x;
    mean /= static_cast<double>(t.size());
    double sxx = 0;
    bool spread = false;
    for (double x : t) {
      sxx += (x - mean) * (x - mean);
      spread = spread || x != mean;
    }
    if (!spread) throw CalcError("the measurement times must not all be equal");
    // squares of differences below about 1e-154 underflow, and the slopes' variance σe²/Sxx would be infinite
    if (!(sxx >= stats::DBL_MIN_)) throw CalcError("the measurement times are too close together to compute with");
    return sxx;
  }
  // m equally spaced times from 0 to D
  const double m = s("m"), d = s("duration");
  const double sxx = d * d * m * (m + 1) / (12 * (m - 1));
  // D² underflows below about 1e-154, and the slopes' variance σe²/Sxx would be infinite
  if (!(sxx >= stats::DBL_MIN_))
    throw CalcError("D (the time from first to last measurement) is too small to compute with");
  return sxx;
}

double slope_sd(const Scenario& s) {
  const double sb = s("sdb"), se = s("sde");
  // in units of the larger SD when both are below about 1e-154, where their squares underflow to 0
  const double u = std::fmax(sb, se) < 1e-150 ? std::fmax(sb, se) : 1.0;
  return u * std::sqrt((sb / u) * (sb / u) + (se / u) * (se / u) / slope_sxx(s));
}

double slopes_power(const Scenario& s, double n1, double n2, double delta) {
  const double se = slope_sd(s) * std::sqrt(1 / n1 + 1 / n2);
  const Alt a = alt(s);
  return s.is("test", "z") ? z_power(delta / se, s("alpha"), a) : t_power(delta / se, n1 + n2 - 2, s("alpha"), a);
}

Procedure slopes() {
  Procedure p;
  p.id = "means.two.slopes";
  p.name = "Two Mean Slopes in a Longitudinal Design (Random Slopes)";
  p.category = "Means";
  p.subcategory = "Repeated Measures";
  p.order = 3;
  p.summary =
      "Sample size and power for comparing the mean rate of change (slope over time) between two groups when each "
      "subject is measured repeatedly, under a linear mixed model with random subject slopes. The test is the "
      "two-sample t-test on the subjects' least-squares slopes, with power from the noncentral t distribution.";
  p.keywords = {"longitudinal", "slope", "rate of change", "mixed model", "random slope", "growth curve",
                "repeated measures", "decline", "trajectory", "progression", "rate of decline"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total Sample Size)"},
                 {"power", "Power", "power", "Power"},
                 {"delta", "Effect Size (δ)", "delta", "Detectable δ"}};
  p.params = {alternative("Two-sided (H1: β1 ≠ β2)", "One-sided (H1: β1 > β2)", "One-sided (H1: β1 < β2)"),
              select("test", "Test Statistic",
                     {{"t", "T-test (variances estimated)", {}}, {"z", "Z-test (variances known)", {}}}, "t")
                  .grp(G_TEST)
                  .hint("The t-test estimates the variability from the data (the usual choice); the z-test "
                        "treats the SDs as known."),
              power_param(), alpha_param()};
  for (auto& q : two_group_params(2)) p.params.push_back(q);
  p.params.push_back(dropout_param());
  for (auto q :
       {num("delta", "δ (Difference in Mean Slopes, β1 − β2)", "0.5")
            .grp(G_EFFECT)
            .sym("δ")
            .sugg({"0.5", "0.25 0.5 1"})
            .hint("Difference between the groups' mean changes per time unit.")
            .show(when(not_solving({"delta"}))),
        num("sdb", "σb (SD of the Subjects' True Slopes)", "0.5")
            .grp(G_EFFECT)
            .sym("σb")
            .at_least(0)
            .sugg({"0.5", "0 0.5 1"})
            .hint("Between-subject standard deviation of the slopes (0 when all subjects share their group's slope)."),
        num("sde", "σe (Residual SD)", "1")
            .grp(G_EFFECT)
            .sym("σe")
            .above(0)
            .sugg({"1", "0.5 1 2"})
            .hint("Within-subject standard deviation of the measurements around the subject's line."),
        select("timing", "Spacing of the Measurements", {{"equal", "Equally spaced", {}}, {"custom", "Custom times", {}}}, "equal")
            .grp("Measurement Schedule")
            .hint("Measure every subject M times, equally spaced from time 0 to D, or list the measurement times "
                  "yourself."),
        num("m", "M (Measurements per Subject)", "5")
            .grp("Measurement Schedule")
            .sym("M")
            .at_least(2)
            .whole()
            .sugg({"5", "3 5 7"})
            .hint("How many times each subject is measured, counting the first (baseline) measurement.")
            .show(when(is("timing", {"equal"}))),
        num("duration", "D (Time from First to Last Measurement)", "4")
            .grp("Measurement Schedule")
            .sym("D")
            .above(0)
            .sugg({"4", "12", "2 4"})
            .hint("Length of follow-up from the first to the last measurement, in the time unit of the slopes "
                  "(e.g. years if δ is a change per year).")
            .show(when(is("timing", {"equal"}))),
        vec("times", "Measurement Times (e.g. 0 1 2 4 8)", "0 1 2 4 8")
            .grp("Measurement Schedule")
            .sym("Times")
            .hint("The times at which every subject is measured, in the time unit of the slopes.")
            .show(when(is("timing", {"custom"})))})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    const double sd = slope_sd(s);
    TwoN n;
    double delta;
    if (s.solve == "n") {
      delta = s("delta");
      if (delta == 0) throw CalcError("δ is zero; power cannot exceed alpha");
      if ((a == Alt::Upper && delta < 0) || (a == Alt::Lower && delta > 0))
        throw CalcError("δ has the sign opposite to the one-sided alternative");
      const double guess = z_factor(s("alpha"), s("power"), a) * 2 * sd * sd / (delta * delta) + 1;
      n = solve_two_group_n(s, [&](double n1, double n2) { return slopes_power(s, n1, n2, delta); }, s("power"), 2,
                            guess, true, ctx);
    } else {
      n = two_group_n(s);
      if (n.n1 < 2 || n.n2 < 2) throw CalcError("each group needs at least 2 subjects");
      if (s.solve == "power") {
        delta = s("delta");
      } else {
        const double sign = a == Alt::Lower ? -1 : 1;
        delta = sign * solve::solve_increasing([&](double d) { return slopes_power(s, n.n1, n.n2, sign * d); },
                                               s("power"), 0, sd, std::fmax(1e9 * sd, solve::BOUND_SPAN), 1e-12,
                                               "the slope difference");
      }
    }
    put_two_group(r, n);
    put_power(s, r, slopes_power(s, n.n1, n.n2, delta));
    r.set("delta", delta);
    r.set("sdb", s("sdb"));
    r.set("sde", s("sde"));
    if (s.is("timing", "custom")) {
      r.set("m", static_cast<double>(s.v("times").size()));
      r.set("times", list_text(s.v("times")));
    } else {
      r.set("m", s("m"));
      r.set("duration", s("duration"));
    }
    r.set("sds", sd);
    r.set("es", delta / sd);
    r.set("alpha", s("alpha"));
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    c.push_back(col("delta", "δ", "num", "Difference in mean slopes β1 − β2"));
    c.push_back(col("sdb", "σb", "num", "SD of the subjects' true slopes"));
    c.push_back(col("sde", "σe", "num", "Residual SD"));
    c.push_back(col("m", "M", "int", "Measurements per subject"));
    if (s.is("timing", "custom"))
      c.push_back(col("times", "Times", "text", "Measurement times"));
    else
      c.push_back(col("duration", "D", "num", "Time from the first to the last measurement"));
    c.push_back(col("sds", "σs", "num", "SD of one subject's estimated slope: √(σb² + σe²/Sxx)"));
    c.push_back(col("es", "δ/σs", "num", "Standardized effect size"));
    c.push_back(col("alpha", "Alpha", "alpha", "Significance level"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string schedule = s.is("timing", "custom")
                                     ? "at times " + text_of(r, "times")
                                     : fint(r.num("m")) + " times over " + fnum(r.num("duration")) + " time unit(s)";
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " (total " +
           fint(r.num("n")) + ") " + achieves(s, r) + " to detect a difference of " + fnum(r.num("delta")) +
           " between the groups' mean slopes when every subject is measured " + schedule +
           ", the SD of the subjects' slopes is " + fnum(r.num("sdb")) + " and the residual SD is " +
           fnum(r.num("sde")) + ", using a " + alt_words(alt(s)) + " " + (s.is("test", "z") ? "z-test" : "t-test") +
           " on the subjects' slopes at significance level α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"N1, N2", "Number of subjects in groups 1 and 2; N = N1 + N2."},
                   {"δ", "Difference β1 − β2 between the groups' mean slopes (change per time unit)."},
                   {"σb", "Standard deviation of the subjects' true slopes within a group."},
                   {"σe", "Residual standard deviation of a measurement around the subject's line."},
                   {"M", "Number of measurements per subject."},
                   {"D", "Time from the first to the last measurement (equally spaced times)."},
                   {"σs", "SD of one subject's least-squares slope, √(σb² + σe²/Sxx), with Sxx = Σ(t − t̄)²; for M "
                          "equally spaced times over D, Sxx = D²·M(M + 1)/(12(M − 1))."}};
  p.notes = {"The calculation assumes that every subject completes every measurement; inflate the sample size for "
             "dropout with the dropout rate, or use the times that most subjects will complete."};
  p.references = {"Fitzmaurice, G. M., Laird, N. M. & Ware, J. H. (2011). Applied Longitudinal Analysis, 2nd ed. "
                  "Wiley.",
                  "Diggle, P. J., Heagerty, P., Liang, K.-Y. & Zeger, S. L. (2002). Analysis of Longitudinal Data, "
                  "2nd ed. Oxford University Press.",
                  "Ahn, C., Heo, M. & Zhang, S. (2015). Sample Size Calculations for Clustered and Longitudinal "
                  "Outcomes in Clinical Research. Chapman & Hall/CRC."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.curve_min = 2;
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"delta", "0.25 0.5"}, {"sdb", "0.5"}, {"sde", "1"}, {"m", "5"}, {"duration", "4"}};
  return p;
}
SSC_REGISTER(slopes);

}  // namespace
}  // namespace ssc::procs
