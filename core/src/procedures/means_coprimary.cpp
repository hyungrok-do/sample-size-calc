// Two groups compared on several co-primary continuous endpoints, all of which must be significant (Sozu et al.
// 2011): the power is the probability that every one-sided test rejects, from the endpoints' correlation.
#include <algorithm>

#include "procedures/common.hpp"
#include "procedures/sim_survival.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

constexpr int KMAX = 4;

std::vector<double> standardized_effects(const Scenario& s) {
  std::vector<double> d;
  if (s.is("es", "raw")) {
    const auto& diff = s.v("diffs");
    const auto& sd = s.v("sds");
    if (diff.size() != sd.size())
      throw CalcError("give one standard deviation per mean difference (" + std::to_string(diff.size()) + ")");
    for (std::size_t k = 0; k < diff.size(); ++k) d.push_back(diff[k] / sd[k]);
  } else {
    d = s.v("deltas");
  }
  if (d.size() < 2 || d.size() > KMAX) throw CalcError("give two to four endpoints (one effect per endpoint)");
  // the tests are one-sided, in favor of group 1, and no effect is negative (by the inputs' limits); with no effect
  // on an endpoint its test, and so the study, succeeds with probability alpha at most, which no sample size raises
  for (double x : d)
    if (!(x > 0) && s.solve == "n")
      throw CalcError("every endpoint's effect must be positive: with no effect on one, the power cannot exceed "
                      "alpha");
  return d;
}

// P(X_i > lo_i for every i) for standard normals with a common correlation rho ≥ 0: with X_i = √ρ·Z0 + √(1 − ρ)·Z_i
// it is the average over Z0 of Π Φ((√ρ·Z0 − lo_i)/√(1 − ρ)), a one-dimensional integral.  Each factor turns from 0
// to 1 around Z0 = lo_i/√ρ over a width of about √((1 − ρ)/ρ): the panels are cut there and kept that narrow.
double equicorrelated_upper(const double* lo, int K, double rho) {
  using namespace stats;
  if (rho == 0) {
    double p = 1;
    for (int i = 0; i < K; ++i) p *= pnorm(-lo[i]);
    return p;
  }
  const double a = std::sqrt(rho), b = std::sqrt(1 - rho);
  std::vector<double> cuts = {-12, 12};
  for (int i = 0; i < K; ++i)
    if (std::fabs(lo[i] / a) < 12) cuts.push_back(lo[i] / a);
  std::sort(cuts.begin(), cuts.end());
  const double width = std::fmax(0.006, std::fmin(1.0, 2 * b / a));
  const GLRule& g = gauss_legendre(32);
  double total = 0;
  for (std::size_t c = 0; c + 1 < cuts.size(); ++c) {
    const int panels = static_cast<int>(std::ceil((cuts[c + 1] - cuts[c]) / width - 1e-9));
    const double h = (cuts[c + 1] - cuts[c]) / std::max(panels, 1);
    for (int q = 0; q < panels; ++q)
      for (std::size_t j = 0; j < g.x.size(); ++j) {
        const double z = cuts[c] + h * (q + (g.x[j] + 1) / 2);
        double f = dnorm(z);
        for (int i = 0; i < K; ++i) f *= pnorm((a * z - lo[i]) / b);
        total += h / 2 * g.w[j] * f;
      }
  }
  return total;
}

// The probability that every test rejects, from the endpoints' one-sided powers: their statistics are normal with
// correlation rho (exactly so with known variances; with t statistics their exact powers are joined through the
// same correlation, a Gaussian copula).
double all_reject(const std::vector<double>& marginal, double rho) {
  const int K = static_cast<int>(marginal.size());
  double lo[KMAX], hi[KMAX], r[KMAX * KMAX];
  for (int i = 0; i < K; ++i) {
    lo[i] = stats::qnorm(marginal[i], false);  // X_i > lo_i has probability marginal[i]
    hi[i] = INFINITY;
    for (int j = 0; j < K; ++j) r[i * K + j] = i == j ? 1 : rho;
  }
  if (K == 2) return stats::clamp01(stats::bvn_upper(lo[0], lo[1], rho));
  if (rho >= 0) return stats::clamp01(equicorrelated_upper(lo, K, rho));
  return stats::clamp01(nph::mvn_probability(K, lo, hi, r));  // a negative common correlation: to about 1.5e-4
}

std::vector<double> marginal_powers(const std::vector<double>& d, double n1, double n2, double alpha, bool t) {
  std::vector<double> p;
  const double scale = std::sqrt(n1 * n2 / (n1 + n2)), df = n1 + n2 - 2;
  for (double x : d)
    p.push_back(t ? t_power(x * scale, df, alpha, Alt::Upper) : z_power(x * scale, alpha, Alt::Upper));
  return p;
}

Procedure coprimary() {
  Procedure p;
  p.id = "means.coprimary";
  p.name = "Two Means with Co-Primary Endpoints (All Must Be Significant)";
  p.category = "Means";
  p.subcategory = "Two Independent Means";
  p.order = 6;
  p.summary =
      "Sample size and power for comparing two groups on two to four correlated continuous endpoints that must all "
      "show an effect (co-primary endpoints), each tested one-sided at the same level without a multiplicity "
      "adjustment (Sozu et al.). The power that all tests succeed is lower than that of any one of them, by less "
      "the more the endpoints are correlated.";
  p.keywords = {"co-primary endpoints", "coprimary", "multiple endpoints", "multiple primary endpoints",
                "all endpoints", "must win", "correlated endpoints", "two means", "intersection-union", "sozu"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      select("test", "Test",
             {{"t", "T-tests (variances estimated)", {}}, {"z", "Z-tests (known variances)", {}}}, "t")
          .grp(G_TEST)
          .hint("The analysis of each endpoint. With known variances the power is exact; with t-tests each "
                "endpoint's power is exact and their joint behavior is approximated."),
      power_param(),
      alpha_param("0.025", "Alpha (One-Sided, Each Endpoint)").range(0, 0.5, true, true)};
  for (auto& q : two_group_params(2, "100")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  for (auto q :
       {select("es", "Effect Size Input",
               {{"std", "Standardized effects δk", {}}, {"raw", "Mean differences and standard deviations", {}}}, "std")
            .grp(G_EFFECT)
            .hint("Give each endpoint's effect as a standardized difference (difference ÷ SD), or as a mean "
                  "difference with its standard deviation."),
        vec("deltas", "δ (Standardized Effects, One per Endpoint)", "0.3 0.35")
            .grp(G_EFFECT)
            .sym("δ")
            .at_least(0)
            .hint("Each endpoint's mean difference between the groups divided by its standard deviation, positive "
                  "in favor of group 1; two to four values.")
            .show(when(is("es", {"std"}))),
        vec("diffs", "Mean Differences μ1 − μ2 (One per Endpoint)", "3 4")
            .grp(G_EFFECT)
            .sym("μ1 − μ2")
            .at_least(0)
            .hint("Each endpoint's expected mean difference, positive in favor of group 1; two to four values.")
            .show(when(is("es", {"raw"}))),
        vec("sds", "Standard Deviations (One per Endpoint)", "10 10")
            .grp(G_EFFECT)
            .sym("σ")
            .above(0)
            .hint("Each endpoint's standard deviation within a group, in the order of the mean differences.")
            .show(when(is("es", {"raw"}))),
        num("rho", "ρ (Correlation between Endpoints)", "0.5")
            .grp(G_EFFECT)
            .sym("ρ")
            .range(-1, 1, true, true)
            .sugg({"0.5", "0 0.3 0.5 0.8"})
            .hint("The correlation between any two endpoints within a subject; 0 if unknown gives the most "
                  "conservative sample size for positively related endpoints.")})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const std::vector<double> d = standardized_effects(s);
    const int K = static_cast<int>(d.size());
    const double rho = s("rho"), alpha = s("alpha");
    if (!(rho > -1.0 / (K - 1)))
      throw CalcError("with " + std::to_string(K) + " endpoints a common correlation must exceed " +
                      fnum(-1.0 / (K - 1), 3));
    const bool t = s.is("test", "t");
    auto power = [&](double n1, double n2) { return all_reject(marginal_powers(d, n1, n2, alpha, t), rho); };
    TwoN n;
    if (s.solve == "n") {
      double dmin = d[0];
      for (double x : d) dmin = std::min(dmin, x);
      n = solve_two_group_n(s, power, s("power"), 2, std::ceil(2 * z_factor(alpha, s("power"), Alt::Upper) / (dmin * dmin)),
                            true, ctx, 1e8);
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    const std::vector<double> each = marginal_powers(d, n.n1, n.n2, alpha, t);
    put_power(s, r, all_reject(each, rho));
    double pmin = 1;
    Table tb;
    tb.title = "Endpoints";
    tb.columns = {{"k", "Endpoint"}, {"delta", "δ"}, {"power", "Power Alone"}};
    tb.formats = {"int", "num", "prob"};
    for (int k = 0; k < K; ++k) {
      pmin = std::min(pmin, each[k]);
      tb.rows.push_back({Cell(static_cast<double>(k + 1)), Cell(d[k]), Cell(each[k])});
    }
    r.tables.push_back(std::move(tb));
    r.set("endpoints", K);
    r.set("pmin", pmin);
    r.set("rho", rho);
    r.set("alpha", alpha);
    r.set("deltas", list_text(d));
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("endpoints", "Endpoints", "int"), col("deltas", "δ", "text", "Standardized effects"),
                   col("rho", "ρ", "num", "Correlation between endpoints"),
                   col("pmin", "Lowest Power Alone", "prob", "The weakest endpoint's power on its own"),
                   col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
           " to show an effect on all " + fint(r.num("endpoints")) + " co-primary endpoints (standardized effects " +
           text_of(r, "deltas") + ", correlation " + fnum(r.num("rho")) + ") with one-sided " +
           (s.is("test", "t") ? "t-tests" : "z-tests") + " at α = " + fnum(r.num("alpha")) +
           " each; the weakest endpoint alone would have " + power_text(r.num("pmin")) + " power.";
  };
  p.definitions = {{"δ", "Standardized effect of an endpoint: (μ1 − μ2)/σ."},
                   {"ρ", "Correlation between two endpoints within a subject (the same for every pair)."},
                   {"Power", "Probability that every endpoint's test is significant."}};
  p.notes = {"No multiplicity adjustment is needed when all endpoints must succeed; the price is that the joint "
             "power is below each endpoint's own.",
             "With t-tests each endpoint's power is exact and they are joined through the correlation ρ (a Gaussian "
             "copula); in simulations with 10 to 30 subjects per group it was within 0.005 of the true joint power, "
             "and closer with more. With known variances the power is exact.",
             "With three or four endpoints and a negative correlation the multivariate normal probability is "
             "computed with a lattice rule, to within about 1.5e-4; otherwise it is exact."};
  p.references = {"Sozu, T., Sugimoto, T. & Hamasaki, T. (2011). Sample size determination in superiority clinical "
                  "trials with multiple co-primary correlated endpoints. Journal of Biopharmaceutical Statistics 21, "
                  "650–668.",
                  "Sugimoto, T., Sozu, T. & Hamasaki, T. (2012). A convenient formula for sample size calculations in "
                  "clinical trials with multiple co-primary continuous endpoints. Pharmaceutical Statistics 11, "
                  "118–128.",
                  "Offen, W. et al. (2007). Multiple co-primary endpoints: medical and statistical solutions. Drug "
                  "Information Journal 41, 31–46."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"deltas", "0.3 0.35"}, {"rho", "0 0.5 0.8"}};
  return p;
}
SSC_REGISTER(coprimary);

}  // namespace
}  // namespace ssc::procs
