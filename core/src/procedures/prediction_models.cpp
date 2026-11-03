// Minimum sample size for developing a clinical prediction model (Riley et al. 2019, 2020): the largest of the
// sizes that keep the expected shrinkage of the predictor effects small, the optimism in the model's R² small and
// the overall risk (or, for a continuous outcome, the residual SD and the intercept) precisely estimated.  The
// calculation follows the pmsampsize software, including its roundings, so that its results can be reproduced.
#include <algorithm>
#include <cctype>

#include "procedures/common.hpp"
#include "procedures/derive.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// The Cox–Snell R² of a model with C-statistic `c` for an outcome with prevalence `phi`, as Riley et al. (2022)
// approximate it: a linear predictor that is normal with variance 1 among those with and without the outcome, its
// means √2·Φ⁻¹(C) apart, to which the logistic model fits exactly.  pmsampsize simulates a million people; this is
// the value in an infinitely large sample, R²CS = 1 − exp(−2I) with I the mutual information of the outcome and the
// predictor (the expected log-likelihood ratio per person).
double cstat_r2cs(double c, double phi) {
  using namespace stats;
  const double mu = std::sqrt(2.0) * qnorm(c), lp = std::log(phi), lq = std::log1p(-phi);
  // ln(e^a + e^b) without overflow
  auto lse = [](double a, double b) { return std::max(a, b) + std::log1p(std::exp(-std::fabs(a - b))); };
  // I = φ·E1[ln f1/f] + (1 − φ)·E0[ln f0/f], f = φ·f1 + (1 − φ)·f0, with ln f1/f0 = μx − μ²/2
  auto ln_f1_over_f = [&](double x) { return -lse(lp, lq - mu * x + mu * mu / 2); };
  auto ln_f0_over_f = [&](double x) { return -lse(lp + mu * x - mu * mu / 2, lq); };
  const GLRule& g = gauss_legendre(64);
  double e1 = 0, e0 = 0;
  for (int panel = 0; panel < 6; ++panel) {  // z from −12 to 12 in panels of 4
    const double a = -12 + 4.0 * panel;
    for (std::size_t i = 0; i < g.x.size(); ++i) {
      const double z = a + 2 * (g.x[i] + 1), w = 2 * g.w[i] * dnorm(z);
      e1 += w * ln_f1_over_f(mu + z);
      e0 += w * ln_f0_over_f(z);
    }
  }
  return -std::expm1(-2 * (phi * e1 + (1 - phi) * e0));
}

double round_to(double x, int decimals) {
  const double f = std::pow(10.0, decimals);
  return std::round(x * f) / f;
}

// Riley et al.'s criterion (i) and (ii) size: the n at which the expected uniform shrinkage is S, for p parameters.
double shrinkage_n(double p, double r2cs, double S) { return std::ceil(p / ((S - 1) * std::log(1 - r2cs / S))); }

// The shrinkage that keeps the optimism in the Nagelkerke R² at δ (criterion (ii) of binary and time-to-event outcomes).
double optimism_shrinkage(double r2cs, double delta, double max_r2) { return r2cs / (r2cs + delta * max_r2); }

// Criterion (iii) of a binary outcome: the overall risk φ within the margin m (1.96, as Riley et al.).
double risk_n(double phi, double m) { return std::ceil(std::pow(1.96 / m, 2) * phi * (1 - phi)); }

// The largest Cox–Snell R² possible for a binary or time-to-event outcome, and the one used to convert a Nagelkerke
// R² (for a time-to-event outcome, also the expected events among the 10,000 people of the first).
struct MaxR2 {
  double max_r2, max_conv, events;
};
MaxR2 max_r2cs(const Scenario& s, bool binary) {
  if (binary) {
    const double phi = s("prev");
    const double m = -std::expm1(2 * (phi * std::log(phi) + (1 - phi) * std::log1p(-phi)));
    return {m, m, NAN};
  }
  // as pmsampsize: the null log-likelihood of an exponential model from 10,000 people with the expected events
  // rounded up (and, to convert a Nagelkerke R², per person)
  const double e = std::ceil(s("rate") * s("fup") * 10000);
  const double e1 = s("rate") * s("fup");
  return {-std::expm1(2 * (e * std::log(e / 10000) - e) / 10000), -std::expm1(2 * (e1 * std::log(e1) - e1)), e};
}

// The anticipated Cox–Snell R²: as pmsampsize, a converted one is rounded, to 3 decimals from a Nagelkerke R² and
// to 4 from a C-statistic.
double anticipated_r2cs(const Scenario& s, double max_conv) {
  if (s.is("r2in", "nag")) return round_to(s("r2nag") * max_conv, 3);
  if (s.is("r2in", "c")) return round_to(cstat_r2cs(s("cstat"), s("prev")), 4);
  return s("r2cs");
}

// Van Houwelingen and Le Cessie's expected shrinkage at m subjects of a linear model with p parameters and adjusted
// R² r2, from the apparent R² that m subjects would give.
double vh_shrinkage(double r2, double P, double m) {
  const double app = (r2 * (m - P - 1) + P) / (m - 1);
  return 1 + (P - 2) / (m * std::log(1 - app));
}

// The two ends of the 95% interval for a residual SD with d degrees of freedom, as multiples of its estimate.
double sd_ratio_low(double d) { return std::sqrt(d / stats::qchisq(0.025, d)); }
double sd_ratio_high(double d) { return std::sqrt(stats::qchisq(0.975, d) / d); }

// The half-width of the 95% interval for the intercept at m subjects.
double intercept_half(double m, double P, double sd, double r2) {
  return stats::qt(0.025, m - P - 1, false) * std::sqrt(sd * sd * (1 - r2) / m);
}

const char* REF_RILEY_I =
    "Riley, R. D., Snell, K. I. E., Ensor, J., Burke, D. L., Harrell, F. E., Moons, K. G. M. & Collins, G. S. "
    "(2019). Minimum sample size for developing a multivariable prediction model: Part I – Continuous outcomes. "
    "Statistics in Medicine 38, 1262–1275.";
const char* REF_RILEY_II =
    "Riley, R. D., Snell, K. I. E., Ensor, J., Burke, D. L., Harrell, F. E., Moons, K. G. M. & Collins, G. S. "
    "(2019). Minimum sample size for developing a multivariable prediction model: Part II – Binary and "
    "time-to-event outcomes. Statistics in Medicine 38, 1276–1296.";
const char* REF_RILEY_BMJ =
    "Riley, R. D., Ensor, J., Snell, K. I. E., Harrell, F. E., Martin, G. P., Reitsma, J. B., Moons, K. G. M., "
    "Collins, G. & van Smeden, M. (2020). Calculating the sample size required for developing a clinical prediction "
    "model. BMJ 368, m441.";
const char* REF_RILEY_C =
    "Riley, R. D., Van Calster, B. & Collins, G. S. (2021). A note on estimating the Cox-Snell R² from a reported C "
    "statistic (AUROC) to inform sample size calculations for developing a prediction model with a binary outcome. "
    "Statistics in Medicine 40, 859–864.";
const char* REF_PMSAMPSIZE =
    "Ensor, J., Martin, E. C. & Riley, R. D. pmsampsize: Calculates the Minimum Sample Size Required for Developing a "
    "Multivariable Prediction Model. R package and Stata module.";

Procedure prediction_model() {
  Procedure p;
  p.id = "reg.prediction";
  p.name = "Developing a Clinical Prediction Model (Minimum Sample Size)";
  p.category = "Regression";
  p.subcategory = "Prediction Models";
  p.order = 1;
  p.summary =
      "Minimum sample size for developing a multivariable prediction model with a binary, time-to-event or "
      "continuous outcome (Riley et al.): the largest of the sizes that keep the expected shrinkage of the predictor "
      "effects small, keep the optimism in R² small and estimate the overall risk — or the residual SD and the "
      "intercept — precisely. The anticipated model performance is given as a Cox–Snell or Nagelkerke R², a "
      "C-statistic or, for a continuous outcome, an adjusted R².";
  p.keywords = {"prediction model", "clinical prediction model", "prognostic model", "risk model", "riley",
                "pmsampsize", "events per variable", "epv", "shrinkage", "overfitting", "logistic regression",
                "cox regression", "model development", "r-squared"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}};
  p.params = {
      select("outcome", "Outcome",
             {{"binary", "Binary (logistic regression)", {}},
              {"survival", "Time to event (Cox or parametric survival)", {}},
              {"continuous", "Continuous (linear regression)", {}}},
             "binary")
          .grp(G_TEST)
          .hint("The type of outcome the model will predict, which decides the criteria and the inputs."),
      num("p", "Number of Candidate Predictor Parameters", "24")
          .grp(G_EFFECT)
          .sym("Parameters")
          .at_least(1)
          .whole()
          .sugg({"24", "10 20 30"})
          .hint("Every regression coefficient the model may estimate before any selection: one per continuous "
                "predictor (more with splines or interactions) and one per category beyond the first."),
      select("r2in", "Anticipated Performance",
             {{"cs", "Cox–Snell R²", {}},
              {"nag", "Nagelkerke R²", {}},
              {"c", "C-statistic (area under the ROC curve)", when(is("outcome", {"binary"}))}},
             "cs")
          .grp(G_EFFECT)
          .hint("How the model's expected performance is given, for example from a previous model in the same "
                "setting.")
          .show(when(is("outcome", {"binary", "survival"}))),
      num("r2cs", "R²CS (Anticipated Cox–Snell R²)", "0.288")
          .grp(G_EFFECT)
          .sym("R²CS")
          .range(0, 1, true, true)
          .sugg({"0.288", "0.05 0.1 0.2"})
          .hint("The Cox–Snell R² the model is expected to reach, below its largest possible value for the outcome "
                "(shown in the results). Riley et al. suggest 15% of that largest value when nothing better is "
                "known.")
          .show(when({is("outcome", {"binary", "survival"}), is("r2in", {"cs"})})),
      num("r2nag", "R²Nag (Anticipated Nagelkerke R²)", "0.15")
          .grp(G_EFFECT)
          .sym("R²Nag")
          .range(0, 1, true, true)
          .sugg({"0.15", "0.1 0.2 0.3"})
          .hint("The Nagelkerke R² the model is expected to reach: the Cox–Snell R² as a share of its largest "
                "possible value.")
          .show(when({is("outcome", {"binary", "survival"}), is("r2in", {"nag"})})),
      num("cstat", "C-Statistic (Anticipated AUC)", "0.8")
          .grp(G_EFFECT)
          .sym("C")
          .range(0.5, 1, true, true)
          .sugg({"0.8", "0.7 0.75 0.8"})
          .hint("The C-statistic (area under the ROC curve) the model is expected to reach; it is converted to a "
                "Cox–Snell R² given the prevalence.")
          .show(when({is("outcome", {"binary"}), is("r2in", {"c"})})),
      num("r2", "R² (Anticipated Adjusted R²)", "0.2")
          .grp(G_EFFECT)
          .sym("R²")
          .range(0, 1, true, true)
          .sugg({"0.2", "0.1 0.2 0.3"})
          .hint("The share of the outcome's variance the model is expected to explain (adjusted R²).")
          .show(when(is("outcome", {"continuous"}))),
      num("prev", "Outcome Prevalence", "0.174")
          .grp(G_EFFECT)
          .sym("Prevalence")
          .range(0, 1, true, true)
          .sugg({"0.174", "0.05 0.1 0.2"})
          .hint("The share of people expected to have the outcome, between 0 and 1.")
          .show(when(is("outcome", {"binary"}))),
      num("rate", "Overall Event Rate", "0.065")
          .grp(G_EFFECT)
          .sym("Rate")
          .above(0)
          .sugg({"0.065", "0.05 0.1"})
          .hint("Events per person per unit of time (the time unit of the follow-up), e.g. 0.065 per year.")
          .show(when(is("outcome", {"survival"}))),
      num("tpt", "Time Point of Interest", "2")
          .grp(G_EFFECT)
          .sym("Time Point")
          .above(0)
          .sugg({"2", "1 5 10"})
          .hint("The time at which the model's predicted risks matter, e.g. 2 years; used to report the precision "
                "of the overall risk.")
          .show(when(is("outcome", {"survival"}))),
      num("fup", "Mean Follow-Up Time", "2.07")
          .grp(G_EFFECT)
          .sym("Mean Follow-Up")
          .above(0)
          .sugg({"2.07", "1 2 5"})
          .hint("The average follow-up per person, in the same unit as the event rate.")
          .show(when(is("outcome", {"survival"}))),
      num("mean", "Mean Outcome (Intercept)", "1.9")
          .grp(G_EFFECT)
          .sym("Mean")
          .sugg({"1.9", "10 50"})
          .hint("The average outcome in the population, which the model's intercept estimates (with centered "
                "predictors).")
          .show(when(is("outcome", {"continuous"}))),
      num("sdy", "SD of the Outcome", "0.6")
          .grp(G_EFFECT)
          .sym("SD")
          .above(0)
          .sugg({"0.6", "1 5"})
          .hint("The standard deviation of the outcome in the population.")
          .show(when(is("outcome", {"continuous"}))),
      num("shrink", "S (Target Expected Shrinkage)", "0.9")
          .grp(G_OPTIONS)
          .sym("S")
          .range(0, 1, true, true)
          .sugg({"0.9", "0.85 0.9 0.95"})
          .hint("The expected uniform shrinkage the predictor effects may need, at least: 0.9 means overfitting "
                "inflates them by no more than about 10%."),
      num("delta", "δ (Largest Optimism in R²)", "0.05")
          .grp(G_OPTIONS)
          .sym("δ")
          .range(0, 1, true, true)
          .sugg({"0.05", "0.02 0.05"})
          .hint("The largest acceptable difference between the apparent and the adjusted R² (for binary and "
                "time-to-event outcomes, the Nagelkerke R²).")
          .show(when(is("outcome", {"binary", "continuous"}))),
      num("margin", "Margin of Error for the Overall Risk", "0.05")
          .grp(G_OPTIONS)
          .sym("Margin")
          .range(0, 0.5, true, false)
          .sugg({"0.05", "0.025 0.05"})
          .hint("The overall risk (the outcome's prevalence) must be estimated within this margin with 95% "
                "confidence.")
          .show(when(is("outcome", {"binary"}))),
      num("mmoe", "Multiplicative Margin of Error", "1.1")
          .grp(G_OPTIONS)
          .sym("MMOE")
          .range(1, 2, true, false)
          .sugg({"1.1", "1.05 1.1"})
          .hint("The residual SD must be estimated within this factor of its value (95% interval), and the "
                "intercept within this factor of the mean outcome.")
          .show(when(is("outcome", {"continuous"})))};
  p.run = [](const Scenario& s, Row& r, RunContext&) {
    const std::string& outcome = s.s("outcome");
    const double P = s("p"), S = s("shrink");
    Table tb;
    tb.title = "Criteria";
    tb.columns = {{"criterion", "Criterion"}, {"n", "N"}, {"shrinkage", "Expected Shrinkage"}};
    tb.formats = {"text", "int", "num"};
    auto add = [&](const std::string& what, double n, double shrinkage) {
      tb.rows.push_back({Cell(what), Cell(n), Cell(shrinkage)});
    };
    double n = 0;
    std::string binding;
    auto bind = [&](double size, const char* what) {
      if (size > n) n = size, binding = what;
    };
    if (outcome == "continuous") {
      const double r2 = s("r2"), sd = s("sdy"), mean = s("mean"), mmoe = s("mmoe"), delta = s("delta");
      // Van Houwelingen and Le Cessie's expected shrinkage at n, from the apparent R² that n subjects would give
      auto svh = [&](double m) { return vh_shrinkage(r2, P, m); };
      // (i) as pmsampsize: from p + 2 subjects on, the first n whose expected shrinkage reaches S.  It first falls
      // with n and then rises, so once below S at p + 2 subjects the condition holds from its first n on.
      double n1 = P + 2;
      if (!(svh(n1) > S))
        n1 = solve::first_whole([&](double m) { return svh(m) >= S; }, P + 3, 1e9,
                                "the target shrinkage cannot be reached with up to a billion subjects");
      // (ii) small optimism in the adjusted R²
      const double n2 = std::ceil(1 + P * (1 - r2) / delta);
      // (iii) the residual SD within the multiplicative margin: both ends of its 95% interval (234 + p for 1.1)
      const double df = solve::first_whole(
          [&](double d) { return sd_ratio_low(d) <= mmoe && sd_ratio_high(d) <= mmoe; },
          1, 1e9, "the residual SD cannot be estimated within this margin");
      const double n3 = df + P + 1;
      // (iv) the intercept within the margin of the mean outcome: from the largest size so far, as pmsampsize
      if (mean == 0) throw CalcError("the margin for the intercept is relative to the mean outcome, which cannot be 0");
      auto half = [&](double m) { return intercept_half(m, P, sd, r2); };
      const double n4 = solve::first_whole([&](double m) { return half(m) <= (mmoe - 1) * std::fabs(mean); },
                                           std::max({n1, n2, n3}), 1e9,
                                           "the intercept cannot be estimated within this margin");
      bind(n1, "Shrinkage");
      bind(n2, "Optimism in R²");
      bind(n3, "Residual SD");
      bind(n4, "Intercept");
      add("(i) Expected shrinkage ≥ " + fnum(S), n1, svh(n1));
      add("(ii) Optimism in R² ≤ " + fnum(delta), n2, svh(n2));
      add("(iii) Residual SD within a factor of " + fnum(mmoe), n3, svh(n3));
      add("(iv) Intercept within a factor of " + fnum(mmoe), n4, svh(n4));
      r.set("n_shrink", n1);
      r.set("n_opt", n2);
      r.set("n_sd", n3);
      r.set("n_int", n4);
      r.set("shrinkage", svh(n));
      r.set("spp", n / P);
      r.set("r2", r2);
      r.set("int_lo", mean - half(n));
      r.set("int_hi", mean + half(n));
    } else {
      const bool binary = outcome == "binary";
      // the largest Cox–Snell R² possible, and the one used to convert a Nagelkerke R²
      const MaxR2 largest = max_r2cs(s, binary);
      const double max_r2 = largest.max_r2;
      if (!(max_r2 > 0))
        throw CalcError("so many events per person leave no room for a Cox–Snell R² (its largest possible value is 0)");
      const double r2cs = anticipated_r2cs(s, largest.max_conv);
      if (!(r2cs > 0)) throw CalcError("the anticipated Cox–Snell R² rounds to 0: the model would explain nothing");
      if (!(r2cs < max_r2))
        throw CalcError("the anticipated Cox–Snell R² (" + fnum(r2cs) + ") must be below its largest possible value "
                        "for this outcome, " + fnum(max_r2, 3));
      if (!(r2cs < S)) throw CalcError("the target shrinkage S must exceed the anticipated Cox–Snell R²");
      // (i) expected shrinkage S; (ii) the shrinkage that keeps the optimism in the Nagelkerke R² at δ (0.05 for
      // time-to-event outcomes, as pmsampsize)
      const double delta = binary ? s("delta") : 0.05;
      const double S2 = optimism_shrinkage(r2cs, delta, max_r2);
      const double n1 = shrinkage_n(P, r2cs, S), n2 = shrinkage_n(P, r2cs, S2);
      bind(n1, "Shrinkage");
      bind(n2, "Optimism in R²");
      add("(i) Expected shrinkage ≥ " + fnum(S), n1, S);
      add("(ii) Optimism in Nagelkerke R² ≤ " + fnum(delta), n2, S2);
      const double final_s = std::max(S, S2);
      if (binary) {
        // (iii) the overall risk within the margin (1.96, as Riley et al.)
        const double phi = s("prev"), m = s("margin");
        const double n3 = risk_n(phi, m);
        bind(n3, "Overall risk");
        add("(iii) Overall risk within ± " + fnum(m), n3, final_s);
        r.set("n_risk", n3);
        r.set("events", n * phi);
        r.set("prev", phi);
      } else {
        // (iii) reported, not imposed (as pmsampsize): the 95% interval for the overall risk at the time point,
        // from the event rate estimated over the follow-up of n people (person-time rounded to 1 decimal)
        const double rate = s("rate"), t = s("tpt"), pt = round_to(s("fup") * n, 1);
        const double se = 1.96 * std::sqrt(rate / pt);
        r.set("risk", -std::expm1(-rate * t));
        r.set("risk_lo", -std::expm1(-(rate - se) * t));
        r.set("risk_hi", -std::expm1(-(rate + se) * t));
        r.set("events", n * rate * s("fup"));
        r.set("rate", rate);
      }
      r.set("n_shrink", n1);
      r.set("n_opt", n2);
      r.set("shrinkage", final_s);
      r.set("epp", r.num("events") / P);
      r.set("r2cs", r2cs);
      r.set("r2max", max_r2);
      r.set("r2nag", r2cs / max_r2);
    }
    r.set("n", n);
    r.set("p", P);
    r.set("binding", binding);
    r.tables.push_back(std::move(tb));
  };
  p.columns = [](const Scenario& s) {
    std::vector<Column> c = {col("n", "N", "int", "Minimum sample size: the largest of the criteria"),
                             col("binding", "Set By", "text", "The criterion that needs the largest sample"),
                             col("p", "Parameters", "int", "Candidate predictor parameters")};
    if (s.is("outcome", "continuous")) {
      for (auto x : {col("n_shrink", "N (i)", "int", "Expected shrinkage at least S"),
                     col("n_opt", "N (ii)", "int", "Small optimism in R²"),
                     col("n_sd", "N (iii)", "int", "Precise residual SD"), col("n_int", "N (iv)", "int", "Precise intercept"),
                     col("shrinkage", "Shrinkage", "num", "Expected shrinkage at N"),
                     col("spp", "SPP", "num", "Subjects per candidate parameter"), col("r2", "R²", "num"),
                     col("int_lo", "Intercept Lower", "num", "95% interval for the intercept at N"),
                     col("int_hi", "Intercept Upper", "num", "95% interval for the intercept at N")})
        c.push_back(x);
      return c;
    }
    for (auto x : {col("n_shrink", "N (i)", "int", "Expected shrinkage at least S"),
                   col("n_opt", "N (ii)", "int", "Small optimism in the Nagelkerke R²")})
      c.push_back(x);
    if (s.is("outcome", "binary")) c.push_back(col("n_risk", "N (iii)", "int", "Precise overall risk"));
    for (auto x : {col("shrinkage", "Shrinkage", "num", "Expected shrinkage at N"),
                   col("events", "Events", "num", "Expected number of outcome events"),
                   col("epp", "EPP", "num", "Events per candidate parameter"), col("r2cs", "R²CS", "num", "Cox–Snell R²"),
                   col("r2max", "Max R²CS", "num", "Largest possible Cox–Snell R²"),
                   col("r2nag", "R²Nag", "num", "Nagelkerke R²")})
      c.push_back(x);
    if (s.is("outcome", "binary")) {
      c.push_back(col("prev", "Prevalence", "num"));
    } else {
      for (auto x : {col("rate", "Event Rate", "num"), col("risk", "Risk at Time Point", "num"),
                     col("risk_lo", "Risk Lower", "num", "95% interval for the overall risk at N"),
                     col("risk_hi", "Risk Upper", "num", "95% interval for the overall risk at N")})
        c.push_back(x);
    }
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const std::string& outcome = s.s("outcome");
    const std::string head = "A minimum of " + fint(r.num("n")) + " participants";
    const std::string set_by = "; the " + text_of(r, "binding") + " criterion sets the size.";
    if (outcome == "continuous")
      return head + " (" + fnum(r.num("spp"), 1) + " per candidate parameter) is needed to develop a linear prediction "
             "model with " + fint(r.num("p")) + " candidate predictor parameters and an anticipated R² of " +
             fnum(r.num("r2")) + ", so that the expected shrinkage is at least " + fnum(s("shrink")) +
             ", the optimism in R² at most " + fnum(s("delta")) + ", and the residual SD and the intercept are "
             "estimated within a factor of " + fnum(s("mmoe")) + set_by;
    const bool binary = outcome == "binary";
    return head + " (" + fint(std::ceil(r.num("events"))) + " events, " + fnum(r.num("epp"), 2) +
           " per candidate parameter) is needed to develop a prediction model for a " +
           (binary ? "binary outcome with prevalence " + fnum(r.num("prev")) : "time-to-event outcome with an event "
                                                                                  "rate of " + fnum(r.num("rate"))) +
           ", " + fint(r.num("p")) + " candidate predictor parameters and an anticipated Cox–Snell R² of " +
           fnum(r.num("r2cs")) + " (Nagelkerke " + fnum(r.num("r2nag"), 3) + "), so that the expected shrinkage is "
           "at least " + fnum(s("shrink")) + " and the optimism in the Nagelkerke R² at most " +
           fnum(binary ? s("delta") : 0.05) +
           (binary ? " and the overall risk is estimated within ± " + fnum(s("margin"))
                   : ", which estimates the risk at time " + fnum(s("tpt")) + " as " + fnum(r.num("risk"), 3) + " (95% "
                     "interval " + fnum(r.num("risk_lo"), 3) + " to " + fnum(r.num("risk_hi"), 3) + ")") +
           set_by;
  };
  // The smallest N that meets every criterion: each criterion's own size, and the largest of them, which decides.
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const std::string& outcome = s.s("outcome");
    const double P = r.num("p"), S = s("shrink"), n = r.num("n");
    const std::string Pp = gvi(P), Ss = gv(S);
    std::string binding = text_of(r, "binding");
    if (!binding.empty()) binding[0] = static_cast<char>(std::tolower(static_cast<unsigned char>(binding[0])));
    const std::string decides = ": the " + binding + " criterion decides";
    Derivation out;
    out.answer = "N = " + dvi(n) + ", set by the " + binding + " criterion";
    if (outcome == "continuous") {
      const double r2 = r.num("r2"), sd = s("sdy"), mean = s("mean"), mmoe = s("mmoe"), delta = s("delta");
      const double n1 = r.num("n_shrink"), n2 = r.num("n_opt"), n3 = r.num("n_sd"), n4 = r.num("n_int");
      const std::string M = gv(mmoe);
      out.find = "the smallest N that meets every criterion: (i) an expected shrinkage of at least S = " + Ss +
                 ", (ii) an optimism in R² of at most δ = " + gv(delta) + ", (iii) the residual SD and (iv) the "
                 "intercept estimated within a factor MMOE = " + M + " (95% intervals)";
      // (i): P + 2 subjects when their expected shrinkage already exceeds S, else the first size from P + 3 that
      // reaches it
      std::string i_numbers;
      if (n1 == P + 2) {
        i_numbers = "S_VH at N = " + dvi(n1) + ": " + dv(vh_shrinkage(r2, P, n1)) + " > " + Ss;
      } else {
        const double below = vh_shrinkage(r2, P, n1 - 1);
        i_numbers = "S_VH at N = " + dvi(n1 - 1) + ": " + dvp(below, S) + (below < S ? " < " : " ≤ ") + Ss +
                    " ≤ S_VH at N = " + dvi(n1) + ": " + dvp(vh_shrinkage(r2, P, n1), S);
      }
      out.steps.push_back({"(i) Expected shrinkage",
                           "N = P + 2 if S_VH(P + 2) > S, else the smallest N from P + 3 with S_VH(N) ≥ S: S_VH(N) = "
                           "1 + (P − 2) / (N ln(1 − R²app)), R²app = (R²(N − P − 1) + P) / (N − 1) the apparent R²",
                           i_numbers});
      out.steps.push_back({"(ii) Optimism in R²", "N = ⌈1 + P(1 − R²) / δ⌉",
                           "⌈1 + " + Pp + " × (1 − " + gv(r2) + ") / " + gv(delta) + "⌉ = " + dvi(n2)});
      // (iii): the first degrees of freedom whose interval for the residual SD lies within the margin at both ends
      const double df = n3 - P - 1;
      const auto ends = [&](double d) { return dvp(sd_ratio_low(d), mmoe) + " and " + dvp(sd_ratio_high(d), mmoe); };
      out.steps.push_back({"(iii) Residual SD",
                           "N = ν + P + 1, ν the smallest whole number with √(ν / χ²(0.025; ν)) ≤ MMOE and "
                           "√(χ²(0.975; ν) / ν) ≤ MMOE, χ²(q; ν) the q quantile",
                           (df > 1 ? "at ν = " + dvi(df - 1) + ": " + ends(df - 1) + ", not both ≤ " + M + "; " : "") +
                               "at ν = " + dvi(df) + ": " + ends(df) + " ≤ " + M + "; N = " + dvi(df) + " + " + Pp +
                               " + 1 = " + dvi(n3)});
      // (iv): from the largest size so far, the first whose interval for the intercept is within the margin
      const double from = std::max({n1, n2, n3}), bound = (mmoe - 1) * std::fabs(mean);
      const auto half = [&](double m) {
        return "t(0.975; " + dvi(m - P - 1) + ") × √(" + gv(sd) + "² × (1 − " + gv(r2) + ") / " + dvi(m) + ") = " +
               dvp(intercept_half(m, P, sd, r2), bound);
      };
      out.steps.push_back({"(iv) Intercept",
                           "N = the smallest size from the largest of (i)–(iii) whose half-width t(0.975; N − P − 1) "
                           "√(SD²(1 − R²) / N) ≤ (MMOE − 1)|mean|",
                           "(MMOE − 1)|mean| = (" + M + " − 1) × |" + gv(mean) + "| = " + dv(bound) + "; " +
                               (n4 > from ? "at N = " + dvi(n4 - 1) + ": " + half(n4 - 1) + " > " + dv(bound) + "; "
                                          : "") +
                               "at N = " + dvi(n4) + ": " + half(n4) + " ≤ " + dv(bound)});
      out.steps.push_back({"Minimum sample size", "N = the largest of the criteria's sizes",
                           "max(" + dvi(n1) + ", " + dvi(n2) + ", " + dvi(n3) + ", " + dvi(n4) + ") = " + dvi(n) +
                               decides});
      out.answer += " (" + dv(r.num("spp"), 1) + " subjects per candidate parameter)";
      return out;
    }
    const bool binary = outcome == "binary";
    const MaxR2 largest = max_r2cs(s, binary);
    const double r2cs = r.num("r2cs"), max_r2 = largest.max_r2;
    // the Cox–Snell R² is given, or converted from the Nagelkerke R² or the C-statistic
    const std::string R2 = s.is("r2in", "cs") ? gv(r2cs) : dv(r2cs);
    // δ is held at 0.05 for a time-to-event outcome (as pmsampsize): a constant then, not an input
    const double delta = binary ? s("delta") : 0.05;
    const std::string D = binary ? gv(delta) : num(delta);
    out.find = "the smallest N that meets every criterion: (i) an expected shrinkage of at least S = " + Ss +
               ", (ii) an optimism in the Nagelkerke R² of at most δ = " + D +
               (binary ? ", (iii) the overall risk estimated within ± " + gv(s("margin")) + " (95% interval)" : "");
    if (binary) {
      const std::string phi = gv(s("prev"));
      out.steps.push_back({"Largest Cox–Snell R²",
                           "max R²CS = 1 − exp(2(φ ln φ + (1 − φ) ln(1 − φ))), φ the prevalence",
                           "1 − exp(2 × (" + phi + " × ln " + phi + " + (1 − " + phi + ") × ln(1 − " + phi +
                               "))) = " + dv(max_r2)});
    } else {
      const double e = largest.events;
      out.steps.push_back({"Largest Cox–Snell R²",
                           "max R²CS = 1 − exp(2(E ln(E / 10000) − E) / 10000), E = ⌈rate × follow-up × 10000⌉ the "
                           "events among 10,000 people (as pmsampsize)",
                           "E = ⌈" + gv(s("rate")) + " × " + gv(s("fup")) + " × 10000⌉ = " + dvi(e) +
                               "; 1 − exp(2 × (" + dvi(e) + " × ln(" + dvi(e) + " / 10000) − " + dvi(e) +
                               ") / 10000) = " + dv(max_r2)});
    }
    if (s.is("r2in", "nag")) {
      if (binary)
        out.steps.push_back({"Anticipated Cox–Snell R²",
                             "R²CS = R²Nag × max R²CS, rounded to 3 decimals (as pmsampsize)",
                             "round(" + gv(s("r2nag")) + " × " + dv(largest.max_conv) + ", 3) = " + dv(r2cs)});
      else
        out.steps.push_back({"Anticipated Cox–Snell R²",
                             "R²CS = R²Nag × (1 − exp(2(e ln e − e))), e = rate × follow-up the events per person, "
                             "rounded to 3 decimals (as pmsampsize)",
                             "e = " + gv(s("rate")) + " × " + gv(s("fup")) + " = " + dv(s("rate") * s("fup")) +
                                 "; round(" + gv(s("r2nag")) + " × " + dv(largest.max_conv) + ", 3) = " + dv(r2cs)});
    } else if (s.is("r2in", "c")) {
      out.steps.push_back({"Anticipated Cox–Snell R²",
                           "R²CS = 1 − exp(−2I), I the mutual information of the outcome and a linear predictor that "
                           "is normal with variance 1 among those with and without it, its means √2 Φ⁻¹(C) apart "
                           "(Riley et al.), rounded to 4 decimals (as pmsampsize)",
                           "√2 × Φ⁻¹(" + gv(s("cstat")) + ") = " + dv(std::sqrt(2.0) * stats::qnorm(s("cstat"))) +
                               "; round(" + dv(cstat_r2cs(s("cstat"), s("prev")), 8) + ", 4) = " + dv(r2cs)});
    }
    const double n1 = r.num("n_shrink"), n2 = r.num("n_opt");
    const double S2 = optimism_shrinkage(r2cs, delta, max_r2);
    const auto shrink_size = [&](const std::string& shrink) {
      return "⌈" + Pp + " / ((" + shrink + " − 1) × ln(1 − " + R2 + " / " + shrink + "))⌉";
    };
    out.steps.push_back({"(i) Expected shrinkage", "N = ⌈P / ((S − 1) ln(1 − R²CS / S))⌉",
                         shrink_size(Ss) + " = " + dvi(n1)});
    out.steps.push_back({"(ii) Optimism in R²",
                         std::string("S₂ = R²CS / (R²CS + δ max R²CS), the shrinkage that keeps the optimism in the "
                                     "Nagelkerke R² at δ") +
                             (binary ? "" : " = 0.05 (as pmsampsize)") + "; N = ⌈P / ((S₂ − 1) ln(1 − R²CS / S₂))⌉",
                         "S₂ = " + R2 + " / (" + R2 + " + " + D + " × " + dv(max_r2) + ") = " + dv(S2) + "; " +
                             shrink_size(dv(S2)) + " = " + dvi(n2)});
    std::string sizes = dvi(n1) + ", " + dvi(n2);
    if (binary) {
      const std::string phi = gv(s("prev"));
      const double n3 = r.num("n_risk");
      out.steps.push_back({"(iii) Overall risk", "N = ⌈(1.96 / margin)² φ(1 − φ)⌉",
                           "⌈(1.96 / " + gv(s("margin")) + ")² × " + phi + " × (1 − " + phi + ")⌉ = " + dvi(n3)});
      sizes += ", " + dvi(n3);
    }
    out.steps.push_back({"Minimum sample size", "N = the largest of the criteria's sizes",
                         "max(" + sizes + ") = " + dvi(n) + decides});
    out.answer += " (" + dvi(std::ceil(r.num("events"))) + " events, " + dv(r.num("epp"), 2) +
                  " per candidate parameter)";
    return out;
  };
  p.definitions = {
      {"Parameters", "Candidate predictor parameters: every coefficient considered, before any selection."},
      {"Shrinkage", "Expected uniform shrinkage factor of the predictor effects (1 = no overfitting)."},
      {"R²CS, R²Nag", "Cox–Snell R², and Nagelkerke R² (the Cox–Snell R² as a share of its largest possible value)."},
      {"EPP, SPP", "Events, or subjects, per candidate predictor parameter."}};
  p.notes = {
      "As in pmsampsize, a Cox–Snell R² converted from a Nagelkerke R² is rounded to three decimals and one converted "
      "from a C-statistic to four. The C-statistic is converted exactly in an infinitely large sample, where "
      "pmsampsize simulates a million people, so its rounded value can occasionally differ in the last decimal.",
      "For time-to-event outcomes the precision of the overall risk is reported, not imposed, and the optimism in "
      "the Nagelkerke R² is held at 0.05, as in pmsampsize.",
      "These are minimum sizes for developing a model; external validation has its own sample size criteria."};
  p.references = {REF_RILEY_BMJ, REF_RILEY_I, REF_RILEY_II, REF_RILEY_C, REF_PMSAMPSIZE,
                  "van Houwelingen, J. C. & Le Cessie, S. (1990). Predictive value of statistical models. Statistics "
                  "in Medicine 9, 1303–1325."};
  p.example_solve = "n";
  p.example = {{"outcome", "binary"}, {"p", "24"}, {"r2cs", "0.288"}, {"prev", "0.174"}};
  return p;
}
SSC_REGISTER(prediction_model);

}  // namespace
}  // namespace ssc::procs
