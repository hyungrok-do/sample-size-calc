// Survival (logrank test with accrual and follow-up), rates (Poisson and
// negative binomial) and variances (chi-square and F tests).
#include "procedures/apart.hpp"
#include "procedures/common.hpp"

namespace ssc::procs {
namespace {

using namespace blocks;

// =====================================================================================
// Logrank test (exponential survival, uniform accrual, exponential dropout)
// =====================================================================================
// Probability that a subject entering uniformly during [0, A] and followed
// until A + F has an event, with event hazard lam and censoring hazard eta.
double event_prob(double lam, double eta, double A, double F) {
  const double h = lam + eta;
  if (!(h > 0) || !(A + F > 0)) return 0;
  // an infinite hazard (a median of about 0): each subject has the event, or drops out, at once
  if (std::isinf(lam) || std::isinf(eta)) return std::isinf(eta) ? (std::isinf(lam) ? NAN : 0) : 1;
  // two finite hazards whose sum overflows: the probability depends on hazard × time and on lam / h only, so it is
  // that of half of each over twice the time
  if (std::isinf(h)) return event_prob(lam / 2, eta / 2, 2 * A, 2 * F);
  const double frac = lam / h;
  if (A <= 0) return frac * -std::expm1(-h * F);
  // 1 − e^{−hF}(1 − e^{−hA})/(hA), written without cancellation for small hazards (a median of millions of
  // time units): −expm1(−hF) + e^{−hF}(x + expm1(−x))/x with x = hA, the last factor by its series when small
  // (and 1 once hA overflows)
  const double x = h * A;
  const double tail = x < 1e-3   ? x / 2 - x * x / 6 + x * x * x / 24 - x * x * x * x / 120
                      : x < 1e300 ? (x + std::expm1(-x)) / x
                                  : 1;
  return frac * (-std::expm1(-h * F) + std::exp(-h * F) * tail);
}

double control_hazard(const Scenario& s) {
  const std::string& t = s.s("survin");
  if (t == "hazard") return s("lam2");
  if (t == "median") return std::log(2.0) / s("med2");
  const double S = s("s2"), t0 = s("t0");
  if (!(S > 0 && S < 1)) throw CalcError("the survival proportion must be between 0 and 1");
  return -std::log(S) / t0;
}

struct LogrankDesign {
  double lam1, lam2, hr, p1, p2;  // hazards and event probabilities
};

LogrankDesign logrank_design(const Scenario& s, double hr) {
  LogrankDesign d;
  d.lam2 = control_hazard(s);
  d.hr = hr;
  d.lam1 = hr * d.lam2;
  const double eta = s("eta"), A = s("accrual"), F = s("follow");
  d.p1 = event_prob(d.lam1, eta, A, F);
  d.p2 = event_prob(d.lam2, eta, A, F);
  return d;
}

// Power from the expected number of events.
double logrank_power(const Scenario& s, double n1, double n2, const LogrankDesign& d) {
  const double alpha = s("alpha");
  const Alt a = alt(s);
  const double events = n1 * d.p1 + n2 * d.p2;
  double drift;
  if (s.is("method", "freedman")) {
    const double phi = n2 / n1;  // control-to-treatment ratio
    drift = std::sqrt(events * phi) * (d.hr - 1) / (d.hr + phi);
  } else {  // Schoenfeld
    const double w1 = n1 / (n1 + n2), w2 = 1 - w1;
    drift = std::sqrt(events * w1 * w2) * std::log(d.hr);
  }
  return z_power(drift, alpha, a);
}

// Help texts of the inputs shared by the two logrank procedures.
const char* const HELP_SURVIN =
    "How you describe survival in the control group, assumed exponential: by its median survival time, its hazard "
    "rate, or the share still event-free at a time T0.";
const char* const HELP_MED2 =
    "Median survival time of the control group (half of the subjects have had the event by then), in the time unit "
    "of accrual and follow-up.";
const char* const HELP_LAM2 =
    "Constant event rate of the control group per unit of time, in the time unit of accrual and follow-up; it equals "
    "ln 2 divided by the median survival.";
const char* const HELP_S2 = "The share of control subjects still event-free at time T0, between 0 and 1.";
const char* const HELP_T0 = "The time at which the control survival S2 is given, in the time unit of the study.";
const char* const HELP_ACCRUAL =
    "Length of the period over which subjects enter the study at an even rate, in the time unit of the study; 0 if "
    "all start together.";
const char* const HELP_FOLLOW =
    "Time from the end of accrual to the final analysis; the last subject to enter is followed for this long.";
const char* const HELP_ETA =
    "Rate per unit of time at which subjects are lost to follow-up for reasons other than the end of the study "
    "(a constant hazard); 0 for none.";

Procedure logrank() {
  Procedure p;
  p.id = "surv.logrank";
  p.name = "Logrank Tests (Exponential Survival, Accrual and Follow-Up)";
  p.category = "Survival";
  p.subcategory = "Two Survival Curves";
  p.order = 1;
  p.summary =
      "Sample size and power for comparing two survival curves with the logrank test, using Schoenfeld's or "
      "Freedman's formula for the number of events and exponential survival with uniform accrual, additional "
      "follow-up and exponential loss to follow-up to convert events into subjects.";
  p.keywords = {"logrank", "survival", "hazard ratio", "time to event", "schoenfeld", "freedman", "accrual",
                "follow-up", "median survival", "oncology", "kaplan-meier", "two survival curves", "exponential survival", "events"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"},
                 {"power", "Power", "power", "Power"},
                 {"hr", "Effect Size (HR)", "hr", "Detectable HR"}};
  p.params = {
      alternative("Two-sided (H1: HR ≠ 1)", "One-sided (H1: HR > 1)", "One-sided (H1: HR < 1)"),
      select("method", "Events Formula",
             {{"schoenfeld", "Schoenfeld (1981)", {}}, {"freedman", "Freedman (1982)", {}}}, "schoenfeld")
          .grp(G_TEST)
          .hint("The formula for the number of events needed: Schoenfeld's, based on the log hazard ratio, is the "
                "usual choice. With equal groups Freedman's asks for a few more events; with unequal groups the two "
                "can differ a lot either way."),
      power_param(),
      alpha_param()};
  for (auto& q : two_group_params(1, "100")) p.params.push_back(q);
  for (auto q :
       {num("hr", "HR (Hazard Ratio, Treatment/Control)", "0.7")
            .grp(G_EFFECT)
            .sym("HR")
            .above(0)
            .sugg({"0.7", "0.6 0.7 0.8", "1.5"})
            .hint("The treatment group's hazard (event rate) divided by the control group's; 0.7 means a 30% lower "
                  "risk of the event at any time.")
            .show(when(not_solving({"hr"}))),
        select("hrdir", "Search Direction for HR", {{"below", "HR < 1", {}}, {"above", "HR > 1", {}}}, "below")
            .grp(G_EFFECT)
            .hint("A two-sided test can detect a hazard ratio below or above 1; choose which of the two to find.")
            .show(when({solving({"hr"}), is("alternative", {"two"})})),
        select("survin", "Control Group Survival Input",
               {{"median", "Median survival time", {}}, {"hazard", "Hazard rate", {}}, {"surv", "Survival proportion at time T0", {}}},
               "median")
            .grp(G_EFFECT)
            .hint(HELP_SURVIN),
        num("med2", "Median Survival (Control)", "12")
            .grp(G_EFFECT)
            .sym("Median2")
            .above(0)
            .sugg({"12", "6 12 24"})
            .hint(HELP_MED2)
            .show(when(is("survin", {"median"}))),
        num("lam2", "λ2 (Control Hazard Rate)", "0.0578")
            .grp(G_EFFECT)
            .sym("λ2")
            .above(0)
            .hint(HELP_LAM2)
            .show(when(is("survin", {"hazard"}))),
        num("s2", "S2 (Control Survival at T0)", "0.5")
            .grp(G_EFFECT)
            .sym("S2")
            .range(0, 1, true, true)
            .hint(HELP_S2)
            .show(when(is("survin", {"surv"}))),
        num("t0", "T0 (Time for S2)", "12")
            .grp(G_EFFECT)
            .sym("T0")
            .above(0)
            .hint(HELP_T0)
            .show(when(is("survin", {"surv"}))),
        num("accrual", "Accrual Time (Uniform Entry)", "12")
            .grp("Study Duration")
            .sym("Accrual")
            .at_least(0)
            .sugg({"12", "24"})
            .hint(HELP_ACCRUAL),
        num("follow", "Additional Follow-Up after Accrual", "12")
            .grp("Study Duration")
            .sym("Follow-Up")
            .at_least(0)
            .sugg({"12", "6 12 24"})
            .hint(HELP_FOLLOW),
        num("eta", "Loss to Follow-Up Hazard Rate", "0")
            .grp("Study Duration")
            .sym("η")
            .at_least(0)
            .sugg({"0", "0.01"})
            .hint(HELP_ETA)})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const Alt a = alt(s);
    double hr;
    TwoN n;
    if (s.solve == "n") {
      hr = s("hr");
      if (hr == 1) throw CalcError("HR = 1: there is no effect to detect");
      if ((a == Alt::Upper && hr < 1) || (a == Alt::Lower && hr > 1))
        throw CalcError("the hazard ratio is on the side opposite to the one-sided alternative");
      const auto d = logrank_design(s, hr);
      if (!(d.p1 > 0 || d.p2 > 0)) throw CalcError("no events are expected with this accrual and follow-up");
      const double ev = z_factor(s("alpha"), s("power"), a) / (0.25 * std::log(hr) * std::log(hr));
      const double guess = ev / ((d.p1 + d.p2) / 2) / 2;
      n = solve_two_group_n(s, [&](double n1, double n2) { return logrank_power(s, n1, n2, d); }, s("power"), 1,
                            std::ceil(guess), true, ctx);
    } else {
      n = two_group_n(s);
      if (s.solve == "power") {
        hr = s("hr");
      } else {
        const bool below = a == Alt::Lower || (a == Alt::Two && !s.is("hrdir", "above"));
        auto f = [&](double x) { return logrank_power(s, n.n1, n.n2, logrank_design(s, x)) - s("power"); };
        if (f(1) >= 0) throw CalcError("the target power must exceed the power at HR = 1 (the significance level)");
        if (below) {
          if (!(f(1e-6) >= 0)) throw CalcError("the target power cannot be reached with these sample sizes");
          hr = numerics::brent(f, 1e-6, 1, 1e-12);
        } else {
          double hi = 2;
          while (f(hi) < 0 && hi < 1e6) hi *= 2;
          if (!(f(hi) >= 0)) throw CalcError("the target power cannot be reached with these sample sizes");
          hr = numerics::brent(f, 1, hi, 1e-12);
        }
      }
    }
    const auto d = logrank_design(s, hr);
    put_two_group(r, n);
    put_power(s, r, logrank_power(s, n.n1, n.n2, d));
    r.set("events", n.n1 * d.p1 + n.n2 * d.p2);
    r.set("e1", n.n1 * d.p1);
    r.set("e2", n.n2 * d.p2);
    r.set("hr", hr);
    r.set("med1", std::log(2.0) / d.lam1);
    r.set("med2", std::log(2.0) / d.lam2);
    r.set("pe1", d.p1);
    r.set("pe2", d.p2);
    r.set("accrual", s("accrual"));
    r.set("follow", s("follow"));
    r.set("alpha", s("alpha"));
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("events", "Events", "num", "Expected total number of events"), col("hr", "HR", "num"),
                   col("med1", "Median1", "num", "Treatment median survival"), col("med2", "Median2", "num", "Control median survival"),
                   col("pe1", "P(Event)1", "num"), col("pe2", "P(Event)2", "num"), col("accrual", "Accrual", "num"),
                   col("follow", "Follow-Up", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sample sizes of " + fint(r.num("n1")) + " (treatment) and " + fint(r.num("n2")) +
           " (control), with about " + fint(std::ceil(r.num("events") - 1e-9)) + " expected events, " + achieves(s, r) +
           " to detect a hazard ratio of " + fnum(r.num("hr")) + " (median survival " + fnum(r.num("med1"), 2) +
           " versus " + fnum(r.num("med2"), 2) + ") with a " + alt_words(alt(s)) + " logrank test at α = " +
           fnum(r.num("alpha")) + " (" + (s.is("method", "freedman") ? "Freedman" : "Schoenfeld") +
           " method), assuming exponential survival, uniform accrual over " + fnum(r.num("accrual")) +
           " time units and " + fnum(r.num("follow")) + " units of additional follow-up.";
  };
  p.definitions = {{"HR", "Hazard ratio λ1/λ2 (treatment relative to control)."},
                   {"Events", "Expected number of events N1·P(Event)1 + N2·P(Event)2."},
                   {"P(Event)", "Probability that a subject has an event before the end of the study."},
                   {"Accrual, Follow-Up", "Length of uniform enrollment and of follow-up after enrollment ends."}};
  p.references = {"Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival "
                  "distributions. Biometrika 68, 316–319.",
                  "Freedman, L. S. (1982). Tables of the number of patients required in clinical trials using the "
                  "logrank test. Statistics in Medicine 1, 121–129.",
                  "Lachin, J. M. & Foulkes, M. A. (1986). Evaluation of sample size and power for analyses of survival "
                  "with allowance for nonuniform patient entry, losses to follow-up, noncompliance, and "
                  "stratification. Biometrics 42, 507–519.",
                  REF_MACHIN};
  p.notes = {"Both formulas are large-sample approximations. With equal groups they agree with simulation. With "
             "unequal groups and a hazard ratio far from 1, Schoenfeld's overstates the power when the larger group "
             "has the higher hazard (HR 0.3 with 1:4 allocation: 0.86, against 0.74 by simulation) and understates "
             "it in the opposite case, while Freedman's is conservative. With fewer than about 30 events both are "
             "rough; the simulation method for two survival curves checks a design directly."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Treatment Group Size)";
  p.curve_min = 1;
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"hr", "0.6 0.7 0.75"}, {"med2", "12"}, {"accrual", "12"}, {"follow", "12"}};
  return p;
}
SSC_REGISTER(logrank);

// =====================================================================================
// Logrank non-inferiority (hazard ratio margin)
// =====================================================================================
// With higher hazard worse (e.g. death), H0: HR ≥ M versus H1: HR < M for a
// margin M > 1; with higher hazard better (e.g. recovery), H0: HR ≤ M versus
// H1: HR > M for M < 1.  The logrank/Cox statistic for ln HR − ln M has drift
// √(D w1 w2)·|ln M − ln HR| (Schoenfeld; Jung et al., 2005).
double logrank_ni_power(double n1, double n2, const LogrankDesign& d, double margin, double alpha) {
  const double events = n1 * d.p1 + n2 * d.p2;
  const double w1 = n1 / (n1 + n2);
  const double drift = std::sqrt(events * w1 * (1 - w1)) * std::fabs(std::log(margin) - std::log(d.hr));
  return stats::pnorm(drift - stats::qnorm(alpha, false));
}

Procedure logrank_ni() {
  Procedure p;
  p.id = "surv.logrank.ni";
  p.name = "Non-Inferiority Logrank Tests (Hazard Ratio Margin)";
  p.category = "Survival";
  p.subcategory = "Two Survival Curves";
  p.order = 2;
  p.summary =
      "Sample size and power for showing that a treatment's hazard ratio relative to control is below a "
      "non-inferiority margin (or above it when a higher hazard is better), with Schoenfeld's formula for the "
      "number of events and exponential survival with uniform accrual and follow-up to convert events into subjects.";
  p.keywords = {"non-inferiority", "logrank", "survival", "hazard ratio margin", "time to event", "active control", "survival non-inferiority"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      select("hazdir", "Direction",
             {{"worse", "Higher hazard is worse (H1: HR < M, margin M > 1)", {}},
              {"better", "Higher hazard is better (H1: HR > M, margin M < 1)", {}}},
             "worse")
          .grp(G_TEST)
          .hint("Whether the event is bad, such as death (a higher hazard is worse and the margin exceeds 1), or "
                "good, such as recovery (the margin is below 1)."),
      power_param(),
      alpha_param("0.025", "Alpha (One-Sided)").range(0, 0.5, true, true)};
  for (auto& q : two_group_params(1, "200")) p.params.push_back(q);
  for (auto q :
       {num("margin", "M (Non-Inferiority Margin for the HR)", "1.3")
            .grp(G_EFFECT)
            .sym("M")
            .above(0)
            .sugg({"1.3", "1.25", "1.2 1.3"})
            .hint("The least favorable hazard ratio (treatment/control) still considered acceptable, e.g. 1.3; "
                  "below 1 when a higher hazard is better."),
        num("hr", "HR (Actual Hazard Ratio, Treatment/Control)", "1")
            .grp(G_EFFECT)
            .sym("HR")
            .above(0)
            .sugg({"1", "0.9 1"})
            .hint("The hazard ratio of treatment to control that you expect to be true, often 1 (no difference); "
                  "it must lie on the acceptable side of the margin."),
        select("survin", "Control Group Survival Input",
               {{"median", "Median survival time", {}}, {"hazard", "Hazard rate", {}}, {"surv", "Survival proportion at time T0", {}}},
               "median")
            .grp(G_EFFECT)
            .hint(HELP_SURVIN),
        num("med2", "Median Survival (Control)", "12")
            .grp(G_EFFECT)
            .sym("Median2")
            .above(0)
            .hint(HELP_MED2)
            .show(when(is("survin", {"median"}))),
        num("lam2", "λ2 (Control Hazard Rate)", "0.0578")
            .grp(G_EFFECT)
            .sym("λ2")
            .above(0)
            .hint(HELP_LAM2)
            .show(when(is("survin", {"hazard"}))),
        num("s2", "S2 (Control Survival at T0)", "0.5")
            .grp(G_EFFECT)
            .sym("S2")
            .range(0, 1, true, true)
            .hint(HELP_S2)
            .show(when(is("survin", {"surv"}))),
        num("t0", "T0 (Time for S2)", "12")
            .grp(G_EFFECT)
            .sym("T0")
            .above(0)
            .hint(HELP_T0)
            .show(when(is("survin", {"surv"}))),
        num("accrual", "Accrual Time (Uniform Entry)", "12")
            .grp("Study Duration")
            .sym("Accrual")
            .at_least(0)
            .hint(HELP_ACCRUAL),
        num("follow", "Additional Follow-Up after Accrual", "12")
            .grp("Study Duration")
            .sym("Follow-Up")
            .at_least(0)
            .hint(HELP_FOLLOW),
        num("eta", "Loss to Follow-Up Hazard Rate", "0").grp("Study Duration").sym("η").at_least(0).hint(HELP_ETA)})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double margin = s("margin"), hr = s("hr"), alpha = s("alpha");
    const bool worse = s.is("hazdir", "worse");
    if (worse && !(margin > 1)) throw CalcError("when a higher hazard is worse, the margin must exceed 1");
    if (!worse && !(margin < 1)) throw CalcError("when a higher hazard is better, the margin must be below 1");
    // on the margin itself the power is alpha, which no sample size raises
    if (s.solve == "n" ? (worse ? !(hr < margin) : !(hr > margin)) : (worse ? hr > margin : hr < margin))
      throw CalcError(s.solve == "n" ? "the actual hazard ratio must lie on the non-inferior side of the margin"
                                     : "the actual hazard ratio cannot lie on the inferior side of the margin");
    const auto d = logrank_design(s, hr);
    if (!(d.p1 > 0 || d.p2 > 0)) throw CalcError("no events are expected with this accrual and follow-up");
    TwoN n;
    if (s.solve == "n") {
      const double z = stats::qnorm(alpha, false) + stats::qnorm(s("power"));
      const double lr = std::log(margin) - std::log(hr);
      const double guess = z * z / (0.25 * lr * lr) / ((d.p1 + d.p2) / 2) / 2;
      n = solve_two_group_n(s, [&](double n1, double n2) { return logrank_ni_power(n1, n2, d, margin, alpha); },
                            s("power"), 1, std::ceil(guess), true, ctx);
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    put_power(s, r, logrank_ni_power(n.n1, n.n2, d, margin, alpha));
    r.set("events", n.n1 * d.p1 + n.n2 * d.p2);
    r.set("margin", margin);
    r.set("hr", hr);
    r.set("med2", std::log(2.0) / d.lam2);
    r.set("pe1", d.p1);
    r.set("pe2", d.p2);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("events", "Events", "num", "Expected total number of events"), col("margin", "Margin M", "num"),
                   col("hr", "HR", "num", "Actual hazard ratio"), col("med2", "Median2", "num"),
                   col("pe1", "P(Event)1", "num"), col("pe2", "P(Event)2", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const bool worse = s.is("hazdir", "worse");
    return "Group sample sizes of " + fint(r.num("n1")) + " (treatment) and " + fint(r.num("n2")) +
           " (control), with about " + fint(std::ceil(r.num("events") - 1e-9)) + " expected events, " + achieves(s, r) +
           " to show that the hazard ratio is " + (worse ? "below" : "above") + " the non-inferiority margin of " +
           fnum(r.num("margin")) + " when the actual hazard ratio is " + fnum(r.num("hr")) +
           ", using a one-sided logrank test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"M", "Non-inferiority margin for the hazard ratio (treatment/control)."},
                   {"HR", "Actual hazard ratio at which power is computed."},
                   {"Events", "Expected number of events N1·P(Event)1 + N2·P(Event)2."}};
  p.references = {"Schoenfeld, D. (1981). The asymptotic properties of nonparametric tests for comparing survival "
                  "distributions. Biometrika 68, 316–319.",
                  "Jung, S.-H., Kang, S. J., McCall, L. M. & Blumenstein, B. (2005). Sample size computation for "
                  "two-sample noninferiority log-rank test. Journal of Biopharmaceutical Statistics 15, 969–979.",
                  REF_MACHIN};
  p.curve_param = "n1";
  p.curve_label = "N1 (Treatment Group Size)";
  p.curve_min = 1;
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"margin", "1.3"}, {"hr", "1"}, {"med2", "12"}, {"accrual", "12"}, {"follow", "12"}};
  return p;
}
SSC_REGISTER(logrank_ni);

// =====================================================================================
// One Poisson rate
// =====================================================================================
// The count X of events in a total exposure of N·t is Poisson with mean
// μ = λ·N·t.  Both tests reject in the tails of X: the exact test where the
// Poisson tail probability under λ0 is at most α (α/2 per side for two-sided
// tests), the z-test where |X − μ0|/√μ0 exceeds the normal critical value.
// Power and the actual significance level are therefore Poisson tail sums.
struct PoisRegion {
  double cl = -1;        // reject if X ≤ cl
  double cu = INFINITY;  // reject if X ≥ cu
};

double pois_at_least(double c, double mu) {
  if (c <= 0) return 1.0;
  if (!std::isfinite(c)) return 0.0;
  return stats::ppois(c - 1, mu, false);
}

double pois_at_most(double c, double mu) { return c < 0 ? 0.0 : stats::ppois(c, mu, true); }

PoisRegion pois_region(bool exact, double mu0, double alpha, Alt a) {
  PoisRegion r;
  const double side = tail_alpha(alpha, a);  // never 0, at which the steps below would never end
  const double z = stats::qnorm(side, false), sd = std::sqrt(mu0);
  if (exact) {
    // start near the normal quantile (with the Poisson skewness term), then step; beyond 10^15 events a step of
    // one no longer changes a count in floating point, and the start is as close as the arithmetic can tell
    const bool step = mu0 < 1e15;
    if (a != Alt::Lower) {  // smallest c with P(X ≥ c) ≤ side
      double c = std::max(0.0, std::ceil(mu0 + z * sd + (z * z - 1) / 6));
      while (step && c > 0 && pois_at_least(c - 1, mu0) <= side) c -= 1;
      while (step && pois_at_least(c, mu0) > side) c += 1;
      r.cu = c;
    }
    if (a != Alt::Upper) {  // largest c with P(X ≤ c) ≤ side
      double c = std::max(-1.0, std::floor(mu0 - z * sd + (z * z - 1) / 6));
      while (step && c >= 0 && pois_at_most(c, mu0) > side) c -= 1;
      while (step && pois_at_most(c + 1, mu0) <= side) c += 1;
      r.cl = c;
    }
    return r;
  }
  if (a != Alt::Lower) r.cu = std::floor(mu0 + z * sd) + 1;                // X > μ0 + z√μ0
  if (a != Alt::Upper) r.cl = std::max(-1.0, std::ceil(mu0 - z * sd) - 1);  // X < μ0 − z√μ0
  return r;
}

double pois_power(const PoisRegion& r, double mu) { return stats::clamp01(pois_at_least(r.cu, mu) + pois_at_most(r.cl, mu)); }

// Beyond 1e15 expected events pois_region no longer steps to the critical counts, and beyond 2^53 the counts that can
// be stored are no longer every whole number; beyond about 1e26 they lie so far apart next to the standard deviation
// that the critical counts, and so the powers, were far off (a two-sided exact test of size 1 at 1e33 events). There
// the critical counts are kept as their distances from μ0, and the tails are the normal ones with the Poisson
// skewness and a continuity correction (Edgeworth; their error is of the order of 1/μ): up and lo are the distances
// of the upper and lower critical counts (NaN: no such tail), dmu = μ − μ0.
constexpr double POIS_FAR = 1e15;

double pois_power_far(double up, double lo, double mu, double dmu) {
  using namespace stats;
  const double sd = std::sqrt(mu);
  // none far out in a tail, where w² overflows next to a density of 0 (μ1 a rounding error above 0)
  auto skew = [&](double w) { return std::fabs(w) < 40 ? dnorm(w) * (w * w - 1) / (6 * sd) : 0.0; };
  double pw = 0;
  if (!std::isnan(up)) {
    const double w = (up - 0.5 - dmu) / sd;  // P(X ≥ μ0 + up)
    pw += pnorm(w, false) + skew(w);
  }
  if (!std::isnan(lo)) {
    const double w = (lo + 0.5 - dmu) / sd;  // P(X ≤ μ0 + lo)
    pw += pnorm(w) - skew(w);
  }
  return clamp01(pw);
}

// The distances from μ0 of the critical counts (NaN: no such tail), placed as pois_region places them; f0, the
// fraction of μ0, makes them whole counts.
std::pair<double, double> pois_region_far(bool exact, double mu0, double alpha, Alt a) {
  const double side = tail_alpha(alpha, a);
  const double z = stats::qnorm(side, false), sd0 = std::sqrt(mu0), f0 = mu0 - std::floor(mu0);
  double up = NAN, lo = NAN;
  if (!exact) {
    if (a != Alt::Lower) up = std::floor(f0 + z * sd0) + 1 - f0;  // X > μ0 + z√μ0
    if (a != Alt::Upper) lo = std::ceil(f0 - z * sd0) - 1 - f0;   // X < μ0 − z√μ0
    return {up, lo};
  }
  // from the normal quantile with the skewness and continuity terms, then stepped by the tails themselves while a
  // step of one still changes the distance (the exact test never exceeds alpha)
  const bool step = z * sd0 < 1e15;
  if (a != Alt::Lower) {  // smallest with P(X ≥ μ0 + up) ≤ side
    up = std::ceil(f0 + z * sd0 + (z * z - 1) / 6 + 0.5) - f0;
    while (step && pois_power_far(up - 1, NAN, mu0, 0) <= side) up -= 1;
    while (step && pois_power_far(up, NAN, mu0, 0) > side) up += 1;
  }
  if (a != Alt::Upper) {  // largest with P(X ≤ μ0 + lo) ≤ side
    lo = std::floor(f0 - z * sd0 + (z * z - 1) / 6 - 0.5) - f0;
    while (step && pois_power_far(NAN, lo, mu0, 0) > side) lo -= 1;
    while (step && pois_power_far(NAN, lo + 1, mu0, 0) <= side) lo += 1;
  }
  return {up, lo};
}

struct OnePoisResult {
  double power, actual_alpha;
  PoisRegion region;
  bool exact;  // power from the Poisson distribution (else the normal approximation)
};

// The power at λ1 with n subjects, and the size.  d, when given, is λ1 − λ0 itself: with many events the difference
// a power needs can be finer than the doubles next to λ0 (8.9e-17 next to λ0 = 1 at 1e33 events, where they lie
// 2.2e-16 apart).  Beyond 1e15 events μ1 − μ0 is taken from the difference of the rates, of which μ1 and μ0 keep few
// digits or none.
OnePoisResult one_pois_eval(const Scenario& s, double n, double lam1, double d = NAN) {
  const double lam0 = s("lam0"), exposure = n * s("time"), mu0 = lam0 * exposure;
  const bool far = mu0 >= POIS_FAR;
  const double mu1 = std::isnan(d) ? lam1 * exposure : mu0 + d * exposure;
  const double dmu = !std::isnan(d) ? d * exposure : far ? (lam1 - lam0) * exposure : mu1 - mu0;
  const double alpha = s("alpha");
  const Alt a = alt(s);
  OnePoisResult res;
  res.exact = s.is("test", "exact") || s.is("method", "exact");
  if (!res.exact) {
    const double z = stats::qnorm(tail_alpha(alpha, a), false);
    const double s0 = std::sqrt(mu0), s1 = std::sqrt(mu1);
    double pw = 0;
    if (a != Alt::Lower) pw += stats::pnorm((dmu - z * s0) / s1);
    if (a != Alt::Upper) pw += stats::pnorm((-dmu - z * s0) / s1);
    res.power = stats::clamp01(pw);
    res.actual_alpha = NAN;
    return res;
  }
  const bool exact_test = s.is("test", "exact");
  if (!far) {
    res.region = pois_region(exact_test, mu0, alpha, a);
    res.power = pois_power(res.region, mu1);
    res.actual_alpha = pois_power(res.region, mu0);
    return res;
  }
  const auto [up, lo] = pois_region_far(exact_test, mu0, alpha, a);
  // the counts shown: exact up to 2^53, beyond it the nearest that can be stored
  res.region.cu = std::isnan(up) ? INFINITY : mu0 + up;
  res.region.cl = std::isnan(lo) ? -1 : mu0 + lo;
  res.power = pois_power_far(up, lo, mu1, dmu);
  res.actual_alpha = pois_power_far(up, lo, mu0, 0);
  return res;
}

// non-breaking spaces keep each bound on one line in reports
std::string pois_region_text(const PoisRegion& r) {
  std::string t;
  if (r.cl >= 0) t += "X\u00a0≤\u00a0" + fint(r.cl);
  if (std::isfinite(r.cu)) t += std::string(t.empty() ? "" : " or ") + "X\u00a0≥\u00a0" + fint(r.cu);
  return t.empty() ? "never rejects" : t;
}

Procedure one_poisson() {
  Procedure p;
  p.id = "rates.one.poisson";
  p.name = "Tests for One Poisson Rate";
  p.category = "Rates and Counts";
  p.subcategory = "One Rate";
  p.order = 0;
  p.summary =
      "Sample size and power for testing an event rate (events per unit of exposure) against a hypothesized value "
      "with the exact Poisson test or a z-test. Power is computed exactly from the Poisson distribution, which also "
      "gives the actual significance level.";
  p.keywords = {"poisson", "event rate", "incidence", "counts", "exact test", "person-time", "single rate", "rate against a reference value"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}, {"lam1", "Effect Size (λ1)", "lam1", "Detectable λ1"}};
  p.params = {
      alternative("Two-sided (H1: λ ≠ λ0)", "One-sided (H1: λ > λ0)", "One-sided (H1: λ < λ0)"),
      select("test", "Test Type", {{"exact", "Exact Poisson test", {}}, {"z", "Z-test using λ0", {}}}, "exact")
          .grp(G_TEST)
          .hint("The test of the analysis. The exact Poisson test never exceeds alpha; two-sided, it puts at "
                "most α/2 in each tail (equal-tailed, unlike R's poisson.test). The z-test compares the event count "
                "with its expected value under λ0 by the normal approximation."),
      select("method", "Power Calculation",
             {{"exact", "Exact (Poisson tail probabilities)", {}}, {"normal", "Normal approximation", {}}}, "exact")
          .grp(G_TEST)
          .hint("The Poisson distribution gives the exact power and the actual alpha of the chosen test.")
          .show(when(is_not("test", {"exact"}))),
      power_param(),
      alpha_param(),
      n_param("n", "N (Number of Subjects)", "50", 1)
          .sym("N")
          .hint("Number of subjects, each observed for the exposure time t; their events are counted together."),
      dropout_param(),
      num("time", "t (Exposure Time per Subject)", "1")
          .grp(G_EFFECT)
          .sym("t")
          .above(0)
          .sugg({"1", "2", "0.5"})
          .hint("How long each subject is observed, in the time unit of the rates (e.g. years for rates per "
                "year)."),
      num("lam0", "λ0 (Event Rate under H0)", "1")
          .grp(G_EFFECT)
          .sym("λ0")
          .above(0)
          .sugg({"1", "0.5", "2"})
          .hint("Events per subject per unit of time assumed by the null hypothesis, such as a historical rate."),
      num("lam1", "λ1 (Event Rate under H1)", "1.5")
          .grp(G_EFFECT)
          .sym("λ1")
          .above(0)
          .sugg({"1.5", "1.25 1.5 2"})
          .hint("The true number of events per subject per unit of time that you expect; the study is sized to "
                "tell it apart from λ0.")
          .show(when(not_solving({"lam1"}))),
      select("l1dir", "Search Direction for λ1", {{"above", "λ1 > λ0", {}}, {"below", "λ1 < λ0", {}}}, "above")
          .grp(G_EFFECT)
          .hint("A two-sided test can detect a λ1 above or below λ0; choose which of the two to find.")
          .show(when({solving({"lam1"}), is("alternative", {"two"})})),
  };
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    const double lam0 = s("lam0"), t = s("time");
    const Alt a = alt(s);
    double n, lam1;
    if (s.solve == "n") {
      lam1 = s("lam1");
      if (lam1 == lam0) throw CalcError("λ1 equals λ0");
      if ((a == Alt::Upper && lam1 < lam0) || (a == Alt::Lower && lam1 > lam0))
        throw CalcError("λ1 lies on the side opposite to the one-sided alternative");
      const double za = stats::qnorm(tail_alpha(s("alpha"), a), false);
      const double zb = stats::qnorm(s("power"));
      const double g = std::pow((za * std::sqrt(lam0) + zb * std::sqrt(lam1)) / (lam1 - lam0), 2) / t;
      solve::NSearch o;
      o.lo = 1;
      o.guess = std::ceil(g);
      o.monotone = !(s.is("test", "exact") || s.is("method", "exact"));
      if (!o.monotone) o.period = 1 / (lam0 * t);
      o.hi = 1e9;
      n = solve::min_n([&](double x) { return one_pois_eval(s, x, lam1).power; }, s("power"), o, &ctx);
    } else {
      n = s("n");
      if (s.solve == "power") {
        lam1 = s("lam1");
      } else {
        const bool up = a == Alt::Upper || (a == Alt::Two && !s.is("l1dir", "below"));
        auto f = [&](double x) { return one_pois_eval(s, n, x).power; };
        const double target = s("power");
        // as far as a rate inside the accepted inputs (1e30) and at least 1e12·λ0
        const double last = std::fmax(1e12 * lam0, 1e30);
        if (!up && f(1e-300) < target) throw CalcError("the target power cannot be reached with this N for any λ1 < λ0");
        if (lam0 * (n * t) >= POIS_FAR) {
          // the difference d = |λ1 − λ0| by itself, which can be finer than the doubles next to λ0 (one_pois_eval),
          // then the double next to λ0 ± d that reaches the target, refused when none keeps d (solve::reach)
          const double sign = up ? 1 : -1;
          auto g = [&](double x) { return one_pois_eval(s, n, lam0 + sign * x, sign * x).power; };
          const double d = solve::solve_increasing(g, target, 0, up ? lam0 : lam0 / 2, up ? last - lam0 : lam0 - 1e-300,
                                                   1e-12, "λ1");
          lam1 = solve::reach(f, lam0 + sign * d, up ? last : 1e-300, up, target, lam0, sign * d, "λ1", "λ0",
                              "the target power cannot be reached for any admissible value of λ1");
        } else if (up) {
          lam1 = solve::solve_increasing(f, target, lam0, 2 * lam0, last, 1e-12, "λ1", "λ0");
        } else {
          lam1 = solve::solve_decreasing(f, target, lam0, lam0 / 2, 1e-300, 1e-12, "λ1", "λ0");
        }
      }
    }
    const auto res = one_pois_eval(s, n, lam1);
    r.set("n", n);
    put_power(s, r, res.power);
    r.set("exposure", n * t);
    r.set("lam0", lam0);
    r.set("lam1", lam1);
    r.set("ratio", lam1 / lam0);
    r.set("mu0", lam0 * n * t);
    r.set("mu1", lam1 * n * t);
    r.set("alpha", s("alpha"));
    if (res.exact) {
      r.set("actual_alpha", res.actual_alpha);
      r.set("region", pois_region_text(res.region));
    }
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("exposure", "N·t", "num", "Total exposure"), col("lam0", "λ0", "num"),
                   col("lam1", "λ1", "num"), col("ratio", "λ1/λ0", "num"),
                   col("mu0", "μ0", "num", "Expected number of events under H0"),
                   col("mu1", "μ1", "num", "Expected number of events under H1"), col("alpha", "Target Alpha", "alpha")})
      c.push_back(x);
    if (s.is("test", "exact") || s.is("method", "exact")) {
      c.push_back(col("actual_alpha", "Actual Alpha", "prob", "Exact size of the test under λ0"));
      c.push_back(col("region", "Reject H0 if", "text", "Rejection region for the number of events X"));
    }
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    const auto [l1, l0] = fnum_apart(r.num("lam1"), r.num("lam0"));
    std::string t = "A sample of " + fint(r.num("n")) + " subjects, each observed for " + fnum(s("time")) +
                    " time unit(s) (total exposure " + fnum(r.num("exposure")) + "), " + achieves_singular(s, r) +
                    " to detect an event rate of " + l1 + " when the rate under H0 is " + l0 + ", using a " +
                    alt_words(alt(s)) + " " +
                    (s.is("test", "exact") ? "exact Poisson test" : "z-test") + " with a target significance level of " +
                    fnum(r.num("alpha"));
    if (r.has("actual_alpha")) t += " (actual significance level " + fprob(r.num("actual_alpha")) + ")";
    return t + ".";
  };
  p.definitions = {{"N", "Number of subjects."},
                   {"t", "Exposure time per subject; N·t is the total exposure."},
                   {"λ0, λ1", "Event rates per unit of exposure under H0 and H1."},
                   {"μ0, μ1", "Expected numbers of events λ·N·t under H0 and H1."},
                   {"Actual Alpha", "Exact probability of rejecting H0 when λ = λ0 (discreteness makes it ≤ target)."},
                   {"Reject H0 if", "Rejection region for the observed number of events X."}};
  p.notes = {"Because the number of events is discrete, exact power does not increase monotonically with N; the "
             "reported N is the smallest sample size whose power reaches the target."};
  p.references = {"Garwood, F. (1936). Fiducial limits for the Poisson distribution. Biometrika 28, 437–442.",
                  "Gu, K., Ng, H. K. T., Tang, M. L. & Schucany, W. R. (2008). Testing the ratio of two Poisson "
                  "rates. Biometrical Journal 50, 283–298.",
                  REF_MACHIN};
  p.curve_param = "n";
  p.curve_label = "N (Number of Subjects)";
  p.curve_min = 1;
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"lam0", "1"}, {"lam1", "1.25 1.5 2"}, {"time", "1"}};
  return p;
}
SSC_REGISTER(one_poisson);

// =====================================================================================
// Two Poisson rates
// =====================================================================================
// A test of two rates finds group 1's rate the higher only with an event in group 1, and group 2's only with one in
// group 2: the power of each side is at most the chance of such an event, 1 − e^log_none. With next to no events
// expected in a group the normal approximation breaks down (the variance under H1 dwarfs the one under H0 and both
// sides tend to 1/2: a two-sided power of 100% with 1e-58 events expected), and so it does where overdispersed
// counts gather a group's events in so few subjects that it likely has none (a dispersion of 1e6 leaves 10
// subjects expecting 100 events a chance of 0.00016 of any); a side's power above that chance by more than the
// precision shown is refused, and an n search goes on to the first sizes that give no side such a power.
bool beyond_events(double side_power, double log_none) { return side_power + std::expm1(log_none) > 5e-5; }

// `kappa`: the dispersion of negative binomial counts (0 for Poisson ones, whose chance of no event is e^−events).
void check_rate_events(int group, double side_power, double log_none, double events, double kappa = 0) {
  if (!beyond_events(side_power, log_none)) return;
  const std::string g = "group " + std::to_string(group), power = power_text(side_power);
  // as many events Poisson-distributed would allow that power: it is the dispersion that leaves so little chance
  if (kappa > 0 && !beyond_events(side_power, -events))
    throw CalcError(g + " expects " + fnum_fine(events) + " events, but with a dispersion of " + fnum(kappa) +
                    " they gather in so few subjects that its chance of any is only " +
                    text::format_g(-std::expm1(log_none), 3) + "; without one no test finds its rate the higher, but "
                    "the normal approximation gives that a power of " + power + ": it fails with so little chance of "
                    "an event");
  throw CalcError(g + " expects only " + fnum_fine(events) + (events == 1 ? " event" : " events") +
                  ", and without one no test finds its rate the higher, but the normal approximation gives that a "
                  "power of " + power + ": it fails with so few events");
}

const char* FEW_EVENTS = "A group expects fewer than 10 events; the normal approximation of the test is rough with so few.";

// Sizes `n` from an n search at which a side's power exceeds the chance of the events it needs (`possible` is false):
// the smallest from which none does, with a note, or `n` itself (refused by check_rate_events) when no sizes up to the
// search's limit qualify.
TwoN sizes_with_events(const Scenario& s, const std::function<double(double, double)>& P,
                       const std::function<bool(double, double)>& possible, const TwoN& n, Row& r,
                       const RunContext& ctx) {
  try {
    const TwoN m = solve_two_group_n(s, [&](double n1, double n2) { return possible(n1, n2) ? P(n1, n2) : 0.0; },
                                     s("power"), 1, n.n1, true, ctx);
    r.notes.push_back("Smaller samples leave a group so little chance of an event that the normal approximation gives "
                      "the test more power than that chance allows; this is the smallest that does not.");
    return m;
  } catch (const CalcError&) {
    return n;
  }
}

Procedure two_poisson() {
  Procedure p;
  p.id = "rates.two.poisson";
  p.name = "Tests for the Ratio of Two Poisson Rates";
  p.category = "Rates and Counts";
  p.subcategory = "Two Rates";
  p.order = 1;
  p.summary =
      "Sample size and power for comparing two Poisson event rates (events per unit of exposure) with a z-test on the "
      "log rate ratio, using the variance under H1 (Wald) or the pooled variance under H0.";
  p.keywords = {"poisson", "rate ratio", "counts", "incidence", "events per person-time", "exposure", "incidence rate ratio", "irr", "two rates"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: λ1/λ2 ≠ 1)", "One-sided (H1: λ1/λ2 > 1)", "One-sided (H1: λ1/λ2 < 1)"),
      select("method", "Variance",
             {{"h0", "Pooled rate under H0 (score-type)", {}}, {"h1", "Rates under H1 (Wald)", {}}}, "h0")
          .grp(G_TEST)
          .hint("Where the test's variance of the log rate ratio comes from: the pooled rate expected under H0 "
                "(score-type), or the two rates under H1 (Wald)."),
      power_param(),
      alpha_param()};
  for (auto& q : two_group_params(1, "50")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  for (auto q : {num("lam1", "λ1 (Event Rate, Group 1)", "1.5")
                     .grp(G_EFFECT)
                     .sym("λ1")
                     .above(0)
                     .sugg({"1.5", "1.2 1.5 2"})
                     .hint("Expected number of events per subject per unit of time in group 1 (e.g. per "
                           "person-year)."),
                 num("lam2", "λ2 (Event Rate, Group 2)", "1")
                     .grp(G_EFFECT)
                     .sym("λ2")
                     .above(0)
                     .sugg({"1"})
                     .hint("Expected number of events per subject per unit of time in group 2; the test compares "
                           "λ1/λ2 with 1."),
                 num("time", "t (Exposure Time per Subject)", "1")
                     .grp(G_EFFECT)
                     .sym("t")
                     .above(0)
                     .sugg({"1", "2"})
                     .hint("How long each subject is observed, in the time unit of the rates.")})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    using namespace stats;
    const double l1 = s("lam1"), l2 = s("lam2"), t = s("time"), alpha = s("alpha");
    const Alt a = alt(s);
    // the variances' reciprocals of exposures below about 5.6e-309 overflow to infinity (and their ratio, the
    // statistic, is then no number)
    if (!std::isfinite(1 / t)) throw CalcError("the exposure time t is too small to compute with");
    if (!std::isfinite(1 / (t * std::fmin(l1, l2))))
      throw CalcError("the expected events per subject (λ × t) are too small to compute with");
    const double lr = std::log(l1 / l2);
    // the chances of rejecting with group 1's rate the higher and with group 2's (0 for a side the test does not have)
    auto sides = [&](double n1, double n2) -> std::pair<double, double> {
      const double v1 = 1 / (n1 * t * l1) + 1 / (n2 * t * l2);
      const double lbar = (n1 * l1 + n2 * l2) / (n1 + n2);
      const double v0 = s.is("method", "h1") ? v1 : (1 / (n1 * t) + 1 / (n2 * t)) / lbar;
      const double zt = qnorm(tail_alpha(alpha, a), false);
      auto tail = [&](double d) { return pnorm((d - zt * std::sqrt(v0)) / std::sqrt(v1)); };
      return {a != Alt::Lower ? tail(lr) : 0.0, a != Alt::Upper ? tail(-lr) : 0.0};
    };
    auto P = [&](double n1, double n2) {
      const auto [up, down] = sides(n1, n2);
      return a == Alt::Two ? clamp01(up + down) : up + down;  // a NaN stays NaN (std::min(1.0, x) gives 1)
    };
    TwoN n;
    if (s.solve == "n") {
      if (l1 == l2) throw CalcError("the rates are equal");
      if ((a == Alt::Upper && (l1 - l2) < 0) || (a == Alt::Lower && (l1 - l2) > 0))
        throw CalcError("the rate ratio λ1/λ2 lies on the side opposite to the one-sided alternative");
      const double g = z_factor(alpha, s("power"), a) * (1 / l1 + 1 / l2) / t / (lr * lr);
      n = solve_two_group_n(s, P, s("power"), 1, std::ceil(g), true, ctx);
      auto possible = [&](double n1, double n2) {
        const auto [up, down] = sides(n1, n2);
        return !beyond_events(up, -n1 * l1 * t) && !beyond_events(down, -n2 * l2 * t);
      };
      if (!possible(n.n1, n.n2)) n = sizes_with_events(s, P, possible, n, r, ctx);
    } else {
      n = two_group_n(s);
    }
    const auto [up, down] = sides(n.n1, n.n2);
    check_rate_events(1, up, -n.n1 * l1 * t, n.n1 * l1 * t);
    check_rate_events(2, down, -n.n2 * l2 * t, n.n2 * l2 * t);
    put_two_group(r, n);
    put_power(s, r, P(n.n1, n.n2));
    if (std::fmin(n.n1 * l1, n.n2 * l2) * t < 10) r.notes.push_back(FEW_EVENTS);
    r.set("lam1", l1);
    r.set("lam2", l2);
    r.set("ratio", l1 / l2);
    r.set("time", t);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("lam1", "λ1", "num"), col("lam2", "λ2", "num"), col("ratio", "λ1/λ2", "num"),
                   col("time", "t", "num", "Exposure per subject"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + ", each followed for " +
           fnum(r.num("time")) + " time unit(s), " + achieves(s, r) + " to detect a rate ratio of " +
           fnum(r.num("ratio")) + " (rates " + fnum(r.num("lam1")) + " and " + fnum(r.num("lam2")) +
           ") with a " + alt_words(alt(s)) + " z-test on the log rate ratio at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"λ1, λ2", "Event rates per unit of exposure."}, {"t", "Exposure time per subject."}};
  p.references = {"Gu, K., Ng, H. K. T., Tang, M. L. & Schucany, W. R. (2008). Testing the ratio of two Poisson rates. "
                  "Biometrical Journal 50, 283–298.",
                  REF_MACHIN};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.curve_min = 1;
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"lam1", "1.25 1.5 2"}, {"lam2", "1"}, {"time", "1"}};
  return p;
}
SSC_REGISTER(two_poisson);

// =====================================================================================
// Two negative binomial rates (Zhu & Lakkis 2014)
// =====================================================================================
Procedure two_negbin() {
  Procedure p;
  p.id = "rates.two.negbin";
  p.name = "Tests for the Ratio of Two Negative Binomial Rates";
  p.category = "Rates and Counts";
  p.subcategory = "Two Rates";
  p.order = 2;
  p.summary =
      "Sample size and power for comparing overdispersed event rates (e.g. exacerbations, relapses) between two "
      "groups with a negative binomial model, following Zhu & Lakkis (2014).";
  p.keywords = {"negative binomial", "overdispersion", "rate ratio", "exacerbation", "recurrent events", "counts", "relapse rate", "count data", "incidence rate ratio"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {
      alternative("Two-sided (H1: R ≠ 1)", "One-sided (H1: R > 1)", "One-sided (H1: R < 1)"),
      select("method", "Variance under H0",
             {{"3", "Method 3: pooled (H0 MLE) rate", {}}, {"2", "Method 2: true rates", {}}, {"1", "Method 1: control rate", {}}},
             "3")
          .grp(G_TEST)
          .hint("Which rates the variance under H0 uses (Zhu and Lakkis): the pooled rate (method 3, the usual "
                "choice), the two true rates (method 2) or the control rate (method 1)."),
      power_param(),
      alpha_param()};
  for (auto& q : two_group_params(1, "100")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  for (auto q :
       {num("r0", "λ2 (Control Event Rate per Unit Time)", "1")
            .grp(G_EFFECT)
            .sym("λ2")
            .above(0)
            .sugg({"0.8", "1", "1.5"})
            .hint("Expected number of events per subject per unit of time in the control group (group 2), such as "
                  "yearly exacerbations."),
        num("rr", "R (Rate Ratio λ1/λ2)", "0.7")
            .grp(G_EFFECT)
            .sym("R")
            .above(0)
            .sugg({"0.7", "0.6 0.7 0.8"})
            .hint("The treatment group's event rate divided by the control group's, to be detected; 0.7 means 30% "
                  "fewer events."),
        num("kappa", "κ (Dispersion: Var = μ + κμ²)", "0.5")
            .grp(G_EFFECT)
            .sym("κ")
            .at_least(0)
            .sugg({"0.5", "0.3 0.5 1"})
            .hint("How much more the counts vary between subjects than Poisson counts would (0 for Poisson "
                  "counts); estimate it from earlier studies of the same outcome."),
        num("time", "Mean Exposure Time per Subject", "1")
            .grp(G_EFFECT)
            .sym("t")
            .above(0)
            .sugg({"1", "0.75 1"})
            .hint("Average time each subject is followed, in the time unit of the rates.")})
    p.params.push_back(q);
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    using namespace stats;
    const double r0 = s("r0"), rr = s("rr"), kappa = s("kappa"), mt = s("time"), alpha = s("alpha");
    const double r1 = r0 * rr;
    const Alt a = alt(s);
    // the variances' reciprocals of exposures below about 5.6e-309 overflow to infinity (and their ratio, the
    // statistic, is then no number)
    if (!std::isfinite(1 / (mt * std::fmin(r0, r1))))
      throw CalcError("the expected events per subject (rate × time) are too small to compute with");
    const std::string& m = s.s("method");
    // the chances of rejecting with the treatment rate the higher and with the control rate (0 for a side the test
    // does not have)
    auto sides = [&](double n1, double n2) -> std::pair<double, double> {
      // n1: treatment (rate r1), n2: control (rate r0)
      const double v1 = (1 / (mt * r1) + kappa) / n1 + (1 / (mt * r0) + kappa) / n2;
      double v0;
      if (m == "1") {
        v0 = (1 / (mt * r0) + kappa) * (1 / n1 + 1 / n2);
      } else if (m == "2") {
        v0 = v1;
      } else {
        const double rt = (n1 * r1 + n2 * r0) / (n1 + n2);
        v0 = (1 / (mt * rt) + kappa) * (1 / n1 + 1 / n2);
      }
      const double lr = std::log(rr);
      const double zt = qnorm(tail_alpha(alpha, a), false);
      auto tail = [&](double d) { return pnorm((d - zt * std::sqrt(v0)) / std::sqrt(v1)); };
      return {a != Alt::Lower ? tail(lr) : 0.0, a != Alt::Upper ? tail(-lr) : 0.0};
    };
    auto P = [&](double n1, double n2) {
      const auto [up, down] = sides(n1, n2);
      return a == Alt::Two ? clamp01(up + down) : up + down;  // a NaN stays NaN (std::min(1.0, x) gives 1)
    };
    // the chance that a subject with mean count μ has no event: (1 + κμ)^(−1/κ), e^(−μ) for κ = 0
    auto log_none = [&](double mu) { return kappa > 0 ? -std::log1p(kappa * mu) / kappa : -mu; };
    TwoN n;
    if (s.solve == "n") {
      if (rr == 1) throw CalcError("a rate ratio of 1 has no effect to detect");
      if ((a == Alt::Upper && (rr - 1) < 0) || (a == Alt::Lower && (rr - 1) > 0))
        throw CalcError("the rate ratio R lies on the side opposite to the one-sided alternative");
      const double v = (1 / (mt * r1) + kappa) + (1 / (mt * r0) + kappa);
      n = solve_two_group_n(s, P, s("power"), 1, std::ceil(z_factor(alpha, s("power"), a) * v / (std::log(rr) * std::log(rr))),
                            true, ctx);
      auto possible = [&](double n1, double n2) {
        const auto [up, down] = sides(n1, n2);
        return !beyond_events(up, n1 * log_none(r1 * mt)) && !beyond_events(down, n2 * log_none(r0 * mt));
      };
      if (!possible(n.n1, n.n2)) n = sizes_with_events(s, P, possible, n, r, ctx);
    } else {
      n = two_group_n(s);
    }
    const auto [up, down] = sides(n.n1, n.n2);
    check_rate_events(1, up, n.n1 * log_none(r1 * mt), n.n1 * r1 * mt, kappa);
    check_rate_events(2, down, n.n2 * log_none(r0 * mt), n.n2 * r0 * mt, kappa);
    put_two_group(r, n);
    put_power(s, r, P(n.n1, n.n2));
    if (std::fmin(n.n1 * r1, n.n2 * r0) * mt < 10) r.notes.push_back(FEW_EVENTS);
    r.set("r1", r1);
    r.set("r0", r0);
    r.set("rr", rr);
    r.set("kappa", kappa);
    r.set("time", mt);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    for (auto x : {col("r1", "λ1", "num", "Treatment rate"), col("r0", "λ2", "num", "Control rate"), col("rr", "R", "num"),
                   col("kappa", "κ", "num"), col("time", "t", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sample sizes of " + fint(r.num("n1")) + " (treatment) and " + fint(r.num("n2")) + " (control) " +
           achieves(s, r) + " to detect a rate ratio of " + fnum(r.num("rr")) + " when the control rate is " +
           fnum(r.num("r0")) + " per unit time, the mean exposure is " + fnum(r.num("time")) +
           " and the dispersion is " + fnum(r.num("kappa")) + ", using a " + alt_words(alt(s)) +
           " test of the log rate ratio at α = " + fnum(r.num("alpha")) + ".";
  };
  p.definitions = {{"λ1, λ2", "Event rates in the treatment and control groups."},
                   {"R", "Rate ratio λ1/λ2."},
                   {"κ", "Negative binomial dispersion (Var = μ + κμ²; κ = 0 is Poisson)."}};
  p.references = {"Zhu, H. & Lakkis, H. (2014). Sample size calculation for comparing two negative binomial rates. "
                  "Statistics in Medicine 33, 376–387.",
                  "Keene, O. N., Jones, M. R. K., Lane, P. W. & Anderson, J. (2007). Analysis of exacerbation rates in "
                  "asthma and chronic obstructive pulmonary disease. Pharmaceutical Statistics 6, 89–97."};
  p.curve_param = "n1";
  p.curve_label = "N1 (Treatment Group Size)";
  p.curve_min = 1;
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"r0", "1"}, {"rr", "0.6 0.7 0.8"}, {"kappa", "0.5"}};
  return p;
}
SSC_REGISTER(two_negbin);

// =====================================================================================
// One variance (chi-square test)
// =====================================================================================
Procedure one_variance() {
  Procedure p;
  p.id = "var.one";
  p.name = "Tests for One Variance (Chi-Square Test)";
  p.category = "Variances";
  p.subcategory = "One Variance";
  p.order = 1;
  p.summary = "Exact power and sample size for the chi-square test of a normal population variance, H0: σ² = σ0².";
  p.keywords = {"variance", "standard deviation", "chi-square", "precision", "variability", "single variance", "standard deviation test"};
  p.solve_for = {{"n", "Sample Size", "n", "N"}, {"power", "Power", "power", "Power"}};
  p.params = {alternative("Two-sided (H1: σ ≠ σ0)", "One-sided (H1: σ > σ0)", "One-sided (H1: σ < σ0)"),
              power_param(),
              alpha_param(),
              n_param("n", "N (Sample Size)", "30", 2).sym("N").hint("Number of subjects whose outcome is measured."),
              dropout_param(),
              num("sd0", "σ0 (SD under H0)", "1")
                  .grp(G_EFFECT)
                  .sym("σ0")
                  .above(0)
                  .sugg({"1"})
                  .hint("The standard deviation assumed by the null hypothesis, such as a specification or a "
                        "historical value."),
              num("sd1", "σ1 (SD under H1)", "1.5")
                  .grp(G_EFFECT)
                  .sym("σ1")
                  .above(0)
                  .sugg({"1.5", "0.5", "1.25 1.5 2"})
                  .hint("The true standard deviation you expect; the study is sized to tell it apart from σ0.")};
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    using namespace stats;
    const double ratio = s("sd1") * s("sd1") / (s("sd0") * s("sd0")), alpha = s("alpha");
    const Alt a = alt(s);
    auto P = [&](double n) {
      const double df = n - 1;
      double pw = 0;
      const double side = tail_alpha(alpha, a);
      if (a != Alt::Lower) pw += pchisq(qchisq(side, df, false) / ratio, df, false);
      if (a != Alt::Upper) pw += pchisq(qchisq(side, df, true) / ratio, df, true);
      return clamp01(pw);
    };
    double n;
    if (s.solve == "n") {
      if (ratio == 1) throw CalcError("σ1 equals σ0");
      if ((a == Alt::Upper && ratio < 1) || (a == Alt::Lower && ratio > 1))
        throw CalcError("σ1 is on the side opposite to the one-sided alternative");
      solve::NSearch o;
      o.lo = 2;
      const double lr = 0.5 * std::log(ratio);
      o.guess = z_factor(alpha, s("power"), a) / (2 * lr * lr) + 2;
      n = solve::min_n(P, s("power"), o, &ctx);
    } else {
      n = s("n");
    }
    r.set("n", n);
    put_power(s, r, P(n));
    r.set("sd0", s("sd0"));
    r.set("sd1", s("sd1"));
    r.set("ratio", ratio);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (auto x : {col("n", "N", "int"), col("sd0", "σ0", "num"), col("sd1", "σ1", "num"),
                   col("ratio", "σ1²/σ0²", "num"), col("alpha", "Alpha", "alpha")})
      c.push_back(x);
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "A sample size of " + fint(r.num("n")) + " " + achieves_singular(s, r) +
           " to detect a standard deviation of " + fnum(r.num("sd1")) + " when the null value is " + fnum(r.num("sd0")) +
           " (variance ratio " + fnum(r.num("ratio")) + ") using a " + alt_words(alt(s)) + " chi-square test at α = " +
           fnum(r.num("alpha")) + ".";
  };
  p.references = {"Ostle, B. & Malone, L. C. (1988). Statistics in Research, 4th ed. Iowa State University Press.", REF_MACHIN};
  p.curve_param = "n";
  p.curve_label = "N (Sample Size)";
  p.dropout_groups = {"n"};
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"sd0", "1"}, {"sd1", "1.25 1.5 2"}};
  return p;
}
SSC_REGISTER(one_variance);

// =====================================================================================
// Two variances (F test)
// =====================================================================================
Procedure two_variances() {
  Procedure p;
  p.id = "var.two";
  p.name = "Tests for Two Variances (F-Test)";
  p.category = "Variances";
  p.subcategory = "Two Variances";
  p.order = 1;
  p.summary = "Exact power and sample size for the F-test comparing two normal variances, H0: σ1² = σ2².";
  p.keywords = {"variance ratio", "f-test", "homogeneity", "standard deviations", "variability", "equality of variances", "two standard deviations", "two variances"};
  p.solve_for = {{"n", "Sample Size", "n", "N (Total)"}, {"power", "Power", "power", "Power"}};
  p.params = {alternative("Two-sided (H1: σ1² ≠ σ2²)", "One-sided (H1: σ1² > σ2²)", "One-sided (H1: σ1² < σ2²)"),
              power_param(), alpha_param()};
  for (auto& q : two_group_params(2, "30")) p.params.push_back(q);
  p.params.push_back(dropout_param());
  p.params.push_back(num("ratio", "σ1²/σ2² (Variance Ratio under H1)", "2")
                         .grp(G_EFFECT)
                         .sym("σ1²/σ2²")
                         .above(0)
                         .sugg({"2", "1.5 2 3", "0.5"})
                         .hint("The variance (squared SD) of group 1 divided by that of group 2 that you expect; 1 "
                               "means equal variability."));
  p.run = [](const Scenario& s, Row& r, RunContext& ctx) {
    using namespace stats;
    const double ratio = s("ratio"), alpha = s("alpha");
    const Alt a = alt(s);
    auto P = [&](double n1, double n2) {
      const double d1 = n1 - 1, d2 = n2 - 1;
      const double side = tail_alpha(alpha, a);
      double pw = 0;
      if (a != Alt::Lower) pw += pf(qf(side, d1, d2, false) / ratio, d1, d2, false);
      if (a != Alt::Upper) pw += pf(qf(side, d1, d2, true) / ratio, d1, d2, true);
      return clamp01(pw);
    };
    TwoN n;
    if (s.solve == "n") {
      if (ratio == 1) throw CalcError("the variance ratio is 1");
      if ((a == Alt::Upper && ratio < 1) || (a == Alt::Lower && ratio > 1))
        throw CalcError("the variance ratio is on the side opposite to the one-sided alternative");
      const double lr = 0.5 * std::log(ratio);
      n = solve_two_group_n(s, P, s("power"), 2, z_factor(alpha, s("power"), a) / (lr * lr) + 2, true, ctx);
    } else {
      n = two_group_n(s);
    }
    put_two_group(r, n);
    put_power(s, r, P(n.n1, n.n2));
    r.set("ratio", ratio);
    r.set("alpha", alpha);
  };
  p.columns = [](const Scenario& s) {
    auto c = power_columns(s);
    for (const auto& x : two_group_columns()) c.push_back(x);
    c.push_back(col("ratio", "σ1²/σ2²", "num"));
    c.push_back(col("alpha", "Alpha", "alpha"));
    return c;
  };
  p.statement = [](const Scenario& s, const Row& r) {
    return "Group sample sizes of " + fint(r.num("n1")) + " and " + fint(r.num("n2")) + " " + achieves(s, r) +
           " to detect a variance ratio of " + fnum(r.num("ratio")) + " using a " + alt_words(alt(s)) +
           " F-test at α = " + fnum(r.num("alpha")) + ".";
  };
  p.references = {"Desu, M. M. & Raghavarao, D. (1990). Sample Size Methodology. Academic Press.", REF_MACHIN};
  p.curve_param = "n1";
  p.curve_label = "N1 (Group 1 Sample Size)";
  p.dropout_groups = {"n1", "n2"};
  p.dropout_total = "n";
  p.example_solve = "n";
  p.example = {{"power", "0.8 0.9"}, {"ratio", "1.5 2 3"}};
  return p;
}
SSC_REGISTER(two_variances);

}  // namespace
}  // namespace ssc::procs
