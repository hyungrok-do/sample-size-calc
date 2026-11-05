// Survival (logrank test with accrual and follow-up), rates (Poisson and
// negative binomial) and variances (chi-square and F tests).
#include <algorithm>
#include <functional>

#include "procedures/apart.hpp"
#include "procedures/common.hpp"
#include "procedures/derive_models.hpp"

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

// The mean of the logrank statistic under H1 (its drift), from the expected number of events.
double logrank_drift(const Scenario& s, double n1, double n2, const LogrankDesign& d) {
  const double events = n1 * d.p1 + n2 * d.p2;
  if (s.is("method", "freedman")) {
    const double phi = n2 / n1;  // control-to-treatment ratio
    return std::sqrt(events * phi) * (d.hr - 1) / (d.hr + phi);
  }
  // Schoenfeld
  const double w1 = n1 / (n1 + n2), w2 = 1 - w1;
  return std::sqrt(events * w1 * w2) * std::log(d.hr);
}

// Power from the expected number of events.
double logrank_power(const Scenario& s, double n1, double n2, const LogrankDesign& d) {
  return z_power(logrank_drift(s, n1, n2, d), s("alpha"), alt(s));
}

// The pieces of the logrank derivations: the hazards, each group's probability of an event (as event_prob computes
// it) and the expected number of events, from which the drift θ of the statistic follows (θ, as λ is a hazard).
// λ2 from the control group's survival input, in symbols and numbers (both empty when it is the hazard itself), and
// λ2 as a step writes it: given when entered as the hazard, else derived.
std::pair<std::string, std::string> control_hazard_text(const Scenario& s, double lam2) {
  using namespace derive;
  const std::string& in = s.s("survin");
  if (in == "median") return {"λ2 = ln 2 / Median2", "λ2 = ln 2 / " + gsig(s("med2")) + " = " + dsig(lam2)};
  if (in == "surv")
    return {"λ2 = −ln S2 / T0", "λ2 = −ln " + gsig(s("s2")) + " / " + gsig(s("t0")) + " = " + dsig(lam2)};
  return {"", ""};
}
std::string lam2_value(const Scenario& s, double lam2) {
  return s.is("survin", "hazard") ? derive::gsig(lam2) : derive::dsig(lam2);
}

Step hazards_step(const Scenario& s, const LogrankDesign& d) {
  using namespace derive;
  const auto [f, t] = control_hazard_text(s, d.lam2);
  const std::string lam1 = "λ1 = " + gv(d.hr) + " × " + lam2_value(s, d.lam2) + " = " + dsig(d.lam1);
  return {"Hazard rates (exponential survival)", (f.empty() ? "" : f + ", ") + "λ1 = HR × λ2",
          (t.empty() ? "" : t + ", ") + lam1};
}

// Pi in symbols (`i` the group, or "i" for both), and with the numbers of group i, its hazard `lam` written
// `lam_text`.
std::string event_prob_formula(const Scenario& s, const std::string& i) {
  const bool loss = s("eta") > 0, entry = s("accrual") > 0;
  const std::string h = loss ? "(λ" + i + " + η)" : "λ" + i;
  std::string t = "1 − e^(−" + h + "F)" + (entry ? " (1 − e^(−" + h + "A)) / (" + h + "A)" : "");
  if (loss) t = "λ" + i + "/(λ" + i + " + η) × [" + t + "]";
  return "P" + i + " = " + t +
         (entry ? ", with uniform entry over the accrual A and follow-up F after it"
                : ", all entering at once and followed F") +
         (loss ? ", and loss to follow-up at the rate η" : "");
}
std::string event_prob_numbers(const Scenario& s, int i, double lam, const std::string& lam_text, double p) {
  using namespace derive;
  const double eta = s("eta"), A = s("accrual"), F = s("follow");
  const bool loss = eta > 0;
  // h × F and h × A, h the hazard of leaving the risk set (the event's alone without loss to follow-up)
  const std::string h = loss ? dsig(lam + eta) : lam_text;
  std::string t = "1 − e^(−" + h + " × " + gsig(F) + ")";
  if (A > 0) t += " (1 − e^(−" + h + " × " + gsig(A) + ")) / (" + h + " × " + gsig(A) + ")";
  if (loss) t = lam_text + "/" + h + " × [" + t + "]";
  return "P" + std::to_string(i) + " = " + t + " = " + dsig(p);
}

Step event_probs_step(const Scenario& s, const LogrankDesign& d) {
  using namespace derive;
  return {"Probability of an event", event_prob_formula(s, "i"),
          event_prob_numbers(s, 1, d.lam1, dsig(d.lam1), d.p1) + "; " +
              event_prob_numbers(s, 2, d.lam2, lam2_value(s, d.lam2), d.p2)};
}

Step events_step(const TwoN& n, bool solved, const LogrankDesign& d) {
  using namespace derive;
  return {"Expected events", "D = N1 P1 + N2 P2",
          size(n.n1, solved) + " × " + dsig(d.p1) + " + " + size(n.n2, solved) + " × " + dsig(d.p2) + " = " +
              dv(n.n1 * d.p1 + n.n2 * d.p2, 2)};
}

// θ from the expected events by Schoenfeld's or Freedman's formula: its expression, the definitions it uses, and its
// numbers (`hr` the hazard ratio as written, given or derived).
struct DriftText {
  std::string expr, defs, numbers;
};
DriftText logrank_drift_text(const Scenario& s, const TwoN& n, const LogrankDesign& d, const std::string& hr) {
  using namespace derive;
  const double events = n.n1 * d.p1 + n.n2 * d.p2, drift = logrank_drift(s, n.n1, n.n2, d);
  if (s.is("method", "freedman")) {
    const std::string phi = dv(n.n2 / n.n1);
    return {"√(D φ) (HR − 1) / (HR + φ)", "φ = N2/N1",
            "√(" + dv(events, 2) + " × " + phi + ") × (" + hr + " − 1) / (" + hr + " + " + phi + ") = " + dv(drift)};
  }
  const double w1 = n.n1 / (n.n1 + n.n2);
  return {"√(D w1 w2) ln HR", "w1 = N1/(N1 + N2), w2 = 1 − w1",
          "√(" + dv(events, 2) + " × " + dv(w1) + " × " + dv(1 - w1) + ") × ln " + hr + " = " + dv(drift)};
}

// The shares of the groups, w1 = N1/N, where the allocation fixes them (equal groups, a ratio, a percentage): what
// a closed form of the events needs; NaN when one group's size is fixed instead.
double logrank_share(const Scenario& s) {
  const std::string& a = s.s("alloc");
  if (a == "equal") return 0.5;
  if (a == "ratio") return 1 / (1 + s("r"));
  if (a == "pct") return s("pct1") / 100;
  return NAN;
}

// The events a logrank test needs by the closed form of its method, and the subjects that many events take: for
// comparison with the search. `effect` is the log hazard ratio's distance from H0 in symbols and with its numbers.
Step logrank_closed_form(const Scenario& s, const LogrankDesign& d, double w1, double alpha, Alt a,
                         const std::string& effect, const std::string& effect_numbers, double effect_value,
                         bool freedman) {
  using namespace derive;
  const double target = s("power"), w2 = 1 - w1, za = z_alpha(alpha, a), zb = stats::qnorm(target);
  const double zz = (za + zb) * (za + zb), phi = w2 / w1;
  const double events = freedman ? zz * (d.hr + phi) * (d.hr + phi) / (phi * (d.hr - 1) * (d.hr - 1))
                                 : zz / (w1 * w2 * effect_value * effect_value);
  const double subjects = events / (w1 * d.p1 + w2 * d.p2);
  const std::string tail = std::string(" and N ≈ D / (w1 P1 + w2 P2), w1 and w2 = 1 − w1 the groups' shares of N") +
                           (freedman ? " and φ = w2/w1" : "") + ", for comparison; the steps below give the answer";
  const std::string per = ", N ≈ " + dv(events, 2) + " / (" + dv(w1) + " × " + dsig(d.p1) + " + " + dv(w2) + " × " +
                          dsig(d.p2) + ") = " + dv(subjects, 2);
  if (freedman)
    return {"Events needed (Freedman)", "D ≈ " + z_sum_formula(a) + " (HR + φ)² / (φ (HR − 1)²)" + tail,
            "D ≈ " + z_sum_numbers(alpha, a, target) + " × (" + gv(d.hr) + " + " + dv(phi) + ")² / (" + dv(phi) +
                " × (" + gv(d.hr) + " − 1)²) = " + dv(events, 2) + per};
  return {"Events needed (Schoenfeld)", "D ≈ " + z_sum_formula(a) + " / (w1 w2 (" + effect + ")²)" + tail,
          "D ≈ " + z_sum_numbers(alpha, a, target) + " / (" + dv(w1) + " × " + dv(w2) + " × (" + effect_numbers +
              ")²) = " + dv(events, 2) + per};
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const Alt a = alt(s);
    const bool solved = s.solve == "n", freedman = s.is("method", "freedman");
    const double alpha = s("alpha");
    const TwoN n{r.num("n1"), r.num("n2")};
    const LogrankDesign d = logrank_design(s, r.num("hr"));
    const double events = n.n1 * d.p1 + n.n2 * d.p2, drift = logrank_drift(s, n.n1, n.n2, d);
    const double power = logrank_power(s, n.n1, n.n2, d);
    const std::string method = freedman ? " (Freedman)" : " (Schoenfeld)";
    Derivation out;
    if (s.solve == "hr") {
      // the target power gives θ; HR solves θ(HR) = θ, where the events change with HR through P1: a root search
      const double target = s("power");
      out.find = find_effect("HR", target, "with " + two_sizes(n, false));
      const auto [f, t] = control_hazard_text(s, d.lam2);
      if (!f.empty()) out.steps.push_back({"Control hazard (exponential survival)", f, t});
      out.steps.push_back({"Probability of an event, control group", event_prob_formula(s, "2"),
                           event_prob_numbers(s, 2, d.lam2, lam2_value(s, d.lam2), d.p2)});
      out.steps.push_back(z_critical(alpha, a));
      out.steps.push_back(renamed(z_ncp_needed(drift, alpha, a, target), "λ", "θ"));
      const DriftText theta = logrank_drift_text(s, n, d, dv(d.hr));
      out.steps.push_back(
          {"Hazard ratio detected" + method,
           "HR solves " + theta.expr + " = θ, " + theta.defs + ", D = N1 P1 + N2 P2 and P1 the treatment group's "
               "probability of an event at λ1 = HR × λ2 (as P2 above): D changes with HR, so HR is found by a root "
               "search",
           "at HR = " + dv(d.hr) + ": λ1 = " + dv(d.hr) + " × " + lam2_value(s, d.lam2) + " = " + dsig(d.lam1) +
               ", P1 = " + dsig(d.p1) + ", D = " + gvi(n.n1) + " × " + dsig(d.p1) + " + " + gvi(n.n2) + " × " +
               dsig(d.p2) + " = " + dv(events, 2) + ", θ = " + theta.numbers});
      out.steps.push_back(effect_search("HR", d.hr, power, target));
      out.answer = "HR = " + dv(d.hr) + ", with " + dv(events, 2) + " expected events";
      return out;
    }
    out.steps.push_back(hazards_step(s, d));
    out.steps.push_back(event_probs_step(s, d));
    if (solved) {
      out.find = find_size(two_group_unknown(s), s("power"));
      const double w1 = logrank_share(s);
      if (std::isfinite(w1))
        out.steps.push_back(
            logrank_closed_form(s, d, w1, alpha, a, "ln HR", "ln " + gv(d.hr), std::log(d.hr), freedman));
    } else {
      out.find = find_power("with " + two_sizes(n, false));
    }
    out.steps.push_back(events_step(n, solved, d));
    const DriftText theta = logrank_drift_text(s, n, d, gv(d.hr));
    out.steps.push_back({"Noncentrality" + method, "θ = " + theta.expr + ", " + theta.defs, theta.numbers});
    out.steps.push_back(z_critical(alpha, a));
    out.steps.push_back(renamed(z_power_step(drift, alpha, a), "λ", "θ"));
    if (solved) {
      out.steps.push_back(two_group_search(
          s, [&](double x1, double x2) { return logrank_power(s, x1, x2, d); }, n, s("power"), 1));
      out.answer = sizes_answer(n, dv(events, 2) + " expected events and power " + dvp(power));
    } else {
      out.answer = "power = " + dvp(power) + ", with " + dv(events, 2) + " expected events";
    }
    return out;
  };
  p.definitions = {{"HR", "Hazard ratio λ1/λ2 (treatment relative to control)."},
                   {"Events", "Expected number of events N1·P(Event)1 + N2·P(Event)2."},
                   {"P(Event)", "Probability that a subject has an event before the end of the study."},
                   {"Accrual, Follow-Up", "Length of uniform enrollment and of follow-up after enrollment ends."}};
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
double logrank_ni_drift(double n1, double n2, const LogrankDesign& d, double margin) {
  const double events = n1 * d.p1 + n2 * d.p2;
  const double w1 = n1 / (n1 + n2);
  return std::sqrt(events * w1 * (1 - w1)) * std::fabs(std::log(margin) - std::log(d.hr));
}

double logrank_ni_power(double n1, double n2, const LogrankDesign& d, double margin, double alpha) {
  return stats::pnorm(logrank_ni_drift(n1, n2, d, margin) - stats::qnorm(alpha, false));
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const bool solved = s.solve == "n";
    const double margin = r.num("margin"), alpha = s("alpha");
    const TwoN n{r.num("n1"), r.num("n2")};
    const LogrankDesign d = logrank_design(s, r.num("hr"));
    const double events = n.n1 * d.p1 + n.n2 * d.p2, drift = logrank_ni_drift(n.n1, n.n2, d, margin);
    const double power = logrank_ni_power(n.n1, n.n2, d, margin, alpha);
    const std::string effect = "|ln " + gv(margin) + " − ln " + gv(d.hr) + "|";
    Derivation out;
    out.steps.push_back(hazards_step(s, d));
    out.steps.push_back(event_probs_step(s, d));
    // a one-sided test (Alt::Upper) that the hazard ratio lies on the acceptable side of the margin
    if (solved) {
      out.find = find_size(two_group_unknown(s), s("power"));
      const double w1 = logrank_share(s);
      if (std::isfinite(w1))
        out.steps.push_back(logrank_closed_form(s, d, w1, alpha, Alt::Upper, "ln M − ln HR",
                                                "ln " + gv(margin) + " − ln " + gv(d.hr),
                                                std::log(margin) - std::log(d.hr), false));
    } else {
      out.find = find_power("with " + two_sizes(n, false));
    }
    out.steps.push_back(events_step(n, solved, d));
    const double w1 = n.n1 / (n.n1 + n.n2);
    out.steps.push_back({"Noncentrality (Schoenfeld)", "θ = √(D w1 w2) |ln M − ln HR|, w1 = N1/(N1 + N2), w2 = 1 − w1",
                         "√(" + dv(events, 2) + " × " + dv(w1) + " × " + dv(1 - w1) + ") × " + effect + " = " +
                             dv(drift)});
    out.steps.push_back(z_critical(alpha, Alt::Upper));
    out.steps.push_back(renamed(z_power_step(drift, alpha, Alt::Upper), "λ", "θ"));
    if (solved) {
      out.steps.push_back(two_group_search(
          s, [&](double x1, double x2) { return logrank_ni_power(x1, x2, d, margin, alpha); }, n, s("power"), 1));
      out.answer = sizes_answer(n, dv(events, 2) + " expected events and power " + dvp(power));
    } else {
      out.answer = "power = " + dvp(power) + ", with " + dv(events, 2) + " expected events";
    }
    return out;
  };
  p.definitions = {{"M", "Non-inferiority margin for the hazard ratio (treatment/control)."},
                   {"HR", "Actual hazard ratio at which power is computed."},
                   {"Events", "Expected number of events N1·P(Event)1 + N2·P(Event)2."}};
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

// The first distance at which ok holds (a tail at most its side), where it holds from some distance on, sought from
// d by whole counts: one at a time for the first 64 steps, then by steps that double, and back by halves.  Next to a
// tail of nearly a half (two-sided alphas next to 1, one-sided ones next to 0.5) z is next to 0, a step moves the
// tail by about 0.4/σ0, and with 1e46 events (σ0 = 1e23) that is a ten-millionth of a double there: one at a time
// the steps went on until their sum got past the tail's rounding, 1e7 times or more, and Stop was not heard.
template <class F>
double pois_first_ok(double d, F ok) {
  const bool down = ok(d);  // the first ok one is d or below it
  double bad = d, good = d;
  for (double step = 1, k = 0;; step = ++k < 64 ? 1 : 2 * step) {
    const double e = down ? good - step : bad + step;
    if (ok(e) == down) {
      (down ? good : bad) = e;
      continue;
    }
    (down ? bad : good) = e;
    break;
  }
  while (good - bad > 1) {
    const double m = bad + std::floor((good - bad) / 2);
    (ok(m) ? good : bad) = m;
  }
  return good;
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
  // step of one still changes the distance (the exact test never exceeds alpha; see pois_first_ok() for the steps)
  const bool step = z * sd0 < 1e15;
  if (a != Alt::Lower) {  // smallest with P(X ≥ μ0 + up) ≤ side
    up = std::ceil(f0 + z * sd0 + (z * z - 1) / 6 + 0.5) - f0;
    if (step) up = pois_first_ok(up, [&](double d) { return pois_power_far(d, NAN, mu0, 0) <= side; });
  }
  if (a != Alt::Upper) {  // largest with P(X ≤ μ0 + lo) ≤ side
    lo = std::floor(f0 - z * sd0 + (z * z - 1) / 6 - 0.5) - f0;
    if (step) lo = -pois_first_ok(-lo, [&](double d) { return pois_power_far(NAN, -d, mu0, 0) <= side; });
  }
  return {up, lo};
}

struct OnePoisResult {
  double power, actual_alpha;
  PoisRegion region;
  bool exact;  // power from the Poisson distribution (else the normal approximation)
  // what the power came from, for the derivation: the expected events and their difference, and beyond POIS_FAR
  // the critical counts' distances from μ0 (NaN: no such tail)
  double mu0 = NAN, mu1 = NAN, dmu = NAN;
  bool far = false;
  double up = NAN, lo = NAN;
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
  res.mu0 = mu0;
  res.mu1 = mu1;
  res.dmu = dmu;
  res.far = far;
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
  res.up = up;
  res.lo = lo;
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

// How N was found where the exact power rises in a saw-tooth (solve::min_n, not monotone): the first size of its
// scan to reach the target, so that the size below it falls short unless the scan began at N.
Step sawtooth_n_search(const std::function<double(double)>& power, double n, double target, double lo) {
  using namespace derive;
  const double at = power(n), below = n - 1 >= lo ? power(n - 1) : NAN;
  if (!(below >= target)) return n_search([&](double x) { return x == n ? at : below; }, n, target, lo);
  return {"Smallest size",
          "N = the first size of the search's scan whose power reaches 1 − β (the exact power rises in a saw-tooth "
          "with N)",
          "N = " + dvi(n) + " is where the scan began, and its power is " + dvp(at, target) + " ≥ " + gv(target) +
              " (the saw-tooth gives N = " + dvi(n - 1) + " a power of " + dvp(below, target) + " too)"};
}

// How the tails are taken beyond POIS_FAR, as a formula says it.
const char* const POIS_FAR_TEXT =
    ", from the normal tails with the Poisson skewness and a continuity correction (Edgeworth), as beyond 10^15 events";

// The rejection region of a one-rate Poisson test as one_pois_eval placed it (critical counts from μ0 alone).
Step pois_region_step(const Scenario& s, const OnePoisResult& res) {
  using namespace derive;
  const Alt a = alt(s);
  const bool upper = a != Alt::Lower, lower = a != Alt::Upper;
  const std::string side = a == Alt::Two ? "α/2" : "α";
  const PoisRegion& g = res.region;
  std::string formula, numbers;
  auto add = [](std::string& to, const std::string& t) { to += (to.empty() ? "" : "; ") + t; };
  if (s.is("test", "exact")) {
    if (upper) {
      add(formula, "c_U = the smallest count with P(X ≥ c_U | μ0) ≤ " + side);
      add(numbers, "c_U = " + dvi(g.cu));
    }
    if (lower) {
      add(formula, "c_L = the largest count with P(X ≤ c_L | μ0) ≤ " + side);
      add(numbers, g.cl >= 0 ? "c_L = " + dvi(g.cl) : "no c_L: even P(X = 0 | μ0) exceeds " + side);
    }
    if (res.far) formula += POIS_FAR_TEXT;
  } else {
    const double z = z_alpha(s("alpha"), a);
    const std::string spread = dv(z) + " × √" + dv(res.mu0);
    if (upper) {
      add(formula, "c_U = ⌊μ0 + z√μ0⌋ + 1");
      add(numbers, "c_U = ⌊" + dv(res.mu0) + " + " + spread + "⌋ + 1 = " + dvi(g.cu));
    }
    if (lower) {
      add(formula, "c_L = ⌈μ0 − z√μ0⌉ − 1");
      add(numbers, g.cl >= 0 ? "c_L = ⌈" + dv(res.mu0) + " − " + spread + "⌉ − 1 = " + dvi(g.cl)
                             : "no c_L: μ0 − z√μ0 = " + dv(res.mu0 - z * std::sqrt(res.mu0)) +
                                   " leaves no count below it");
    }
    formula += ", the counts beyond μ0 ± z√μ0";
  }
  return {"Rejection region", formula, numbers};
}

// The chance that a one-rate Poisson test rejects at the mean count `at` (`sym` in symbols, at − μ0 = dmu), as
// one_pois_eval sums it (Poisson tails, or beyond POIS_FAR their Edgeworth approximation): in symbols, and with its
// terms (ending before the total) and the total, told apart from `against`.
struct PoisReject {
  std::string formula, terms, total;
};
PoisReject pois_reject(const Scenario& s, const OnePoisResult& res, const std::string& sym, double at, double dmu,
                       double total, double against = NAN) {
  using namespace derive;
  const Alt a = alt(s);
  const bool upper = a != Alt::Lower, lower = a != Alt::Upper;
  const PoisRegion& g = res.region;
  const std::string far = res.far ? POIS_FAR_TEXT : "";
  const auto [hi, lo] = res.far ? std::pair<double, double>{upper ? pois_power_far(res.up, NAN, at, dmu) : 0.0,
                                                            lower ? pois_power_far(NAN, res.lo, at, dmu) : 0.0}
                                : std::pair<double, double>{pois_at_least(g.cu, at), pois_at_most(g.cl, at)};
  // beyond POIS_FAR the means differ by less than the digits shown: μ1 as μ0 and the difference
  const std::string mean =
      res.far && dmu != 0 ? dv(res.mu0) + (dmu < 0 ? " − " : " + ") + dv(std::fabs(dmu)) : dv(at);
  // a lower tail without a critical count adds nothing
  const bool hi_tail = upper && std::isfinite(g.cu), lo_tail = lower && g.cl >= 0;
  if (!hi_tail && !lo_tail)
    return {"P(reject | " + sym + "), 0 as no count is rare enough to reject", "P(reject | " + mean + ")",
            dvp(total, against)};
  std::string f, t;
  if (hi_tail) {
    f = "P(X ≥ c_U | " + sym + ")";
    t = "P(X ≥ " + dvi(g.cu) + " | " + mean + ")";
  }
  if (lo_tail) {
    f += std::string(hi_tail ? " + " : "") + "P(X ≤ c_L | " + sym + ")";
    t += std::string(hi_tail ? " + " : "") + "P(X ≤ " + dvi(g.cl) + " | " + mean + ")";
  }
  if (hi_tail && lo_tail) t += " = " + dvp(hi) + " + " + dvp(lo);
  return {f + far, t, dvp(total, against)};
}

// The power of the z-test of one rate by the normal approximation of one_pois_eval, X ~ N(μ1, μ1) against the
// critical values μ0 ± z√μ0: in symbols, its terms and the total.
PoisReject pois_normal_power(const Scenario& s, const OnePoisResult& res, double against = NAN) {
  using namespace derive;
  const Alt a = alt(s);
  const double z = z_alpha(s("alpha"), a), s0 = std::sqrt(res.mu0), s1 = std::sqrt(res.mu1);
  auto term = [&](double d) {
    return "Φ((" + dv(d) + " − " + dv(z) + " × √" + dv(res.mu0) + ") / √" + dv(res.mu1) + ")";
  };
  const std::string up = "Φ((μ1 − μ0 − z√μ0) / √μ1)", down = "Φ((μ0 − μ1 − z√μ0) / √μ1)";
  if (a == Alt::Two)
    return {up + " + " + down,
            term(res.dmu) + " + " + term(-res.dmu) + " = " + dvp(stats::pnorm((res.dmu - z * s0) / s1)) + " + " +
                dvp(stats::pnorm((-res.dmu - z * s0) / s1)),
            dvp(res.power, against)};
  if (a == Alt::Upper) return {up, term(res.dmu), dvp(res.power, against)};
  return {down, term(-res.dmu), dvp(res.power, against)};
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const double n = r.num("n"), t = s("time"), lam0 = r.num("lam0"), lam1 = r.num("lam1"), alpha = s("alpha");
    const Alt a = alt(s);
    const bool solved = s.solve == "n", effect = s.solve == "lam1";
    const OnePoisResult res = one_pois_eval(s, n, lam1);
    // λ1 with the digits that tell it from λ0 (which a rate found next to a large λ0 can need)
    const auto texts = apart(lam1, lam0);
    const std::string l1 = effect ? derived(texts.first) : given(texts.first), l0 = given(texts.second);
    const std::string N = size(n, solved), T = gsig(t);
    const std::string mu0 = "μ0 = " + l0 + " × " + N + " × " + T + " = " + dv(res.mu0);
    std::string mu1 = "μ1 = " + l1 + " × " + N + " × " + T + " = " + dv(res.mu1);
    if (res.far) mu1 += ", μ1 − μ0 = " + dv(res.dmu);  // which one_pois_eval takes from the rates' difference
    const std::string mu1_formula = res.far ? "μ1 = λ1 N t, μ1 − μ0 = (λ1 − λ0) N t" : "μ1 = λ1 N t";
    // the power in symbols and with its terms at μ1 (Poisson tails, or the normal approximation)
    const double target = s.solve == "power" ? NAN : s("power");
    const PoisReject pw = res.exact ? pois_reject(s, res, "μ1", res.mu1, res.dmu, res.power, target)
                                    : pois_normal_power(s, res, target);
    Derivation out;
    // the critical value or counts, and the actual size of an exact power: all from μ0 alone
    const auto region = [&] {
      if (!s.is("test", "exact")) out.steps.push_back(z_critical(alpha, a));
      if (!res.exact) return;
      out.steps.push_back(pois_region_step(s, res));
      const PoisReject size = pois_reject(s, res, "μ0", res.mu0, 0, res.actual_alpha);
      out.steps.push_back({"Actual significance level", "α* = " + size.formula, size.terms + " = " + size.total});
    };
    if (effect) {  // λ1 is where the power at μ1 = λ1 N t reaches the target
      out.find = find_effect("λ1", target, "with N = " + gvi(n));
      out.steps.push_back({"Expected events under H0", "μ0 = λ0 N t", mu0});
      region();
      out.steps.push_back({"Expected events at λ1", mu1_formula, mu1});
      out.steps.push_back(effect_found("λ1", l1, pw.formula + ", μ1 = λ1 N t", pw.terms, res.power, target));
      out.answer = "λ1 = " + l1 + " (λ1/λ0 = " + derived(apart(lam1 / lam0, 1).first) + ")";
      return out;
    }
    if (solved) {
      // the closed form of the normal approximation, beside the answer the steps below give
      out.find = find_size("N", target);
      const double za = z_alpha(alpha, a), zb = stats::qnorm(target);
      const double g = std::pow((za * std::sqrt(lam0) + zb * std::sqrt(lam1)) / (lam1 - lam0), 2) / t;
      out.steps.push_back({"Normal approximation",
                           std::string("N ≈ (z(") + level_symbol(a) +
                               ") √λ0 + z(1 − β) √λ1)² / (t (λ1 − λ0)²), for comparison; the steps below give the "
                               "answer",
                           "(" + dv(za) + " × √" + l0 + " + " + dv(zb) + " × √" + l1 + ")² / (" + T + " × (" + l1 +
                               " − " + l0 + ")²) = " + dv(g)});
    } else {
      out.find = find_power("with N = " + gvi(n));
    }
    out.steps.push_back({"Expected events", "μ0 = λ0 N t, " + mu1_formula, mu0 + ", " + mu1});
    region();
    out.steps.push_back({"Power", pw.formula, pw.terms + " = " + pw.total});
    if (solved) {
      const auto power = [&](double x) { return one_pois_eval(s, x, lam1).power; };
      out.steps.push_back(res.exact ? sawtooth_n_search(power, n, target, 1) : n_search(power, n, target, 1));
      out.answer = size_answer(n, res.power);
    } else {
      out.answer = "power = " + dvp(res.power);
    }
    return out;
  };
  p.definitions = {{"N", "Number of subjects."},
                   {"t", "Exposure time per subject; N·t is the total exposure."},
                   {"λ0, λ1", "Event rates per unit of exposure under H0 and H1."},
                   {"μ0, μ1", "Expected numbers of events λ·N·t under H0 and H1."},
                   {"Actual Alpha", "Exact probability of rejecting H0 when λ = λ0 (discreteness makes it ≤ target)."},
                   {"Reject H0 if", "Rejection region for the observed number of events X."}};
  p.notes = {"Because the number of events is discrete, exact power does not increase monotonically with N; the "
             "reported N is the smallest sample size whose power reaches the target."};
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
const char* SIZES_WITH_EVENTS =
    "Smaller samples leave a group so little chance of an event that the normal approximation gives the test more "
    "power than that chance allows; this is the smallest that does not.";

// Sizes `n` from an n search at which a side's power exceeds the chance of the events it needs (`possible` is false):
// the smallest from which none does, with a note, or `n` itself (refused by check_rate_events) when no sizes up to the
// search's limit qualify.
TwoN sizes_with_events(const Scenario& s, const std::function<double(double, double)>& P,
                       const std::function<bool(double, double)>& possible, const TwoN& n, Row& r,
                       const RunContext& ctx) {
  try {
    const TwoN m = solve_two_group_n(s, [&](double n1, double n2) { return possible(n1, n2) ? P(n1, n2) : 0.0; },
                                     s("power"), 1, n.n1, true, ctx);
    r.notes.push_back(SIZES_WITH_EVENTS);
    return m;
  } catch (const CalcError&) {
    return n;
  }
}

// The power of a two-rate z-test of the log ratio `lr` (written `sym`) from its variances under H0 and H1, with the
// chances `up` and `down` of rejecting on each side as the procedure's sides computed them, and their total.
Step log_ratio_power_step(const std::string& sym, double lr, double v0, double v1, std::pair<double, double> sides,
                          double power, double alpha, Alt a) {
  using namespace derive;
  const double z = z_alpha(alpha, a);
  const std::string up = "Φ((" + sym + " − z√V0) / √V1)", down = "Φ((−" + sym + " − z√V0) / √V1)";
  auto term = [&](double d) {
    return "Φ((" + dv(d) + " − " + dv(z) + " × √" + dsig(v0) + ") / √" + dsig(v1) + ")";
  };
  if (a == Alt::Two)
    return {"Power", up + " + " + down,
            term(lr) + " + " + term(-lr) + " = " + dvp(sides.first) + " + " + dvp(sides.second) + " = " + dvp(power)};
  if (a == Alt::Upper) return {"Power", up, term(lr) + " = " + dvp(power)};
  return {"Power", down, term(-lr) + " = " + dvp(power)};
}

// The closed form of a two-rate test's group sizes where group 2 follows group 1 (equal groups or a ratio): with
// v0 = N1 V0 and v1 = N1 V1, which do not change with N1 there, N1 ≈ (z √v0 + z(1 − β) √v1)² / (log ratio)², for
// comparison with the search (the form of Zhu and Lakkis's closed formula). M is the test (TwoPois or TwoNegbin).
template <class M>
void rate_closed_form(const Scenario& s, const M& m, double lr, const std::string& lr_sym, Derivation& out) {
  using namespace derive;
  const std::string& alloc = s.s("alloc");
  if (alloc != "equal" && alloc != "ratio") return;
  const double k = alloc == "ratio" ? s("r") : 1, v0 = m.v0(1, k), v1 = m.v1(1, k);
  const double za = z_alpha(m.alpha, m.a), zb = stats::qnorm(s("power"));
  const double g = std::pow(za * std::sqrt(v0) + zb * std::sqrt(v1), 2) / (lr * lr);
  out.steps.push_back(
      {"Normal approximation",
       std::string(alloc == "equal" ? "N1 = N2" : "N1") + " ≈ (z(" + level_symbol(m.a) + ") √v0 + z(1 − β) √v1)² / (" +
           lr_sym + ")², v0 = N1 V0 and v1 = N1 V1 from the variances below at N1 = 1, N2 = " +
           (alloc == "equal" ? "1" : "r") + ", for comparison; the steps below give the answer",
       "(" + dv(za) + " × √" + dsig(v0) + " + " + dv(zb) + " × √" + dsig(v1) + ")² / " + dvs(lr) + "² = " + dv(g)});
}

// How a two-rate test's group sizes were found: as two_group_search shows it, or, where sizes_with_events moved them
// on, the smallest whose power reaches the target without a side's power above its group's chance of an event.  M is
// the test (TwoPois or TwoNegbin).
template <class M>
Step rate_sizes_step(const Scenario& s, const Row& r, const M& m) {
  using namespace derive;
  const TwoN n{r.num("n1"), r.num("n2")};
  const double target = s("power");
  const auto P = [&](double n1, double n2) { return m.power(n1, n2); };
  if (std::find(r.notes.begin(), r.notes.end(), SIZES_WITH_EVENTS) == r.notes.end())
    return two_group_search(s, P, n, target, 1);
  Step st = two_group_search(
      s, [&](double n1, double n2) { return m.possible(n1, n2) ? P(n1, n2) : 0.0; }, n, target, 1);
  st.formula += ", where no side's power exceeds its group's chance of an event (the normal approximation fails there)";
  const TwoN less = smaller_two_group(s, n);
  if (less.n1 >= 1 && less.n2 >= 1 && !m.possible(less.n1, less.n2)) {
    // the side that the sizes below fail on, and the chance of an event in the group it needs one in
    const auto [up, down] = m.sides(less.n1, less.n2);
    const auto [none1, none2] = m.log_none_groups(less.n1, less.n2);
    const bool first = beyond_events(up, none1);
    const double side = first ? up : down, chance = -std::expm1(first ? none1 : none2);
    st.numbers = "at " + two_sizes(less, true) + " the normal approximation gives group " + (first ? "1" : "2") +
                 "'s side a power of " + dvp(side, chance) + ", above that group's chance of an event, " +
                 dvp(chance, side) + "; power at " + two_sizes(n, true) + ": " + dvp(P(n.n1, n.n2), target) + " ≥ " +
                 gv(target);
  }
  return st;
}

// The z-test of the log ratio of two Poisson rates: its variances under H1 and H0, and the chances of rejecting with
// group 1's rate the higher and with group 2's (0 for a side the test does not have).
struct TwoPois {
  double l1, l2, t, alpha, lr;
  Alt a;
  bool wald;  // the variance under H0 is the one under H1
  explicit TwoPois(const Scenario& s)
      : l1(s("lam1")), l2(s("lam2")), t(s("time")), alpha(s("alpha")), lr(std::log(l1 / l2)), a(alt(s)),
        wald(s.is("method", "h1")) {}
  double v1(double n1, double n2) const { return 1 / (n1 * t * l1) + 1 / (n2 * t * l2); }
  double pooled(double n1, double n2) const { return (n1 * l1 + n2 * l2) / (n1 + n2); }
  double v0(double n1, double n2) const {
    return wald ? v1(n1, n2) : (1 / (n1 * t) + 1 / (n2 * t)) / pooled(n1, n2);
  }
  std::pair<double, double> sides(double n1, double n2) const {
    const double var1 = v1(n1, n2), var0 = v0(n1, n2);
    const double zt = stats::qnorm(tail_alpha(alpha, a), false);
    auto tail = [&](double d) { return stats::pnorm((d - zt * std::sqrt(var0)) / std::sqrt(var1)); };
    return {a != Alt::Lower ? tail(lr) : 0.0, a != Alt::Upper ? tail(-lr) : 0.0};
  }
  double power(double n1, double n2) const {
    const auto [up, down] = sides(n1, n2);
    return a == Alt::Two ? stats::clamp01(up + down) : up + down;  // a NaN stays NaN (std::min(1.0, x) gives 1)
  }
  // the logs of each group's chance of no event: e^−events for Poisson counts
  std::pair<double, double> log_none_groups(double n1, double n2) const { return {-n1 * l1 * t, -n2 * l2 * t}; }
  // no side's power above its group's chance of an event (see beyond_events)
  bool possible(double n1, double n2) const {
    const auto [up, down] = sides(n1, n2);
    const auto [none1, none2] = log_none_groups(n1, n2);
    return !beyond_events(up, none1) && !beyond_events(down, none2);
  }
};

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
    const double l1 = s("lam1"), l2 = s("lam2"), t = s("time"), alpha = s("alpha");
    const Alt a = alt(s);
    // the variances' reciprocals of exposures below about 5.6e-309 overflow to infinity (and their ratio, the
    // statistic, is then no number)
    if (!std::isfinite(1 / t)) throw CalcError("the exposure time t is too small to compute with");
    if (!std::isfinite(1 / (t * std::fmin(l1, l2))))
      throw CalcError("the expected events per subject (λ × t) are too small to compute with");
    const TwoPois m(s);
    const double lr = m.lr;
    auto sides = [&](double n1, double n2) { return m.sides(n1, n2); };
    auto P = [&](double n1, double n2) { return m.power(n1, n2); };
    TwoN n;
    if (s.solve == "n") {
      if (l1 == l2) throw CalcError("the rates are equal");
      if ((a == Alt::Upper && (l1 - l2) < 0) || (a == Alt::Lower && (l1 - l2) > 0))
        throw CalcError("the rate ratio λ1/λ2 lies on the side opposite to the one-sided alternative");
      const double g = z_factor(alpha, s("power"), a) * (1 / l1 + 1 / l2) / t / (lr * lr);
      n = solve_two_group_n(s, P, s("power"), 1, std::ceil(g), true, ctx);
      auto possible = [&](double n1, double n2) { return m.possible(n1, n2); };
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const TwoPois m(s);
    const bool solved = s.solve == "n";
    const TwoN n{r.num("n1"), r.num("n2")};
    const double v1 = m.v1(n.n1, n.n2), v0 = m.v0(n.n1, n.n2), power = m.power(n.n1, n.n2);
    const std::string N1 = size(n.n1, solved), N2 = size(n.n2, solved), T = gsig(m.t), L1 = gsig(m.l1),
                      L2 = gsig(m.l2);
    Derivation out;
    out.steps.push_back({"Log rate ratio", "ln(λ1/λ2)", "ln(" + L1 + " / " + L2 + ") = " + dv(m.lr)});
    if (solved) {
      out.find = find_size(two_group_unknown(s), s("power"));
      rate_closed_form(s, m, m.lr, "ln(λ1/λ2)", out);
    } else {
      out.find = find_power("with " + two_sizes(n, false));
    }
    out.steps.push_back({"Variance under H1", "V1 = 1/(N1 t λ1) + 1/(N2 t λ2), one over each group's expected events",
                         "1/(" + N1 + " × " + T + " × " + L1 + ") + 1/(" + N2 + " × " + T + " × " + L2 + ") = " +
                             dsig(v1)});
    if (m.wald) {
      out.steps.push_back({"Variance under H0", "V0 = V1 (Wald)", "V0 = " + dsig(v0)});
    } else {
      const double pooled = m.pooled(n.n1, n.n2);
      out.steps.push_back({"Variance under H0",
                           "V0 = (1/(N1 t) + 1/(N2 t)) / λ̄, λ̄ = (N1 λ1 + N2 λ2) / (N1 + N2), the pooled rate",
                           "λ̄ = (" + N1 + " × " + L1 + " + " + N2 + " × " + L2 + ") / (" + N1 + " + " + N2 + ") = " +
                               dsig(pooled) + ", V0 = (1/(" + N1 + " × " + T + ") + 1/(" + N2 + " × " + T + ")) / " +
                               dsig(pooled) + " = " + dsig(v0)});
    }
    out.steps.push_back(z_critical(m.alpha, m.a));
    out.steps.push_back(log_ratio_power_step("ln(λ1/λ2)", m.lr, v0, v1, m.sides(n.n1, n.n2), power, m.alpha, m.a));
    if (solved) {
      out.steps.push_back(rate_sizes_step(s, r, m));
      out.answer = sizes_answer(n, "power " + dvp(power));
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
  };
  p.definitions = {{"λ1, λ2", "Event rates per unit of exposure."}, {"t", "Exposure time per subject."}};
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
// The test of the log rate ratio: its variances under H1 and H0 (by Zhu and Lakkis's method 1, 2 or 3), and the
// chances of rejecting with the treatment rate the higher and with the control rate (0 for a side the test does not
// have).  Group 1 is treatment (rate r1), group 2 control (rate r0).
struct TwoNegbin {
  double r0, rr, r1, kappa, mt, alpha;
  Alt a;
  std::string method;
  explicit TwoNegbin(const Scenario& s)
      : r0(s("r0")), rr(s("rr")), r1(r0 * rr), kappa(s("kappa")), mt(s("time")), alpha(s("alpha")), a(alt(s)),
        method(s.s("method")) {}
  double v1(double n1, double n2) const { return (1 / (mt * r1) + kappa) / n1 + (1 / (mt * r0) + kappa) / n2; }
  double pooled(double n1, double n2) const { return (n1 * r1 + n2 * r0) / (n1 + n2); }
  double v0(double n1, double n2) const {
    if (method == "1") return (1 / (mt * r0) + kappa) * (1 / n1 + 1 / n2);
    if (method == "2") return v1(n1, n2);
    return (1 / (mt * pooled(n1, n2)) + kappa) * (1 / n1 + 1 / n2);
  }
  std::pair<double, double> sides(double n1, double n2) const {
    const double var1 = v1(n1, n2), var0 = v0(n1, n2);
    const double lr = std::log(rr);
    const double zt = stats::qnorm(tail_alpha(alpha, a), false);
    auto tail = [&](double d) { return stats::pnorm((d - zt * std::sqrt(var0)) / std::sqrt(var1)); };
    return {a != Alt::Lower ? tail(lr) : 0.0, a != Alt::Upper ? tail(-lr) : 0.0};
  }
  double power(double n1, double n2) const {
    const auto [up, down] = sides(n1, n2);
    return a == Alt::Two ? stats::clamp01(up + down) : up + down;  // a NaN stays NaN (std::min(1.0, x) gives 1)
  }
  // the chance that a subject with mean count μ has no event: (1 + κμ)^(−1/κ), e^(−μ) for κ = 0
  double log_none(double mu) const { return kappa > 0 ? -std::log1p(kappa * mu) / kappa : -mu; }
  std::pair<double, double> log_none_groups(double n1, double n2) const {
    return {n1 * log_none(r1 * mt), n2 * log_none(r0 * mt)};
  }
  // no side's power above its group's chance of an event (see beyond_events)
  bool possible(double n1, double n2) const {
    const auto [up, down] = sides(n1, n2);
    const auto [none1, none2] = log_none_groups(n1, n2);
    return !beyond_events(up, none1) && !beyond_events(down, none2);
  }
};

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
    const double r0 = s("r0"), rr = s("rr"), kappa = s("kappa"), mt = s("time"), alpha = s("alpha");
    const double r1 = r0 * rr;
    const Alt a = alt(s);
    // the variances' reciprocals of exposures below about 5.6e-309 overflow to infinity (and their ratio, the
    // statistic, is then no number)
    if (!std::isfinite(1 / (mt * std::fmin(r0, r1))))
      throw CalcError("the expected events per subject (rate × time) are too small to compute with");
    const TwoNegbin m(s);
    auto sides = [&](double n1, double n2) { return m.sides(n1, n2); };
    auto P = [&](double n1, double n2) { return m.power(n1, n2); };
    auto log_none = [&](double mu) { return m.log_none(mu); };
    TwoN n;
    if (s.solve == "n") {
      if (rr == 1) throw CalcError("a rate ratio of 1 has no effect to detect");
      if ((a == Alt::Upper && (rr - 1) < 0) || (a == Alt::Lower && (rr - 1) > 0))
        throw CalcError("the rate ratio R lies on the side opposite to the one-sided alternative");
      const double v = (1 / (mt * r1) + kappa) + (1 / (mt * r0) + kappa);
      n = solve_two_group_n(s, P, s("power"), 1, std::ceil(z_factor(alpha, s("power"), a) * v / (std::log(rr) * std::log(rr))),
                            true, ctx);
      auto possible = [&](double n1, double n2) { return m.possible(n1, n2); };
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const TwoNegbin m(s);
    const bool solved = s.solve == "n";
    const TwoN n{r.num("n1"), r.num("n2")};
    const double v1 = m.v1(n.n1, n.n2), v0 = m.v0(n.n1, n.n2), power = m.power(n.n1, n.n2), lr = std::log(m.rr);
    const std::string N1 = size(n.n1, solved), N2 = size(n.n2, solved);
    // 1/(t λ) + κ, the variance of a subject's log rate in a group, written out (`rate` as written)
    auto per = [&](const std::string& rate) {
      return "(1/(" + gsig(m.mt) + " × " + rate + ") + " + gsig(m.kappa) + ")";
    };
    const std::string inv = "(1/" + N1 + " + 1/" + N2 + ")";
    Derivation out;
    out.steps.push_back({"Treatment rate", "λ1 = R × λ2, and ln R the log rate ratio tested",
                         "λ1 = " + gv(m.rr) + " × " + gsig(m.r0) + " = " + dsig(m.r1) + ", ln R = ln " + gv(m.rr) +
                             " = " + dv(lr)});
    if (solved) {
      out.find = find_size(two_group_unknown(s), s("power"));
      rate_closed_form(s, m, lr, "ln R", out);
    } else {
      out.find = find_power("with " + two_sizes(n, false));
    }
    out.steps.push_back({"Variance under H1", "V1 = (1/(t λ1) + κ)/N1 + (1/(t λ2) + κ)/N2",
                         per(dsig(m.r1)) + "/" + N1 + " + " + per(gsig(m.r0)) + "/" + N2 + " = " + dsig(v1)});
    if (m.method == "1") {
      out.steps.push_back({"Variance under H0 (method 1)", "V0 = (1/(t λ2) + κ)(1/N1 + 1/N2), from the control rate",
                           per(gsig(m.r0)) + " × " + inv + " = " + dsig(v0)});
    } else if (m.method == "2") {
      out.steps.push_back({"Variance under H0 (method 2)", "V0 = V1, from the true rates", "V0 = " + dsig(v0)});
    } else {
      const double pooled = m.pooled(n.n1, n.n2);
      out.steps.push_back({"Variance under H0 (method 3)",
                           "V0 = (1/(t λ̄) + κ)(1/N1 + 1/N2), λ̄ = (N1 λ1 + N2 λ2) / (N1 + N2), the pooled rate",
                           "λ̄ = (" + N1 + " × " + dsig(m.r1) + " + " + N2 + " × " + gsig(m.r0) + ") / (" + N1 +
                               " + " + N2 + ") = " + dsig(pooled) + ", V0 = " + per(dsig(pooled)) + " × " + inv +
                               " = " + dsig(v0)});
    }
    out.steps.push_back(z_critical(m.alpha, m.a));
    out.steps.push_back(log_ratio_power_step("ln R", lr, v0, v1, m.sides(n.n1, n.n2), power, m.alpha, m.a));
    if (solved) {
      out.steps.push_back(rate_sizes_step(s, r, m));
      out.answer = sizes_answer(n, "power " + dvp(power));
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
  };
  p.definitions = {{"λ1, λ2", "Event rates in the treatment and control groups."},
                   {"R", "Rate ratio λ1/λ2."},
                   {"κ", "Negative binomial dispersion (Var = μ + κμ²; κ = 0 is Poisson)."}};
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
// The exact power with n subjects at the variance ratio σ1²/σ0²: (N − 1)S²/σ0² is σ1²/σ0² times a chi-square
// variable with N − 1 degrees of freedom.
double one_var_power(double n, double ratio, double alpha, Alt a) {
  using namespace stats;
  const double df = n - 1;
  double pw = 0;
  const double side = tail_alpha(alpha, a);
  if (a != Alt::Lower) pw += pchisq(qchisq(side, df, false) / ratio, df, false);
  if (a != Alt::Upper) pw += pchisq(qchisq(side, df, true) / ratio, df, true);
  return clamp01(pw);
}

// The critical values and the power of a test of variances whose statistic, divided by the variance ratio
// `ratio_sym`, has the law `q` (χ² or F, its degrees of freedom `df`, defined by `df_def`) under H1: c_U and c_L are
// its quantiles (quantile(lower, p)), and the power the chances beyond them once divided by the ratio (tail(lower, x),
// as the power function computes them).
// `ratio_text` is the ratio as written (given, or derived from the SDs), `df_numbers` the degrees of freedom.
void variance_test_steps(const std::string& q, const std::string& df, const std::string& df_def,
                         const std::string& df_numbers, const std::string& ratio_sym, const std::string& ratio_text,
                         double ratio, double alpha, Alt a, double power,
                         const std::function<double(bool, double)>& quantile,
                         const std::function<double(bool, double)>& tail, derive::Steps& out) {
  using namespace derive;
  const double side = tail_alpha(alpha, a);
  const std::string level = a == Alt::Two ? "α/2" : "α", level_lo = gv(alpha) + (a == Alt::Two ? "/2" : "");
  const bool hi = a != Alt::Lower, lo = a != Alt::Upper;
  const double cu = hi ? quantile(false, side) : NAN, cl = lo ? quantile(true, side) : NAN;
  std::string f, t, pf, pt;
  if (hi) {
    f = "c_U = " + q + "(1 − " + level + "; " + df + ")";
    t = "c_U = " + q + "(" + level_numbers(alpha, a) + "; " + df_numbers + ") = " + dv(cu);
    pf = "P(" + q + " > c_U / " + ratio_sym + ")";
    pt = "P(" + q + " > " + dv(cu) + " / " + ratio_text + ")";
  }
  if (lo) {
    const std::string join = hi ? ", " : "", plus = hi ? " + " : "";
    f += join + "c_L = " + q + "(" + level + "; " + df + ")";
    t += join + "c_L = " + q + "(" + level_lo + "; " + df_numbers + ") = " + dv(cl);
    pf += plus + "P(" + q + " < c_L / " + ratio_sym + ")";
    pt += plus + "P(" + q + " < " + dv(cl) + " / " + ratio_text + ")";
  }
  out.push_back({"Critical values", f + ", the " + q + " quantiles with " + df_def, t});
  pt += " = ";
  if (hi && lo) pt += dvp(tail(false, cu / ratio)) + " + " + dvp(tail(true, cl / ratio)) + " = ";
  out.push_back({"Power", pf + ", " + q + " with " + df + " degrees of freedom", pt + dvp(power)});
}

// The normal approximation of a test of variances on the log scale, where ln S² has a variance of about 2/(N − 1):
// the size whose log ratio ln(`ratio`) the test detects with the target power, N − 1 ≈ c (z + z(1 − β))² / (ln
// ratio)², for comparison (`c` 2 for one variance, 4 for two equal groups, 2 (1 + 1/r) for a ratio, written
// `c_formula` and `c_numbers`).
Step log_variance_approx(const std::string& lhs, const std::string& ratio_sym, const std::string& ratio_text,
                         double ratio, const std::string& c_formula, const std::string& c_numbers, double c,
                         double alpha, Alt a, double target) {
  using namespace derive;
  const double za = z_alpha(alpha, a), zb = stats::qnorm(target), lr = std::log(ratio);
  return {"Normal approximation (log scale)",
          lhs + " ≈ 1 + " + c_formula + " " + z_sum_formula(a) + " / (ln " + ratio_sym +
              ")², as ln S² has a variance of about 2/(N − 1), for comparison; the steps below give the answer",
          "1 + " + c_numbers + " × " + z_sum_numbers(alpha, a, target) + " / (ln " + ratio_text + ")² = " +
              dv(1 + c * (za + zb) * (za + zb) / (lr * lr))};
}

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
    const double ratio = s("sd1") * s("sd1") / (s("sd0") * s("sd0")), alpha = s("alpha");
    const Alt a = alt(s);
    auto P = [&](double n) { return one_var_power(n, ratio, alpha, a); };
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const double n = r.num("n"), ratio = r.num("ratio"), alpha = s("alpha"), df = n - 1;
    const Alt a = alt(s);
    const bool solved = s.solve == "n";
    const double power = one_var_power(n, ratio, alpha, a);
    Derivation out;
    out.steps.push_back(
        {"Variance ratio", "σ1²/σ0²", gsig(s("sd1")) + "² / " + gsig(s("sd0")) + "² = " + dsig(ratio)});
    if (solved) {
      out.find = find_size("N", s("power"));
      out.steps.push_back(
          log_variance_approx("N", "(σ1²/σ0²)", dsig(ratio), ratio, "2", "2", 2, alpha, a, s("power")));
    } else {
      out.find = find_power("with N = " + gvi(n));
    }
    out.steps.push_back({"Degrees of freedom", "ν = N − 1", size(n, solved) + " − 1 = " + dvi(df)});
    variance_test_steps(
        "χ²", "ν", "ν = N − 1", dvi(df), "(σ1²/σ0²)", dsig(ratio), ratio, alpha, a, power,
        [&](bool lower, double p) { return stats::qchisq(p, df, lower); },
        [&](bool lower, double x) { return stats::pchisq(x, df, lower); }, out.steps);
    if (solved) {
      out.steps.push_back(n_search([&](double x) { return one_var_power(x, ratio, alpha, a); }, n, s("power"), 2));
      out.answer = size_answer(n, power);
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
  };
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
// The exact power with group sizes n1 and n2 at the variance ratio σ1²/σ2²: S1²/S2² is σ1²/σ2² times an F variable
// with N1 − 1 and N2 − 1 degrees of freedom.
double two_var_power(double n1, double n2, double ratio, double alpha, Alt a) {
  using namespace stats;
  const double d1 = n1 - 1, d2 = n2 - 1;
  const double side = tail_alpha(alpha, a);
  double pw = 0;
  if (a != Alt::Lower) pw += pf(qf(side, d1, d2, false) / ratio, d1, d2, false);
  if (a != Alt::Upper) pw += pf(qf(side, d1, d2, true) / ratio, d1, d2, true);
  return clamp01(pw);
}

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
    const double ratio = s("ratio"), alpha = s("alpha");
    const Alt a = alt(s);
    auto P = [&](double n1, double n2) { return two_var_power(n1, n2, ratio, alpha, a); };
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
  p.derivation = [](const Scenario& s, const Row& r) {
    using namespace derive;
    const TwoN n{r.num("n1"), r.num("n2")};
    const double ratio = r.num("ratio"), alpha = s("alpha"), d1 = n.n1 - 1, d2 = n.n2 - 1;
    const Alt a = alt(s);
    const bool solved = s.solve == "n";
    const double power = two_var_power(n.n1, n.n2, ratio, alpha, a);
    Derivation out;
    if (solved) {
      out.find = find_size(two_group_unknown(s), s("power"));
      // ln(S1²/S2²) has a variance of about 2/ν1 + 2/ν2: 4/ν1 with equal groups, 2 (1 + 1/r)/ν1 with N2 = r N1
      const std::string& alloc = s.s("alloc");
      if (alloc == "equal")
        out.steps.push_back(
            log_variance_approx("N1 = N2", "(σ1²/σ2²)", gsig(ratio), ratio, "4", "4", 4, alpha, a, s("power")));
      else if (alloc == "ratio")
        out.steps.push_back(log_variance_approx("N1", "(σ1²/σ2²)", gsig(ratio), ratio, "2 (1 + 1/r)",
                                                "2 × (1 + 1/" + gv(s("r")) + ")", 2 * (1 + 1 / s("r")), alpha, a,
                                                s("power")));
    } else {
      out.find = find_power("with " + two_sizes(n, false));
    }
    out.steps.push_back({"Degrees of freedom", "ν1 = N1 − 1, ν2 = N2 − 1",
                         "ν1 = " + size(n.n1, solved) + " − 1 = " + dvi(d1) + ", ν2 = " + size(n.n2, solved) +
                             " − 1 = " + dvi(d2)});
    variance_test_steps(
        "F", "ν1, ν2", "ν1 = N1 − 1 and ν2 = N2 − 1", dvi(d1) + ", " + dvi(d2), "(σ1²/σ2²)", gsig(ratio), ratio,
        alpha, a, power, [&](bool lower, double p) { return stats::qf(p, d1, d2, lower); },
        [&](bool lower, double x) { return stats::pf(x, d1, d2, lower); }, out.steps);
    if (solved) {
      out.steps.push_back(two_group_search(
          s, [&](double x1, double x2) { return two_var_power(x1, x2, ratio, alpha, a); }, n, s("power"), 2));
      out.answer = sizes_answer(n, "power " + dvp(power));
    } else {
      out.answer = "power = " + dvp(power);
    }
    return out;
  };
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
