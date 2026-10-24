// Simulation procedures for further designs: two means with non-normal data,
// several treatments against a control, ordinal outcomes, the win ratio,
// repeated measures with dropout (MMRM) and stepped-wedge trials.  Simulated
// power is checked against exact power where one exists, against R (the
// values quoted were simulated with R's own functions, as stated), and for
// its type I error, reproducibility and the sample-size search.
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "engine/runner.hpp"
#include "procedures/sim_mmrm.hpp"
#include "procedures/sim_ordinal.hpp"
#include "procedures/sim_winratio.hpp"
#include "sim/random.hpp"
#include "stats/stats.hpp"
#include "testing.hpp"
#include "util/json.hpp"
#include "util/numtext.hpp"

using namespace ssc;

namespace {

json::Value run(const std::string& request) {
  Token token;
  return run_request(json::parse(request), &token);
}

json::Value first_row(const std::string& request) {
  const json::Value res = run(request);
  CHECK(res.get("ok").as_bool(false));
  return res.get("rows").elements().at(0);
}

double num(const json::Value& row, const char* key) { return row.get(key).as_number(NAN); }

std::string power_request(const std::string& procedure, const std::string& inputs, int threads = 0) {
  return R"({"procedure":")" + procedure + R"(","solve_for":"power","options":{"power_curve":false,"threads":)" +
         std::to_string(threads) + R"(},"inputs":{)" + inputs + "}}";
}

/// 4 standard errors of the difference between two simulated rates p (R1 and R2 studies).
double tolerance(double p, double r1, double r2) { return 4 * std::sqrt(p * (1 - p) * (1 / r1 + 1 / r2)); }

}  // namespace

// ------------------------------------------------------------------ two means with non-normal data

TEST_CASE("robust: normal data reproduce the exact t-test power and P(X1 > X2)") {
  const json::Value r = first_row(power_request(
      "sim.means.robust", R"("dist":"normal","test":"t","n1":"25","mu1":"0.7","sd1":"1","sd2":"1","reps":"20000")"));
  const json::Value e = first_row(
      R"({"procedure":"means.two.t","solve_for":"power","inputs":{"n1":"25","mu1":"0.7","mu2":"0","sigma":"1"}})");
  CHECK_CLOSE(num(r, "power"), num(e, "power"), 0, 4 * num(r, "power_se"));
  CHECK_CLOSE(num(r, "approx_power"), num(e, "power"), 1e-9, 0);
  CHECK_CLOSE(num(r, "pgt"), stats::pnorm(0.7 / std::sqrt(2.0)), 1e-8, 0);
  CHECK_CLOSE(num(r, "mean_delta"), 0.7, 0, 0.01);
}

TEST_CASE("robust: simulated power agrees with R's t.test and wilcox.test") {
  // R: 4,000 studies of 30 + 30 subjects, means 1 and 0, SDs 2, each analyzed with t.test (pooled and Welch)
  // and wilcox.test(exact = FALSE, correct = TRUE)
  struct Case {
    const char* inputs;
    double t, welch, mw;
  };
  const Case cases[] = {{R"("dist":"lognormal","skew":"2")", 0.507, 0.506, 0.678},
                        {R"("dist":"gamma","skew":"-1.5")", 0.491, 0.491, 0.626},
                        {R"("dist":"t","df":"3")", 0.564, 0.563, 0.732},
                        {R"("dist":"contaminated","eps":"0.1","kout":"5")", 0.538, 0.536, 0.823},
                        {R"("dist":"uniform")", 0.466, 0.466, 0.424}};
  for (const Case& c : cases) {
    const json::Value r = first_row(power_request(
        "sim.means.robust", std::string(c.inputs) + R"(,"n1":"30","mu1":"1","sd1":"2","sd2":"2","reps":"20000")"));
    CHECK_CLOSE(num(r, "power_t"), c.t, 0, tolerance(c.t, 4000, 20000));
    CHECK_CLOSE(num(r, "power_welch"), c.welch, 0, tolerance(c.welch, 4000, 20000));
    CHECK_CLOSE(num(r, "power_mw"), c.mw, 0, tolerance(c.mw, 4000, 20000));
    CHECK_EQ(num(r, "power"), num(r, "power_welch"));  // Welch is the default test
  }
}

TEST_CASE("robust: type I errors with skewed data, unequal SDs and unequal groups") {
  // R (20,000 studies): pooled t 0.0076, Welch 0.0505, Mann–Whitney 0.1734 (P(X1 > X2) ≠ 1/2 here)
  const json::Value r = first_row(power_request(
      "sim.means.robust",
      R"("dist":"lognormal","skew":"2","n1":"20","r":"3","mu1":"0","mu2":"0","sd1":"1","sd2":"3","reps":"20000")"));
  CHECK_CLOSE(num(r, "power_welch"), 0.05, 0, 0.0065);
  CHECK(num(r, "power_t") < 0.015);
  CHECK_CLOSE(num(r, "power_mw"), 0.1734, 0, tolerance(0.1734, 20000, 20000));
  CHECK(num(r, "pgt") > 0.5);  // group 2 is more spread out, so its median lies further below its mean
}

TEST_CASE("robust: reproducible across threads, and the sample-size search is minimal") {
  const std::string inputs = R"("dist":"t","df":"4","test":"mw","n1":"18","mu1":"1","reps":"600","seed":"5")";
  const json::Value one = first_row(power_request("sim.means.robust", inputs, 1));
  const json::Value four = first_row(power_request("sim.means.robust", inputs, 4));
  CHECK_EQ(num(one, "power"), num(four, "power"));
  CHECK_EQ(num(one, "mean_delta"), num(four, "mean_delta"));
  const json::Value n = first_row(
      R"({"procedure":"sim.means.robust","solve_for":"n","options":{"power_curve":false},)"
      R"("inputs":{"dist":"gamma","skew":"3","test":"mw","power":"0.85","reps":"1000","seed":"8"}})");
  CHECK(num(n, "power") >= 0.85);
  const json::Value below = first_row(power_request(
      "sim.means.robust", R"("dist":"gamma","skew":"3","test":"mw","reps":"1000","seed":"8","n1":")" +
                              text::format_g(num(n, "n1") - 1, 10) + "\""));
  CHECK(num(below, "power") < 0.85);
}

TEST_CASE("robust: designs without an effect to detect give clear messages") {
  const json::Value res = run(
      R"({"procedure":"sim.means.robust","solve_for":"n","inputs":{"dist":"normal","mu1":"0","mu2":"0"}})");
  const std::string err = res.get("rows").elements().at(0).get("_error").to_text();
  CHECK(err.find("μ1 equals μ2") != std::string::npos);
}

// ------------------------------------------------------------------ several treatments against a control

TEST_CASE("multiarm: Dunnett's critical values give a familywise error of exactly alpha") {
  // critical values from R's mvtnorm::qmvt (its quasi-Monte Carlo error is about 1e-4); mvtnorm::pmvt with
  // an absolute error of 1e-9 gives 0.95000 coverage at the engine's values
  struct Case {
    const char* inputs;
    double crit;
  };
  const Case cases[] = {{R"("mu":"0 0 0","n":"40")", 2.372068},
                        {R"("mu":"0 0 0","n":"20")", 2.396905},
                        {R"("mu":"0 0 0 0","n":"25","r0":"2","alternative":"upper")", 2.218782},
                        {R"("mu":"0 0 0 0 0 0","n":"12","r0":"2.5")", 2.663576}};
  for (const Case& c : cases) {
    const json::Value r = first_row(
        R"({"procedure":"sim.multiarm","solve_for":"power","options":{"approximate":true},"inputs":{)" +
        std::string(c.inputs) + "}}");
    CHECK_CLOSE(num(r, "crit"), c.crit, 0, 5e-4);
    CHECK_CLOSE(num(r, "power"), 0.05, 0, 1e-8);  // no effective treatment: the familywise error
  }
}

TEST_CASE("multiarm: Dunnett's critical values in the far tail, where 1 − alpha rounds to 1") {
  // references: the familywise error integrated in double precision over the common normal component and the
  // SD's chi distribution, with 1 − (1 − p)^k as −expm1(k·log1p(−p)) and the nodes following the tail (relative
  // error about 1e-12); before, the value was the unadjusted one below an alpha of about 1e-17, and lost digits
  // from about 1e-10
  struct Case {
    const char* inputs;
    double crit;
  };
  const Case cases[] = {
      {R"("mu":"0 0 0","n":"40","alpha":"1e-12")", 7.9592321202},              // Bonferroni's 7.9592764
      {R"("mu":"0 0 0","n":"804","alpha":"1e-30")", 11.742200475},             // Bonferroni's to 2e-13
      {R"("mu":"0 0 0","n":"800","r0":"0.1","alpha":"1e-30")", 11.7784991977},  // Bonferroni's 11.779508
      {R"("mu":"0 0 0","n":"2","alpha":"1e-30")", 61141209.83},  // 4 df: 6% below Bonferroni's, by the shared SD
      {R"("mu":"0 0 0 0 0 0 0 0","n":"3","r0":"10","alpha":"1e-30")", 30.3825182176}};
  for (const Case& c : cases) {
    const json::Value r = first_row(
        R"({"procedure":"sim.multiarm","solve_for":"power","options":{"approximate":true},"inputs":{)" +
        std::string(c.inputs) + "}}");
    CHECK_CLOSE(num(r, "crit"), c.crit, 1e-8, 0);
  }
  // at 1e-30 with three arms Dunnett's value is Bonferroni's to 13 digits, and so is the sample size
  auto n_for = [](const char* adjust) {
    return num(first_row(R"({"procedure":"sim.multiarm","solve_for":"n","options":{"approximate":true,)"
                         R"("power_curve":false},"inputs":{"alpha":"1e-30","adjust":")" +
                         std::string(adjust) + "\"}}"),
               "n");
  };
  CHECK_EQ(n_for("dunnett"), 817.0);
  CHECK_EQ(n_for("bonferroni"), 817.0);
  // with no effective treatment the power is the familywise error, also this far out (it was 0.29·alpha)
  for (const char* adjust : {"dunnett", "bonferroni"}) {
    const json::Value r = first_row(
        R"({"procedure":"sim.multiarm","solve_for":"power","options":{"approximate":true},)"
        R"("inputs":{"mu":"0 0 0","n":"817","alpha":"1e-30","adjust":")" + std::string(adjust) + "\"}}");
    CHECK_CLOSE(num(r, "power"), 1e-30, 1e-9, 0);
  }
}

TEST_CASE("multiarm: a small control arm and few degrees of freedom") {
  // references: the familywise error and the powers integrated over z on panels as fine as 1/8 of the width of an
  // arm's rejection, √(1 − λ²), and over the SD's chi distribution in log s (relative error about 1e-11).  With a
  // control arm a tenth of a treatment arm (λ = 0.95), the 64 nodes in z were too coarse for that width: the
  // critical value was 2.333920 (a familywise error of 0.05014) and the power that all eight arms reject 0.10833
  const json::Value a = first_row(
      R"({"procedure":"sim.multiarm","solve_for":"power","options":{"approximate":true},)"
      R"("inputs":{"mu":"0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3","n":"300","r0":"0.1","goal":"all"}})");
  CHECK_CLOSE(num(a, "crit"), 2.33508721495, 1e-10, 0);
  CHECK_CLOSE(num(a, "power_all"), 0.108874948216, 1e-8, 0);
  CHECK_CLOSE(num(a, "power_any"), 0.363652587996, 1e-8, 0);
  const json::Value z = first_row(
      R"({"procedure":"sim.multiarm","solve_for":"power","options":{"approximate":true},)"
      R"("inputs":{"mu":"0 0 0 0 0 0 0 0","n":"300","r0":"0.1"}})");
  CHECK_CLOSE(num(z, "power"), 0.05, 1e-10, 0);
  // 3 degrees of freedom (2 arms of 2): 26.453136 before
  const json::Value b = first_row(
      R"({"procedure":"sim.multiarm","solve_for":"power","options":{"approximate":true},)"
      R"("inputs":{"mu":"0 0","n":"2","alternative":"upper","alpha":"1e-4"}})");
  CHECK_CLOSE(num(b, "crit"), 26.4507712356, 1e-8, 0);
}

TEST_CASE("multiarm: powers of tiny effects at a tiny alpha, and with few degrees of freedom") {
  // references: the arm's power and that all arms reject, by nested adaptive quadrature (scipy's quad, relative
  // error about 1e-12).  Before, single_step's nodes stopped at z = 8.5 and at the SD's own range, short of where
  // such rejections lie: the power that all eight reject was 1.1e-65 for 8.2e-41, an arm's 3.6% of its value with
  // r0 = 0.1, and the power to find any of three arms 0.29·α at 1e-30 (a jump from α with no effect at all)
  auto request = [](const std::string& inputs) {
    return R"({"procedure":"sim.multiarm","solve_for":"power","options":{"approximate":true},"inputs":{)" + inputs +
           "}}";
  };
  auto arm = [](const json::Value& row) {
    return row.get("_tables").elements().at(0).get("rows").elements().at(0).elements().at(4).as_number(NAN);
  };
  const std::string tiny8 = R"("mu":"1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12")";
  const json::Value a = first_row(request(tiny8 + R"(,"n":"817","alternative":"upper","alpha":"1e-20")"));
  CHECK_CLOSE(num(a, "power_all"), 8.219122489484e-41, 1e-8, 0);
  CHECK_CLOSE(arm(a), 1.250000163088e-21, 1e-9, 0);
  const json::Value b = first_row(request(tiny8 + R"(,"n":"300","r0":"0.1","alpha":"1e-20")"));
  CHECK_CLOSE(num(b, "power_all"), 7.649290908681e-25, 1e-9, 0);
  CHECK_CLOSE(arm(b), 1.375460309792e-21, 1e-9, 0);
  // with effects of 1e-12 the power to find any arm is that with none, the familywise error
  for (const char* alpha : {"1e-16", "1e-30"}) {
    const json::Value near = first_row(request(R"("mu":"1e-12 1e-12 1e-12","n":"817","alpha":")" + std::string(alpha) +
                                               "\""));
    CHECK_CLOSE(num(near, "power_any"), std::stod(alpha), 1e-9, 0);
  }
  // one arm and 2 or 4 degrees of freedom: the noncentral t probability integrated with mpmath over S (48 nodes over
  // the SD's whole range in log s gave 0.0621924 for 0.0622364, 1.1% too much at 0.001 and 0.6% at 1e-8)
  const json::Value c = first_row(request(R"("mu":"1.5","n":"2","alternative":"upper","alpha":"0.01")"));
  CHECK_CLOSE(num(c, "power"), 6.223643770602e-02, 1e-10, 0);
  const json::Value e = first_row(request(R"("mu":"3","n":"2","alpha":"0.001")"));
  CHECK_CLOSE(num(e, "power"), 9.946206373842e-03, 1e-10, 0);
  const json::Value f = first_row(request(R"("mu":"3","n":"3","alpha":"1e-8")"));
  CHECK_CLOSE(num(f, "power"), 8.866865848107e-07, 1e-10, 0);
}

TEST_CASE("multiarm: exact powers of Dunnett's test agree with mvtnorm") {
  // R: 1 − pmvt(−c, c, delta, df = 156, corr 0.5) and the sum over the sign patterns of the tails
  const json::Value a = first_row(
      R"({"procedure":"sim.multiarm","solve_for":"power","options":{"approximate":true},"inputs":{"n":"40"}})");
  CHECK_CLOSE(num(a, "power_any"), 0.7834828, 0, 5e-5);
  CHECK_CLOSE(num(a, "power_all"), 0.1252485, 0, 5e-5);
  const json::Value b = first_row(
      R"({"procedure":"sim.multiarm","solve_for":"power","options":{"approximate":true},)"
      R"("inputs":{"n":"25","r0":"2","mu":"0.5 0.4 0 0.7","alternative":"upper"}})");
  CHECK_CLOSE(num(b, "power_any"), 0.8296233, 0, 5e-5);
}

TEST_CASE("multiarm: simulated powers agree with the exact ones and with R's step-wise procedures") {
  const json::Value d = first_row(power_request("sim.multiarm", R"("n":"40","reps":"20000")"));
  CHECK_CLOSE(num(d, "power_any"), 0.7834828, 0, 4 * std::sqrt(0.78 * 0.22 / 20000));
  CHECK_CLOSE(num(d, "power_all"), 0.1252485, 0, 4 * std::sqrt(0.125 * 0.875 / 20000));
  // R: 4,000 studies, lm() on the raw data and p.adjust(..., "holm" / "hochberg")
  const json::Value h = first_row(power_request("sim.multiarm", R"("n":"40","adjust":"holm","reps":"20000")"));
  CHECK_CLOSE(num(h, "power_any"), 0.7708, 0, tolerance(0.77, 4000, 20000));
  CHECK_CLOSE(num(h, "power_all"), 0.2127, 0, tolerance(0.21, 4000, 20000));
  const json::Value g = first_row(power_request("sim.multiarm", R"("n":"40","adjust":"hochberg","reps":"20000")"));
  CHECK_CLOSE(num(g, "power_any"), 0.7750, 0, tolerance(0.77, 4000, 20000));
  CHECK_CLOSE(num(g, "power_all"), 0.2253, 0, tolerance(0.22, 4000, 20000));
}

TEST_CASE("multiarm: on the same studies Hochberg finds at least what Holm finds, and Holm what Bonferroni finds") {
  double any[3], all[3];
  const char* methods[] = {"bonferroni", "holm", "hochberg"};
  for (int i = 0; i < 3; ++i) {
    const json::Value r = first_row(power_request(
        "sim.multiarm", R"("n":"30","mu":"0.5 0.4 0.3 0","reps":"4000","adjust":")" + std::string(methods[i]) + "\""));
    any[i] = num(r, "power_any");
    all[i] = num(r, "power_all");
    CHECK(num(r, "fwer") <= 0.05 + 4 * std::sqrt(0.05 * 0.95 / 4000));
  }
  CHECK(any[0] <= any[1] && any[1] <= any[2]);
  CHECK(all[0] <= all[1] && all[1] <= all[2]);
  CHECK(all[2] > all[0]);
}

TEST_CASE("multiarm: the familywise error of every procedure is at most alpha") {
  for (const char* adj : {"dunnett", "holm", "hochberg", "bonferroni"}) {
    const json::Value r = first_row(power_request(
        "sim.multiarm", R"("n":"15","mu":"0 0 0 0 0","reps":"40000","adjust":")" + std::string(adj) + "\""));
    CHECK(num(r, "power") <= 0.05 + 4 * std::sqrt(0.05 * 0.95 / 40000));
    if (std::string(adj) == "dunnett") CHECK_CLOSE(num(r, "power"), 0.05, 0, 4 * std::sqrt(0.05 * 0.95 / 40000));
  }
  const json::Value none = first_row(power_request("sim.multiarm", R"("n":"15","mu":"0 0 0 0 0","adjust":"none","reps":"4000")"));
  CHECK(num(none, "power") > 0.15);  // five unadjusted tests
}

TEST_CASE("multiarm: reproducible across threads, and the sample-size search is minimal") {
  const std::string inputs = R"("n":"25","adjust":"holm","goal":"all","reps":"800","seed":"3")";
  const json::Value one = first_row(power_request("sim.multiarm", inputs, 1));
  const json::Value four = first_row(power_request("sim.multiarm", inputs, 4));
  CHECK_EQ(num(one, "power"), num(four, "power"));
  const json::Value n = first_row(
      R"({"procedure":"sim.multiarm","solve_for":"n","options":{"power_curve":false},)"
      R"("inputs":{"adjust":"hochberg","goal":"each","mu":"0.5 0.6","power":"0.8","reps":"1000","seed":"4"}})");
  CHECK(num(n, "power") >= 0.8);
  const json::Value below = first_row(power_request(
      "sim.multiarm", R"("adjust":"hochberg","goal":"each","mu":"0.5 0.6","reps":"1000","seed":"4","n":")" +
                          text::format_g(num(n, "n") - 1, 10) + "\""));
  CHECK(num(below, "power") < 0.8);
}

// ------------------------------------------------------------------ ordinal outcomes

TEST_CASE("ordinal: the proportional odds fit equals R's clm and polr") {
  // R: ordinal::clm(y ~ g) (Newton's method, gradient tolerance 1e-12) on 40 + 40 subjects, g = 1 for the
  // first table; polr with reltol = 1e-15 agrees to 9 digits
  const double g1[4] = {5, 9, 14, 12}, g2[4] = {10, 15, 12, 3};
  const procs::ordinal::PoFit f = procs::ordinal::po_fit(g1, g2, 4);
  CHECK(f.ok && !f.separated);
  CHECK_CLOSE(f.beta, 1.1816474269, 0, 1e-9);
  CHECK_CLOSE(f.se, 0.4235270380, 0, 1e-9);
  CHECK_CLOSE(f.loglik, -104.281081407, 0, 1e-9);
  CHECK_CLOSE(f.loglik0, -108.3368, 0, 1e-4);
  // an empty category changes nothing
  const double e1[5] = {5, 9, 0, 14, 12}, e2[5] = {10, 15, 0, 12, 3};
  const procs::ordinal::PoFit e = procs::ordinal::po_fit(e1, e2, 5);
  CHECK_CLOSE(e.beta, f.beta, 1e-10, 0);
  CHECK_CLOSE(e.se, f.se, 1e-10, 0);
  // complete separation: an infinite estimate, and the likelihood of the two groups' own distributions
  const double s1[3] = {0, 0, 7}, s2[3] = {4, 3, 0};
  const procs::ordinal::PoFit s = procs::ordinal::po_fit(s1, s2, 3);
  CHECK(s.separated && std::isinf(s.beta) && s.beta > 0);
  CHECK_CLOSE(s.loglik, 4 * std::log(4.0 / 7) + 3 * std::log(3.0 / 7), 1e-12, 0);
}

TEST_CASE("ordinal: simulated power agrees with R's polr and wilcox.test") {
  // R: 4,000 studies each, polr Wald (pnorm of the t value) and likelihood-ratio tests, wilcox.test
  struct Case {
    const char* inputs;
    double wald, lr, mw;
  };
  const Case cases[] = {{R"("n1":"100","or":"1.8")", 0.6455, 0.6475, 0.6435},
                        {R"("n1":"100","effect":"probs","p1":"0.1 0.2 0.3 0.25 0.15")", 0.6332, 0.6368, 0.6305},
                        {R"("n1":"150","r":"2","p2":"0.05 0.1 0.13 0.17 0.2 0.15 0.2","or":"0.7")", 0.5270, 0.5288,
                         0.5255}};
  for (const Case& c : cases) {
    const json::Value r = first_row(power_request("sim.ordinal", std::string(c.inputs) + R"(,"reps":"20000")"));
    CHECK_CLOSE(num(r, "power_wald"), c.wald, 0, tolerance(c.wald, 4000, 20000));
    CHECK_CLOSE(num(r, "power_lr"), c.lr, 0, tolerance(c.lr, 4000, 20000));
    CHECK_CLOSE(num(r, "power_mw"), c.mw, 0, tolerance(c.mw, 4000, 20000));
  }
}

TEST_CASE("ordinal: type I error, Whitehead's approximation and the search") {
  const json::Value null = first_row(power_request("sim.ordinal", R"("n1":"40","or":"1","reps":"40000")"));
  for (const char* key : {"power_wald", "power_lr", "power_mw"})
    CHECK_CLOSE(num(null, key), 0.05, 0, 4 * std::sqrt(0.05 * 0.95 / 40000) + 0.003);
  // Whitehead: N = 3 (z_α/2 + z_β)² / ((log OR)² (1 − Σ p̄³)) for equal groups
  const json::Value a = first_row(
      R"({"procedure":"sim.ordinal","solve_for":"power","options":{"approximate":true},"inputs":{"n1":"120"}})");
  const double p2[5] = {0.2, 0.25, 0.25, 0.2, 0.1};
  double below = 0, cum = 0, cube = 0;
  for (int j = 0; j < 5; ++j) {
    cum += p2[j];
    const double up = j == 4 ? 1 : 1 / (1 + std::exp(-(std::log(cum / (1 - cum)) - std::log(1.8))));
    const double avg = (p2[j] + (up - below)) / 2;
    cube += avg * avg * avg;
    below = up;
  }
  const double ncp = std::log(1.8) / std::sqrt(3 / (240 * 0.25 * (1 - cube)));
  CHECK_CLOSE(num(a, "power"), stats::pnorm(ncp - 1.959963984540054) + stats::pnorm(-ncp - 1.959963984540054),
              1e-10, 0);
  const json::Value n = first_row(
      R"({"procedure":"sim.ordinal","solve_for":"n","options":{"power_curve":false},)"
      R"("inputs":{"test":"lr","power":"0.9","reps":"1000","seed":"6"}})");
  CHECK(num(n, "power") >= 0.9);
  const json::Value below_n = first_row(power_request(
      "sim.ordinal", R"("test":"lr","reps":"1000","seed":"6","n1":")" + text::format_g(num(n, "n1") - 1, 10) + "\""));
  CHECK(num(below_n, "power") < 0.9);
}

TEST_CASE("ordinal: the approximation gives the expected share of separated studies for the Wald test") {
  // Whitehead's formula gave 100% where nearly every study is separated and the Wald test cannot reject; the share
  // is P(max Y1 ≤ min Y2) + P(max Y2 ≤ min Y1) less twice P(every outcome in one category), here against a sum over
  // every pair of tables
  auto approx = [](const std::string& inputs) {
    return first_row(R"({"procedure":"sim.ordinal","solve_for":"power","options":{"approximate":true},"inputs":{)" +
                     inputs + "}}");
  };
  CHECK_CLOSE(num(approx(R"("p2":"0.5 0.5","n1":"5")"), "separated"), 0.16726510438359024, 1e-12, 0);
  CHECK_CLOSE(num(approx(R"("p2":"0.2 0.5 0.3","or":"4","n1":"4","r":"1.5")"), "separated"), 0.23775891079458425,
              1e-12, 0);
  const json::Value far = approx(R"("or":"1e10")");
  CHECK(num(far, "power") == 1);
  CHECK(num(far, "separated") > 0.9999);
  CHECK(far.get("_statement").as_string().find("In nearly all such studies one group's outcomes would all lie at or "
                                               "above the other's, where the Wald test cannot reject") !=
        std::string::npos);
  CHECK(approx(R"("or":"1e5")").get("_statement").as_string().find("In an expected 99.1% of such studies") !=
        std::string::npos);
  // not below 1%, nor for the tests that can reject a separated study
  CHECK(approx("").get("separated").is_null());
  CHECK(approx(R"("or":"1e10","test":"lr")").get("separated").is_null());
}

// ------------------------------------------------------------------ repeated measures with dropout (MMRM)

namespace {

/// The monotone data set also fitted in R: 15 + 15 patients, 4 visits, patient i observed at visits
/// 1 … 1 + (7i mod 4).
procs::mmrm::Stats mmrm_data() {
  procs::mmrm::Stats st;
  st.visits = 4;
  for (int i = 1; i <= 30; ++i) {
    const int g = i <= 15 ? 1 : 0, k = 1 + (i * 7) % 4;
    double y[4];
    for (int t = 1; t <= k; ++t)
      y[t - 1] = std::sin(2.3 * i) + 0.25 * g * t + 0.8 * std::sin(1.3 * i + 0.7 * t + 0.4 * i * t);
    st.add(g == 1 ? 0 : 1, k, y);
  }
  return st;
}

}  // namespace

TEST_CASE("mmrm: the maximum likelihood fit equals nlme::gls") {
  // R: gls(y ~ 0 + group:visit, correlation = corSymm(~ visit | id), weights = varIdent(~ 1 | visit), "ML",
  // tolerance 1e-12); gls's optimizer stops within about 1e-5 of the closed form.  gls scales the
  // covariance of the coefficients by N/(N − p) (77 observations, 8 means) even for ML.
  const procs::mmrm::Fit f = procs::mmrm::fit(mmrm_data(), false);
  CHECK(f.ok);
  CHECK_CLOSE(f.delta, 0.7524530571, 0, 5e-6);
  CHECK_CLOSE(f.se * std::sqrt(77.0 / 69), 0.4899937569, 2e-5, 0);
  CHECK_CLOSE(f.mu[0][3], 0.8790832595, 0, 5e-6);
  CHECK_CLOSE(f.mu[1][3], 0.1266302024, 0, 5e-6);
  const double sigma[16] = {0.8064073, 0.5156544, 0.4452442, 0.4021781, 0.5156544, 0.8608824, 0.3684058, 0.6330725,
                            0.4452442, 0.3684058, 0.7402863, 0.4879956, 0.4021781, 0.6330725, 0.4879956, 0.7951023};
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b) CHECK_CLOSE(f.sigma[a * procs::mmrm::TMAX + b], sigma[a * 4 + b], 0, 2e-5);
  CHECK_EQ(f.observed_last, 8.0);
}

TEST_CASE("mmrm: the REML fit equals nlme::gls with REML, and is exact with complete data") {
  // R (REML): delta 0.7524551490, se 0.5143913535, variances 0.8640051 0.9350190 0.8278465 0.9273541
  const procs::mmrm::Fit f = procs::mmrm::fit(mmrm_data(), true);
  CHECK_CLOSE(f.delta, 0.7524551490, 0, 5e-6);
  CHECK_CLOSE(f.se, 0.5143913535, 1e-5, 0);
  const double var[4] = {0.8640051, 0.9350190, 0.8278465, 0.9273541};
  for (int a = 0; a < 4; ++a) CHECK_CLOSE(f.sigma[a * procs::mmrm::TMAX + a], var[a], 1e-5, 0);
  // complete data: REML Σ̂ = E'E/(n − 2), and the difference at the last visit is the difference of means
  procs::mmrm::Stats st;
  st.visits = 3;
  double s1 = 0, s2 = 0, ss = 0;
  std::vector<std::array<double, 3>> rows;
  for (int i = 0; i < 20; ++i) {
    double y[3];
    for (int t = 0; t < 3; ++t) y[t] = std::sin(1.1 * i + 0.9 * t) + 0.5 * std::cos(0.7 * i * (t + 1)) + (i < 10 ? 0.3 : 0);
    st.add(i < 10 ? 0 : 1, 3, y);
    (i < 10 ? s1 : s2) += y[2];
    rows.push_back({y[0], y[1], y[2]});
  }
  const double m1 = s1 / 10, m2 = s2 / 10;
  for (int i = 0; i < 20; ++i) ss += std::pow(rows[i][2] - (i < 10 ? m1 : m2), 2);
  const procs::mmrm::Fit c = procs::mmrm::fit(st, true);
  CHECK_CLOSE(c.delta, m1 - m2, 1e-12, 1e-12);
  CHECK_CLOSE(c.se, std::sqrt(ss / 18 * (0.1 + 0.1)), 1e-10, 0);
}

TEST_CASE("mmrm: the fit is the same in combinations of the visits, and says why it fails") {
  // y = Le with L lower triangular: the fit of e asked for the combination in L's last row is the fit of y at the
  // last visit, since each visit's regression on the earlier ones changes with them.  With visits correlated
  // within 1e-13 of 1, y's sums of squares keep too few digits of the visits' differences and its fit fails.
  for (double rho : {0.6, 1 - 1e-13}) {
    const int T = 4;
    double l[T][T] = {};
    for (int j = 0; j < T; ++j) {  // Cholesky factor of the exchangeable correlation
      double v = 1;
      for (int k = 0; k < j; ++k) v -= l[j][k] * l[j][k];
      l[j][j] = std::sqrt(v);
      for (int i = j + 1; i < T; ++i) {
        double w = rho;
        for (int k = 0; k < j; ++k) w -= l[i][k] * l[j][k];
        l[i][j] = w / l[j][j];
      }
    }
    procs::mmrm::Stats se, sy;
    se.visits = sy.visits = T;
    for (int i = 1; i <= 40; ++i) {
      const int g = i <= 20 ? 0 : 1, k = 1 + (i * 7) % 4;
      double e[T], y[T] = {};
      for (int t = 0; t < T; ++t) e[t] = std::sin(2.3 * i + 0.9 * t * t) + 0.3 * std::cos(1.7 * i * (t + 1));
      for (int t = 0; t < T; ++t)
        for (int c = 0; c <= t; ++c) y[t] += l[t][c] * e[c];
      se.add(g, k, e);
      sy.add(g, k, y);
    }
    const procs::mmrm::Fit fe = procs::mmrm::fit(se, true, l[T - 1]), fy = procs::mmrm::fit(sy, true);
    CHECK(fe.ok);
    if (rho == 0.6) {
      CHECK(fy.ok);
      CHECK_CLOSE(fe.delta, fy.delta, 1e-12, 1e-14);
      CHECK_CLOSE(fe.se, fy.se, 1e-12, 0);
      CHECK_EQ(fe.observed_last, fy.observed_last);
    } else {
      CHECK(!fy.ok && fy.failure == procs::mmrm::Failure::Arithmetic);
    }
  }
  // a group with nobody left at the last visit, and 5 patients there where the regression on the group and 3
  // earlier visits needs 6
  procs::mmrm::Stats empty, few;
  empty.visits = few.visits = 4;
  for (int i = 0; i < 20; ++i) {
    double y[4];
    for (int t = 0; t < 4; ++t) y[t] = std::sin(1.9 * i + 1.3 * t);
    empty.add(i % 2, i % 2 == 0 ? 4 : 3, y);
    few.add(i % 2, i < 5 ? 4 : 3, y);
  }
  const procs::mmrm::Fit a = procs::mmrm::fit(empty, true), b = procs::mmrm::fit(few, true);
  CHECK(!a.ok && a.failure == procs::mmrm::Failure::Empty && a.visit == 4);
  CHECK(!b.ok && b.failure == procs::mmrm::Failure::Few && b.visit == 4);
  // one patient of each group at the first visit leaves its variance nothing to go on: no earlier visits are
  // involved, and the statement says so
  procs::mmrm::Stats first;
  first.visits = 2;
  for (int g = 0; g < 2; ++g) {
    const double y[2] = {0.3 + g, -0.4 * g};
    first.add(g, 2, y);
  }
  const procs::mmrm::Fit c = procs::mmrm::fit(first, true);
  CHECK(!c.ok && c.failure == procs::mmrm::Failure::Few && c.visit == 1);
  const json::Value one = first_row(power_request("sim.mmrm", R"("visits":"1","n1":"3","drop1":"0.5","drop2":"0.5",)"
                                                              R"("reps":"200")"));
  const std::string why = one.get("failures_why").as_string();
  CHECK(why.find("only one patient of each group was measured at the first visit, too few to estimate its "
                 "variance") != std::string::npos);
  CHECK(why.find("earlier visits") == std::string::npos);
}

TEST_CASE("mmrm: the statement shows a correlation or a dropout near an end of its range as it is") {
  // ρ = 0.9999999999993 read "exchangeable correlation 1", outside the open interval allowed, and a dropout of
  // 0.4% read 0%
  auto statement = [](const std::string& inputs) {
    return first_row(R"({"procedure":"sim.mmrm","solve_for":"power","options":{"approximate":true},"inputs":{)" +
                     inputs + "}}")
        .get("_statement")
        .as_string();
  };
  CHECK(statement(R"("rho":"0.9999999999993")").find("exchangeable correlation 0.9999999999993)") !=
        std::string::npos);
  CHECK(statement(R"("rho":"-0.9999999999999999","corr":"ar1")").find("AR(1) correlation -0.9999999999999999)") !=
        std::string::npos);
  CHECK(statement(R"("rho":"-0.2499999")").find("exchangeable correlation -0.2499999)") != std::string::npos);
  CHECK(statement(R"("rho":"0.6")").find("exchangeable correlation 0.6)") != std::string::npos);
  CHECK(statement(R"("drop1":"0.004","drop2":"0.9999")").find("with 0.4% of group 1 and 99.99% of group 2 lost") !=
        std::string::npos);
  CHECK(statement(R"("drop1":"0.2","drop2":"0.2")").find("with 20% of each group lost") != std::string::npos);
  // past 12 decimals, up to the largest share below 1
  CHECK(statement(R"("drop1":"0.999999999999999","analysis":"completers")")
            .find("with 99.9999999999999% of group 1 and 20% of group 2 lost") != std::string::npos);
  CHECK(statement(R"("drop1":"0.9999999999999999","analysis":"completers")")
            .find("with 99.99999999999999% of group 1 and 20% of group 2 lost") != std::string::npos);
}

TEST_CASE("mmrm: visits correlated within 1e-12 of 1 or -1 keep their power") {
  // The fits of the outcomes failed there (4% to 17% where the power is near 80%, or none analyzed at all); the
  // trials are now fitted in independent combinations of the visits, the same simulated errors whatever ρ.  The
  // approximation, which inverted the visits' information, drifted from 1 − 1e-14 on and gave α at 1 − 1.1e-16.
  const json::Value near = first_row(power_request("sim.mmrm", R"("rho":"0.99999999999","reps":"200")"));
  for (const char* inputs :
       {R"("rho":"0.9999999999993")", R"("rho":"0.9999999999999999")", R"("rho":"0.9999999999999999","corr":"ar1")",
        R"("rho":"-0.9999999999999999","corr":"ar1")"}) {
    const json::Value r = first_row(power_request("sim.mmrm", std::string(inputs) + R"(,"reps":"200")"));
    CHECK_EQ(num(r, "failures"), 0.0);
    CHECK_EQ(num(r, "power"), num(near, "power"));
    CHECK_CLOSE(num(r, "approx_power"), num(near, "approx_power"), 1e-9, 0);
  }
}

TEST_CASE("mmrm: without dropout the MMRM is the t-test; with dropout it agrees with the expected GLS variance") {
  const json::Value c = first_row(power_request("sim.mmrm", R"("drop1":"0","drop2":"0","n1":"60","reps":"4000")"));
  CHECK_EQ(num(c, "power_mmrm"), num(c, "power_completers"));
  const json::Value t = first_row(
      R"({"procedure":"means.two.t","solve_for":"power","inputs":{"n1":"60","mu1":"0.4","mu2":"0","sigma":"1"}})");
  CHECK_CLOSE(num(c, "approx_power"), num(t, "power"), 1e-9, 0);
  struct Case {
    const char* inputs;
  };
  for (const char* inputs : {R"("n1":"100","drop1":"0.2","drop2":"0.2")",
                             R"("n1":"80","visits":"12","drop1":"0.4","drop2":"0.3","rho":"0.5")",
                             R"("n1":"70","corr":"ar1","rho":"0.9","drop1":"0.4","drop2":"0.2")"}) {
    const json::Value r = first_row(power_request("sim.mmrm", std::string(inputs) + R"(,"reps":"20000")"));
    CHECK_CLOSE(num(r, "power"), num(r, "approx_power"), 0, 4 * num(r, "power_se") + 0.01);
    CHECK_CLOSE(num(r, "mean_delta"), 0.4, 0, 0.01);
    CHECK(num(r, "power_mmrm") > num(r, "power_completers"));  // the earlier visits of patients who leave help
  }
}

TEST_CASE("mmrm: the type I error is close to alpha") {
  for (const char* inputs : {R"("n1":"30","drop1":"0.3","drop2":"0.1")", R"("n1":"60","visits":"8","corr":"ar1")"}) {
    const json::Value r =
        first_row(power_request("sim.mmrm", std::string(inputs) + R"(,"delta":"0","reps":"40000")"));
    CHECK_CLOSE(num(r, "power_mmrm"), 0.05, 0, 4 * std::sqrt(0.05 * 0.95 / 40000) + 0.002);
    CHECK_CLOSE(num(r, "power_completers"), 0.05, 0, 4 * std::sqrt(0.05 * 0.95 / 40000));
  }
}

TEST_CASE("mmrm: reproducible across threads, and the sample-size search is minimal") {
  const std::string inputs = R"("n1":"40","drop1":"0.25","visits":"6","reps":"700","seed":"12")";
  const json::Value one = first_row(power_request("sim.mmrm", inputs, 1));
  const json::Value four = first_row(power_request("sim.mmrm", inputs, 4));
  CHECK_EQ(num(one, "power"), num(four, "power"));
  CHECK_EQ(num(one, "mean_delta"), num(four, "mean_delta"));
  const json::Value n = first_row(
      R"({"procedure":"sim.mmrm","solve_for":"n","options":{"power_curve":false},)"
      R"("inputs":{"power":"0.85","delta":"0.5","reps":"1000","seed":"13"}})");
  CHECK(num(n, "power") >= 0.85);
  const json::Value below = first_row(power_request(
      "sim.mmrm", R"("delta":"0.5","reps":"1000","seed":"13","n1":")" + text::format_g(num(n, "n1") - 1, 10) + "\""));
  CHECK(num(below, "power") < 0.85);
}

// ------------------------------------------------------------------ win ratio

TEST_CASE("winratio: win, loss and tie probabilities and power agree with an independent R implementation") {
  // R: 400,000 pairs (wins, losses, ties) and 1,500 trials analyzed with its own pairwise comparisons, the
  // Finkelstein–Schoenfeld test and the win ratio's U-statistic Wald test
  const json::Value base = first_row(power_request("sim.winratio", R"("reps":"5000")"));
  CHECK_CLOSE(num(base, "pwin"), 0.4196, 0, 0.003);
  CHECK_CLOSE(num(base, "ploss"), 0.3221, 0, 0.003);
  CHECK_CLOSE(num(base, "ptie"), 0.2583, 0, 0.003);
  CHECK_CLOSE(num(base, "power_fs"), 0.5947, 0, tolerance(0.6, 1500, 5000));
  CHECK_CLOSE(num(base, "power_wr"), 0.5920, 0, tolerance(0.6, 1500, 5000));
  const json::Value score = first_row(power_request("sim.winratio", R"("tier3":"score","delta_s":"0.2","reps":"5000")"));
  CHECK_CLOSE(num(score, "pwin"), 0.5491, 0, 0.003);
  CHECK_CLOSE(num(score, "ptie"), 0.0227, 0, 0.002);
  CHECK_CLOSE(num(score, "power_fs"), 0.7600, 0, tolerance(0.76, 1500, 5000));
  CHECK_CLOSE(num(score, "power_wr"), 0.7587, 0, tolerance(0.76, 1500, 5000));
  const json::Value other = first_row(power_request(
      "sim.winratio", R"("tier2":"none","tier3":"score","delta_s":"0.3","margin_s":"0.2","accrual":"12",)"
                      R"("follow":"24","med_d":"36","hr_d":"0.7","n1":"150","r":"2","reps":"5000")"));
  CHECK_CLOSE(num(other, "power_fs"), 0.8560, 0, tolerance(0.85, 1500, 5000));
  CHECK_CLOSE(num(other, "power_wr"), 0.8480, 0, tolerance(0.85, 1500, 5000));
}

TEST_CASE("winratio: the type I error is close to alpha, and results do not depend on threads") {
  const json::Value null = first_row(
      power_request("sim.winratio", R"("hr_d":"1","hr_h":"1","tier3":"score","delta_s":"0","n1":"80","reps":"20000")"));
  CHECK_CLOSE(num(null, "power_fs"), 0.05, 0, 4 * std::sqrt(0.05 * 0.95 / 20000));
  CHECK_CLOSE(num(null, "power_wr"), 0.05, 0, 4 * std::sqrt(0.05 * 0.95 / 20000) + 0.002);
  const std::string inputs = R"("n1":"60","hr_d":"0.6","reps":"300","seed":"21")";
  const json::Value one = first_row(power_request("sim.winratio", inputs, 1));
  const json::Value four = first_row(power_request("sim.winratio", inputs, 4));
  CHECK_EQ(num(one, "power"), num(four, "power"));
  CHECK_EQ(num(one, "mean_nb"), num(four, "mean_nb"));
  const json::Value n = first_row(
      R"({"procedure":"sim.winratio","solve_for":"n","options":{"power_curve":false},)"
      R"("inputs":{"power":"0.8","hr_d":"0.6","hr_h":"0.6","reps":"400","seed":"22"}})");
  CHECK(num(n, "power") >= 0.8);
  const json::Value below = first_row(power_request(
      "sim.winratio", R"("hr_d":"0.6","hr_h":"0.6","reps":"400","seed":"22","n1":")" +
                          text::format_g(num(n, "n1") - 1, 10) + "\""));
  CHECK(num(below, "power") < 0.8);
}

// ------------------------------------------------------------------ stepped wedge with a mixed model

TEST_CASE("sw: the quick approximation is Hussey and Hughes's power when the cluster effect is constant") {
  const json::Value a = first_row(
      R"({"procedure":"sim.sw","solve_for":"power","options":{"approximate":true},"inputs":{"df":"z"}})");
  const json::Value e = first_row(
      R"({"procedure":"cluster.sw","solve_for":"power","inputs":{"steps":"4","g":"2","m":"20","icc":"0.05","delta":"0.3"}})");
  CHECK_CLOSE(num(a, "power"), num(e, "power"), 1e-9, 0);
  // with cluster-period effects (CAC < 1) the variance grows, so the power falls
  const json::Value c = first_row(
      R"({"procedure":"sim.sw","solve_for":"power","options":{"approximate":true},"inputs":{"df":"z","cac":"0.7"}})");
  CHECK(num(c, "power") < num(a, "power"));
}

TEST_CASE("sw: simulated fits agree with lmerTest") {
  // R: lmer(y ~ period + trt + (1 | cluster)) with Satterthwaite's t-test on 8 clusters (4 steps of 2), 20
  // subjects per cluster and period: power 0.714 (1,500 trials), type I error 0.0553 (3,000 trials), and
  // mean df 382.5; with a cluster autocorrelation of 0.7 and the model lmer(… + (1 | cluster) + (1 | cluster:period)):
  // power 0.5945 (2,000), type I error 0.0626 (5,000; the cluster-period variance at 0 in 19.8% of the trials,
  // mean df 91.1); the random-intercept model then rejects 9.1% (2,000)
  const json::Value hh = first_row(power_request("sim.sw", R"("reps":"20000")"));
  CHECK_CLOSE(num(hh, "power"), 0.714, 0, tolerance(0.714, 1500, 20000));
  CHECK_CLOSE(num(hh, "mean_delta"), 0.3, 0, 0.004);
  CHECK_CLOSE(num(hh, "mean_df"), 382.5, 0.05, 0);
  const json::Value hh0 = first_row(power_request("sim.sw", R"("delta":"0","reps":"20000")"));
  CHECK_CLOSE(num(hh0, "power"), 0.0553, 0, tolerance(0.055, 3000, 20000));
  const json::Value ho = first_row(power_request("sim.sw", R"("cac":"0.7","analysis":"hooper","reps":"20000")"));
  CHECK_CLOSE(num(ho, "power"), 0.5945, 0, tolerance(0.6, 2000, 20000));
  const json::Value ho0 = first_row(power_request("sim.sw", R"("cac":"0.7","analysis":"hooper","delta":"0","reps":"20000")"));
  CHECK_CLOSE(num(ho0, "power"), 0.0626, 0, tolerance(0.06, 5000, 20000));
  CHECK_CLOSE(num(ho0, "mean_df"), 91.1, 0.08, 0);
  const json::Value mis = first_row(power_request("sim.sw", R"("cac":"0.7","delta":"0","reps":"20000")"));
  CHECK_CLOSE(num(mis, "power"), 0.0910, 0, tolerance(0.09, 2000, 20000));
}

TEST_CASE("sw: reproducible across threads, and the search for clusters per step is minimal") {
  const std::string inputs = R"("steps":"5","g":"2","m":"10","cac":"0.8","analysis":"hooper","reps":"600","seed":"31")";
  const json::Value one = first_row(power_request("sim.sw", inputs, 1));
  const json::Value four = first_row(power_request("sim.sw", inputs, 4));
  CHECK_EQ(num(one, "power"), num(four, "power"));
  CHECK_EQ(num(one, "mean_df"), num(four, "mean_df"));
  const json::Value g = first_row(
      R"({"procedure":"sim.sw","solve_for":"g","options":{"power_curve":false},)"
      R"("inputs":{"power":"0.9","steps":"3","m":"15","delta":"0.25","reps":"1000","seed":"32"}})");
  CHECK(num(g, "power") >= 0.9);
  if (num(g, "g") > 1) {
    const json::Value below = first_row(power_request(
        "sim.sw", R"("steps":"3","m":"15","delta":"0.25","reps":"1000","seed":"32","g":")" +
                      text::format_g(num(g, "g") - 1, 10) + "\""));
    CHECK(num(below, "power") < 0.9);
  }
}

TEST_CASE("winratio: the counts by sorting equal the counts over every pair, ties included") {
  // patients with times rounded to whole units (many exact ties of death, event and follow-up times),
  // equal follow-up for some, values rounded to tenths and margins that hit them exactly
  procs::winratio::Counts fast, direct;
  for (int trial = 0; trial < 60; ++trial) {
    sim::Rng rng(77, static_cast<std::uint64_t>(trial));
    const bool event = trial % 3 != 0, score = trial % 2 == 0, round_times = trial % 4 < 2;
    const double margin = trial % 5 == 0 ? 0.0 : 0.1 * (trial % 5);
    const long n1 = 5 + static_cast<long>(rng.uniform() * 40), n2 = 3 + static_cast<long>(rng.uniform() * 40);
    std::vector<procs::winratio::Patient> p;
    for (long i = 0; i < n1 + n2; ++i) {
      auto t = [&](double x) { return round_times ? std::round(x) : x; };
      const double end = trial % 7 == 0 ? 10.0 : t(5 + 10 * rng.uniform());
      const double death = t(rng.exponential() * 20), ev = t(rng.exponential() * 8);
      procs::winratio::Patient q;
      q.death = death <= end ? death : INFINITY;
      q.alive = std::fmin(death, end);
      q.event = event && ev <= end && ev < death ? ev : INFINITY;
      q.free = event ? std::fmin(ev, q.alive) : INFINITY;
      q.value = score && death > end ? std::round(10 * rng.normal()) / 10 : NAN;
      p.push_back(q);
    }
    procs::winratio::pair_counts(p, n1, margin, fast);
    procs::winratio::pair_counts_direct(p, n1, margin, direct);
    for (int g = 0; g < 2; ++g)
      for (std::size_t k = 0; k < p.size(); ++k) {
        CHECK_EQ(fast.win[g][k], direct.win[g][k]);
        CHECK_EQ(fast.loss[g][k], direct.loss[g][k]);
      }
  }
}
