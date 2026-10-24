// What a user can type at the extremes (a ratio of 1e-30 or 1e-300, a noncentrality of 1e300, a count beyond a
// billion) ends at once with an answer or a message, Stop ends what takes long, and a long simulation keeps only
// a block of studies in memory with the same result as keeping them all.
#include <chrono>
#include <cmath>
#include <string>
#include <thread>
#include <vector>

#include "engine/runner.hpp"
#include "engine/solve.hpp"
#include "procedures/common.hpp"
#include "procedures/sim_common.hpp"
#include "sim/simulate.hpp"
#include "stats/dispatch.hpp"
#include "stats/numerics.hpp"
#include "stats/stats.hpp"
#include "testing.hpp"
#include "util/json.hpp"

using namespace ssc;

namespace {

double seconds_since(std::chrono::steady_clock::time_point t) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

json::Value run(const std::string& procedure, const std::string& solve, const std::string& inputs) {
  Token token;
  return run_request(json::parse(R"({"procedure":")" + procedure + R"(","solve_for":")" + solve +
                                 R"(","options":{"power_curve":true,"approximate":true},"inputs":)" + inputs + "}"),
                     &token);
}

// The message of a result: its error, an input's error or the first row's error.
std::string message(const json::Value& res) {
  if (!res.get("ok").as_bool(false)) {
    std::string m = res.get("error").as_string();
    for (const auto& kv : res.get("input_errors").items()) m += " " + kv.second.as_string();
    return m;
  }
  return res.get("rows").elements().at(0).get("_error").as_string();
}

}  // namespace

TEST_CASE("extremes: noncentral distributions at enormous and infinite noncentrality") {
  const auto t0 = std::chrono::steady_clock::now();
  // δ²/2 overflows beyond 1.3e154: the mass sits far above any finite t
  CHECK_EQ(stats::pnt(2.0, 10, 1e300, true), 0.0);
  CHECK_EQ(stats::pnt(2.0, 10, 1e300, false), 1.0);
  CHECK(stats::pnt(1e9, 10, 1e300, true) < 1e-300);
  CHECK_EQ(stats::pnchisq(10, 3, INFINITY, true), 0.0);
  CHECK_EQ(stats::pnchisq(10, 3, INFINITY, false), 1.0);
  CHECK_EQ(stats::pnf(2, 3, 20, INFINITY, true), 0.0);
  CHECK_EQ(stats::pnf(2, 3, 20, INFINITY, false), 1.0);
  CHECK(seconds_since(t0) < 1);
}

TEST_CASE("extremes: the smallest group 1 that gives group 2 enough subjects") {
  for (double r : {1.0, 0.7, 0.5, 0.3, 0.013, 2.5}) {
    auto n2 = [&](double n1) { return solve::round_up(r * n1); };
    for (double need : {1.0, 2.0, 3.0}) {
      double counted = need;
      while (n2(counted) < need) counted += 1;
      CHECK_EQ(solve::first_whole([&](double x) { return n2(x) >= need; }, need, 1e8, "none"), counted);
    }
  }
  bool thrown = false;
  try {
    solve::first_whole([](double x) { return std::ceil(1e-300 * x - 1e-9) >= 2; }, 2, 1e8, "too small");
  } catch (const CalcError& e) {
    thrown = std::string(e.what()) == "too small";
  }
  CHECK(thrown);
}

TEST_CASE("extremes: a group from a tiny ratio has 1 subject, n × the ratio rounded up") {
  // the slack for the rounding of the product, 1e-9 absolute, left a group of 0 where n × the ratio was 1e-9 or less
  CHECK_EQ(solve::round_up(1e-300), 1.0);
  CHECK_EQ(solve::round_up(1e-9), 1.0);
  CHECK_EQ(solve::round_up(1.0000000000000002), 1.0);
  CHECK_EQ(solve::round_up(20 * 1.1), 22.0);
  CHECK_EQ(solve::round_up(0), 0.0);
  CHECK(message(run("means.two.t", "power", R"({"alloc":"ratio","n1":"10","r":"5e-11"})")) ==
        "group 2 would have N2 = 1 (R × N1, rounded up); each group needs at least 2");
  // where a group may have 1 subject, the answer is that of N2 = 1 typed in (it was refused for an empty group 2)
  for (const char* pid : {"surv.logrank", "rates.two.poisson", "gs.logrank"}) {
    const json::Value a = run(pid, "power", R"({"alloc":"ratio","n1":"10","r":"1e-12"})");
    const json::Value b = run(pid, "power", R"({"alloc":"sep","n1":"10","n2":"1"})");
    const json::Value& row = a.get("rows").elements().at(0);
    CHECK_EQ(row.get("n2").as_number(0), 1.0);
    CHECK_EQ(row.get("power").as_number(0), b.get("rows").elements().at(0).get("power").as_number(1));
  }
}

TEST_CASE("extremes: a level reads as itself, never as a level next to it") {
  // a confidence of 0.5000000000000001 read as 50%, the level a one-sided interval is refused at
  CHECK(procs::level_text(0.5000000000000001) == "50.00000000000001%");
  CHECK(procs::level_text(0.4999999999999999) == "49.99999999999999%");
  CHECK(procs::level_text(0.95000000001) == "95.000000001%");
  CHECK(procs::level_text(0.9999999999999999) == "99.99999999999999%");
  CHECK(procs::level_text(0.95) == "95%");
  CHECK(procs::level_text(0.975) == "97.5%");
  CHECK(procs::level_text(0.00012) == "0.012%");
  CHECK(procs::level_text(1e-7) == "1e-05%");
  CHECK(procs::level_text(3.25e-30) == "3.25e-28%");
  const std::string stmt = run("means.one.ci", "n", R"({"sides":"one","conf":"0.5000000000000001","e":"0.1"})")
                               .get("rows").elements().at(0).get("_statement").as_string();
  CHECK(stmt.find(" 50.00000000000001% ") != std::string::npos);
}

TEST_CASE("extremes: a simulation in blocks summarizes as all its studies at once") {
  auto study = [](long i) {
    sim::Outcome o;
    o.reject = i % 3 == 0;
    o.failed = i % 97 == 0;
    o.estimate = std::sin(0.001 * static_cast<double>(i));
    o.se = 1 + 0.5 * std::cos(0.002 * static_cast<double>(i));
    o.df = 10 + static_cast<double>(i % 7);
    o.extra = static_cast<double>(i % 11);
    o.also = static_cast<unsigned>(i % 4);
    return o;
  };
  const long reps = 2 * sim::BLOCK + 123;
  std::vector<sim::Outcome> all;
  for (long i = 0; i < reps; ++i) all.push_back(study(i));
  const sim::Summary whole = sim::summarize(all);
  Token token;
  for (int threads : {1, 3, 8}) {
    const RunContext ctx{&token, threads, false, nullptr};
    const sim::Summary blocks = sim::simulate(reps, study, ctx);
    CHECK_EQ(blocks.reps, whole.reps);
    CHECK_EQ(blocks.rejections, whole.rejections);
    CHECK_EQ(blocks.failures, whole.failures);
    CHECK_EQ(blocks.power, whole.power);
    CHECK_EQ(blocks.mean_estimate, whole.mean_estimate);
    CHECK_EQ(blocks.sd_estimate, whole.sd_estimate);
    CHECK_EQ(blocks.mean_se, whole.mean_se);
    CHECK_EQ(blocks.mean_df, whole.mean_df);
    CHECK_EQ(blocks.mean_extra, whole.mean_extra);
    for (int b = 0; b < sim::ALSO_MAX; ++b) CHECK_EQ(blocks.also[b], whole.also[b]);
  }
}

TEST_CASE("extremes: typed inputs that used to run forever end at once with a message") {
  struct Case {
    const char* procedure;
    const char* solve;
    const char* inputs;
    const char* says;  // part of the message, or "" for an answer
  };
  // at the smallest and largest sizes an input may have: the values beyond them that once ran forever (1e-300,
  // 1e300, subnormal ones) are refused, as the next test checks
  const Case cases[] = {
      {"means.two.cluster", "k", R"({"kr":"1e-30"})", "cluster ratio is too small"},
      // group 2 has 1 cluster at any k1, which limits the power
      {"prop.two.cluster", "k", R"({"kr":"1e-30"})", "group 2 has 1 cluster at every k1"},
      {"diag.auc", "n", R"({"ratio":"1e-30"})", "ratio of controls to cases is too small"},
      {"sim.lmm.cluster", "k", R"({"kr":"1e-30"})", "cluster ratio is too small"},
      {"sim.glmm.cluster", "k", R"({"kr":"1e-30"})", "cluster ratio is too small"},
      {"sim.lmm.longitudinal", "n", R"({"r":"1e-30"})", "allocation ratio is too small"},
      {"sim.mmrm", "n", R"({"r":"1e-30"})", "allocation ratio is too small"},
      {"sim.means.robust", "n", R"({"r":"1e-30"})", "allocation ratio is too small"},
      {"sim.ordinal", "n", R"({"r":"1e-30"})", "allocation ratio is too small"},
      {"sim.winratio", "n", R"({"r":"1e-30"})", "allocation ratio is too small"},
      {"prop.one", "power", R"({"n":"1e10"})", "at most 1,000,000,000"},
      {"corr.one", "power", R"({"n":"1e300"})", "at most 1,000,000,000"},
      {"prop.one.simon", "power", R"({"r1":"1","n1":"999999999","r":"5","n":"1000000000"})", "up to 100,000"},
      {"prop.paired.mcnemar", "power", R"({"n":"10000000"})", "up to 1,000,000 pairs"},
      {"means.two.ancova", "mu1", R"({"ncov":"1000"})", "number of covariates"},
      {"means.anova2", "power", R"({"a":"20","b":"214748365"})", "at most 100"},  // A × B wrapped around to 4
      {"means.two.welch", "power", R"({"sd1":"1e30"})", ""},
      {"means.two.t", "mu1", R"({"sigma":"1e30"})", ""},
      {"sim.multiarm", "n", R"({"alpha":"1e-30"})", ""},  // 1 − α rounds to 1: Dunnett's value from the tail
      {"rates.one.poisson", "power", R"({"lam0":"1e30"})", ""},
      {"rates.one.poisson", "n", R"({"time":"1e30"})", ""},
      {"means.crossover", "n", R"({"delta":"1e30"})", ""},
      {"means.anova", "power", R"({"sigma":"1e-30"})", ""},
      {"means.rm.anova", "n", R"({"sigma":"1e-30"})", ""},
      // loops that a value rounding to 0 or infinity stopped from advancing (and that Stop could not end), at the
      // nearest sizes allowed
      {"sim.surv.nph", "power", R"({"med2":"1e-30"})", ""},
      {"sim.surv.nph", "n", R"({"med2":"1e-30"})", "cannot be reached"},  // every event before T*, where HR1 = 1
      {"rates.one.poisson", "power", R"({"alpha":"1e-30","test":"exact"})", ""},
      {"sim.means.robust", "power", R"({"dist":"lognormal","skew":"1e-9"})", ""},
      {"sim.means.robust", "n", R"({"dist":"gamma","skew":"-1e-8"})", ""},
      // 1e6 / prevalence is 1e36 subjects, past any int
      {"diag.paired", "n", R"({"method":"exact","goal":"se","prev":"1e-30"})", "cannot be reached"},
  };
  {  // a full simulation whose every study would be enormous (Stop acts between studies)
    Token token;
    const auto t0 = std::chrono::steady_clock::now();
    const json::Value res = run_request(
        json::parse(R"({"procedure":"sim.lmm.cluster","solve_for":"power","options":{"power_curve":false},)"
                    R"("inputs":{"kr":"1e10","reps":"100"}})"),
        &token);
    CHECK(message(res).find("limited to 1,000,000 clusters") != std::string::npos);
    CHECK(seconds_since(t0) < 5);
  }
  for (const Case& c : cases) {
    const auto t0 = std::chrono::steady_clock::now();
    const json::Value res = run(c.procedure, c.solve, c.inputs);
    const std::string m = message(res);
    if (seconds_since(t0) > 5 || (*c.says ? m.find(c.says) == std::string::npos : !m.empty()))
      ::ssc::test::report_failure(std::string(c.procedure) + " " + c.inputs + ": '" + m + "'", __FILE__, __LINE__);
    ::ssc::test::count_check();
  }
}

TEST_CASE("extremes: values beyond 1e-30 and 1e30 in size are refused, in the terms of each input's limits") {
  struct Case {
    const char* procedure;
    const char* solve;
    const char* inputs;
    const char* id;    // the input refused
    const char* says;  // its message
  };
  const Case cases[] = {
      // values that once ran forever, or were computed with as 0 or infinity
      {"means.two.welch", "power", R"({"sd1":"1e300"})", "sd1", "values must be between 1e-30 and 1e30"},
      {"means.two.t", "mu1", R"({"sigma":"1e300"})", "sigma", "values must be between 1e-30 and 1e30"},
      {"sim.surv.nph", "power", R"({"med2":"5e-324"})", "med2", "values must be between 1e-30 and 1e30"},
      {"sim.surv.nph", "n", R"({"med2":"1e-320"})", "med2", "values must be between 1e-30 and 1e30"},
      {"rates.one.poisson", "power", R"({"alpha":"5e-324","test":"exact"})", "alpha", "values must be at least 1e-30"},
      {"reg.cox", "n", R"({"alpha":"5e-324"})", "alpha", "values must be at least 1e-30"},
      {"means.be.ratio", "power", R"({"cv":"5e-324"})", "cv", "values must be between 1e-30 and 1e30"},
      {"cluster.sw", "power", R"({"sigma":"1e-300"})", "sigma", "values must be between 1e-30 and 1e30"},
      {"means.two.slopes", "delta", R"({"duration":"1e-160"})", "duration", "values must be between 1e-30 and 1e30"},
      {"prop.two.ci", "n", R"({"conf":"1e-300"})", "conf", "values must be at least 1e-30"},
      {"prop.cc", "n", R"({"p0":"5e-324"})", "p0", "values must be at least 1e-30"},
      {"rates.two.negbin", "power", R"({"r0":"5e-324"})", "r0", "values must be between 1e-30 and 1e30"},
      {"gs.logrank", "power", R"({"alpha":"5e-324","boundary":"wt_obf"})", "alpha", "values must be at least 1e-30"},
      {"sim.winratio", "power", R"({"med_d":"5e-324"})", "med_d", "values must be between 1e-30 and 1e30"},
      {"means.two.t", "power", R"({"alloc":"ratio","r":"1e308"})", "r", "values must be between 1e-30 and 1e30"},
      // 0 allowed or not, either sign, one end or both; every value of a list, a range and a list of lists
      {"means.two.t", "power", R"({"mu1":"-1e31"})", "mu1", "values must be 0 or between 1e-30 and 1e30 in size"},
      {"means.two.t", "power", R"({"mu1":"0.5 1e-31"})", "mu1", "values must be 0 or between 1e-30 and 1e30 in size"},
      {"means.two.t", "n", R"({"dropout":"0 1e-31"})", "dropout", "values must be 0 or at least 1e-30"},
      {"means.two.t", "power", R"({"sigma":"1e-31 to 1 by 0.5"})", "sigma", "values must be between 1e-30 and 1e30"},
      // the ends of a range as typed, as each would be alone: one the values do not reach, and a start kept as typed
      // (rounded to 12 digits, it would be 1e-30)
      {"means.two.t", "power", R"({"mu1":"0 to 1.5e30 by 1e30"})", "mu1",
       "values must be 0 or between 1e-30 and 1e30 in size"},
      {"means.two.t", "power", R"({"mu1":"9.9999999999999e-31 to 9.9999999999999e-31"})", "mu1",
       "values must be 0 or between 1e-30 and 1e30 in size"},
      {"corr.one", "n", R"({"rho1":"-1e-40"})", "rho1", "values must be 0 or at least 1e-30 in size"},
      {"prop.two.equiv", "n", R"({"dl":"-1e-31"})", "dl", "values must be at least 1e-30 in size"},
      {"prop.two.equiv", "n", R"({"scale":"or","oru":"1e31"})", "oru", "values must be at most 1e30"},
      {"means.anova", "n", R"({"means":"0 0.5 1; 0 1e31 2"})", "means",
       "values must be 0 or between 1e-30 and 1e30 in size"},
      {"means.two.slopes", "delta", R"({"timing":"custom","times":"0 1e-160 2e-160"})", "times",
       "values must be 0 or between 1e-30 and 1e30 in size"},
  };
  for (const Case& c : cases) {
    const json::Value res = run(c.procedure, c.solve, c.inputs);
    const json::Value& errors = res.get("input_errors");
    if (res.get("ok").as_bool(true) || errors.size() != 1 || errors.get(c.id).as_string() != c.says)
      ::ssc::test::report_failure(std::string(c.procedure) + " " + c.inputs + ": '" + message(res) + "'", __FILE__,
                                  __LINE__);
    ::ssc::test::count_check();
  }
  // the sizes themselves, 0 (of either sign) and counts (which have their own limit) are allowed
  const json::Value res = run("means.two.t", "power", R"({"sigma":"1e-30 1e30","mu1":"0 -0 -1e30 1e-30"})");
  CHECK(res.get("ok").as_bool(false));
  CHECK_EQ(res.get("rows").size(), std::size_t{8});
  CHECK(message(run("corr.one", "power", R"({"n":"1000000000"})")).empty());
  // so is a range that crosses 0 at a small scale, whose value there is 0 rather than the rounding noise of the sum
  // (3.9e-31 and 3.5e-46 here)
  for (const char* inputs :
       {R"({"mu1":"-3e-15 to 3e-15 by 1e-15","sigma":"2e-15"})", R"({"mu1":"-3e-30 to 3e-30 by 1e-30"})"}) {
    const json::Value range = run("means.two.t", "power", inputs);
    CHECK(range.get("ok").as_bool(false));
    CHECK_EQ(range.get("rows").at(3).get("mu1").as_number(NAN), 0.0);
  }
}

TEST_CASE("extremes: a share of the smallest alpha stays positive") {
  CHECK_EQ(procs::half_alpha(0.05), 0.025);
  CHECK_EQ(procs::half_alpha(1e-323), 5e-324);
  CHECK_EQ(procs::half_alpha(5e-324), 5e-324);  // half of it rounds to 0
  CHECK_EQ(procs::alpha_share(5e-324, 3), 5e-324);
  CHECK_EQ(procs::tail_alpha(0.05, blocks::Alt::Two), 0.025);
  CHECK_EQ(procs::tail_alpha(0.05, blocks::Alt::Upper), 0.05);
}

TEST_CASE("extremes: the smallest inputs allowed give an answer, or a message that says why not") {
  struct Case {
    const char* procedure;
    const char* solve;
    const char* inputs;
    const char* says;  // part of the message, or "" for an answer
  };
  // at subnormal values (now refused) each used to say that the target "cannot be reached", or gave an answer
  // computed from a 0 or an infinity; at 1e-30, the smallest size an input may have, each is computed
  const Case cases[] = {
      // the smallest alpha, two-sided
      {"reg.cox", "n", R"({"alpha":"1e-30"})", ""},
      {"rates.two.poisson", "n", R"({"alpha":"1e-30"})", ""},
      {"sim.lmm.cluster", "k", R"({"alpha":"1e-30"})", ""},
      {"sim.multiarm", "n", R"({"alpha":"1e-30","adjust":"bonferroni"})", ""},
      // squares of tiny inputs: CV², σ², D²
      {"means.be.ratio", "n", R"({"cv":"1e-30"})", ""},
      {"means.be.ratio", "power", R"({"cv":"1e-30"})", ""},
      {"cluster.sw", "power", R"({"sigma":"1e-30"})", ""},
      {"means.two.t", "mu1", R"({"sigma":"1e-30"})", ""},
      {"means.two.slopes", "delta", R"({"duration":"1e-30"})", ""},
      // SDs squared together: the Welch degrees of freedom, the paired SD, the slopes' SD
      {"means.two.welch", "mu1", R"({"sd1":"1e-30","sd2":"1e-30"})", ""},
      {"means.paired.t", "mud", R"({"sdtype":"corr","sd1":"1e-30","sd2":"1e-30"})", ""},
      {"means.two.slopes", "delta", R"({"sdb":"1e-30","sde":"1e-30"})", ""},
      // intervals at a confidence near 0 close on the estimate
      {"prop.two.ci", "n", R"({"conf":"1e-30"})", ""},
      {"prop.one.ci", "n", R"({"conf":"1e-15","method":"wilson"})", ""},
      // proportions, rates and alphas this small: a study of about 1e30 subjects, or an answer
      {"prop.cc", "n", R"({"p0":"1e-30"})", "cannot be reached"},
      {"prop.two", "n", R"({"es":"or","p2":"1e-30"})", "cannot be reached"},
      {"rates.two.poisson", "n", R"({"time":"1e-30"})", "cannot be reached"},
      // (a normal approximation that gives a side of the test more power than the chance of the events it needs,
      // 7e-29 in group 1, is refused)
      {"rates.two.negbin", "power", R"({"r0":"1e-30"})", "it fails with so few events"},
      // (and so is one where a dispersion of 1e6 gathers 100 expected events in so few subjects that the chance of
      // any is 0.00016)
      {"rates.two.negbin", "power", R"({"kappa":"1e6","r0":"10","rr":"1","n1":"10","n2":"10"})",
       "it fails with so little chance of an event"},
      {"gs.means", "n", R"({"alpha":"1e-30"})", ""},
      {"gs.logrank", "power", R"({"alpha":"1e-30","boundary":"wt_obf"})", ""},
      {"sim.winratio", "power", R"({"med_d":"1e-30"})", ""},
  };
  for (const Case& c : cases) {
    const auto t0 = std::chrono::steady_clock::now();
    const json::Value res = run(c.procedure, c.solve, c.inputs);
    const std::string m = message(res);
    if (seconds_since(t0) > 5 || (*c.says ? m.find(c.says) == std::string::npos : !m.empty()))
      ::ssc::test::report_failure(std::string(c.procedure) + " " + c.inputs + ": '" + m + "'", __FILE__, __LINE__);
    ::ssc::test::count_check();
  }
}

TEST_CASE("extremes: an effect on a tiny scale is solved to the tolerance relative to it") {
  // a tolerance of 1e-12 absolute took the first value tried (the scale itself) for the root 2.5 × scale
  for (double scale : {1e-8, 1e-15, 1e-160}) {
    auto power = [&](double x) { return stats::pnorm(x / scale - 2.5); };
    CHECK_CLOSE(solve::solve_increasing(power, 0.5, 0, scale, 1e9 * scale, 1e-12), 2.5 * scale, 1e-9, 0);
  }
  // a tolerance from the bracket (here 1e4 times the root) left the fourth digit wrong; relative to the root, the
  // answer is the same at every scale
  for (double scale : {1e-10, 1.0, 1e10}) {
    auto power = [&](double x) { return stats::pnorm(x / scale - 2.5); };
    CHECK_CLOSE(solve::solve_increasing(power, 0.5, 0, 1e4 * scale, 1e9 * scale, 1e-12) / scale, 2.5, 1e-11, 0);
  }
  // from a null value far from 0 the tolerance is relative to the distance from it
  auto power = [](double x) { return stats::pnorm((x - 1e6) / 1e-3 - 2.5); };
  CHECK_CLOSE(solve::solve_increasing(power, 0.5, 1e6, 1e6 + 1e-3, 1e7, 1e-12) - 1e6, 2.5e-3, 1e-6, 0);
  CHECK_CLOSE(solve::solve_decreasing([](double x) { return stats::pnorm((1e6 - x) / 1e-3 - 2.5); }, 0.5, 1e6,
                                      1e6 - 1e-3, 0, 1e-12) - 1e6, -2.5e-3, 1e-6, 0);
}

TEST_CASE("extremes: an effect that a large value next to it rounds away is refused") {
  CHECK_EQ(solve::offset(1e6, 0.5, "μ1", "μ2"), 1e6 + 0.5);
  CHECK_EQ(solve::offset(0, 1e-300, "μ1", "μ2"), 1e-300);
  // next to 1e10 the difference keeps six digits (0.89999962), more than the four shown
  CHECK_EQ(solve::offset(1e10, 0.9, "μ1", "μ2", true), 1e10 + 0.9);
  auto refusal = [](auto&& f) {
    try {
      f();
    } catch (const CalcError& e) {
      return std::string(e.what());
    }
    return std::string();
  };
  // nothing of it is kept next to 1e16 (0 or 2), and three digits next to 1e13 (0.90039)
  for (double x0 : {1e16, -1e17, 1e13}) {
    const std::string m = refusal([&] { solve::offset(x0, 0.9, "μ1", "μ2", true); });
    CHECK(m.find("μ1 cannot be told apart from μ2 = ") == 0);
    CHECK((m.find("is lost to rounding") != std::string::npos) == (x0 != 1e13));
    CHECK(m.find("measure the values from a nearby reference instead (μ2 = 0, for example)") != std::string::npos);
  }
  // next to 1e11 it keeps five digits (0.89999390), but a power that steep is far from the target there
  auto steep = [](double d) { return stats::pnorm((d - 0.9) / 1e-6); };
  auto gentle = [](double d) { return stats::pnorm(d - 0.9); };
  const std::string m = refusal([&] { solve::offset(1e11, 0.9, steep, 0.5, "μ1", "μ2"); });
  CHECK(m.find("μ1 cannot be placed next to μ2") == 0);
  CHECK_EQ(solve::offset(1e11, 0.9, gentle, 0.5, "μ1", "μ2"), 1e11 + 0.9);
  // the search itself, from a null value whose neighbors are coarser than the effect
  auto power = [](double x) { return stats::pnorm((x - 0.9999999999999) / 1e-17 - 2.5); };
  bool refused = false;
  try {
    solve::solve_increasing(power, 0.5, 0.9999999999999, 0.99999999999995, 1 - 1e-16, 1e-12, "ρ1", "ρ0");
  } catch (const CalcError& e) {
    refused = std::string(e.what()).find("ρ1 cannot be told apart from ρ0 = 0.9999999999999") == 0;
  }
  CHECK(refused);
}

TEST_CASE("extremes: a search for an effect never reaches past its last admissible value") {
  // ρ0 + (last − ρ0) rounds to 1 for many ρ0, where the power was taken as 1 (a false bracket: ρ1 = 1 at 100%)
  const double last = std::nextafter(1.0, 0.0);
  bool beyond = false;
  auto power = [&](double x) {
    beyond = beyond || x > last;
    return x < 1 ? 0.05 + 0.5 * (x - 0.3) : 1.0;
  };
  bool unreachable = false;
  try {
    solve::solve_increasing(power, 0.6, 0.3, 0.65, last, 1e-11, "ρ1", "ρ0");
  } catch (const CalcError& e) {
    unreachable = std::string(e.what()) == "the target power cannot be reached for any admissible value of ρ1";
  }
  CHECK(unreachable && !beyond);
  // and below: 0.3 − (0.3 − first) lands on the singular point, where 1 + 3ρ = 0
  const double first = std::nextafter(-1.0 / 3, 0.0);
  auto lower = [&](double x) {
    beyond = beyond || x < first;
    return 1 + 3 * x > 0 ? 0.05 + 0.1 * (0.3 - x) : 1.0;
  };
  unreachable = false;
  try {
    solve::solve_decreasing(lower, 0.5, 0.3, 0.0, first, 1e-11, "ρ1", "ρ0");
  } catch (const CalcError& e) {
    unreachable = std::string(e.what()) == "the target power cannot be reached for any admissible value of ρ1";
  }
  CHECK(unreachable && !beyond);
}

TEST_CASE("extremes: an answer on doubles too coarse for the root is the first that reaches the target") {
  // a power that moves by 1e-6 from one double to the next above 0.75 (whose doubles lie 1.1e-16 apart)
  const double x0 = 0.75, step = std::nextafter(x0, 1.0) - x0;
  auto power = [&](double x) { return 0.8 + 1e-6 * ((x - x0) / step - 0.4); };
  const double last = std::nextafter(1.0, 0.0);
  // from a few doubles to either side, the first double whose power reaches 0.8 (the target lies 0.4 of a step past
  // x0: x0 itself falls short by 4e-7, which a power shown to four decimals would hide, but the next one reaches it)
  // (the effect is measured from 0.5, so that every double here keeps it)
  const double d = x0 + 0.4 * step - 0.5;
  for (int off : {-5, -1, 0, 1, 5}) {
    double x = x0;
    for (int i = 0; i < std::abs(off); ++i) x = std::nextafter(x, off < 0 ? 0.0 : 1.0);
    CHECK_EQ(solve::reach(power, x, last, true, 0.8, 0.5, d, "ρ1", "ρ0", "unreachable"), std::nextafter(x0, 1.0));
  }
  // a target between two doubles whose powers both miss it by more than the precision shown
  auto steep = [&](double x) { return 0.8 + 1e-3 * ((x - x0) / step - 0.5); };
  auto message = [&](auto&& pw, double target, double from, double dd) {
    try {
      solve::reach(pw, x0, last, true, target, from, dd, "ρ1", "ρ0", "unreachable");
    } catch (const CalcError& e) {
      return std::string(e.what());
    }
    return std::string();
  };
  CHECK(message(steep, 0.8, 0.5, d).find("the power cannot be brought to the target at any ρ1 that can be stored") ==
        0);
  // a power that moves unevenly from one double to the next is not blamed on the doubles
  auto uneven = [&](double x) { return 0.8 + (std::lround((x - x0) / step) % 2 == 0 ? 1e-3 : -1e-3); };
  CHECK(message(uneven, 0.8, 0.5, d).find("the power cannot be computed to within 5e-05 of the target") == 0);
  // the double next to the root keeps fewer than four digits of a difference 3,000 doubles wide: refused as such,
  // though the next one reaches the target to the precision shown
  auto flat = [&](double x) { return 0.8 + 1e-9 * ((x - x0) / step - 0.5); };
  CHECK(message(flat, 0.8, x0 - 3000 * step, 3000.5 * step).find("keeps fewer than four significant digits") !=
        std::string::npos);
  // and one that keeps it to four digits, short of the target by less than the precision shown, is the answer when
  // the next overshoots by more
  auto coarse = [&](double x) { return 0.8 + 1e-4 * ((x - x0) / step - 0.3); };
  CHECK_EQ(solve::reach(coarse, x0, last, true, 0.8, 0.5, x0 - 0.5, "ρ1", "ρ0", "unreachable"), x0);
}

TEST_CASE("extremes: an answer near a target of 100% meets it to the finer precision shown there") {
  const double x0 = 0.75, step = std::nextafter(x0, 1.0) - x0, last = std::nextafter(1.0, 0.0);
  for (double target : {1 - 1e-11, 1 - 1e-10}) {
    const double slack = solve::power_slack(target);
    // a power that moves by 4e-12 per double, overshooting at x0 by more than the precision shown
    auto rising = [&](double x) { return std::fmin(1.0, target + 6e-12 + 4e-12 * ((x - x0) / step)); };
    const double x = solve::reach(rising, x0, last, true, target, 0.5, x0 - 0.5, "ρ1", "ρ0", "unreachable");
    CHECK(std::fabs(rising(x) - target) <= slack);
    // and one short by 8e-10 at x0, which the 1e-9 that counts as reaching an ordinary target would have accepted
    auto flat = [&](double y) { return target - 8e-10 + 1e-10 * ((y - x0) / step); };
    const double y = solve::reach(flat, x0, last, true, target, 0.5, x0 - 0.5, "ρ1", "ρ0", "unreachable");
    CHECK(std::fabs(flat(y) - target) <= slack);
  }
}

TEST_CASE("extremes: a power within the precision shown over many doubles is refused for the difference it needs") {
  // next to 1 − 1e-10 the power moves by 1e-13 from one double to the next, so that all the doubles tried count as
  // reaching the target: when the double nearest the root does not keep the difference, none does (the refusal was
  // "the answer cannot be computed for these inputs")
  const double target = 1 - 1e-10, x0 = 0.75, step = std::nextafter(x0, 1.0) - x0, last = std::nextafter(1.0, 0.0);
  auto power = [&](double x) { return target + 1e-13 * ((x - x0) / step); };
  std::string message;
  try {
    solve::reach(power, x0, last, true, target, x0 - 3000.5 * step, 3000.5 * step, "ρ1", "ρ0", "unreachable");
  } catch (const CalcError& e) {
    message = e.what();
  }
  CHECK(message.find("ρ1 cannot be told apart from ρ0 = ") == 0);
  CHECK(message.find("keeps fewer than four significant digits") != std::string::npos);
  // and when it does keep it, it is the answer, though the root converted to a value lay 20 doubles off
  const double from = x0 - 30000.25 * step;
  double x = x0;
  for (int i = 0; i < 20; ++i) x = std::nextafter(x, 1.0);
  CHECK_EQ(solve::reach(power, x, last, true, target, from, 30000.25 * step, "ρ1", "ρ0", "unreachable"), x0);
}

TEST_CASE("extremes: n times a ratio rounded up is a whole product however large") {
  // the absolute slack of 1e-9 was less than the rounding of a product above 2^24: 1e8 × 1.1 is
  // 110000000.00000001, and a group of n × 1.1 got one subject more than the rule gives
  CHECK_EQ(solve::round_up(1e8 * 1.1), 110000000.0);
  CHECK_EQ(solve::round_up(3e7 * 1.1), 33000000.0);
  CHECK_EQ(solve::round_up(1e12 * 0.7), 700000000000.0);
  CHECK_EQ(solve::round_up(110000000.5), 110000001.0);
  // and a whole product stays whole past 2^50, where four units in the last place are a whole unit or more: 1e16
  // became 9999999999999992
  CHECK_EQ(solve::round_up(1e8 * 1e8), 1e16);
  CHECK_EQ(solve::round_up(0x1p50), 0x1p50);
  CHECK_EQ(solve::round_up(0x1p52 + 1), 0x1p52 + 1);
  CHECK_EQ(solve::round_up(1e30), 1e30);
}

TEST_CASE("extremes: the correlation's distribution on Fisher's scale") {
  for (double r : {-0.9, 0.0, 0.3, 0.95})
    for (double rho : {-0.5, 0.0, 0.6})
      CHECK_EQ(stats::pcorr_z(std::atanh(r), 20, std::atanh(rho), false), stats::pcorr(r, 20, rho, false));
  // near 1, where the doubles next to the critical values are too coarse, the tail on Fisher's scale depends on
  // the distance from ζ alone, as ρ = tanh ζ nears 1 (the same at ζ = 10 and ζ = 14)
  CHECK_CLOSE(stats::pcorr_z(14.5, 50, 14, false), stats::pcorr_z(10.5, 50, 10, false), 1e-7, 0);
  CHECK_EQ(stats::pcorr_z(INFINITY, 50, 14, false), 0.0);
  CHECK_EQ(stats::pcorr_z(-INFINITY, 50, 14, false), 1.0);
}

TEST_CASE("extremes: cluster sizes at a CV whose square underflows are all the mean") {
  // 1/CV² was infinite, and every size 10⁹
  const procs::simc::ClusterSizes sizes(20, 1e-160);
  sim::Rng rng(1, 2);
  CHECK_EQ(sizes.draw(rng), 20L);
  CHECK(sizes.grid() == std::vector<double>{20});
}

TEST_CASE("extremes: odds-ratio limits far from 1 never lower the power, and an empty group 2 is refused") {
  auto power = [](const std::string& procedure, const std::string& inputs) {
    return run(procedure, "power", "{" + inputs + "}").get("rows").elements().at(0).get("power").as_number(NAN);
  };
  // the restricted estimates took the root of the wrong sign once a limit was some 10^13 from 1, and 1 − P1 lost
  // its digits: an upper limit of 10^15 gave 0% power, though widening a limit can never lower it
  for (const char* method : {"normal", "exact"}) {
    const std::string common =
        std::string(R"("scale":"or","p1":"0.6","p2":"0.6","n1":"100","method":")") + method + "\"";
    const double upper = power("prop.two.equiv", common + R"(,"orl":"0.67","oru":"1e6")");
    const double lower = power("prop.two.equiv", common + R"(,"orl":"1e-6","oru":"1.5")");
    CHECK(upper > 0.4 && lower > 0.4);
    for (const char* limit : {"1e13", "1e15", "1e30"})
      CHECK_CLOSE(power("prop.two.equiv", common + R"(,"orl":"0.67","oru":")" + limit + "\""), upper, 1e-9, 0);
    for (const char* limit : {"1e-20", "1e-30"})
      CHECK_CLOSE(power("prop.two.equiv", common + R"(,"orl":")" + limit + R"(","oru":"1.5")"), lower, 1e-9, 0);
    // a margin whose boundary P1 rounds to 1 was refused
    const double ni = power("prop.two.ni", common + R"(,"h1":"lower","or0":"1e8")");
    CHECK(ni > 0.99);
    CHECK_CLOSE(power("prop.two.ni", common + R"(,"h1":"lower","or0":"1e30")"), ni, 1e-9, 0);
  }
  // a cluster ratio that left group 2 no clusters gave a power of 100%; rounded up, it leaves one
  const auto tiny = run("prop.two.cluster", "power", R"({"k1":"10","kr":"1e-12"})").get("rows").elements().at(0);
  CHECK_EQ(tiny.get("k2").as_number(NAN), 1.0);
  CHECK(tiny.get("power").as_number(NAN) < 1);
}

TEST_CASE("extremes: difference margins next to 0 and ratio margins far from 1 keep the tables' variances") {
  auto row = [](const std::string& inputs) {
    return run("prop.two.ni", "power", "{" + inputs + "}").get("rows").elements().at(0);
  };
  // in the table where every subject responded, 1 − P1 = −D0 was lost to rounding below about 5.6e-17: no variance,
  // and an infinite statistic from the numerator −D0
  const std::string diff = R"("scale":"diff","method":"exact","p1":"0.95","p2":"0.95","n1":"20","alpha":"0.025")";
  const double near = row(diff + R"(,"h1":"upper","d0":"-1e-16")").get("power").as_number(NAN);
  CHECK(near > 0.005 && near < 0.007);
  for (const char* d0 : {"-5e-17", "-1e-20", "-1e-30"})
    CHECK_CLOSE(row(diff + R"(,"h1":"upper","d0":")" + d0 + "\"").get("power").as_number(NAN), near, 1e-12, 0);
  CHECK_CLOSE(row(diff + R"(,"h1":"lower","d0":"1e-20")").get("power").as_number(NAN), near, 1e-12, 0);
  // the smaller root of the ratio's quadratic cancelled to 0 below a margin of about 1e-16, with the same effect
  const json::Value ratio =
      row(R"("scale":"ratio","method":"exact","p1":"0.001","p2":"0.2","n1":"50","alpha":"0.025","h1":"lower",)"
          R"("r0":"1e-17")");
  CHECK(ratio.get("power").as_number(NAN) < 1e-12);
  CHECK(ratio.get("actual_alpha").as_number(NAN) < 1e-12);
  // the exact size at the boundary itself, not at a P1 kept 1e-9 from 1 (150-digit enumeration)
  const json::Value far =
      row(R"("scale":"or","method":"exact","p1":"0.6","p2":"0.6","n1":"4","alpha":"0.025","h1":"lower","or0":"1e30")");
  CHECK_CLOSE(far.get("actual_alpha").as_number(NAN), 2.5984e-30, 1e-9, 0);
}

TEST_CASE("extremes: requests that took minutes or gigabytes end at once with an answer or a message") {
  struct Case {
    const char* request;
    const char* says;  // part of the message, or "" for an answer
  };
  const std::string lists = [] {
    std::string t;
    for (int i = 0; i < 50; ++i) t += (i ? ";" : "") + std::string("1 to 10000");
    return t;
  }();
  const std::string shorter = [] {
    std::string t;
    for (int i = 0; i < 10; ++i) t += (i ? ";" : "") + std::string("1 to 1000");
    return t;
  }();
  const std::string many = R"({"procedure":"means.two.t","solve_for":"power","inputs":{"mu1":"0.1 to 2.1 by 0.01",)"
                           R"("mu2":"0 to 1.4 by 0.01"},"options":{"power_curve":false,"max_scenarios":)";
  const std::vector<std::pair<std::string, std::string>> cases = {
      // Fisher's test with one proportion near 0: few table totals, each with a wide null distribution (without
      // the power curve, whose smaller designs are enumerated: each is bounded, and Stop ends them)
      {R"({"procedure":"prop.two","solve_for":"power","options":{"power_curve":false},)"
       R"("inputs":{"test":"fisher","es":"or","or1":"1e-9","n1":"10000000"}})", ""},
      {R"({"procedure":"prop.cc","solve_for":"or","options":{"power_curve":false},)"
       R"("inputs":{"test":"fisher","ncase":"1000000","m":"100"}})", ""},
      // a group above 2^31 − 1 is not enumerated (its count is an int)
      {R"({"procedure":"prop.two","solve_for":"power","inputs":{"alloc":"ratio","n1":"1000000000","r":"3"}})", ""},
      {R"({"procedure":"prop.two.equiv","solve_for":"power","inputs":{"alloc":"ratio","n1":"1000000000","r":"3"}})", ""},
      {R"({"procedure":"means.two.repeated","solve_for":"power","inputs":{"struct":"ar1","m":"1000000000"}})", ""},
      // lists: 10,000 values in all the alternatives, and 1,000,000 over all the scenarios
      {R"({"procedure":"means.anova","solve_for":"n","inputs":{"means":")" + lists + R"("}})",
       "too many values in all the lists together"},
      {R"({"procedure":"means.anova","solve_for":"n","inputs":{"means":")" + shorter +
           R"(","alpha":"0.01 to 0.2 by 0.001"}})",
       "Too many values in the lists over all the combinations"},
      // the option lowers the limit of scenarios and never raises it
      {many + "1e12}}", "limit 20000"},
      {many + "-1}}", "limit 1)"},
      // with the power curve (on here): its points use the normal approximation when enumerating them all would
      // take long (74 s), and the paired diagnostic search computes each number of pairs once (3 s)
      {R"({"procedure":"prop.cc","solve_for":"n","inputs":{"test":"fisher","m":"100","or":"1.05"}})", ""},
      {R"({"procedure":"diag.paired","solve_for":"n","inputs":{"method":"exact","goal":"se","se1":"0.9","se2":"0.8",)"
       R"("prev":"0.001"}})", ""},
      // Fisher's test costs more a table than the cells it enumerates: its exact search took 20 s
      {R"({"procedure":"prop.two","solve_for":"n","options":{"power_curve":false},)"
       R"("inputs":{"test":"fisher","method":"exact","p1":"0.428","p2":"0.4"}})", ""},
      // and its search for P1, some 40 exact evaluations at nearly a million a group: 58 s
      {R"({"procedure":"prop.two","solve_for":"p1","options":{"power_curve":false},)"
       R"("inputs":{"test":"fisher","method":"exact","n1":"918515","p2":"0.7814"}})", ""},
      // a noncentral F with an enormous noncentrality: Poisson terms with a shape near 1e25, where the continued
      // fraction of the incomplete beta cannot converge (it ran to its cap, 2 s a term, for minutes)
      {R"({"procedure":"means.anova","solve_for":"power","options":{"power_curve":false},)"
       R"("inputs":{"estype":"f","f":"2e12","k":"2","alloc":"ratio","ratios":"1 0.5","n":"2","alpha":"1e-9"}})", ""},
      {R"({"procedure":"means.anova","solve_for":"n","options":{"power_curve":false},)"
       R"("inputs":{"means":"0 1","sigma":"1e-13","alloc":"ratio","ratios":"1 0.5","alpha":"1e-9"}})", ""},
      // 20,000 group-sequential designs with 20 analyses: about 25 minutes of work, solving for the sample size or
      // for the power (each alpha, and each power with β-spending futility, is a design of its own)
      {R"({"procedure":"gs.means","solve_for":"n","inputs":{"looks":"20","alpha":"0.0001 to 0.05 by 0.0000049995",)"
       R"("power":"0.8 0.9"}})",
       "designs this long to compute"},
      {R"({"procedure":"gs.means","solve_for":"power","inputs":{"looks":"20","alternative":"upper","futility":"beta",)"
       R"("binding":"binding","alpha":"0.0001 to 0.05 by 0.0000049995","power":"0.8 0.9","n1":"50"}})",
       "designs this long to compute"},
  };
  for (const auto& c : cases) {
    Token token;
    const auto t0 = std::chrono::steady_clock::now();
    const json::Value res = run_request(json::parse(c.first), &token);
    const std::string m = message(res);
    if (seconds_since(t0) > 5 || (c.second.empty() ? !m.empty() : m.find(c.second) == std::string::npos))
      ::ssc::test::report_failure(c.first.substr(0, 160) + ": '" + m + "'", __FILE__, __LINE__);
    ::ssc::test::count_check();
  }
  {  // typed text is quoted briefly in a message (a page lays out a long unbroken word for seconds)
    Token token;
    const json::Value res = run_request(json::parse(R"({"procedure":"means.two.t","solve_for":"n","inputs":{"power":")" +
                                                    std::string(100000, 'W') + R"("}})"),
                                        &token);
    CHECK(res.get("input_errors").get("power").as_string().size() < 100);
  }
  {  // detail tables come with the first 1,000 rows only (20,000 rows of boundaries took a gigabyte)
    Token token;
    const json::Value res = run_request(
        json::parse(R"({"procedure":"gs.means","solve_for":"power","options":{"power_curve":false},)"
                    R"("inputs":{"looks":"3","mu1":"0.002 to 2.4 by 0.002","n1":"50"}})"),
        &token);
    const auto& rows = res.get("rows").elements();
    CHECK_EQ(rows.size(), static_cast<std::size_t>(1200));
    CHECK(rows.at(999).get("_tables").size() > 0);
    CHECK(rows.at(1000).find("_tables") == nullptr);
  }
  {  // too large to enumerate: the normal approximation, and a note that says so (it used to be a power of 0)
    Token token;
    const json::Value res = run_request(json::parse(cases[2].first), &token);
    const json::Value& row = res.get("rows").elements().at(0);
    CHECK(row.get("power").as_number(0) > 0.99);
    CHECK(row.get("_notes").size() == 1);
  }
}

TEST_CASE("extremes: Stop ends an exact enumeration on one thread at once") {
  // Fisher's test with many table totals, and an equivalence test with 65 rows of 380,000 tables each
  for (const char* inputs :
       {R"({"alternative":"upper","test":"fisher","es":"or","or1":"1e-9","n1":"2000000"})",
        R"({"method":"exact","alloc":"ratio","n1":"64","r":"32812500","p1":"0.5","p2":"0.5","dl":"-0.4","du":"0.4"})"}) {
    Token token;
    std::thread stop([&] {
      std::this_thread::sleep_for(std::chrono::milliseconds(300));
      token.cancel();
    });
    const std::string procedure = std::string(inputs).find("fisher") != std::string::npos ? "prop.two" : "prop.two.equiv";
    const auto t0 = std::chrono::steady_clock::now();
    const json::Value res = run_request(json::parse(R"({"procedure":")" + procedure +
                                                    R"(","solve_for":"power","options":{"threads":1,"power_curve":false},"inputs":)" +
                                                    inputs + "}"),
                                        &token);
    const double took = seconds_since(t0);
    stop.join();
    CHECK(res.get("cancelled").as_bool(false));
    CHECK(took < 1.0);
  }
}

TEST_CASE("extremes: a request object with many keys is read in linear time") {
  std::string text = "{";
  for (int i = 0; i < 200000; ++i) text += "\"k" + std::to_string(i) + "\":1,";
  text += R"("a":1,"b":2,"a":3})";
  const auto t0 = std::chrono::steady_clock::now();
  const json::Value v = json::parse(text);
  CHECK(seconds_since(t0) < 2);  // about 0.1 s; it was minutes
  CHECK_EQ(v.size(), static_cast<std::size_t>(200002));
  CHECK_EQ(v.get("a").as_number(0), 3.0);  // a repeated key keeps its place and takes the last value
  CHECK_EQ(v.items()[200000].first, std::string("a"));
}

TEST_CASE("extremes: counts beyond a 32-bit integer give the same answers on every platform") {
  auto request = [](const std::string& text) {
    Token token;
    return run_request(json::parse(text), &token);
  };
  {  // a dropout rate of 1e-12: the last visit, as an int, overflowed and dropped every patient (power 0)
    const json::Value res = request(
        R"({"procedure":"sim.mmrm","solve_for":"power","options":{"power_curve":false},)"
        R"("inputs":{"n1":"100","drop1":"1e-12","drop2":"1e-12","reps":"200"}})");
    CHECK(res.get("rows").elements().at(0).get("power").as_number(0) > 0.6);
  }
  {  // 1e9 subjects per cluster-period: the observations per cluster fit no 32-bit long (WebAssembly, Windows)
    const json::Value res = request(
        R"({"procedure":"sim.sw","solve_for":"power","options":{"power_curve":false},)"
        R"("inputs":{"m":"1000000000","g":"2","delta":"4.242640687119285e-05","reps":"100"}})");
    CHECK(res.get("rows").elements().at(0).get("mean_df").as_number(0) > 2e10);
  }
  {  // a power curve stays within the input's range (whole numbers up to 1e9), where the integer searches hold
    const json::Value res = request(R"({"procedure":"prop.one","solve_for":"power","inputs":{"n":"800000000"}})");
    const json::Value& series = res.get("plots").elements().at(0).get("series").elements().at(0);
    const auto& xs = series.get("x").elements();
    const auto& ys = series.get("y").elements();
    bool ok = !xs.empty();
    for (std::size_t i = 0; i < xs.size(); ++i)
      if (xs[i].as_number(0) > 1e9 || (xs[i].as_number(0) > 1e6 && ys[i].as_number(0) < 0.5)) ok = false;
    CHECK(ok);
  }
}

TEST_CASE("extremes: distribution functions outside their domain give NaN at once") {
  // a negative shape or an infinite count made the internal routines loop for ever (ssc_eval passes them on)
  const auto t0 = std::chrono::steady_clock::now();
  double out = 0;
  const double bad[][4] = {{-INFINITY, 0, 0, 0}, {-1e300, 0, 0, 0}};
  for (const auto& a : bad) {
    CHECK(stats::eval("stirlerr", a, 1, out));
    CHECK(std::isnan(out));
  }
  const double dbeta_args[3] = {0.5, -1e300, 2}, pbeta_args[4] = {0.5, -1e300, 2, 1}, qbinom_args[3] = {0.5, 1e300, NAN};
  CHECK(stats::eval("dbeta", dbeta_args, 3, out) && std::isnan(out));
  CHECK(stats::eval("pbeta", pbeta_args, 4, out) && std::isnan(out));
  CHECK(stats::eval("qbinom", qbinom_args, 3, out) && std::isnan(out));
  const double ok_args[4] = {0.5, 2, 3, 1};
  CHECK(stats::eval("pbeta", ok_args, 4, out) && out > 0 && out < 1);
  CHECK(seconds_since(t0) < 1);
}

TEST_CASE("extremes: the bivariate normal at an infinite limit is that of the other variable") {
  // an endpoint whose power is 1 to double precision has Φ⁻¹ = −∞, and the series met 0·∞ (NaN)
  for (double r : {-0.95, -0.5, 0.0, 0.5, 0.95})
    for (double k : {-2.0, 0.0, 1.5}) {
      CHECK(stats::bvn_upper(-INFINITY, k, r) == stats::pnorm(-k));
      CHECK(stats::bvn_upper(k, -INFINITY, r) == stats::pnorm(-k));
      CHECK(stats::bvn_upper(INFINITY, k, r) == 0);
      CHECK(stats::bvn_upper(k, INFINITY, r) == 0);
    }
  CHECK(stats::bvn_upper(-INFINITY, -INFINITY, 0.3) == 1);
  CHECK(std::isnan(stats::bvn_upper(NAN, 0, 0.3)));
}

TEST_CASE("extremes: a NaN integrand ends the adaptive integration at once") {
  // a NaN error estimate never passed the test, so every panel was split down to max_depth (6e10 calls)
  long calls = 0;
  const double nan = numerics::integrate([&](double) { ++calls; return NAN; }, 0, 1);
  CHECK(std::isnan(nan));
  CHECK(calls < 2000);
  // a 0/0 at the midpoint, a node of the first panel, still clears
  const double si = numerics::integrate([](double x) { return std::sin(x) / x; }, -1, 1);
  CHECK(std::fabs(si - 1.8921661407343660) < 1e-13);
}

TEST_CASE("extremes: a search reads a long text in part and lists each matched term once") {
  // a megabyte of description took 1.5 s and 700 MB, and its answer held 3,000,000 matched terms (26 MB)
  std::string q;
  for (int i = 0; i < 100000; ++i) q += "odds ratio ";
  const auto t0 = std::chrono::steady_clock::now();
  const json::Value res = search_json(q, 20);
  CHECK(seconds_since(t0) < 1);
  CHECK(res.get("query").as_string().size() <= 2000);
  CHECK(res.get("hits").size() == 20);
  CHECK(res.get("hits").at(0).get("terms").size() <= 2);
  CHECK(res.dump().size() < 20000);
}

TEST_CASE("extremes: the incomplete beta with shapes beyond what its series resolve") {
  // both shapes large, far from the mean: a tail far below the smallest double (and 1e300 gave 1, not 0)
  CHECK(stats::pbeta(0.4, 1e300, 1e300, true) == 0);
  CHECK(stats::pbeta(0.6, 1e300, 1e300, true) == 1);
  CHECK(stats::pbeta(0.4, 1e7, 1e7, false) == 1);
  // one shape astronomically larger than the other: the gamma limit (mpmath, 50 digits: 8.0800219368998048e-9)
  CHECK_CLOSE(stats::pbeta(3e-13, 3315, 1e16, true), 8.0800219368998048e-9, 1e-12, 0);
  // and the other way round, given the complement exactly (1 − x in doubles would lose its digits)
  CHECK_CLOSE(stats::pbeta_raw(1 - 3e-13, 1e16, 3315, false, 3e-13), 8.0800219368998048e-9, 1e-12, 0);
}

TEST_CASE("extremes: binomial, F and beta probabilities at a billion trials agree with R") {
  // shapes beyond 1e8 took the series or continued fraction near the mean, which needs about √a terms
  // (0.5008 for 0.5 at a billion trials; the reference values are R's, from TOMS 708)
  CHECK_CLOSE(stats::pbinom(5e8, 1e9, 0.5, true), 0.50001261566260691, 1e-12, 0);
  CHECK_CLOSE(stats::pbinom(5e8 + std::sqrt(2.5e8), 1e9, 0.5, true), 0.8413464554586505, 1e-12, 0);
  CHECK_CLOSE(stats::pbinom(5e8 - 3 * std::sqrt(2.5e8), 1e9, 0.5, true), 0.0013498040997050719, 1e-11, 0);
  CHECK_CLOSE(stats::pf(1, 1e9, 1e9, true), 0.5, 1e-12, 0);
  // a − (a + b)·x is rounded once where the compiler fuses it (Apple arm64) and twice elsewhere (as in R): one
  // rounding of (a + b)·x ≈ 1e10 moves the probability by about 1e-11, so the check is to the problem's conditioning
  CHECK_CLOSE(stats::pbeta(1.0 / 3 + std::sqrt(2e10 / (9e20 * (3e10 + 1))), 1e10, 2e10, true), 0.50000453231113418,
              2e-11, 0);
}
