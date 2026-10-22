// What a user can type at the extremes (a ratio of 1e-300, a noncentrality of 1e300, a count beyond a
// billion) ends at once with an answer or a message, Stop ends what takes long, and a long simulation keeps only
// a block of studies in memory with the same result as keeping them all.
#include <chrono>
#include <cmath>
#include <string>
#include <thread>
#include <vector>

#include "engine/runner.hpp"
#include "engine/solve.hpp"
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
    auto n2 = [&](double n1) { return std::ceil(r * n1 - 1e-9); };
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
  const Case cases[] = {
      {"means.two.cluster", "k", R"({"kr":"1e-300"})", "cluster ratio is too small"},
      {"prop.two.cluster", "k", R"({"kr":"1e-300"})", "cluster ratio is too small"},
      {"diag.auc", "n", R"({"ratio":"1e-300"})", "ratio of controls to cases is too small"},
      {"sim.lmm.cluster", "k", R"({"kr":"1e-300"})", "cluster ratio is too small"},
      {"sim.glmm.cluster", "k", R"({"kr":"1e-300"})", "cluster ratio is too small"},
      {"sim.lmm.longitudinal", "n", R"({"r":"1e-300"})", "allocation ratio is too small"},
      {"sim.mmrm", "n", R"({"r":"1e-300"})", "allocation ratio is too small"},
      {"sim.means.robust", "n", R"({"r":"1e-300"})", "allocation ratio is too small"},
      {"sim.ordinal", "n", R"({"r":"1e-300"})", "allocation ratio is too small"},
      {"sim.winratio", "n", R"({"r":"1e-300"})", "allocation ratio is too small"},
      {"prop.one", "power", R"({"n":"1e10"})", "at most 1,000,000,000"},
      {"corr.one", "power", R"({"n":"1e300"})", "at most 1,000,000,000"},
      {"prop.one.simon", "power", R"({"r1":"1","n1":"999999999","r":"5","n":"1000000000"})", "up to 100,000"},
      {"prop.paired.mcnemar", "power", R"({"n":"10000000"})", "up to 1,000,000 pairs"},
      {"means.two.ancova", "mu1", R"({"ncov":"1000"})", "number of covariates"},
      {"means.anova2", "power", R"({"a":"20","b":"214748365"})", "at most 100"},  // A × B wrapped around to 4
      {"means.two.welch", "power", R"({"sd1":"1e300"})", "cannot be computed"},
      {"means.two.t", "mu1", R"({"sigma":"1e300"})", ""},  // the search bound 1e9·σ overflowed to infinity
      {"sim.multiarm", "n", R"({"alpha":"1e-300"})", ""},  // 1 − α rounds to 1: Dunnett's value is Bonferroni's
      {"rates.one.poisson", "power", R"({"lam0":"1e300"})", ""},
      {"rates.one.poisson", "n", R"({"time":"1e300"})", ""},
      {"means.crossover", "n", R"({"delta":"1e300"})", ""},
      {"means.anova", "power", R"({"sigma":"1e-300"})", ""},
      {"means.rm.anova", "n", R"({"sigma":"1e-300"})", ""},
      // loops that a value rounding to 0 or infinity stopped from advancing (and that Stop could not end)
      {"sim.surv.nph", "power", R"({"med2":"5e-324"})", "too short to compute with"},
      {"sim.surv.nph", "n", R"({"med2":"1e-320"})", "too short to compute with"},
      {"rates.one.poisson", "power", R"({"alpha":"5e-324","test":"exact"})", "too small to compute with"},
      {"sim.means.robust", "power", R"({"dist":"lognormal","skew":"1e-9"})", ""},
      {"sim.means.robust", "n", R"({"dist":"gamma","skew":"-1e-8"})", ""},
      {"diag.paired", "n", R"({"method":"exact","goal":"se","prev":"1e-303"})", "cannot be reached"},
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
