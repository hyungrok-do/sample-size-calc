// Two survival curves with non-proportional hazards by simulation
// (sim.surv.nph): the test statistics on a classic data set, the multivariate
// normal probability behind MaxCombo, the quick approximation and simulated
// power against the logrank formula under proportional hazards, the type I
// error of every test, the ranking of the tests under a delayed effect,
// reproducibility and the sample-size search.
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "engine/runner.hpp"
#include "procedures/sim_survival.hpp"
#include "stats/stats.hpp"
#include "testing.hpp"
#include "util/json.hpp"
#include "util/numtext.hpp"

using namespace ssc;
namespace nph = ssc::procs::nph;

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

// A request for `procedure`, with the inputs and options given as the bodies of JSON objects.
std::string request(const char* procedure, const char* solve, const std::string& inputs,
                    const std::string& options = "") {
  return std::string(R"({"procedure":")") + procedure + R"(","solve_for":")" + solve + R"(","inputs":{)" + inputs +
         R"(},"options":{)" + options + "}}";
}

std::string nph_request(const char* solve, const std::string& inputs, const std::string& options = "") {
  return request("sim.surv.nph", solve, inputs, options);
}

}  // namespace

// ------------------------------------------------------------------ statistics of one data set

TEST_CASE("nph: weighted logrank and RMST statistics of the 6-MP leukemia data") {
  // remission times in weeks (Freireich et al. 1963; negative = censored), with many ties; the reference
  // values come from a direct implementation, and the logrank chi-square 16.79 is the classic result
  const double mp[21] = {6, 6, 6, -6, 7, -9, 10, -10, -11, 13, 16, -17, -19, -20, 22, 23, -25, -32, -32, -34, -35};
  const double placebo[21] = {1, 1, 2, 2, 3, 4, 4, 5, 5, 8, 8, 8, 8, 11, 11, 12, 12, 15, 17, 22, 23};
  std::vector<nph::Obs> obs;
  for (double t : mp) obs.push_back({std::fabs(t), 0, t > 0 ? 1 : 0});
  for (double t : placebo) obs.push_back({t, 1, 1});
  std::stable_sort(obs.begin(), obs.end(), [](const nph::Obs& a, const nph::Obs& b) { return a.time < b.time; });
  const double rho[4] = {0, 0, 1, 1}, gamma[4] = {0, 1, 0, 1};
  const nph::WeightedLogrank w = nph::weighted_logrank(obs, 21, 21, 4, rho, gamma);
  CHECK_CLOSE(w.u[0], -10.25050094803113, 1e-12, 0);  // observed − expected events on 6-MP
  CHECK_CLOSE(w.v[0], 6.2569605736755, 1e-12, 0);
  CHECK_CLOSE(w.u[0] * w.u[0] / w.v[0], 16.792940989216532, 1e-12, 0);
  CHECK_CLOSE(w.u[1], -3.3734559104592137, 1e-12, 0);  // FH(0, 1)
  CHECK_CLOSE(w.v[5], 0.8721500239394372, 1e-12, 0);
  CHECK_CLOSE(w.u[2], -6.877045037571917, 1e-12, 0);  // FH(1, 0)
  CHECK_CLOSE(w.v[10], 3.2713049093714295, 1e-12, 0);
  CHECK_CLOSE(w.u[3], -1.6646510205668954, 1e-12, 0);  // FH(1, 1)
  CHECK_CLOSE(w.v[15], 0.2174833381921807, 1e-12, 0);
  CHECK_CLOSE(w.v[1], 1.928902844121754, 1e-12, 0);
  CHECK_CLOSE(w.v[4], w.v[1], 1e-15, 0);
  CHECK_CLOSE(w.v[7], 0.413761666303133, 1e-12, 0);
  // the weights S and 1 − S add up to 1, so FH(1, 0) = FH(0, 0) − FH(0, 1) and MaxCombo's matrix is singular
  CHECK_CLOSE(w.u[2], w.u[0] - w.u[1], 0, 1e-12);
  // the logrank test alone (k = 1) gives the same score
  const double zero = 0;
  CHECK_CLOSE(nph::weighted_logrank(obs, 21, 21, 1, &zero, &zero).u[0], w.u[0], 1e-14, 0);

  const nph::Rmst r10 = nph::restricted_mean(obs, 21, 21, 10);
  CHECK_CLOSE(r10.rmst[0], 9.277310924369749, 1e-12, 0);
  CHECK_CLOSE(r10.var[0], 0.10677805477109285, 1e-10, 0);
  CHECK_CLOSE(r10.rmst[1], 6.6190476190476195, 1e-12, 0);
  CHECK_CLOSE(r10.var[1], 0.5373069862865782, 1e-10, 0);
  const nph::Rmst r23 = nph::restricted_mean(obs, 21, 21, 23);  // every placebo subject has relapsed by 23
  CHECK_CLOSE(r23.rmst[0], 17.909243697478995, 1e-12, 0);
  CHECK_CLOSE(r23.var[0], 2.4123991085312047, 1e-10, 0);
  CHECK_CLOSE(r23.rmst[1], 8.666666666666668, 1e-12, 0);
  CHECK_CLOSE(r23.var[1], 1.8972033257747551, 1e-10, 0);
}

TEST_CASE("nph: multivariate normal probabilities") {
  const double inf = stats::INF;
  // two dimensions: the bivariate normal distribution function
  for (double r : {-0.9, -0.3, 0.0, 0.5, 0.95})
    for (const auto& hk : {std::make_pair(0.0, 0.0), std::make_pair(1.2, -0.4), std::make_pair(-1.5, 2.0)}) {
      const double corr[4] = {1, r, r, 1}, lo[2] = {-inf, -inf}, hi[2] = {hk.first, hk.second};
      CHECK_CLOSE(nph::mvn_probability(2, lo, hi, corr), stats::pbvnorm(hk.first, hk.second, r), 0, 1e-5);
    }
  {  // a rectangle, by inclusion and exclusion
    const double r = 0.7, corr[4] = {1, r, r, 1}, lo[2] = {-1, -0.5}, hi[2] = {1.5, 2};
    const double want = stats::pbvnorm(1.5, 2, r) - stats::pbvnorm(-1, 2, r) - stats::pbvnorm(1.5, -0.5, r) +
                        stats::pbvnorm(-1, -0.5, r);
    CHECK_CLOSE(nph::mvn_probability(2, lo, hi, corr), want, 0, 1e-5);
  }
  {  // independence: products of normal probabilities
    double corr[16] = {};
    for (int i = 0; i < 4; ++i) corr[i * 5] = 1;
    const double lo[4] = {-1, -inf, -2, -0.3}, hi[4] = {1, 0.4, inf, 1.7};
    double want = 1;
    for (int i = 0; i < 4; ++i) want *= stats::pnorm(hi[i]) - stats::pnorm(lo[i]);
    CHECK_CLOSE(nph::mvn_probability(4, lo, hi, corr), want, 1e-12, 0);
  }
  // equicorrelated with ρ = 1/2: P(all X_i < 0) = 1/(m + 1), e.g. 1/4 in three dimensions
  for (int m = 2; m <= 4; ++m) {
    double corr[16], lo[4], hi[4];
    for (int i = 0; i < m * m; ++i) corr[i] = i % (m + 1) == 0 ? 1 : 0.5;
    for (int i = 0; i < m; ++i) lo[i] = -inf, hi[i] = 0;
    CHECK_CLOSE(nph::mvn_probability(m, lo, hi, corr), 1.0 / (m + 1), 0, 1e-5);
  }
  {  // singular matrices: X3 = (X1 ± X2)/√2 with X1, X2 independent
    const double s = std::sqrt(0.5), lo[4] = {-inf, -inf, -inf, -inf}, hi[4] = {0, 0, 0, 0};
    const double plus[9] = {1, 0, s, 0, 1, s, s, s, 1}, minus[9] = {1, 0, s, 0, 1, -s, s, -s, 1};
    CHECK_CLOSE(nph::mvn_probability(3, lo, hi, plus), 0.25, 0, 1e-5);    // X1, X2 < 0 imply X3 < 0
    CHECK_CLOSE(nph::mvn_probability(3, lo, hi, minus), 0.125, 0, 1e-5);  // X1 < X2 < 0
    const double four[16] = {1, 0, 0, s, 0, 1, 0, s, 0, 0, 1, 0, s, s, 0, 1};  // X4 = (X1 + X2)/√2, X3 apart
    CHECK_CLOSE(nph::mvn_probability(4, lo, hi, four), 0.125, 0, 1e-5);
  }
  {  // mvtnorm::pmvnorm with Miwa's algorithm
    const double r3[9] = {1, .3, .6, .3, 1, .2, .6, .2, 1}, lo3[3] = {-1, -inf, -0.5}, hi3[3] = {1.5, 0.8, 2};
    CHECK_CLOSE(nph::mvn_probability(3, lo3, hi3, r3), 0.434091956878, 0, 1e-5);
    const double r4[16] = {1, .5, .3, .2, .5, 1, .4, .3, .3, .4, 1, .6, .2, .3, .6, 1};
    const double lo4[4] = {-inf, -inf, -inf, -inf}, hi4[4] = {1, 0.5, 1.5, 0.2}, lo4b[4] = {-2, -1.5, -2.5, -1};
    CHECK_CLOSE(nph::mvn_probability(4, lo4, hi4, r4), 0.405747776315, 0, 1e-5);
    CHECK_CLOSE(nph::mvn_probability(4, lo4b, hi4, r4), 0.243429535391, 0, 1e-5);
    const double rn[16] = {1, -.4, .2, .1, -.4, 1, -.3, .25, .2, -.3, 1, -.5, .1, .25, -.5, 1};
    const double lon[4] = {-1, -1, -1, -1}, hin[4] = {2, 1, 0.5, 1.5};
    CHECK_CLOSE(nph::mvn_probability(4, lon, hin, rn), 0.263962518849, 0, 1e-5);
  }
}

TEST_CASE("nph: the probability of leaving a box, also far below the lattice's error") {
  // P(max X_j ≥ c) for four equicorrelated standard normals, from the one-dimensional integral of
  // −expm1(4 ln Φ((c − √ρ z)/√(1 − ρ))) (scipy, relative error 1e-13); 1 − mvn_probability was 9% off at
  // ρ = 0.7 and c = 4, and 0 from c = 8.5
  struct Case {
    double rho, c, want;
  };
  const Case cases[] = {{0.7, 4, 1.1432734108004246e-4},   {0.7, 8, 2.4859854186757747e-15},
                        {0.7, 11, 7.6426013306743705e-28}, {0.95, 11, 6.934087999642074e-28},
                        {0.3, 6, 3.9463097543928375e-9}};
  for (const Case& t : cases) {
    double r[16], lo[4], hi[4];
    for (int i = 0; i < 16; ++i) r[i] = i % 5 == 0 ? 1 : t.rho;
    for (int j = 0; j < 4; ++j) {
      lo[j] = -stats::INF;
      hi[j] = t.c;
    }
    CHECK_CLOSE(nph::mvn_outside(4, lo, hi, r), t.want, 5e-5, 0);
    for (int j = 0; j < 4; ++j) lo[j] = -t.c;  // both tails: twice as likely (an overlap of opposite tails is nil)
    CHECK_CLOSE(nph::mvn_outside(4, lo, hi, r), 2 * t.want, 5e-5, 0);
  }
}

// ------------------------------------------------------------------ approximation and simulated power

TEST_CASE("nph: the quick approximation for proportional hazards is the logrank formula") {
  // HR1 = HR2 (or T* = 0) and the logrank test: Schoenfeld's formula with the expected events, as in
  // surv.logrank, with accrual, follow-up and loss to follow-up
  struct Case {
    const char* both;  // inputs common to the two procedures
    const char* hr;
    const char* nph_alloc;
    const char* lr_alloc;
  };
  const Case cases[] = {
      {R"("n1":"200","med2":"12","accrual":"12","follow":"12","eta":"0")", "0.75", R"("r":"1")",
       R"("alloc":"equal")"},
      {R"("n1":"80","med2":"8","accrual":"0","follow":"18","eta":"0.02")", "0.6", R"("r":"1")", R"("alloc":"equal")"},
      {R"("n1":"150","med2":"20","accrual":"18","follow":"6","eta":"0.01","alternative":"upper")", "1.3",
       R"("r":"2")", R"("alloc":"ratio","r":"2")"}};
  for (const Case& c : cases) {
    const std::string hr = c.hr;
    const json::Value a = first_row(nph_request(
        "power", std::string(c.both) + "," + c.nph_alloc + R"(,"hr1":")" + hr + R"(","hr2":")" + hr + R"(")",
        R"("approximate":true)"));
    const json::Value b = first_row(
        request("surv.logrank", "power", std::string(c.both) + "," + c.lr_alloc + R"(,"hr":")" + hr + R"(")"));
    CHECK(a.get("method").to_text() == "approximation");
    CHECK_CLOSE(num(a, "power"), num(b, "power"), 0, 1e-8);
    CHECK_CLOSE(num(a, "events"), num(b, "events"), 1e-8, 0);
    CHECK_CLOSE(num(a, "ahr"), std::stod(hr), 1e-12, 0);
  }
  // the sample size, with T* = 0 (proportional hazards with HR2)
  const json::Value n = first_row(nph_request(
      "n", R"("power":"0.9","hr1":"3","hr2":"0.7","tchange":"0","med2":"10","accrual":"12","follow":"12")",
      R"("approximate":true)"));
  const json::Value m = first_row(
      request("surv.logrank", "n", R"("power":"0.9","hr":"0.7","med2":"10","accrual":"12","follow":"12")"));
  CHECK_EQ(num(n, "n1"), num(m, "n1"));
  CHECK_CLOSE(num(n, "power"), num(m, "power"), 0, 1e-8);
}

TEST_CASE("nph: simulated power agrees with the logrank formula under proportional hazards") {
  // about 290 events: Schoenfeld's approximation is accurate, so the simulation must agree closely
  const std::string design = R"("n1":"250","hr1":"0.75","hr2":"0.75","med2":"12","accrual":"12","follow":"12")";
  const json::Value s = first_row(nph_request("power", design + R"(,"reps":"20000","seed":"11")"));
  const json::Value a = first_row(
      request("surv.logrank", "power", R"("n1":"250","hr":"0.75","med2":"12","accrual":"12","follow":"12")"));
  CHECK_CLOSE(num(s, "power"), num(a, "power"), 0, 0.025);
  CHECK_CLOSE(num(s, "events"), num(a, "events"), 0.005, 0);  // the average number of events
  CHECK_CLOSE(num(s, "approx_power"), num(a, "power"), 0, 1e-8);
  CHECK(num(s, "power_lo") < num(s, "power") && num(s, "power") < num(s, "power_hi"));
}

TEST_CASE("nph: the type I error is close to alpha for every test") {
  const std::string null = R"("n1":"200","hr1":"1","hr2":"1","tau":"12","reps":"20000","seed":"5")";
  for (const char* test : {R"("test":"logrank")", R"("test":"fh")", R"("test":"maxcombo")", R"("test":"rmst")",
                           R"("test":"maxcombo","alternative":"lower")", R"("test":"fh","alternative":"upper")"}) {
    const json::Value r = first_row(nph_request("power", null + "," + test));
    CHECK_CLOSE(num(r, "power"), 0.05, 0, 0.0062);  // 4 Monte Carlo standard errors
  }
}

TEST_CASE("nph: the MaxCombo approximation agrees with the simulation") {
  // the four statistics jointly normal with the null correlation of the expected events: the critical
  // value keeps the type I error at alpha exactly, and power is close to the simulated power
  for (const char* alt : {"two", "lower"}) {
    const std::string base = std::string(R"("test":"maxcombo","n1":"150","alternative":")") + alt + "\",";
    const json::Value null = first_row(nph_request("power", base + R"("hr1":"1","hr2":"1")", R"("approximate":true)"));
    CHECK_CLOSE(num(null, "power"), 0.05, 0, 2e-4);
    const std::string delayed = base + R"("hr1":"1","hr2":"0.6","tchange":"3","reps":"4000","seed":"77")";
    const json::Value approx = first_row(nph_request("power", delayed, R"("approximate":true)"));
    const json::Value sim = first_row(nph_request("power", delayed));
    CHECK_CLOSE(num(approx, "power"), num(sim, "power"), 0, 0.03);
  }
}

TEST_CASE("nph: weighted tests beat the logrank test under a delayed effect") {
  // no difference for 6 months, then HR 0.6 (control median 12, accrual 12, follow-up 12, 150 per group):
  // FH(0, 1) and MaxCombo weigh the late differences, and the RMST up to τ = 12 is dominated by the early
  // months without an effect (powers about 0.59, 0.51, 0.39 and 0.08)
  const std::string design =
      R"("n1":"150","hr1":"1","hr2":"0.6","tchange":"6","med2":"12","accrual":"12","follow":"12","tau":"12",)"
      R"("reps":"2000","seed":"3")";
  double power[4];
  const char* tests[4] = {"logrank", "fh", "maxcombo", "rmst"};
  for (int i = 0; i < 4; ++i)
    power[i] = num(first_row(nph_request("power", design + R"(,"test":")" + tests[i] + R"(")")), "power");
  CHECK(power[1] > power[0] + 0.1);
  CHECK(power[2] > power[0] + 0.05);
  CHECK(power[0] > power[3] + 0.1);
}

TEST_CASE("nph: results depend on the seed only, not on the number of threads") {
  for (const char* test : {"maxcombo", "rmst"}) {
    const std::string body = R"("n1":"60","test":")" + std::string(test) + R"(","reps":"400","seed":"77")";
    const json::Value one = first_row(nph_request("power", body, R"("threads":1)"));
    const json::Value four = first_row(nph_request("power", body, R"("threads":4)"));
    CHECK_EQ(num(one, "power"), num(four, "power"));
    CHECK_EQ(num(one, "events"), num(four, "events"));
    if (std::string(test) == "rmst") CHECK_EQ(num(one, "mean_drmst"), num(four, "mean_drmst"));
  }
}

TEST_CASE("nph: the sample-size search returns the smallest size reaching the target") {
  // a small and a large design (approximation for the second: N1 = 778)
  for (const char* design : {R"("test":"fh","hr2":"0.65","reps":"1000","seed":"9")",
                             R"("test":"logrank","hr2":"0.78","reps":"1000","seed":"4")"}) {
    auto power_at = [&](double n1) {
      return first_row(nph_request("power", R"("n1":")" + text::format_g(n1, 10) + R"(",)" + design));
    };
    const json::Value r = first_row(nph_request("n", std::string(R"("power":"0.8",)") + design));
    const double n1 = num(r, "n1");
    CHECK(num(r, "power") >= 0.8);
    // with the same seed, one subject fewer falls short (common random numbers)
    CHECK(num(power_at(n1 - 1), "power") < 0.8);
    CHECK_EQ(num(power_at(n1), "power"), num(r, "power"));
  }
}

TEST_CASE("nph: designs without an effect to detect give clear messages") {
  auto error_of = [](const std::string& req) {
    const json::Value res = run(req);
    const json::Value* e = res.get("rows").elements().at(0).find("_error");
    return e ? e->as_string() : std::string();
  };
  CHECK(error_of(nph_request("n", R"("hr1":"1","hr2":"1","reps":"200")")).find("no effect") != std::string::npos);
  CHECK(error_of(nph_request("n", R"("hr1":"1","tchange":"30","reps":"200")")).find("no effect") != std::string::npos);
  CHECK(error_of(nph_request("n", R"("alternative":"upper","reps":"200")")).find("opposite") != std::string::npos);
  CHECK(error_of(nph_request("power", R"("test":"rmst","tau":"30","reps":"200")")).find("τ") != std::string::npos);
  // the statement describes the effect and the test
  const json::Value ex = first_row(nph_request("power", R"("test":"maxcombo","reps":"200")"));
  const std::string st = ex.get("_statement").as_string();
  CHECK(st.find("delayed treatment effect") != std::string::npos && st.find("MaxCombo") != std::string::npos);
}
