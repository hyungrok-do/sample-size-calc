#include "dag/model.hpp"
#include "dag/statistics.hpp"
#define SSC_STATIC  // Tests link the core directly; do not import a Windows DLL.
#include "ssc/ssc.h"
#include "testing.hpp"
#include <fstream>
#include <sstream>

using namespace ssc::dag;
namespace {
J reference() {
  std::ifstream input(std::string(SSC_FIXTURE_DIR) + "/dag_reference.json");
  if (!input)
    throw std::runtime_error("DAG reference fixture missing");
  std::stringstream text;
  text << input.rdbuf();
  return ssc::json::parse(text.str());
}
J target(const J &item) {
  J e = item.get("estimand");
  e.set("kind", "marginal");
  e.set("measure", "mean_difference");
  e.set("alternative", .5);
  return e;
}
} // namespace
TEST_CASE(
    "dag: C++ estimators agree with independent NumPy/SciPy reference fits") {
  const J fixtures = reference();
  for (const auto &item : fixtures.get("cases").elements()) {
    const auto result = run_request(object(
        {{"op", "analyze"},
         {"graph", item.get("graph")},
         {"estimand", item.get("estimand")},
         {"data",
          object({{"n", item.get("n")}, {"columns", item.get("columns")}})}}));
    const auto &expected = item.get("reference");
    CHECK_CLOSE(result.get("estimate").as_number(),
                expected.get("coefficient").as_number(), 1e-6, 1e-7);
    CHECK_CLOSE(result.get("se").as_number(), expected.get("se").as_number(),
                1e-6, 1e-8);
    CHECK_CLOSE(result.get("pValue").as_number(), expected.get("p").as_number(),
                1e-5, 1e-9);
  }
}
TEST_CASE("dag: KM survival, RMST, Greenwood variance and tail support") {
  const auto k = kaplan_meier({1, 2, 3, 4}, {1, 0, 1, 0}, 4);
  CHECK_CLOSE(k.survival, .375, 0, 1e-14);
  CHECK_CLOSE(k.survivalVariance, .08203125, 0, 1e-14);
  CHECK_CLOSE(k.rmst, 2.875, 0, 1e-14);
  CHECK_CLOSE(k.rmstVariance, .36328125, 0, 1e-14);
  bool failed = false;
  try {
    kaplan_meier({1, 2, 3, 4}, {1, 0, 1, 0}, 5);
  } catch (const Error &e) {
    failed = e.code == "insufficient_followup";
  }
  CHECK(failed);
  CHECK_CLOSE(kaplan_meier({1, 2, 3, 4}, {1, 1, 1, 1}, 3).rmst, 2.25, 0, 1e-14);
  CHECK_CLOSE(p_value(2.262157162798205, 9, "two-sided"), .05, 1e-11, 1e-14);
  CHECK_CLOSE(p_value(-2.262157162798205, 9, "less"), .025, 1e-11, 1e-14);
  CHECK(p_value(10, INFINITY, "greater") > 0);
  CHECK(p_value(10, INFINITY, "greater") < 1e-20);
}
TEST_CASE("dag: singular and constant regression terms fail explicitly") {
  Data d;
  d.n = 6;
  d.values = {{1, 2, 3, 4, 5, 6}, {1, 2, 3, 4, 5, 6}};
  bool singular = false;
  try {
    fit_linear(model_matrix(d, {{0, -1, 0}, {1, -1, 0}}, true), d.values[0], 1);
  } catch (const Error &e) {
    singular = e.code == "singular_model";
  }
  CHECK(singular);
  d.values[0] = {1, 1, 1, 1, 1, 1};
  bool constant = false;
  try {
    model_matrix(d, {{0, -1, 0}}, true);
  } catch (const Error &e) {
    constant = e.code == "constant_predictor";
  }
  CHECK(constant);
}
TEST_CASE("dag: power is reproducible across the calculator thread pool") {
  const J fixture = reference().get("cases").at(0);
  const J e = target(fixture);
  const J config = object({{"sampleSizes", array({80, 300})},
                           {"repetitions", 150},
                           {"populationSize", 1000}});
  J req = object({{"op", "power"},
                  {"graph", fixture.get("graph")},
                  {"estimand", e},
                  {"settings", config},
                  {"threads", 1}});
  const J single = run_request(req);
  req.set("threads", 4);
  const J multi = run_request(req);
  CHECK_EQ(single.get("results").dump(), multi.get("results").dump());
  CHECK_EQ(single.get("alternative").dump(), multi.get("alternative").dump());
  CHECK_EQ(single.get("nullScenario").dump(), multi.get("nullScenario").dump());
  CHECK(multi.get("results")
            .at(1)
            .get("alternative")
            .get("probability")
            .as_number() > multi.get("results")
                               .at(0)
                               .get("alternative")
                               .get("probability")
                               .as_number());
  CHECK_CLOSE(
      single.get("alternative").get("truth").get("estimate").as_number(), .5, 0,
      .025);
}
TEST_CASE("dag: native request cancellation and invalid input are enforced") {
  const J item = reference().get("cases").at(0);
  const J req = object({{"op", "simulate"}, {"graph", item.get("graph")}});
  ssc::Token token;
  token.cancel();
  bool cancelled = false;
  try {
    run_request(req, &token);
  } catch (const ssc::Cancelled &) {
    cancelled = true;
  }
  CHECK(cancelled);
  char *text = ssc_dag_json("{\"op\":\"simulate\",\"graph\":null}", nullptr);
  CHECK(text != nullptr);
  if (text) {
    CHECK(!ssc::json::parse(text).get("ok").as_bool());
    ssc_free(text);
  }
  J bad = req;
  bad.set("options", object({{"n", 2.5}}));
  bool invalid = false;
  try {
    run_request(bad);
  } catch (const Error &) {
    invalid = true;
  }
  CHECK(invalid);
}
