#include "dag/model.hpp"
#include "dag/statistics.hpp"
#define SSC_STATIC  // Tests link the core directly; do not import a Windows DLL.
#include "ssc/ssc.h"
#include "testing.hpp"
#include <fstream>
#include <sstream>

using namespace ssc::dag;
namespace {
J reference(const std::string &name = "dag_reference.json") {
  std::ifstream input(std::string(SSC_FIXTURE_DIR) + "/" + name);
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
J variable(const std::string &id, const std::string &kind,
           const std::string &family, double intercept = 0) {
  return object({{"id", id},
                 {"label", id},
                 {"kind", kind},
                 {"family", family},
                 {"intercept", intercept},
                 {"sd", 1},
                 {"width", 2},
                 {"shape", 2},
                 {"scale", 10},
                 {"censoring", object({{"mode", "none"},
                                       {"time", 20},
                                       {"variable", ""}})},
                 {"interactions", J::array()},
                 {"x", 0},
                 {"y", 0}});
}
J arrow(const std::string &source, const std::string &target, double value) {
  return object(
      {{"source", source}, {"target", target}, {"coefficient", value}});
}
// X and W confound A -> Y; M mediates part of the effect of A.
J confounded_graph(double xa = 1) {
  J nodes = J::array(), edges = J::array();
  for (const auto &v :
       {variable("X", "continuous", "gaussian"),
        variable("W", "binary", "bernoulli"),
        variable("A", "binary", "bernoulli", -.2),
        variable("M", "continuous", "gaussian"),
        variable("Y", "continuous", "gaussian")})
    nodes.push(v);
  for (const auto &e : {arrow("X", "A", xa), arrow("W", "A", .7),
                        arrow("X", "Y", 1), arrow("W", "Y", .5),
                        arrow("A", "M", .5), arrow("M", "Y", .4),
                        arrow("A", "Y", .3)})
    edges.push(e);
  return object({{"version", 1}, {"nodes", nodes}, {"edges", edges}});
}
J marginal(const std::string &outcome = "Y") {
  return object({{"version", 1},
                 {"kind", "marginal"},
                 {"population", "Graph population"},
                 {"definition", ""},
                 {"exposure", "A"},
                 {"outcome", outcome},
                 {"measure", "mean_difference"},
                 {"reference", 0},
                 {"active", 1},
                 {"horizon", J()},
                 {"nullValue", 0},
                 {"alternative", J()}});
}
std::string refusal(const J &graph, const J &estimand) {
  try {
    run_request(
        object({{"op", "plan"}, {"graph", graph}, {"estimand", estimand}}));
  } catch (const Error &e) {
    return e.what();
  }
  return "";
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
TEST_CASE("dag: AIPW for a confounded exposure agrees with a NumPy reference") {
  const J fixtures = reference("dag_aipw_reference.json");
  for (const auto &item : fixtures.get("cases").elements()) {
    const auto result = run_request(object(
        {{"op", "analyze"},
         {"graph", item.get("graph")},
         {"estimand", item.get("estimand")},
         {"data",
          object({{"n", item.get("n")}, {"columns", item.get("columns")}})}}));
    const auto &expected = item.get("reference");
    CHECK_CLOSE(result.get("estimate").as_number(),
                expected.get("estimate").as_number(), 1e-6, 1e-9);
    CHECK_CLOSE(result.get("se").as_number(), expected.get("se").as_number(),
                1e-6, 1e-10);
    CHECK_CLOSE(result.get("pValue").as_number(), expected.get("p").as_number(),
                1e-5, 1e-12);
  }
  // The mediator M stays out of the adjustment set; the outcome's other
  // parent Z joins it for precision; the propensity model is A's own DGP.
  const J &item = fixtures.get("cases").at(0);
  const J plan = run_request(object({{"op", "plan"},
                                     {"graph", item.get("graph")},
                                     {"estimand", item.get("estimand")}}));
  CHECK_EQ(plan.get("method").as_string(), std::string("aipw"));
  CHECK_EQ(plan.get("terms").dump(),
           std::string(R"([{"left":"X"},{"left":"W"},{"left":"Z"},)"
                       R"({"left":"X","right":"W"}])"));
  CHECK_EQ(plan.get("propensityTerms").dump(),
           std::string(R"([{"left":"X"},{"left":"W"},)"
                       R"({"left":"X","right":"W"}])"));
}
TEST_CASE("dag: AIPW power under confounding is unbiased and holds its size") {
  const J graph = confounded_graph(), e = marginal();
  // The total effect is 0.3 + 0.5 * 0.4 = 0.5; the crude difference in means
  // is far larger because X and W raise both A and Y.
  const J truth = run_request(object({{"op", "population"},
                                      {"graph", graph},
                                      {"estimand", e},
                                      {"options", object({{"n", 2000}})}}));
  CHECK_CLOSE(truth.get("estimate").as_number(), .5, 0, 1e-12);
  J req = object({{"op", "power"},
                  {"graph", graph},
                  {"estimand", e},
                  {"settings", object({{"sampleSizes", array({150, 600})},
                                       {"repetitions", 800},
                                       {"populationSize", 2000}})},
                  {"threads", 1}});
  const J single = run_request(req);
  req.set("threads", 4);
  const J multi = run_request(req);
  CHECK_EQ(single.get("results").dump(), multi.get("results").dump());
  CHECK_EQ(multi.get("analysis").get("method").as_string(),
           std::string("aipw"));
  double previous = 0;
  for (const auto &row : multi.get("results").elements()) {
    const auto &alt = row.get("alternative"), &null = row.get("null");
    CHECK_EQ(alt.get("failures").as_number(), 0.0);
    // Mean estimate within 4 Monte Carlo SEs of the true total effect.
    CHECK_CLOSE(alt.get("meanEstimate").as_number(), .5, 0,
                4 * alt.get("estimateSD").as_number() / std::sqrt(800.0));
    CHECK_CLOSE(null.get("meanEstimate").as_number(), 0, 0,
                4 * null.get("estimateSD").as_number() / std::sqrt(800.0));
    // The model-based SE tracks the spread of the estimates.
    CHECK_CLOSE(alt.get("meanSE").as_number(),
                alt.get("estimateSD").as_number(), .12, 0);
    // Type I error near alpha: 0.05 +/- 4 binomial SEs at 800 studies.
    CHECK_CLOSE(null.get("probability").as_number(), .05, 0, .031);
    CHECK(alt.get("probability").as_number() > previous);
    previous = alt.get("probability").as_number();
  }
  CHECK(previous > .95);
}
TEST_CASE("dag: AIPW refuses targets it cannot identify or estimate") {
  J graph = confounded_graph();
  J nodes = graph.get("nodes");
  J survival = variable("T", "survival", "weibull");
  nodes.push(survival);
  J edges = graph.get("edges");
  edges.push(arrow("A", "T", .4));
  graph.set("nodes", nodes);
  graph.set("edges", edges);
  J rmst = marginal("T");
  rmst.set("measure", "rmst_difference");
  rmst.set("horizon", 5);
  CHECK(refusal(graph, rmst).find("censoring-adjusted") != std::string::npos);
  // A censored event time cannot be a confounder: its latent value is unseen.
  survival.set("censoring", object({{"mode", "administrative"},
                                    {"time", 5},
                                    {"variable", ""}}));
  J censored = confounded_graph();
  J cnodes = censored.get("nodes"), cedges = censored.get("edges");
  cnodes.push(survival);
  cedges.push(arrow("T", "A", .1));
  censored.set("nodes", cnodes);
  censored.set("edges", cedges);
  CHECK(refusal(censored, marginal()).find("censored event time") !=
        std::string::npos);
  J continuous = marginal();
  continuous.set("exposure", "X");
  continuous.set("active", 1);
  CHECK(refusal(confounded_graph(), continuous).find("binary comparison") !=
        std::string::npos);
  // A study whose propensity fit separates the arms fails, never a clipped
  // weight: it counts as a non-rejection.
  const J separated = confounded_graph(40);
  const J result = run_request(
      object({{"op", "power"},
              {"graph", separated},
              {"estimand", marginal()},
              {"settings", object({{"sampleSizes", array({60})},
                                   {"repetitions", 20},
                                   {"populationSize", 1000},
                                   {"checkNull", false}})}}));
  CHECK(result.get("results")
            .at(0)
            .get("alternative")
            .get("failures")
            .as_number() > 0);
}
