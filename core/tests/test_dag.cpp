#include "dag/model.hpp"
#include "dag/statistics.hpp"
#include "dag/learners.hpp"
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
TEST_CASE("dag: targets AIPW does not cover run with the g-formula or a note") {
  const auto plan = [](const J &graph, const J &estimand) {
    return run_request(
        object({{"op", "plan"}, {"graph", graph}, {"estimand", estimand}}));
  };
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
  CHECK_EQ(plan(graph, rmst).get("method").as_string(), std::string("gformula"));
  // A censored confounder is adjusted for by its latent value, with a note.
  survival.set("censoring", object({{"mode", "administrative"},
                                    {"time", 5},
                                    {"variable", ""}}));
  J censored = confounded_graph();
  J cnodes = censored.get("nodes"), cedges = censored.get("edges");
  cnodes.push(survival);
  cedges.push(arrow("T", "A", .1));
  censored.set("nodes", cnodes);
  censored.set("edges", cedges);
  const J noted = plan(censored, marginal());
  CHECK_EQ(noted.get("method").as_string(), std::string("aipw"));
  CHECK(noted.get("notes").at(0).as_string().find("latent values") !=
        std::string::npos);
  J continuous = marginal();
  continuous.set("exposure", "X");
  continuous.set("active", 1);
  const J refit = plan(confounded_graph(), continuous);
  CHECK_EQ(refit.get("method").as_string(), std::string("gformula"));
  CHECK_EQ(refit.get("formula").as_string(),
           std::string("Y: X = 1 versus 0; refitted A ~ X + W, M ~ A, "
                       "Y ~ X + W + M + A; other variables as observed"));
  CHECK_EQ(refusal(confounded_graph(), continuous), std::string());
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
namespace {
// A confounded continuous exposure with a mediator; the total effect of one
// unit is 0.3 + 0.8 * 0.6 = 0.78.
J continuous_graph() {
  J nodes = J::array(), edges = J::array();
  for (const auto &v : {variable("X", "continuous", "gaussian"),
                        variable("A", "continuous", "gaussian"),
                        variable("M", "continuous", "gaussian"),
                        variable("Y", "continuous", "gaussian")})
    nodes.push(v);
  for (const auto &e : {arrow("X", "A", .7), arrow("X", "M", .4),
                        arrow("X", "Y", .4), arrow("A", "M", .8),
                        arrow("M", "Y", .6), arrow("A", "Y", .3)})
    edges.push(e);
  return object({{"version", 1}, {"nodes", nodes}, {"edges", edges}});
}
// A confounded binary exposure and an event time whose censoring shares the
// confounder: Kaplan-Meier and AIPW identify neither.
J confounded_survival() {
  J t = variable("T", "survival", "weibull"), c = variable("C", "survival",
                                                           "exponential");
  t.set("shape", 1.5);
  t.set("censoring",
        object({{"mode", "variable"}, {"time", 20}, {"variable", "C"}}));
  c.set("scale", 18);
  J nodes = J::array(), edges = J::array();
  for (const auto &v : {variable("X", "continuous", "gaussian"),
                        variable("A", "binary", "bernoulli"), c, t})
    nodes.push(v);
  for (const auto &e : {arrow("X", "A", .8), arrow("X", "T", .5),
                        arrow("A", "T", -.5), arrow("X", "C", .6)})
    edges.push(e);
  return object({{"version", 1}, {"nodes", nodes}, {"edges", edges}});
}
J power_of(const J &graph, const J &estimand, int n, int repetitions) {
  return run_request(object(
      {{"op", "power"},
       {"graph", graph},
       {"estimand", estimand},
       {"settings", object({{"sampleSizes", array({static_cast<double>(n)})},
                            {"repetitions", repetitions},
                            {"populationSize", 20000}})},
       {"threads", 4}}));
}
// Mean estimate near the truth, model SE near the estimates' spread, and the
// null rejected at about alpha.
void check_calibrated(const J &result, double scale = 1) {
  const auto &row = result.get("results").at(0);
  const auto &alt = row.get("alternative"), &null = row.get("null");
  const double reps = alt.get("repetitions").as_number(),
               truth =
                   result.get("alternative").get("truth").get("estimate").as_number();
  CHECK_EQ(alt.get("failures").as_number(), 0.0);
  CHECK_CLOSE(alt.get("meanEstimate").as_number(), truth, 0,
              4.5 * alt.get("estimateSD").as_number() / std::sqrt(reps));
  CHECK_CLOSE(null.get("meanEstimate").as_number(), scale == 1 ? 0 : 1, 0,
              4.5 * null.get("estimateSD").as_number() / std::sqrt(reps));
  CHECK_CLOSE(alt.get("meanSE").as_number() * scale,
              alt.get("estimateSD").as_number(), .15, 0);
  CHECK_CLOSE(null.get("probability").as_number(), .05, 0, .035);
}
} // namespace
TEST_CASE("dag: the g-formula recovers continuous, count and event-time targets") {
  J e = marginal();
  check_calibrated(power_of(continuous_graph(), e, 300, 500));
  J rmst = marginal("T");
  rmst.set("measure", "rmst_difference");
  rmst.set("horizon", 10);
  const J survival = power_of(confounded_survival(), rmst, 400, 400);
  CHECK_EQ(survival.get("analysis").get("method").as_string(),
           std::string("gformula"));
  check_calibrated(survival);
  J share = rmst;
  share.set("measure", "survival_difference");
  check_calibrated(power_of(confounded_survival(), share, 400, 400));
  // a count exposure compared at 2 versus 0, Poisson outcome
  J nodes = J::array(), edges = J::array();
  nodes.push(variable("X", "continuous", "gaussian"));
  nodes.push(variable("A", "count", "poisson", .3));
  nodes.push(variable("Y", "count", "poisson", .5));
  for (const auto &a : {arrow("X", "A", .3), arrow("X", "Y", .3),
                        arrow("A", "Y", .15)})
    edges.push(a);
  J count = marginal();
  count.set("active", 2);
  check_calibrated(power_of(object({{"version", 1}, {"nodes", nodes},
                                    {"edges", edges}}),
                            count, 400, 400));
}
TEST_CASE("dag: the g-formula extrapolates past follow-up and notes the "
          "assumptions it cannot check") {
  J t = variable("T", "survival", "weibull");
  t.set("shape", 1.3);
  t.set("censoring",
        object({{"mode", "administrative"}, {"time", 8}, {"variable", ""}}));
  J nodes = J::array(), edges = J::array();
  nodes.push(variable("A", "binary", "bernoulli"));
  nodes.push(t);
  edges.push(arrow("A", "T", -.6));
  const J graph = object({{"version", 1}, {"nodes", nodes}, {"edges", edges}});
  J e = marginal("T");
  e.set("measure", "survival_difference");
  e.set("horizon", 12);
  const J result = power_of(graph, e, 400, 400);
  CHECK_EQ(result.get("analysis").get("method").as_string(),
           std::string("gformula"));
  CHECK(result.get("warnings").at(0).as_string().find("extrapolates") !=
        std::string::npos);
  check_calibrated(result);
  // Censoring caused by the event time itself: the run goes ahead and says so.
  J informative = graph;
  J inodes = J::array(), iedges = J::array();
  J it = variable("T", "survival", "weibull");
  it.set("censoring",
         object({{"mode", "variable"}, {"time", 20}, {"variable", "C"}}));
  inodes.push(variable("A", "binary", "bernoulli"));
  inodes.push(it);
  inodes.push(variable("C", "continuous", "lognormal", 1.5));
  iedges.push(arrow("A", "T", -.5));
  iedges.push(arrow("T", "C", .05));
  informative.set("nodes", inodes);
  informative.set("edges", iedges);
  J rmst = marginal("T");
  rmst.set("measure", "rmst_difference");
  rmst.set("horizon", 8);
  const J plan = run_request(
      object({{"op", "plan"}, {"graph", informative}, {"estimand", rmst}}));
  CHECK(plan.get("notes").at(0).as_string().find("depends on T itself") !=
        std::string::npos);
  CHECK_EQ(plan.get("models").at(0).get("node").as_string(), std::string("T"));
}
TEST_CASE("dag: the censored Weibull fit recovers its shape and coefficient") {
  J t = variable("T", "survival", "weibull");
  t.set("shape", 1.7);
  t.set("scale", 5);
  t.set("censoring",
        object({{"mode", "administrative"}, {"time", 6}, {"variable", ""}}));
  J nodes = J::array(), edges = J::array();
  nodes.push(variable("A", "binary", "bernoulli"));
  nodes.push(t);
  edges.push(arrow("A", "T", -.7));
  const auto g =
      parse_graph(object({{"version", 1}, {"nodes", nodes}, {"edges", edges}}));
  const auto d = simulate(g, 40000, 11);
  const auto m = model_matrix(d, g.nodes[1].terms, true);
  const auto fit = fit_weibull(m, d.times[1], d.events[1], false);
  CHECK_CLOSE(std::exp(fit.beta[2]), 1.7, 0, .03);
  CHECK_CLOSE(fit.beta[1] / m.scales[1], -.7, 0, .03);
  // H(t | A) = (t / 5)^1.7 exp(-0.7 A): the intercept at A = 0 is -1.7 log 5
  CHECK_CLOSE(fit.beta[0] - fit.beta[1] * m.means[1] / m.scales[1],
              -1.7 * std::log(5.0), 0, .06);
}
TEST_CASE("dag: a comparison without a causal path still runs") {
  const auto plan = [](const J &estimand) {
    return run_request(object({{"op", "plan"},
                               {"graph", confounded_graph()},
                               {"estimand", estimand}}));
  };
  // X causes M through A: M cannot change X, in any study.
  J reverse = marginal("X");
  reverse.set("exposure", "M");
  const J fixed = plan(reverse);
  CHECK_EQ(fixed.get("method").as_string(), std::string("gformula"));
  CHECK_EQ(fixed.get("models").size(), std::size_t(0));
  CHECK(fixed.get("notes").at(0).as_string().find("X lies upstream of M") !=
        std::string::npos);
  const auto row = power_of(confounded_graph(), reverse, 100, 50)
                       .get("results")
                       .at(0)
                       .get("alternative");
  CHECK_EQ(row.get("failures").as_number(), 0.0);
  CHECK_EQ(row.get("probability").as_number(), 0.0);
  CHECK_EQ(row.get("meanEstimate").as_number(), 0.0);
  // A confounded binary comparison whose parent is the outcome: not AIPW,
  // which would adjust for the outcome itself.
  J parent = marginal("W");
  parent.set("measure", "risk_difference");
  CHECK_EQ(plan(parent).get("method").as_string(), std::string("gformula"));
  // No path between W and X: the test runs at its false-positive rate, and a
  // requested effect that no arrow can carry keeps the drawn model.
  J apart = marginal("X");
  apart.set("exposure", "W");
  apart.set("alternative", .3);
  CHECK(plan(apart).get("notes").at(0).as_string().find("No arrow path") !=
        std::string::npos);
  const J result = power_of(confounded_graph(), apart, 200, 400);
  bool kept = false;
  for (const auto &w : result.get("warnings").elements())
    kept |= w.as_string().find("planning effect 0.3 cannot be set") !=
            std::string::npos;
  CHECK(kept);
  CHECK_CLOSE(result.get("results").at(0).get("alternative").get("probability")
                  .as_number(),
              .05, 0, .035);
}
TEST_CASE("dag: Student t, Beta and negative binomial variables have their moments") {
  J nodes = J::array();
  J t = variable("T", "continuous", "student", 1);
  t.set("sd", 2); t.set("df", 5);
  J b = variable("B", "continuous", "beta", std::log(.3 / .7));
  b.set("precision", 4);
  J c = variable("C", "count", "negbin", std::log(3.0));
  c.set("overdispersion", .5);
  for (const auto &v : {t, b, c}) nodes.push(v);
  const J graph = object({{"version", 1}, {"nodes", nodes}, {"edges", J::array()}});
  const auto data = run_request(object({{"op", "simulate"}, {"graph", graph},
      {"options", object({{"n", 40000}, {"seed", 11}})}}));
  const auto moments_of = [&](const char *id) {
    Vec values;
    for (const auto &v : data.get("columns").get(id).elements()) values.push_back(v.as_number());
    return moments(values);
  };
  // t: mean 1, variance sd^2 df / (df - 2) = 4 * 5 / 3
  const auto mt = moments_of("T"), mb = moments_of("B"), mc = moments_of("C");
  CHECK_CLOSE(mt.mean, 1, 0, .04);
  CHECK_CLOSE(mt.variance, 20.0 / 3, .06, 0);
  // Beta: mean .3, variance .3 * .7 / (1 + 4)
  CHECK_CLOSE(mb.mean, .3, 0, .005);
  CHECK_CLOSE(mb.variance, .21 / 5, .04, 0);
  // negative binomial: mean 3, variance 3 + .5 * 9
  CHECK_CLOSE(mc.mean, 3, 0, .04);
  CHECK_CLOSE(mc.variance, 7.5, .05, 0);
  for (const auto &v : data.get("columns").get("B").elements())
    CHECK(v.as_number() > 0 && v.as_number() < 1);
  // a Student t variable needs a mean
  J bad = graph;
  J badNodes = J::array();
  J t1 = t; t1.set("df", 1);
  badNodes.push(t1);
  bad.set("nodes", badNodes);
  bool refused = false;
  try { parse_graph(bad); } catch (const Error &e) { refused = std::string(e.what()).find("more than 1 degree") != std::string::npos; }
  CHECK(refused);
}
TEST_CASE("dag: unmeasured variables are simulated but never analyzed") {
  // U confounds A -> Y and is unmeasured; X is measured.
  J nodes = J::array(), edges = J::array();
  J u = variable("U", "continuous", "gaussian");
  u.set("measured", false);
  for (const auto &v : {u, variable("X", "continuous", "gaussian"),
                        variable("A", "binary", "bernoulli"), variable("Y", "continuous", "gaussian")})
    nodes.push(v);
  for (const auto &e : {arrow("U", "A", .8), arrow("U", "Y", .5), arrow("X", "A", .5),
                        arrow("X", "Y", .4), arrow("A", "Y", .3)})
    edges.push(e);
  const J graph = object({{"version", 1}, {"nodes", nodes}, {"edges", edges}});
  const J plan = run_request(object({{"op", "plan"}, {"graph", graph}, {"estimand", marginal()}}));
  CHECK_EQ(plan.get("method").as_string(), std::string("aipw"));
  for (const char *key : {"terms", "propensityTerms"})
    for (const auto &t : plan.get(key).elements())
      CHECK(t.get("left").as_string() != "U");
  CHECK(plan.get("notes").at(0).as_string().find("U is unmeasured") != std::string::npos);
  // the data still carry U, and the analysis is biased by it
  const auto result = power_of(graph, marginal(), 400, 200);
  const auto row = result.get("results").at(0).get("alternative");
  CHECK(row.get("bias").as_number() > .2);
  // a model term on an unmeasured parent becomes its measured causes
  const auto g = parse_graph(graph);
  const auto terms = measured_terms(g, g.nodes[g.index("A")].terms);
  CHECK_EQ(terms.size(), std::size_t(1));
  CHECK_EQ(g.nodes[terms[0].left].id, std::string("X"));
  // the compared variable and the outcome must be measured
  J hidden = marginal();
  hidden.set("exposure", "U");
  CHECK(refusal(graph, hidden).find("must be measured") != std::string::npos);
}
TEST_CASE("dag: log terms, quasi-likelihood tests and the interaction test") {
  // Y ~ Poisson-gamma counts over follow-up F with log(F) as an offset-like term.
  J nodes = J::array(), edges = J::array();
  J f = variable("F", "continuous", "uniform", 1);
  f.set("width", 1);
  J y = variable("Y", "count", "negbin", std::log(1.5));
  y.set("overdispersion", .6);
  J functions = J::array();
  functions.push(object({{"left", "F"}, {"transform", object({{"kind", "log"}})}, {"coefficient", 1}}));
  y.set("functions", functions);
  J terms = J::array();
  terms.push(object({{"left", "A"}, {"right", "F"}, {"coefficient", 0}}));
  y.set("interactions", terms);
  for (const auto &v : {variable("A", "binary", "bernoulli"), f, y}) nodes.push(v);
  for (const auto &e : {arrow("A", "Y", std::log(.7)), arrow("F", "Y", 0)}) edges.push(e);
  const J graph = object({{"version", 1}, {"nodes", nodes}, {"edges", edges}});
  J coefficient = marginal();
  coefficient.set("kind", "coefficient");
  coefficient.set("measure", "exponentiated_coefficient");
  coefficient.set("nullValue", 1);
  const J plan = run_request(object({{"op", "plan"}, {"graph", graph}, {"estimand", coefficient}}));
  CHECK(plan.get("label").as_string().find("sandwich") != std::string::npos);
  CHECK(plan.get("formula").as_string().find("log(F)") != std::string::npos);
  const J result = power_of(graph, coefficient, 400, 200);
  const auto row = result.get("results").at(0);
  // with A x F in the model, A's coefficient is the rate ratio at F = 0: its
  // intervals still cover, and the sandwich test holds its size
  CHECK_CLOSE(row.get("alternative").get("coverage").as_number(), .95, 0, .04);
  CHECK_CLOSE(row.get("null").get("probability").as_number(), .05, 0, .035);
  // the zero interaction is the true null of the interaction test
  J interaction = coefficient;
  interaction.set("kind", "joint");
  interaction.set("measure", "interaction");
  interaction.set("nullValue", 0);
  const J joint = run_request(object({{"op", "plan"}, {"graph", graph}, {"estimand", interaction}}));
  CHECK_EQ(joint.get("label").as_string(), std::string("Wald test of the interaction (1 df)"));
  CHECK_CLOSE(power_of(graph, interaction, 400, 300).get("results").at(0).get("alternative")
                  .get("probability").as_number(), .05, 0, .035);
  // a log term of a variable that reaches zero stops the simulation
  J zero = graph, zn = J::array();
  for (J v : graph.get("nodes").elements()) {
    if (v.get("id").as_string() == "F") { v.set("intercept", 0); }
    zn.push(v);
  }
  zero.set("nodes", zn);
  bool stopped = false;
  try { run_request(object({{"op", "simulate"}, {"graph", zero}, {"options", object({{"n", 200}})}})); }
  catch (const Error &e) { stopped = std::string(e.what()).find("log(x)") != std::string::npos; }
  CHECK(stopped);
}
TEST_CASE("dag: relationships of unknown form and the Super Learner") {
  // Y is almost exactly a smooth unknown function of X1 and X2, of SD 0.8.
  const auto unknown = [](const J &inputs, double size, double smoothness, int seed) {
    return object({{"inputs", inputs}, {"size", size}, {"smoothness", smoothness}, {"seed", seed}});
  };
  J inputs = J::array(); inputs.push("X1"); inputs.push("X2");
  J nodes = J::array(), edges = J::array();
  J y = variable("Y", "continuous", "gaussian", 2);
  y.set("sd", .01); y.set("unknown", unknown(inputs, .8, 1.5, 4));
  for (const auto &v : {variable("X1", "continuous", "gaussian"), variable("X2", "continuous", "gaussian"), y}) nodes.push(v);
  for (const auto &e : {arrow("X1", "Y", 0), arrow("X2", "Y", 0)}) edges.push(e);
  const J graph = object({{"version", 1}, {"nodes", nodes}, {"edges", edges}});
  const auto data = run_request(object({{"op", "simulate"}, {"graph", graph}, {"options", object({{"n", 20000}, {"seed", 9}})}}));
  Vec values;
  for (const auto &v : data.get("columns").get("Y").elements()) values.push_back(v.as_number());
  const auto m = moments(values);
  CHECK_CLOSE(m.mean, 2, 0, .05);                  // centred on the intercept
  CHECK_CLOSE(std::sqrt(m.variance), .8, .08, 0);  // its SD over the population is the size
  // the same seed draws the same function; another seed another one
  const auto again = run_request(object({{"op", "simulate"}, {"graph", graph}, {"options", object({{"n", 50}, {"seed", 9}})}}));
  CHECK_EQ(again.get("columns").get("Y").at(3).as_number(), data.get("columns").get("Y").at(3).as_number());
  J other = graph, otherNodes = J::array();
  for (J v : graph.get("nodes").elements()) { if (v.get("id").as_string() == "Y") v.set("unknown", unknown(inputs, .8, 1.5, 5)); otherNodes.push(v); }
  other.set("nodes", otherNodes);
  const auto changed = run_request(object({{"op", "simulate"}, {"graph", other}, {"options", object({{"n", 50}, {"seed", 9}})}}));
  CHECK(changed.get("columns").get("Y").at(3).as_number() != data.get("columns").get("Y").at(3).as_number());
  // The Super Learner fits the shape it was not told: its predictions beat a
  // main-terms regression, and its weights form a convex combination.
  const auto g = parse_graph(graph);
  const auto d = simulate(g, 600, 21);
  std::vector<int> rows(600);
  for (int i = 0; i < 600; ++i) rows[i] = i;
  const auto learned = super_learner(d, {0, 1}, d.values[2], rows, "gaussian");
  const auto linear = fitted_means(model_matrix(d, {Term{0}, Term{1}}, true), d.values[2], rows, "gaussian");
  double learnedError = 0, linearError = 0, total = 0;
  for (int i = 0; i < 600; ++i) {
    learnedError += (learned.prediction[i] - d.values[2][i]) * (learned.prediction[i] - d.values[2][i]);
    linearError += (linear[i] - d.values[2][i]) * (linear[i] - d.values[2][i]);
  }
  for (double w : learned.weights) { CHECK(w >= 0); total += w; }
  CHECK_CLOSE(total, 1, 0, 1e-9);
  CHECK(learnedError < .5 * linearError);
  // The weights do not depend on where the outcome sits: shifted by 1000 (an
  // outcome far from zero relative to its spread), the predictions shift and
  // the weights stay. Normalizing non-negative least squares on the raw
  // predictions instead put nearly all weight on the training mean.
  Vec shifted(d.values[2]);
  for (double &v : shifted) v += 1000;
  const auto moved = super_learner(d, {0, 1}, shifted, rows, "gaussian");
  for (std::size_t l = 0; l < LEARNERS.size(); ++l) CHECK_CLOSE(moved.weights[l], learned.weights[l], 0, 1e-6);
  for (int i = 0; i < 600; i += 60) CHECK_CLOSE(moved.prediction[i] - 1000, learned.prediction[i], 0, 1e-6);
  CHECK(learned.weights[0] < .2);  // the training mean is not the answer
  // Confounding of unknown form: the automatic analysis is TMLE with a Super
  // Learner and cross-fitting, and recovers the effect.
  J a = variable("A", "binary", "bernoulli");
  a.set("unknown", unknown(inputs, 1, 2, 3));
  J y2 = variable("Y", "continuous", "gaussian");
  y2.set("unknown", unknown(inputs, 1, 2, 5));
  J nodes2 = J::array(), edges2 = J::array();
  for (const auto &v : {variable("X1", "continuous", "gaussian"), variable("X2", "continuous", "gaussian"), a, y2}) nodes2.push(v);
  for (const auto &e : {arrow("X1", "A", .3), arrow("X2", "A", 0), arrow("X1", "Y", .3), arrow("X2", "Y", .3), arrow("A", "Y", .3)}) edges2.push(e);
  const J confounded = object({{"version", 1}, {"nodes", nodes2}, {"edges", edges2}});
  const J plan = run_request(object({{"op", "plan"}, {"graph", confounded}, {"estimand", marginal()}}));
  CHECK_EQ(plan.get("method").as_string(), std::string("tmle"));
  CHECK_EQ(plan.get("learner").as_string(), std::string("superlearner"));
  CHECK_EQ(plan.get("folds").as_number(), 3.0);
  // a Super Learner asked for without folds is cross-fitted as well
  const J asked = run_request(object({{"op", "plan"}, {"graph", confounded}, {"estimand", marginal()},
                                      {"analysis", object({{"method", "tmle"}, {"learner", "superlearner"}})}}));
  CHECK_EQ(asked.get("folds").as_number(), 3.0);
  bool told = false;
  for (const auto &n : asked.get("notes").elements()) told |= n.as_string().find("cross-fitted in 3 folds") != std::string::npos;
  CHECK(told);
  // In a study of large clusters, a cluster-level input is standardized from
  // its own spread (SD 3), not from the few clusters that fit 4,000 rows.
  J w = variable("W", "continuous", "gaussian"), ca = variable("A", "binary", "bernoulli"), cy = variable("Y", "continuous", "gaussian");
  w.set("sd", 3); w.set("scope", "cluster"); ca.set("scope", "cluster"); cy.set("scope", "subject");
  J wOnly = J::array(); wOnly.push("W");
  cy.set("unknown", unknown(wOnly, 1, 1, 7));
  J studyNodes = J::array(), studyEdges = J::array(), times = J::array(), sizes = J::array();
  for (const auto &v : {w, ca, cy}) studyNodes.push(v);
  for (const auto &e : {arrow("W", "Y", 0), arrow("A", "Y", .5)}) studyEdges.push(e);
  for (int t = 0; t < 10; ++t) times.push(t);
  sizes.push(400);
  const J study = object({{"version", 1}, {"nodes", studyNodes}, {"edges", studyEdges},
                          {"study", object({{"times", times}, {"clusterSizes", sizes}})}});
  const auto standardized = parse_graph(study).nodes[2].unknown;
  CHECK_CLOSE(standardized.scale[0], 3, .1, 0);
  CHECK_CLOSE(standardized.center[0], 0, 0, .3);
  const auto row = power_of(confounded, marginal(), 400, 60).get("results").at(0).get("alternative");
  CHECK_CLOSE(row.get("bias").as_number(), 0, 0, .1);
  CHECK(row.get("targetedDiagnostics").get("learners").get("outcome").size() == LEARNERS.size());
  // no coefficient describes an effect of unknown form
  J heterogeneous = confounded, hn = J::array();
  J withA = inputs; withA.push("A");
  for (J v : confounded.get("nodes").elements()) { if (v.get("id").as_string() == "Y") v.set("unknown", unknown(withA, .5, 2, 5)); hn.push(v); }
  heterogeneous.set("nodes", hn);
  J coefficient = marginal();
  coefficient.set("kind", "coefficient"); coefficient.set("measure", "coefficient");
  CHECK(refusal(heterogeneous, coefficient).find("unknown form") != std::string::npos);
}
