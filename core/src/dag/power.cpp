#include "model.hpp"
#include "sim/random.hpp"
#include "sim/simulate.hpp"
#include "statistics.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <set>

namespace ssc::dag {
namespace {
bool positive(const std::string &measure) {
  return measure == "risk_ratio" || measure == "odds_ratio" ||
         measure == "exponentiated_coefficient";
}
std::uint64_t seed_value(const J &o, const std::string &key = "seed") {
  if (!o.find(key))
    return 2026;
  const double v = number(o, key);
  if (v < 0 || v > 4294967295.0 || v != std::floor(v))
    fail("Seed must be a whole number from 0 to 4294967295.");
  return static_cast<std::uint64_t>(v);
}
void validate_estimand(const Graph &g, const J &e) {
  if (e.get("version").as_number() != 1)
    fail("Choose an estimand type.");
  const auto kind = e.get("kind").as_string(),
             measure = e.get("measure").as_string();
  if (kind == "custom")
    fail("A custom text definition needs an executable estimator and test. "
         "Choose a supported target to run power.");
  if (kind != "marginal" && kind != "coefficient")
    fail("Choose an estimand type.");
  if (e.get("population").as_string().empty())
    fail("Describe the target population.");
  const int a = g.index(e.get("exposure").as_string()),
            y = g.index(e.get("outcome").as_string());
  if (a == y)
    fail("The outcome and comparison variable must differ.");
  const auto &outcome = g.nodes[y];
  const auto &exposure = g.nodes[a];
  const double lo = number(e, "reference"), hi = number(e, "active");
  bool compatible = false;
  if (kind == "coefficient") {
    bool arrow = false;
    for (const auto &term : outcome.terms)
      if (term.left == a && term.right < 0)
        arrow = true;
    if (!arrow)
      fail("A conditional coefficient requires an arrow from the selected "
           "predictor to the outcome.");
    compatible = measure == "coefficient" ||
                 (measure == "exponentiated_coefficient" &&
                  outcome.family != "gaussian" && outcome.family != "uniform");
  } else {
    if (lo == hi)
      fail("Enter two different comparison levels.");
    for (double level : {lo, hi}) {
      if (exposure.kind == "binary" && level != 0 && level != 1)
        fail("Binary comparison levels must be 0 or 1.");
      if (exposure.kind == "count" && (level < 0 || level != std::floor(level)))
        fail("Count comparison levels must be nonnegative integers.");
      if (exposure.kind == "survival" && level <= 0)
        fail("Survival-time intervention levels must be positive.");
    }
    compatible =
        (measure == "mean_difference" && outcome.kind != "survival") ||
        ((measure == "risk_difference" || measure == "risk_ratio" ||
          measure == "odds_ratio") &&
         outcome.kind == "binary") ||
        ((measure == "survival_difference" || measure == "rmst_difference") &&
         outcome.kind == "survival");
  }
  if (!compatible)
    fail("Choose an effect measure compatible with the outcome distribution.");
  double horizon = INFINITY;
  if (measure == "survival_difference" || measure == "rmst_difference") {
    horizon = number(e, "horizon");
    if (!(horizon > 0))
      fail("Enter a positive time horizon.");
  }
  for (const auto &key : {"nullValue", "alternative"}) {
    if (std::string(key) == "alternative" && e.get(key).is_null())
      continue;
    const double v = number(e, key);
    if (positive(measure) && !(v > 0))
      fail("Null and alternative values must be positive on a ratio scale.");
    if ((measure == "risk_difference" || measure == "survival_difference" ||
         (kind == "marginal" && measure == "mean_difference" &&
          outcome.kind == "binary")) &&
        std::fabs(v) > 1)
      fail("A probability difference must lie between -1 and 1.");
    if (measure == "rmst_difference" && std::fabs(v) > horizon)
      fail("An RMST difference cannot exceed its horizon in magnitude.");
  }
  if (!e.get("alternative").is_null() &&
      e.get("alternative").as_number() == e.get("nullValue").as_number())
    fail("The planning alternative must differ from the null value.");
}
struct Plan {
  std::string method, label, formula, testScale, family;
  int exposure = 0, outcome = 0, targetIndex = -1;
  std::vector<Term> terms;
  std::vector<int> edges;
};
Plan compile_plan(const Graph &g, const J &e) {
  validate_estimand(g, e);
  Plan p;
  p.exposure = g.index(e.get("exposure").as_string());
  p.outcome = g.index(e.get("outcome").as_string());
  const auto &a = g.nodes[p.exposure];
  const auto &y = g.nodes[p.outcome];
  const auto kind = e.get("kind").as_string(),
             measure = e.get("measure").as_string();
  p.family = y.family;
  p.edges = calibration_edges(g, e);
  p.testScale = positive(measure) ? "log" : "original";
  if (kind == "marginal") {
    if (a.kind != "binary" || !a.terms.empty())
      fail("Marginal arm-comparison power currently requires a binary "
           "comparison variable with no parents (independent randomization). A "
           "confounded or continuous comparison needs a different estimator.");
    const double allocation = 1 / (1 + std::exp(-a.intercept));
    if (!(allocation > 0 && allocation < 1))
      fail("The comparison variable gives an empty arm. Use a nondegenerate "
           "allocation probability.");
    if (measure == "mean_difference" && y.kind != "binary") {
      p.method = "welch";
      p.label = "Welch t test of the marginal mean difference";
    } else if (measure == "risk_difference" || measure == "mean_difference") {
      const bool score = e.get("nullValue").as_number() == 0;
      p.method = score ? "risk_score" : "risk_wald";
      p.label = score ? "Two-arm proportion score test of the risk difference"
                      : "Wald test of the risk difference";
    } else if (measure == "risk_ratio") {
      p.method = "risk_ratio";
      p.label = "Wald test of the marginal log risk ratio";
    } else if (measure == "odds_ratio") {
      p.method = "odds_ratio";
      p.label = "Wald test of the marginal log odds ratio";
    } else if (measure == "survival_difference") {
      p.method = "km_survival";
      p.label = "Kaplan–Meier / Greenwood Wald test at time " +
                e.get("horizon").to_text();
    } else if (measure == "rmst_difference") {
      p.method = "km_rmst";
      p.label = "Kaplan–Meier RMST Wald test through time " +
                e.get("horizon").to_text();
    }
    if (y.kind == "survival") {
      if (y.censor == "administrative" &&
          e.get("horizon").as_number() > y.cutoff)
        fail("The estimand horizon exceeds the administrative follow-up "
             "cutoff. It cannot be estimated by this KM analysis.");
      if (y.censor == "variable") {
        const auto t = ancestors(g, p.outcome),
                   c = ancestors(g, y.censor_index);
        std::string shared;
        for (int j = 0; j < static_cast<int>(g.nodes.size()); ++j)
          if (j != p.exposure && t[j] && c[j])
            shared += (shared.empty() ? "" : ", ") + g.nodes[j].id;
        if (!shared.empty())
          fail("The event and censoring mechanisms share " + shared +
               " beyond the randomized arm. Unadjusted KM does not identify "
               "this marginal target; change the censoring DGP or use a "
               "censoring-adjusted estimator.");
      }
    }
    p.formula = y.id + ": " + a.id + " = " + e.get("active").to_text() +
                " versus " + e.get("reference").to_text();
  } else {
    p.terms = y.terms;
    for (int j = 0; j < static_cast<int>(p.terms.size()); ++j) {
      const auto &t = p.terms[j];
      if (t.right >= 0)
        continue;
      if (t.left == p.exposure)
        p.targetIndex = j;
      const auto &parent = g.nodes[t.left];
      if (parent.kind == "survival" && parent.censor != "none")
        fail("The planned regression includes a censored event-time predictor. "
             "It cannot use its latent values as observed covariates.");
    }
    if (y.family == "gaussian" || y.family == "uniform" ||
        y.family == "lognormal") {
      p.method = y.family == "lognormal" ? "ols_log" : "ols";
      p.label = y.family == "lognormal"
                    ? "Log-response linear regression coefficient t test"
                    : "Linear regression coefficient t test";
      if (y.family == "uniform")
        p.label += " (approximate for Uniform errors)";
    } else if (y.family == "bernoulli" || y.family == "poisson" ||
               y.family == "gamma") {
      p.method = "glm";
      p.label = std::string(y.family == "bernoulli" ? "Logistic"
                            : y.family == "poisson" ? "Poisson log-link"
                                                    : "Gamma log-link") +
                " regression coefficient Wald test";
    } else {
      p.method = "cox";
      p.label = "Cox PH coefficient Wald test (Breslow ties)";
      if (y.censor == "variable" && ancestors(g, y.censor_index)[p.outcome])
        fail("The censoring variable is a descendant of the event time. The "
             "proposed Cox analysis requires independent censoring given its "
             "fitted parents.");
    }
    p.formula = (p.method == "ols_log" ? "log(" + y.id + ")" : y.id) + " ~ ";
    for (std::size_t j = 0; j < p.terms.size(); ++j) {
      const auto &t = p.terms[j];
      if (j)
        p.formula += " + ";
      p.formula += g.nodes[t.left].id;
      if (t.right >= 0)
        p.formula += " × " + g.nodes[t.right].id;
    }
  }
  return p;
}
J plan_json(const Graph &g, const Plan &p) {
  J terms = J::array(), edges = J::array();
  for (const auto &t : p.terms) {
    J x = object({{"left", g.nodes[t.left].id}});
    if (t.right >= 0)
      x.set("right", g.nodes[t.right].id);
    terms.push(x);
  }
  for (int i : p.edges)
    edges.push(g.json.get("edges").at(i));
  return object({{"method", p.method},
                 {"label", p.label},
                 {"formula", p.formula},
                 {"testScale", p.testScale},
                 {"terms", terms},
                 {"calibrationEdges", edges},
                 {"targetIndex", p.targetIndex},
                 {"outcome", g.nodes[p.outcome].id},
                 {"exposure", g.nodes[p.exposure].id},
                 {"family", p.family}});
}
struct Analysis {
  double estimate = 0, testEstimate = 0, se = 0, df = INFINITY, statistic = 0;
};
Analysis analyze(const Graph &g, const J &e, const Data &d, const Plan &p) {
  Analysis out;
  const auto measure = e.get("measure").as_string();
  if (e.get("kind").as_string() == "coefficient") {
    const auto m = model_matrix(d, p.terms, p.method != "cox");
    const int target = p.targetIndex + (p.method == "cox" ? 0 : 1);
    Fit fit;
    if (p.method == "cox")
      fit = fit_cox(m, d.times[p.outcome], d.events[p.outcome], target);
    else if (p.method == "glm")
      fit = fit_glm(m, d.values[p.outcome], target, p.family);
    else {
      Vec y = d.values[p.outcome];
      if (p.method == "ols_log")
        for (double &v : y)
          v = std::log(v);
      fit = fit_linear(m, y, target);
    }
    out.testEstimate = fit.coefficient;
    out.se = fit.se;
    out.df = fit.df;
    out.estimate = measure == "exponentiated_coefficient"
                       ? std::exp(fit.coefficient)
                       : fit.coefficient;
  } else {
    Vec values[2], events[2];
    const bool km = p.method == "km_survival" || p.method == "km_rmst";
    for (int i = 0; i < d.n; ++i) {
      const double a = d.values[p.exposure][i];
      const int arm = a == e.get("reference").as_number() ? 0
                      : a == e.get("active").as_number()  ? 1
                                                          : -1;
      if (arm >= 0) {
        values[arm].push_back(km ? d.times[p.outcome][i]
                                 : d.values[p.outcome][i]);
        if (km)
          events[arm].push_back(d.events[p.outcome][i]);
      }
    }
    if (values[0].size() < 2 || values[1].size() < 2)
      throw Error("empty_arm",
                  "A randomized arm has fewer than two observations.");
    if (km) {
      const auto lo = kaplan_meier(values[0], events[0],
                                   e.get("horizon").as_number()),
                 hi = kaplan_meier(values[1], events[1],
                                   e.get("horizon").as_number());
      const bool rmst = p.method == "km_rmst";
      out.estimate = rmst ? hi.rmst - lo.rmst : hi.survival - lo.survival;
      out.se = std::sqrt(rmst ? lo.rmstVariance + hi.rmstVariance
                              : lo.survivalVariance + hi.survivalVariance);
      out.testEstimate = out.estimate;
    } else {
      const auto s0 = moments(values[0]), s1 = moments(values[1]);
      const double p0 = s0.mean, p1 = s1.mean;
      if (p.method == "welch") {
        const double v0 = s0.variance / s0.n, v1 = s1.variance / s1.n;
        out.estimate = p1 - p0;
        out.se = std::sqrt(v0 + v1);
        out.df = (v0 + v1) * (v0 + v1) /
                 (v0 * v0 / (s0.n - 1) + v1 * v1 / (s1.n - 1));
        out.testEstimate = out.estimate;
      } else if (p.method == "risk_score" || p.method == "risk_wald") {
        out.estimate = p1 - p0;
        const double pooled = (p0 * s0.n + p1 * s1.n) / (s0.n + s1.n);
        out.se =
            std::sqrt(p.method == "risk_score"
                          ? pooled * (1 - pooled) * (1.0 / s0.n + 1.0 / s1.n)
                          : p0 * (1 - p0) / s0.n + p1 * (1 - p1) / s1.n);
        out.testEstimate = out.estimate;
      } else {
        if (!(p0 > 0 && p1 > 0) ||
            (p.method == "odds_ratio" && !(p0 < 1 && p1 < 1)))
          throw Error("boundary_response",
                      "The ratio analysis has a zero or boundary cell.");
        if (p.method == "risk_ratio") {
          out.estimate = p1 / p0;
          out.se = std::sqrt((1 - p0) / (s0.n * p0) + (1 - p1) / (s1.n * p1));
        } else {
          out.estimate = (p1 / (1 - p1)) / (p0 / (1 - p0));
          out.se = std::sqrt(1 / (s0.n * p0 * (1 - p0)) +
                             1 / (s1.n * p1 * (1 - p1)));
        }
        out.testEstimate = std::log(out.estimate);
      }
    }
  }
  if (!std::isfinite(out.estimate) || !std::isfinite(out.testEstimate) ||
      !(out.se > 0) || !std::isfinite(out.se) || std::isnan(out.df))
    throw Error("zero_variance",
                "The selected test has no finite positive standard error.");
  const double null = e.get("nullValue").as_number();
  out.statistic =
      (out.testEstimate - (positive(measure) ? std::log(null) : null)) / out.se;
  return out;
}
J analysis_json(const Analysis &a) {
  return object({{"estimate", a.estimate},
                 {"testEstimate", a.testEstimate},
                 {"se", a.se},
                 {"df", std::isfinite(a.df) ? J(a.df) : J()},
                 {"statistic", a.statistic}});
}
J population(const Graph &g, const J &e, int n, std::uint64_t seed,
             const Token *token) {
  validate_estimand(g, e);
  const auto kind = e.get("kind").as_string(),
             measure = e.get("measure").as_string();
  const int a = g.index(e.get("exposure").as_string()),
            y = g.index(e.get("outcome").as_string());
  if (kind == "coefficient") {
    double beta = 0;
    for (const auto &term : g.nodes[y].terms)
      if (term.left == a && term.right < 0)
        beta = term.coefficient;
    const double effect = positive(measure) ? std::exp(beta) : beta;
    if (!std::isfinite(effect))
      fail("The DGP coefficient ratio is not finite.");
    return object({{"estimate", effect},
                   {"mcSE", 0},
                   {"interval", J()},
                   {"n", 0},
                   {"exact", true}});
  }
  Vec arms[2];
  int arm = 0;
  for (const auto &key : {"reference", "active"}) {
    auto data =
        simulate(g, n, seed, {{a, e.get(key).as_number()}}, false, token);
    arms[arm] = std::move(data.values[y]);
    for (double &v : arms[arm]) {
      if (measure == "survival_difference")
        v = v > e.get("horizon").as_number() ? 1 : 0;
      else if (measure == "rmst_difference")
        v = std::min(v, e.get("horizon").as_number());
    }
    ++arm;
  }
  const double m0 = moments(arms[0]).mean, m1 = moments(arms[1]).mean;
  Vec influence(n);
  double estimate;
  if (!positive(measure)) {
    estimate = m1 - m0;
    for (int i = 0; i < n; ++i)
      influence[i] = arms[1][i] - arms[0][i];
  } else {
    if (!(m0 > 0 && m1 > 0) || (measure == "odds_ratio" && !(m0 < 1 && m1 < 1)))
      fail("The population simulation has boundary risks; increase its sample "
           "size or revise the DGP.");
    const bool rr = measure == "risk_ratio";
    estimate = rr ? m1 / m0 : (m1 / (1 - m1)) / (m0 / (1 - m0));
    const double d0 = rr ? m0 : m0 * (1 - m0), d1 = rr ? m1 : m1 * (1 - m1);
    for (int i = 0; i < n; ++i)
      influence[i] = (arms[1][i] - m1) / d1 - (arms[0][i] - m0) / d0;
  }
  const double se = std::sqrt(moments(influence).variance / n),
               mcse = positive(measure) ? estimate * se : se;
  if (!std::isfinite(estimate) || !std::isfinite(mcse))
    fail("The population estimand overflowed. Review the DGP parameters.");
  const Vec interval = positive(measure)
                           ? Vec{std::exp(std::log(estimate) - 1.96 * se),
                                 std::exp(std::log(estimate) + 1.96 * se)}
                           : Vec{estimate - 1.96 * se, estimate + 1.96 * se};
  return object({{"estimate", estimate},
                 {"mcSE", mcse},
                 {"interval", array(interval)},
                 {"n", n},
                 {"exact", false}});
}
struct Settings {
  J json;
  std::vector<int> sizes;
  int repetitions = 500, populationSize = 20000;
  std::uint64_t seed = 2026;
  double alpha = .05, targetPower = .8, lower = -3, upper = 3;
  std::string side = "two-sided", edge;
  bool checkNull = true;
};
Settings settings(const J &raw) {
  Settings s;
  const auto &sizes = raw.get("sampleSizes");
  if (sizes.is_null())
    s.sizes = {100, 200, 400};
  else {
    if (!sizes.is_array() || sizes.size() < 1 || sizes.size() > 10)
      fail("Choose 1–10 total sample sizes.");
    for (const auto &n : sizes.elements()) {
      const double v = n.as_number(-1);
      if (v < 4 || v > 50000 || v != std::floor(v))
        fail("Sample sizes must be whole numbers from 4 to 50,000.");
      s.sizes.push_back(static_cast<int>(v));
    }
    std::sort(s.sizes.begin(), s.sizes.end());
    s.sizes.erase(std::unique(s.sizes.begin(), s.sizes.end()), s.sizes.end());
  }
  s.repetitions = integer(raw, "repetitions", 10, 10000, 500);
  s.populationSize = integer(raw, "populationSize", 1000, 50000, 20000);
  s.seed = seed_value(raw);
  if (raw.find("alpha"))
    s.alpha = number(raw, "alpha");
  if (raw.find("targetPower"))
    s.targetPower = number(raw, "targetPower");
  if (!(s.alpha > 0 && s.alpha < 1 && s.targetPower > 0 && s.targetPower < 1))
    fail("Alpha and target power must be strictly between 0 and 1.");
  if (raw.find("side"))
    s.side = raw.get("side").as_string();
  if (s.side != "two-sided" && s.side != "less" && s.side != "greater")
    fail("Choose a valid rejection direction.");
  if (raw.find("calibrationLower"))
    s.lower = number(raw, "calibrationLower");
  if (raw.find("calibrationUpper"))
    s.upper = number(raw, "calibrationUpper");
  if (s.lower >= s.upper)
    fail("The lower calibration bound must be smaller than the upper bound.");
  if (raw.find("calibrationEdge")) {
    if (!raw.get("calibrationEdge").is_string())
      fail("Choose a calibration arrow.");
    s.edge = raw.get("calibrationEdge").as_string();
  }
  if (raw.find("checkNull")) {
    if (!raw.get("checkNull").is_bool())
      fail("The null-check setting must be true or false.");
    s.checkNull = raw.get("checkNull").as_bool();
  }
  J ns = J::array();
  for (int n : s.sizes)
    ns.push(n);
  s.json = object({{"sampleSizes", ns},
                   {"repetitions", s.repetitions},
                   {"alpha", s.alpha},
                   {"side", s.side},
                   {"seed", static_cast<unsigned long long>(s.seed)},
                   {"populationSize", s.populationSize},
                   {"calibrationEdge", s.edge},
                   {"calibrationLower", s.lower},
                   {"calibrationUpper", s.upper},
                   {"checkNull", s.checkNull},
                   {"targetPower", s.targetPower}});
  return s;
}
std::string edge_key(const J &edge) {
  return edge.get("source").as_string() + ">" + edge.get("target").as_string();
}
J calibrate(const Graph &g, const J &e, double target, const Settings &s,
            std::uint64_t seed, const Token *token, const Progress &progress) {
  const auto candidates = calibration_edges(g, e);
  int chosen = -1;
  for (int i : candidates)
    if (s.edge.empty() || edge_key(g.json.get("edges").at(i)) == s.edge) {
      chosen = i;
      break;
    }
  const auto effect = [&](const Graph &trial, std::uint64_t stream) {
    return population(trial, e, s.populationSize, stream, token);
  };
  const auto answer = [&](const J &graph, const J &truth, const J &parameter,
                          double tolerance) {
    return object({{"graph", graph},
                   {"truth", truth},
                   {"parameter", parameter},
                   {"target", target},
                   {"tolerance", tolerance}});
  };
  const auto parameter = [&](double value) {
    return object({{"edge", edge_key(g.json.get("edges").at(chosen))},
                   {"drawn", g.json.get("edges").at(chosen).get("coefficient")},
                   {"used", value}});
  };
  if (e.get("kind").as_string() == "coefficient") {
    if (chosen < 0)
      fail("The target coefficient arrow is missing.");
    const double value =
        positive(e.get("measure").as_string()) ? std::log(target) : target;
    const auto trial = parse_graph(change_arrow(g, chosen, value));
    return answer(trial.json, effect(trial, seed), parameter(value), 0);
  }
  const auto initial = effect(g, seed);
  const double tolerance = std::max(
      1e-5,
      .01 * std::max(std::fabs(target - e.get("nullValue").as_number()), .01));
  const auto verification = sim::stream_seed(seed, 0x564552494659ULL);
  if (std::fabs(initial.get("estimate").as_number() - target) <= tolerance)
    return answer(g.json, effect(g, verification), J(), tolerance);
  if (chosen < 0)
    fail("The requested planning alternative differs from the DGP, but no "
         "tunable arrow lies on a path from the comparison variable to the "
         "outcome.");
  if (progress)
    progress(object(
        {{"stage", "calibration"},
         {"message", "Calibrating " + edge_key(g.json.get("edges").at(chosen)) +
                         " to the selected estimand in C++/WebAssembly…"}}));
  struct Point {
    double value, error;
  };
  std::map<double, std::shared_ptr<Point>> cache;
  const auto evaluate = [&](double value) -> std::shared_ptr<Point> {
    if (token)
      token->check();
    if (cache.count(value))
      return cache.at(value);
    try {
      const auto truth =
          effect(parse_graph(change_arrow(g, chosen, value)), seed);
      auto point = std::make_shared<Point>(
          Point{value, truth.get("estimate").as_number() - target});
      cache[value] = point;
      return point;
    } catch (const Error &) {
      cache[value] = nullptr;
      return nullptr;
    }
  };
  auto a = evaluate(s.lower), b = evaluate(s.upper);
  if (!(a && b && a->error * b->error <= 0)) {
    std::shared_ptr<Point> previous;
    bool bracket = false;
    for (int i = 0; i <= 16; ++i) {
      auto current = evaluate(s.lower + (s.upper - s.lower) * i / 16);
      if (!current)
        continue;
      if (previous && previous->error * current->error <= 0) {
        a = previous;
        b = current;
        bracket = true;
        break;
      }
      previous = current;
    }
    if (!bracket)
      fail(
          "The planning effect cannot be reached within the calibration bounds "
          "using this arrow. Revise the bounds, tunable arrow, or DGP.");
  }
  auto best = std::fabs(a->error) < std::fabs(b->error) ? a : b;
  for (int i = 0; i < 40 && std::fabs(best->error) > tolerance; ++i) {
    auto middle = evaluate((a->value + b->value) / 2);
    if (!middle)
      fail("A numerical failure interrupted calibration. Revise the parameter "
           "bounds.");
    if (std::fabs(middle->error) < std::fabs(best->error))
      best = middle;
    if (a->error * middle->error <= 0)
      b = middle;
    else
      a = middle;
  }
  if (std::fabs(best->error) > tolerance)
    fail("The effect calibration did not meet its tolerance. Increase the "
         "population simulation size or revise the model.");
  const auto trial = parse_graph(change_arrow(g, chosen, best->value));
  return answer(trial.json, effect(trial, verification), parameter(best->value),
                tolerance);
}
J power(const Graph &g, const J &e, const Settings &s, int requested,
        Token *token, const Progress &progress) {
  const auto plan = compile_plan(g, e);
  if (!s.edge.empty() &&
      std::none_of(plan.edges.begin(), plan.edges.end(), [&](int i) {
        return edge_key(g.json.get("edges").at(i)) == s.edge;
      }))
    fail("The selected calibration arrow no longer belongs to this estimand. "
         "Choose it again.");
  const double draws = std::accumulate(s.sizes.begin(), s.sizes.end(), 0.0) *
                       s.repetitions * g.nodes.size() * (s.checkNull ? 2 : 1);
  if (draws > 200000000)
    fail("This run exceeds the browser work limit. Reduce sample sizes, "
         "repetitions, or the number of variables.");
  // Bound concurrent workspaces, especially for wide regressions. The
  // calculator's pool is still shared.
  const double terms = plan.terms.size() + 1,
               bytes = 8 * (s.sizes.back() * (g.nodes.size() * 3 + terms + 10) +
                            8 * terms * terms);
  const int memoryThreads = static_cast<int>(
      std::min(128.0, std::floor(256e6 / std::max(1.0, bytes))));
  if (memoryThreads < 1)
    fail("One fitted study exceeds the engine memory budget. Reduce the sample "
         "size or the number of analysis terms.");
  int threads = std::max(
      1, std::min({resolve_threads(requested), s.repetitions, memoryThreads}));
#if defined(__EMSCRIPTEN__) && !defined(__EMSCRIPTEN_PTHREADS__)
  threads = 1;
#endif
  const auto truthSeed = sim::stream_seed(s.seed, 0x5452555448ULL);
  if (progress)
    progress(object({{"stage", "calibration"},
                     {"message", "Linking the saved estimand to the "
                                 "alternative DGP in C++/WebAssembly…"}}));
  J alternative = e.get("alternative").is_null()
                      ? object({{"graph", g.json},
                                {"truth", population(g, e, s.populationSize,
                                                     truthSeed, token)},
                                {"parameter", J()},
                                {"target", J()},
                                {"tolerance", 0}})
                      : calibrate(g, e, e.get("alternative").as_number(), s,
                                  truthSeed, token, progress);
  J nullScenario;
  if (s.checkNull) {
    if (e.get("kind").as_string() == "marginal" &&
        e.get("nullValue").as_number() ==
            (positive(e.get("measure").as_string()) ? 1 : 0)) {
      const auto ng = parse_graph(null_graph(g, e));
      nullScenario = object(
          {{"graph", ng.json},
           {"truth", population(ng, e, s.populationSize, truthSeed, token)},
           {"exactNull", true},
           {"construction",
            "All effects of the comparison variable on causal paths to the "
            "outcome, including its interactions, are set to zero. Other paths "
            "and observation rules are retained."}});
    } else {
      nullScenario =
          calibrate(g, e, e.get("nullValue").as_number(), s,
                    sim::stream_seed(s.seed, 0x4e554c4cULL), token, progress);
      const bool exact = e.get("kind").as_string() == "coefficient";
      nullScenario.set("exactNull", exact);
      nullScenario.set(
          "construction",
          exact
              ? "The selected conditional coefficient is set to its null value."
              : "The reference DGP is Monte Carlo calibrated to the shifted "
                "null; its effect uncertainty must be considered when "
                "interpreting rejection.");
    }
  }
  J warnings = J::array(), results = J::array(), first;
  if (!alternative.get("target").is_null() &&
      std::fabs(alternative.get("truth").get("estimate").as_number() -
                alternative.get("target").as_number()) >
          alternative.get("tolerance").as_number() +
              3 * alternative.get("truth").get("mcSE").as_number())
    warnings.push(
        "Independent population verification differs from the planning effect "
        "beyond its numerical tolerance and Monte Carlo uncertainty. Increase "
        "population size and recalibrate.");
  const int scenarios = s.checkNull ? 2 : 1,
            total =
                static_cast<int>(s.sizes.size()) * s.repetitions * scenarios;
  int completed = 0;
  if (token)
    token->set_total(total);
  for (int n : s.sizes) {
    J row = object({{"n", n}});
    for (int scenario = 0; scenario < scenarios; ++scenario) {
      const auto runGraph =
          parse_graph((scenario ? nullScenario : alternative).get("graph"));
      sim::Accumulator accumulator;
      std::map<std::string, int> reasons;
      // Batches keep progress responsive and aggregation in replicate order,
      // independent of thread scheduling.
      constexpr int batchSize = 64;
      for (int begin = 0; begin < s.repetitions; begin += batchSize) {
        const int count = std::min(batchSize, s.repetitions - begin);
        std::vector<sim::Outcome> outcomes(count);
        std::vector<std::string> failures(count);
        parallel_for(
            count, threads,
            [&](std::size_t i) {
              auto &out = outcomes[i];
              try {
                const auto data =
                    simulate(runGraph, n,
                             sim::stream_seed(
                                 s.seed, static_cast<std::uint64_t>(begin + i),
                                 static_cast<std::uint64_t>(scenario)),
                             {}, true, token);
                const auto fit = analyze(runGraph, e, data, plan);
                out.reject = p_value(fit.statistic, fit.df, s.side) < s.alpha;
                out.estimate = fit.estimate;
                out.se = fit.se;
                out.df = fit.df;
              } catch (const Error &error) {
                out.failed = true;
                failures[i] = error.code;
              }
              if (token)
                token->add_done();
            },
            token);
        for (int i = 0; i < count; ++i) {
          accumulator.add(outcomes[i]);
          if (!failures[i].empty())
            ++reasons[failures[i]];
        }
        completed += count;
        if (progress)
          progress(object({{"stage", "simulation"},
                           {"completed", completed},
                           {"total", total},
                           {"n", n},
                           {"scenario", scenario ? "null" : "alternative"},
                           {"threads", threads}}));
      }
      const auto summary = accumulator.finish();
      J failures = J::object();
      for (const auto &kv : reasons)
        failures.set(kv.first, kv.second);
      row.set(
          scenario ? "null" : "alternative",
          object(
              {{"rejected", summary.rejections},
               {"failures", summary.failures},
               {"valid", summary.reps - summary.failures},
               {"repetitions", summary.reps},
               {"probability", summary.power},
               {"mcSE", summary.se},
               {"interval",
                array({summary.rejections ? summary.lo : 0,
                       summary.rejections == summary.reps ? 1 : summary.hi})},
               {"meanEstimate", summary.mean_estimate},
               {"estimateSD", summary.sd_estimate},
               {"meanSE", summary.mean_se},
               {"failureReasons", failures}}));
    }
    if (first.is_null() &&
        row.get("alternative").get("probability").as_number() >= s.targetPower)
      first = n;
    results.push(row);
  }
  return object(
      {{"version", 1},
       {"graph", g.json},
       {"estimand", e},
       {"power", s.json},
       {"analysis", object({{"method", plan.method},
                            {"label", plan.label},
                            {"formula", plan.formula},
                            {"testScale", plan.testScale}})},
       {"alternative", alternative},
       {"nullScenario", nullScenario},
       {"results", results},
       {"firstTestedN", first},
       {"warnings", warnings},
       {"engine", object({{"backend", "C++"},
                          {"threads", threads},
                          {"rng", "xoshiro256++ / SplitMix64"}})},
       {"denominator",
        "Every simulated study is counted. Generation or analysis failures "
        "count as non-rejections and are reported separately."}});
}
} // namespace
J run_request(const J &request, Token *token, const Progress &progress) {
  const auto graph = parse_graph(request.get("graph"));
  const auto op = request.get("op").as_string();
  const auto &e = request.get("estimand");
  if (token)
    token->check();
  if (op == "simulate") {
    const auto &options = request.get("options");
    const int n = integer(options, "n", 1, 50000, 500);
    const auto seed = seed_value(options);
    std::map<int, double> interventions;
    const auto &raw = options.get("interventions");
    if (!raw.is_null() && !raw.is_object())
      fail("Interventions must map variable IDs to fixed numeric values.");
    for (const auto &kv : raw.items())
      interventions[graph.index(kv.first)] = number(raw, kv.first);
    return data_json(graph,
                     simulate(graph, n, seed, interventions,
                              options.get("observe").as_bool(true), token),
                     seed, options.get("summarize").as_bool(true));
  }
  if (op == "plan")
    return plan_json(graph, compile_plan(graph, e));
  if (op == "analyze") {
    const auto fit = analyze(graph, e, read_data(graph, request.get("data")),
                             compile_plan(graph, e));
    J out = analysis_json(fit);
    const auto side = request.get("side").as_string();
    out.set("pValue",
            p_value(fit.statistic, fit.df, side.empty() ? "two-sided" : side));
    return out;
  }
  if (op == "population")
    return population(graph, e,
                      integer(request.get("options"), "n", 1, 50000, 20000),
                      seed_value(request.get("options")), token);
  if (op == "nullGraph") {
    validate_estimand(graph, e);
    return null_graph(graph, e);
  }
  const auto s = settings(request.get("settings"));
  if (op == "calibrate") {
    validate_estimand(graph, e);
    const double target = number(request, "target");
    if (positive(e.get("measure").as_string()) && !(target > 0))
      fail("A ratio target must be positive.");
    return calibrate(graph, e, target, s, seed_value(request), token, progress);
  }
  if (op == "power")
    return power(graph, e, s, integer(request, "threads", 0, 128, 0), token,
                 progress);
  fail("Unknown DAG engine operation: " + op);
}
} // namespace ssc::dag
