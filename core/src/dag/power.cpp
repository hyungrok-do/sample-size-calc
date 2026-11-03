#include "gformula.hpp"
#include "learners.hpp"
#include "model.hpp"
#include "study.hpp"
#include "sim/random.hpp"
#include "sim/simulate.hpp"
#include "statistics.hpp"
#include "stats/stats.hpp"
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
// The regression an analysis fits for a variable of this distribution:
// least squares, logistic (a Beta mean by quasi-likelihood), Poisson (a
// negative binomial mean by quasi-likelihood) or Gamma.
std::string model_family(const std::string &family) {
  return family == "bernoulli" || family == "beta"    ? "bernoulli"
         : family == "poisson" || family == "negbin"  ? "poisson"
         : family == "gamma"                          ? "gamma"
                                                      : "gaussian";
}
// Quasi-likelihood fits whose model-based variance is wrong for the
// distribution: their tests use a sandwich covariance.
bool quasi(const std::string &family) {
  return family == "beta" || family == "negbin";
}
// The terms a joint test of `measure` covers for predictor `a`.
bool joint_term(const Term &t, int a, const std::string &measure) {
  if (measure == "nonlinearity") return nonlinear_in(t, a);
  if (measure == "interaction") return t.right >= 0 && (t.left == a || t.right == a);
  return t.left == a || t.right == a;
}
bool unknown_form(const Graph &g, int v) { return !g.nodes[v].unknown.inputs.empty(); }
std::string unknown_note(const Graph &g, int v) {
  std::string inputs;
  for (int k : g.nodes[v].unknown.inputs) inputs += (inputs.empty() ? "" : ", ") + g.nodes[k].id;
  return g.nodes[v].id + " depends on " + inputs + " through a relationship of unknown form that the planned parametric models do not represent; inspect bias and coverage.";
}
// fitted and judged on the same rows, flexible learners overfit: in simulated
// confounded designs, 13-21% of studies rejected a true null without folds
const char *const CROSS_FIT_NOTE = "The Super Learner's nuisance models are cross-fitted in 3 folds. Fitted and evaluated on the same rows, flexible learners overfit and the test rejects a true null too often.";
const char *const SUPER_LEARNER_NOTE = "Nuisance models are fitted by a Super Learner: a 5-fold cross-validated, non-negative weighted ensemble of the training mean, main-terms and pairwise-product GLMs, additive restricted cubic splines, boosted depth-2 trees and kernel ridge regression. It uses the variables of the planned terms, not their functional forms.";
std::string unmeasured_note(const Graph &g, int a, int y) {
  std::string names;
  const auto before = ancestors(g, y);
  for (int v = 0; v < static_cast<int>(g.nodes.size()); ++v)
    if (!g.nodes[v].measured && before[v] && v != a)
      names += (names.empty() ? "" : ", ") + g.nodes[v].id;
  return names.empty() ? std::string()
                       : names + (names.find(',') == std::string::npos ? " is" : " are") +
                             " unmeasured, so the analysis cannot use " +
                             (names.find(',') == std::string::npos ? "it" : "them") +
                             "; the results show any bias that leaves.";
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
  if (kind != "marginal" && kind != "coefficient" && kind != "joint")
    fail("Choose an estimand type.");
  if (e.get("population").as_string().empty())
    fail("Describe the target population.");
  const int a = g.index(e.get("exposure").as_string()),
            y = g.index(e.get("outcome").as_string());
  if (a == y)
    fail("The outcome and comparison variable must differ.");
  const auto &outcome = g.nodes[y];
  const auto &exposure = g.nodes[a];
  if (!outcome.measured || !exposure.measured)
    fail("The compared variable and the outcome must be measured.");
  const double lo = number(e, "reference"), hi = number(e, "active");
  bool compatible = false;
  if (kind == "coefficient" || kind == "joint") {
    bool arrow = false;
    for (const auto &term : outcome.terms)
      if (term.left == a && term.right < 0 && term.transform.kind == "identity")
        arrow = true;
    if (!arrow)
      fail("A conditional coefficient requires an arrow from the selected "
           "predictor to the outcome.");
    for (int k : outcome.unknown.inputs)
      if (k == a)
        fail("The effect of " + exposure.id + " on " + outcome.id + " has an unknown form, so no coefficient describes it. Compare groups or doses instead.");
    compatible = measure == "coefficient" ||
                 (measure == "exponentiated_coefficient" &&
                  outcome.family != "gaussian" && outcome.family != "uniform" &&
                  outcome.family != "student");
    if (kind == "joint") {
      compatible = measure == "overall_effect" || measure == "nonlinearity" ||
                   measure == "interaction";
      if (number(e,"nullValue") != 0 || !e.get("alternative").is_null())
        fail("Joint tests use a zero coefficient-block null and the graph as drawn.");
      if (measure == "nonlinearity" && exposure.kind == "binary")
        fail("A binary predictor has no nonlinear dose-response component.");
    }
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
        !g.study.enabled && std::fabs(v) > 1)
      fail("A probability difference must lie between -1 and 1.");
    if (measure == "rmst_difference" && std::fabs(v) > horizon)
      fail("An RMST difference cannot exceed its horizon in magnitude.");
  }
  if (!e.get("alternative").is_null() &&
      e.get("alternative").as_number() == e.get("nullValue").as_number())
    fail("The planning alternative must differ from the null value.");
}
struct Plan {
  std::string method, label, formula, testScale, family, outcomeModel, learner = "glm";
  int exposure = 0, outcome = 0, targetIndex = -1, folds = 1;
  bool bounded = false;
  double lower = 0, upper = 1;
  std::vector<int> tested;
  J specification;
  // For AIPW, `terms` are the outcome model's and `propensity` the exposure
  // model's.
  std::vector<Term> terms, propensity;
  std::vector<int> edges;
  // For the g-formula, the refitted models in generation order.
  std::vector<ChainModel> chain;
  // What the analysis assumes that the drawn DGP does not grant: reported
  // with the results instead of refusing the run.
  std::vector<std::string> notes;
};
void note(Plan &p, const std::string &text) {
  if (std::find(p.notes.begin(), p.notes.end(), text) == p.notes.end())
    p.notes.push_back(text);
}
bool censored(const Node &node) {
  return node.kind == "survival" && node.censor != "none";
}
std::string latent_note(const Node &node) {
  return node.id + " is a censored event time; the analysis uses its latent "
                   "values as if they were measured.";
}
std::string informative_note(const Graph &g, const Node &node) {
  const auto &c = g.nodes[node.censor_index].id;
  return c + " censors " + node.id + " but depends on " + node.id +
         " itself, so no analysis of the observed data is unbiased; the "
         "results show the bias.";
}
std::string reverse_note(const Node &a, const Node &y) {
  return y.id + " lies upstream of " + a.id + ", so " + a.id +
         " cannot change " + y.id +
         ": the effect is zero by construction and no study can detect one. "
         "Swap the two roles to study the effect of " + y.id + " on " + a.id +
         ".";
}
std::string no_path_note(const Node &a, const Node &y) {
  return "No arrow path leads from " + a.id + " to " + y.id + ", so " + a.id +
         " has no effect on " + y.id +
         " here; the rejection rate is the test's false-positive rate.";
}
std::string term_text(const Graph &g, const std::vector<Term> &terms) {
  std::string out;
  for (std::size_t j = 0; j < terms.size(); ++j) {
    const auto &t = terms[j];
    if (j)
      out += " + ";
    out += transform_text(g.nodes[t.left].id, t.transform);
    if (t.right >= 0)
      out += " × " + transform_text(g.nodes[t.right].id, t.rightTransform);
  }
  return out;
}
// A confounded binary exposure: the doubly robust (AIPW) estimator, with the
// exposure's own DGP terms as the propensity model and, within each arm, an
// outcome model on the exposure's parents plus the outcome's parents that are
// not affected by the exposure. Conditioning on the exposure's parents
// identifies the marginal target in this DAG, and the added outcome parents
// only improve precision.
void plan_aipw(const Graph &g, const J &e, Plan &p) {
  p.terms.clear(); p.propensity.clear(); p.chain.clear();
  const auto &a = g.nodes[p.exposure];
  const auto &y = g.nodes[p.outcome];
  const auto aterms = measured_terms(g, a.terms), yterms = measured_terms(g, y.terms);
  const auto affected = descendants(g, p.exposure);
  const int count = static_cast<int>(g.nodes.size());
  std::vector<bool> chosen(count);
  for (const auto &t : aterms)
    for (int v : {t.left, t.right})
      if (v >= 0) {
        const auto &node = g.nodes[v];
        if (censored(node))
          note(p, latent_note(node));
        chosen[v] = true;
      }
  for (const auto &t : yterms) {
    if (t.right >= 0 || affected[t.left])
      continue;
    const auto &node = g.nodes[t.left];
    if (!(node.kind == "survival" && node.censor != "none") &&
        (!g.study.enabled || node.scope=="subject" || node.scope=="cluster"))
      chosen[t.left] = true;
  }
  for (int v = 0; v < count; ++v)
    if (chosen[v])
      p.terms.push_back({v, -1, 0});
  // Within an exposure arm its factors are constants: A * f(X) contributes
  // f(X) to the arm-specific model. Keep that basis even without a main f(X).
  for (const auto *source : {&aterms, &yterms}) for (auto t : *source) {
    if (t.left == p.exposure && t.right >= 0) {
      t.left=t.right; t.transform=t.rightTransform; t.right=-1; t.rightTransform=Transform{};
    } else if (t.right == p.exposure) { t.right=-1; t.rightTransform=Transform{}; }
    if (chosen[t.left] && (t.right < 0 || chosen[t.right]) &&
        std::none_of(p.terms.begin(), p.terms.end(), [&](const Term &u) {
          return term_key(u) == term_key(t);
        })) p.terms.push_back(t);
  }
  for (const auto &t : aterms)
    p.propensity.push_back(t);
  p.outcomeModel = y.family == "lognormal" ? "gamma" : model_family(y.family);
  if (const auto text = unmeasured_note(g, p.exposure, p.outcome); !text.empty())
    note(p, text);
  const auto measure = e.get("measure").as_string();
  p.method = "aipw";
  p.label = "Doubly robust AIPW Wald test of the marginal " +
            std::string(measure == "risk_ratio"   ? "log risk ratio"
                        : measure == "odds_ratio" ? "log odds ratio"
                        : measure == "risk_difference" ||
                                y.kind == "binary"
                            ? "risk difference"
                            : "mean difference") +
            (p.terms.empty() ? std::string(", adjusting for nothing measured")
                             : ", adjusting for " + term_text(g, p.terms));
  const std::string model =
      p.outcomeModel == "gaussian"    ? "linear regression"
      : p.outcomeModel == "bernoulli" ? "logistic regression"
      : p.outcomeModel == "poisson"   ? "Poisson log-link regression"
                                      : "Gamma log-link regression";
  p.formula = y.id + ": " + a.id + " = " + e.get("active").to_text() +
              " versus " + e.get("reference").to_text() + "; propensity " +
              a.id + " ~ " + term_text(g, p.propensity) +
              " (logistic); outcome " + y.id + " ~ " + term_text(g, p.terms) +
              " within each arm (" + model + ")";
}
std::string measure_scale(const J &e, const Node &y) {
  const auto measure = e.get("measure").as_string();
  return measure == "risk_ratio"            ? std::string("log risk ratio")
         : measure == "odds_ratio"          ? std::string("log odds ratio")
         : measure == "survival_difference" ? "survival difference at time " +
                                                  e.get("horizon").to_text()
         : measure == "rmst_difference"
             ? "RMST difference through time " + e.get("horizon").to_text()
         : measure == "risk_difference" || y.kind == "binary"
             ? std::string("risk difference")
             : std::string("mean difference");
}
// Any other marginal target: the parametric g-formula. The analysis refits
// the DGP's own conditional model for the outcome and for every variable the
// comparison affects on the way to it (event times by censored Weibull
// likelihood), then averages the outcome over the study's rows with the
// comparison set to each level. Conditioning each model on its parents
// identifies the target in the drawn DAG, whatever the comparison variable's
// type or parents, and censoring that shares causes with the event is
// ignorable given the event time's parents.
void plan_gformula(const Graph &g, const J &e, Plan &p) {
  p.terms.clear(); p.propensity.clear(); p.chain.clear();
  const auto &a = g.nodes[p.exposure];
  const auto &y = g.nodes[p.outcome];
  const auto after = descendants(g, p.exposure),
             before = ancestors(g, p.outcome);
  if (descendants(g, p.outcome)[p.exposure]) {
    // The outcome causes the comparison variable: nothing to refit, and the
    // contrast is zero in every study.
    p.method = "gformula";
    p.label = "Parametric g-formula Wald test of the marginal " +
              std::string(measure_scale(e, y));
    p.formula = y.id + ": " + a.id + " = " + e.get("active").to_text() +
                " versus " + e.get("reference").to_text() +
                "; zero by construction";
    return;
  }
  // Unmeasured variables are not refitted: each model uses the measured
  // variables in their place.
  for (int v : g.order)
    if (v != p.exposure && after[v] && before[v] && g.nodes[v].measured)
      p.chain.push_back({v, measured_terms(g, g.nodes[v].terms)});
  if (const auto text = unmeasured_note(g, p.exposure, p.outcome); !text.empty())
    note(p, text);
  if (!after[p.outcome]) {
    // No path: the analysis still estimates a comparison, which is null.
    auto terms = measured_terms(g, y.terms);
    terms.insert(terms.begin(), Term{p.exposure, -1, 0});
    p.chain.push_back({p.outcome, terms});
  }
  for (const auto &c : p.chain)
    for (const auto &t : c.terms)
      for (int v : {t.left, t.right})
        if (v >= 0 && v != p.exposure && censored(g.nodes[v]))
          note(p, latent_note(g.nodes[v]));
  for (const auto &c : p.chain) {
    const auto &node = g.nodes[c.node];
    if (node.censor == "variable" && descendants(g, c.node)[node.censor_index])
      note(p, informative_note(g, node));
  }
  for (const auto &c : p.chain)
    if (unknown_form(g, c.node)) note(p, unknown_note(g, c.node));
  if (y.censor == "administrative" && e.get("horizon").as_number() > y.cutoff)
    note(p, "The horizon " + e.get("horizon").to_text() +
                " lies past the follow-up cutoff " + J(y.cutoff).to_text() +
                "; the Weibull model extrapolates beyond the data.");
  p.terms = p.chain.back().terms;
  p.method = "gformula";
  p.label = "Parametric g-formula Wald test of the marginal " +
            measure_scale(e, y);
  std::string models;
  for (const auto &c : p.chain)
    models += (models.empty() ? "" : ", ") + g.nodes[c.node].id + " ~ " +
              term_text(g, c.terms);
  p.formula = y.id + ": " + a.id + " = " + e.get("active").to_text() +
              " versus " + e.get("reference").to_text() + "; refitted " +
              models + "; other variables as observed";
}
Plan default_plan(const Graph &g, const J &e) {
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
    const bool upstream = descendants(g, p.outcome)[p.exposure];
    if (upstream)
      note(p, reverse_note(a, y));
    else if (!descendants(g, p.exposure)[p.outcome])
      note(p, no_path_note(a, y));
    // A two-level comparison of a confounded binary variable uses AIPW; any
    // other comparison variable, or a confounded one with an event-time
    // outcome, the g-formula.
    // Parents the analysis cannot measure do not make it adjust; they leave
    // the bias that the unmeasured note reports.
    const bool confounded = !measured_terms(g, a.terms).empty();
    if (upstream || a.kind != "binary" || (confounded && y.kind == "survival")) {
      plan_gformula(g, e, p);
      return p;
    }
    if (confounded) {
      plan_aipw(g, e, p);
      return p;
    }
    if (const auto text = unmeasured_note(g, p.exposure, p.outcome);
        !text.empty() && !a.terms.empty())
      note(p, text);
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
    // Kaplan–Meier needs follow-up to the horizon and censoring independent
    // of the event time within arms; otherwise the g-formula.
    if (y.kind == "survival") {
      bool shared = false;
      if (y.censor == "variable") {
        const auto t = ancestors(g, p.outcome),
                   c = ancestors(g, y.censor_index);
        for (int j = 0; j < static_cast<int>(g.nodes.size()); ++j)
          shared |= j != p.exposure && t[j] && c[j];
      }
      if (shared || (y.censor == "administrative" &&
                     e.get("horizon").as_number() > y.cutoff)) {
        plan_gformula(g, e, p);
        return p;
      }
    }
    p.formula = y.id + ": " + a.id + " = " + e.get("active").to_text() +
                " versus " + e.get("reference").to_text();
  } else {
    p.terms = measured_terms(g, y.terms);
    if (const auto text = unmeasured_note(g, p.exposure, p.outcome); !text.empty())
      note(p, text);
    if (unknown_form(g, p.outcome)) note(p, unknown_note(g, p.outcome));
    for (int j = 0; j < static_cast<int>(p.terms.size()); ++j) {
      const auto &t = p.terms[j];
      if (t.right >= 0 || t.transform.kind != "identity")
        continue;
      if (t.left == p.exposure)
        p.targetIndex = j;
      if (censored(g.nodes[t.left]))
        note(p, latent_note(g.nodes[t.left]));
    }
    if (y.family == "gaussian" || y.family == "uniform" ||
        y.family == "lognormal" || y.family == "student") {
      p.method = y.family == "lognormal" ? "ols_log" : "ols";
      p.label = y.family == "lognormal"
                    ? "Log-response linear regression coefficient t test"
                    : "Linear regression coefficient t test";
      if (y.family == "uniform")
        p.label += " (approximate for Uniform errors)";
      if (y.family == "student")
        p.label += " (approximate for Student t errors)";
    } else if (y.kind != "survival") {
      p.method = "glm";
      p.label = std::string(y.family == "bernoulli" ? "Logistic"
                            : y.family == "beta"    ? "Quasi-binomial logit-mean"
                            : y.family == "poisson" ? "Poisson log-link"
                            : y.family == "negbin"  ? "Quasi-Poisson log-link"
                                                    : "Gamma log-link") +
                " regression coefficient Wald test" +
                (quasi(y.family) ? " with a sandwich standard error" : "");
    } else {
      p.method = "cox";
      p.label = "Cox PH coefficient Wald test (Breslow ties)";
      if (y.censor == "variable" && ancestors(g, y.censor_index)[p.outcome])
        note(p, informative_note(g, y));
    }
    p.formula = (p.method == "ols_log" ? "log(" + y.id + ")" : y.id) +
                " ~ " + term_text(g, p.terms);
  }
  return p;
}
// Configure the point-treatment targeted estimator. Repeated outcomes are
// first reduced to the prespecified subject-level contrast, not treated as
// independent rows and not conditioned on post-treatment visit variables.
void configure_targeted(const Graph &g, const J &e, Plan &p) {
  const auto &spec=p.specification;
  p.terms=analysis_terms(g,p.terms,spec,p.outcome);
  p.propensity=analysis_terms(g,p.propensity,spec,p.exposure,true);
  const auto affected=descendants(g,p.exposure);
  for(const auto *terms:{&p.terms,&p.propensity}) for(const auto &t:*terms)
    for(int v:{t.left,t.right}) if(v>=0) {
      if(affected[v]) fail("AIPW/TMLE adjustment cannot include the exposure or its descendants.");
      if(g.study.enabled && g.nodes[v].scope!="subject" && g.nodes[v].scope!="cluster")
        fail("Trajectory TMLE adjustments must be measured before treatment at subject or cluster level.");
    }
  if(p.method=="tmle") {
    const bool binary=g.nodes[p.outcome].kind=="binary";
    p.bounded=binary || spec.find("outcomeBounds");
    if(binary && spec.find("outcomeBounds")) fail("Binary outcome bounds are determined by the requested contrast.");
    if(binary && g.study.enabled) {
      p.lower=p.upper=0;
      for(double w:study_weights(g,e)) { p.lower+=std::min(0.0,w); p.upper+=std::max(0.0,w); }
    } else if(spec.find("outcomeBounds")) {
      p.lower=spec.get("outcomeBounds").at(0).as_number(); p.upper=spec.get("outcomeBounds").at(1).as_number();
    }
    // A Gaussian initial model also permits signed trajectory summaries. For
    // bounded continuous summaries it is rescaled and clipped before the
    // logistic fluctuation, without clipping observed responses.
    p.outcomeModel=binary && !g.study.enabled ? "bernoulli" : "gaussian";
    note(p,p.bounded ? "TMLE uses a logistic fluctuation within known outcome-summary bounds. Initial scaled predictions are clipped to [0.000001, 0.999999]; observed outcomes outside the prespecified bounds fail the analysis." : "TMLE uses a least-squares initial regression and an identity-link fluctuation for an unbounded outcome mean. It does not constrain predictions to a probability range.");
    if(p.learner!="superlearner") note(p,"TMLE is a semiparametric estimator, but these nuisance regressions use the chosen finite-dimensional bases; choose the Super Learner to fit them without a functional form. Power still depends on the simulated data-generating model.");
  }
  if(p.learner=="superlearner") note(p,SUPER_LEARNER_NOTE);
  else for(int v:{p.exposure,p.outcome}) if(unknown_form(g,v)) note(p,unknown_note(g,v));
  if(spec.get("propensityBound").as_number()>0) note(p,"Propensity truncation is enabled. The reported truncation count and largest inverse weight show its use; truncation can introduce bias and does not repair lack of overlap.");
  note(p,"Influence-function Wald inference is asymptotic. Rare outcomes, weak overlap and few independent units can give poor coverage. Inspect null rejection, coverage and failed fits.");
  p.label=(p.method=="tmle" ? "TMLE" : "Doubly robust AIPW")+std::string(g.study.enabled ? " of the complete trajectory contrast" : " Wald test of the marginal ")+ (g.study.enabled ? "" : measure_scale(e,g.nodes[p.outcome]))+
    (p.folds>1 ? " with "+std::to_string(p.folds)+"-fold cross-fitting" : "")+(p.learner=="superlearner" ? std::string(p.folds>1 ? " and" : " with")+" Super Learner nuisance models" : "");
  if(p.learner=="superlearner") {
    const auto variables=[&](const std::vector<Term> &terms) {
      std::vector<bool> used(g.nodes.size()); std::string out;
      for(const auto &t:terms) { used[t.left]=true; if(t.right>=0) used[t.right]=true; }
      for(int v=0;v<static_cast<int>(g.nodes.size());++v) if(used[v]) out+=(out.empty()?"":", ")+g.nodes[v].id;
      return "Super Learner("+(out.empty() ? std::string("no covariates") : out)+")";
    };
    p.formula="Propensity "+g.nodes[p.exposure].id+" ~ "+variables(p.propensity)+"; arm-specific outcome "+g.nodes[p.outcome].id+" ~ "+variables(p.terms);
  } else
  p.formula="Propensity "+g.nodes[p.exposure].id+" ~ "+term_text(g,p.propensity)+" (logistic); arm-specific outcome "+g.nodes[p.outcome].id+" ~ "+term_text(g,p.terms)+" ("+p.outcomeModel+")";
}
Plan compile_plan(const Graph &g, const J &e, const J &raw = J()) {
  const J spec = parse_analysis(raw);
  for(const auto &kv:spec.get("models").items()) g.index(kv.first);
  if(g.study.enabled) {
    validate_estimand(g,e);
    const auto p0=study_plan(g,e,spec);
    Plan p; p.method=p0.get("method").as_string(); p.label=p0.get("label").as_string(); p.formula=p0.get("formula").as_string();
    p.exposure=g.index(e.get("exposure").as_string()); p.outcome=g.index(e.get("outcome").as_string());
    p.family=g.nodes[p.outcome].family; p.testScale="identity"; p.specification=spec; p.edges=calibration_edges(g,e);
    for(const auto &n:p0.get("notes").elements()) p.notes.push_back(n.as_string());
    if(p.method=="tmle") {
      if(descendants(g,p.outcome)[p.exposure]) fail("The outcome cannot be upstream of treatment for TMLE.");
      plan_aipw(g,e,p); p.method="tmle";
      p.folds=integer(spec,"folds",1,10,1); p.learner=spec.find("learner") ? spec.get("learner").as_string() : "glm";
      if(p.learner=="superlearner" && p.folds<2) { p.folds=3; note(p,CROSS_FIT_NOTE); }
      configure_targeted(g,e,p);
    } else {
      // structured estimators refuse folds in study_plan; a learner would be ignored
      if(spec.find("learner") && spec.get("learner").as_string()=="superlearner")
        fail("The Super Learner fits the nuisance models of AIPW and TMLE. Select one of those methods.");
    }
    return p;
  }
  if(spec.get("method").as_string()=="mmrm" || spec.get("method").as_string()=="lmm" || spec.get("method").as_string()=="trajectory") fail("A trajectory analysis requires study structure.");
  Plan p = default_plan(g,e);
  const auto method=spec.get("method").as_string(), kind=e.get("kind").as_string();
  const auto &a=g.nodes[p.exposure], &y=g.nodes[p.outcome];
  if(kind != "marginal" && method != "auto") fail("Conditional and joint targets use a regression test.");
  if(method=="gformula") plan_gformula(g,e,p);
  if(method=="aipw" || method=="tmle") {
    if(a.kind!="binary" || y.kind=="survival" || descendants(g,p.outcome)[p.exposure])
      fail("AIPW/TMLE require a binary exposure, a nonsurvival outcome, and an outcome downstream or unrelated to exposure.");
    plan_aipw(g,e,p);
    p.method=method;
  }
  p.specification=spec;
  p.folds=integer(spec,"folds",1,10,1);
  p.learner=spec.find("learner") ? spec.get("learner").as_string() : "glm";
  // The analyst knows which variables matter but not how: TMLE with a Super
  // Learner and cross-fitting needs no functional form for them.
  if(method=="auto" && p.method=="aipw" && (unknown_form(g,p.exposure) || unknown_form(g,p.outcome))) {
    p.method="tmle"; p.learner="superlearner"; if(p.folds==1) p.folds=3;
    note(p,"Automatic analysis chose TMLE with a Super Learner and cross-fitting: a nuisance model has a relationship of unknown form.");
  }
  if(p.learner=="superlearner" && p.method!="aipw" && p.method!="tmle")
    fail("The Super Learner fits the nuisance models of AIPW and TMLE. Select one of those methods.");
  if(p.folds>1 && p.method!="aipw" && p.method!="tmle")
    fail("Cross-fitting is available for AIPW and TMLE. Select one of those methods.");
  if(p.learner=="superlearner" && p.folds<2) { p.folds=3; note(p,CROSS_FIT_NOTE); }
  if(p.method=="aipw" || p.method=="tmle") configure_targeted(g,e,p);
  else p.terms=analysis_terms(g,p.terms,spec,p.outcome);
  if(p.method=="gformula") {
    for(auto &c:p.chain) {
      c.terms=analysis_terms(g,c.terms,spec,c.node);
      const auto later=descendants(g,c.node);
      for(const auto &t:c.terms) if(later[t.left] || (t.right>=0 && later[t.right]))
        fail("A g-computation model cannot condition on itself or its descendants.");
    }
    if(!p.chain.empty()) p.terms=p.chain.back().terms;
    std::string models;
    for(const auto &c:p.chain) models+=(models.empty()?"":", ")+g.nodes[c.node].id+" ~ "+term_text(g,c.terms);
    if(!p.chain.empty()) p.formula=y.id+": "+a.id+" = "+e.get("active").to_text()+" versus "+e.get("reference").to_text()+"; refitted "+models+"; other variables as observed";
  }
  if(kind=="coefficient" || kind=="joint") {
    p.targetIndex=-1;
    for(int j=0;j<static_cast<int>(p.terms.size());++j) {
      const auto &t=p.terms[j];
      if(t.left==p.exposure && t.right<0 && t.transform.kind=="identity") p.targetIndex=j;
      if(kind=="joint" && joint_term(t,p.exposure,e.get("measure").as_string())) p.tested.push_back(j+1);
    }
    if(kind=="coefficient" && p.targetIndex<0) fail("The planned model omits the selected linear coefficient.");
    if(kind=="joint") {
      if(p.tested.empty()) fail("The planned model has no terms for the selected joint test.");
      p.method="joint";
      const auto measure=e.get("measure").as_string();
      p.label=std::string(p.tested.size()==1 && measure=="interaction" ? "Wald test of the interaction" : std::string("Joint Wald test of ")+(measure=="nonlinearity" ? "nonlinearity" : measure=="interaction" ? "the interactions" : "all predictor terms"))+" ("+std::to_string(p.tested.size())+" df)";
      if(y.kind=="survival") note(p,"The joint survival test fits a parametric Weibull proportional-hazards model.");
      else note(p,"The joint Wald test uses a sandwich covariance and a large-sample chi-square reference distribution.");
    }
    p.formula=(y.family=="lognormal" ? "log("+y.id+")" : y.id)+" ~ "+term_text(g,p.terms);
  }
  if(spec.get("basis").as_string()!="match" || spec.find("propensity") || spec.find("models"))
    note(p,"Analysis terms are specified separately from the generating model. Inspect bias, interval coverage and null rejection when they differ.");
  if(kind=="marginal" && p.method!="gformula" && p.method!="aipw" && p.method!="tmle" && unknown_form(g,p.outcome))
    note(p,y.id+" depends on its inputs through a relationship of unknown form; the two-arm test does not model it, so it stays valid.");
  if(kind=="marginal" && p.method!="gformula" && p.method!="aipw" && p.method!="tmle" &&
      (spec.get("basis").as_string()!="match" || spec.find("propensity") || spec.find("models")))
    note(p,"The automatic two-arm test does not fit covariate models. Select g-computation, AIPW or TMLE to use the specified regression terms.");
  return p;
}
struct Analysis {
  double estimate = 0, testEstimate = 0, se = 0, df = INFINITY, statistic = 0, jointP = -1;
  J diagnostics;
};
J terms_json(const Graph &g, const std::vector<Term> &list) {
  J terms = J::array();
  for (const auto &t : list) {
    terms.push(term_json(g,t));
  }
  return terms;
}
J plan_json(const Graph &g, const Plan &p) {
  J terms = terms_json(g, p.terms), edges = J::array();
  for (int i : p.edges)
    edges.push(g.json.get("edges").at(i));
  J out = object({{"method", p.method},
                 {"label", p.label},
                 {"formula", p.formula},
                 {"testScale", p.testScale},
                 {"terms", terms},
                 {"calibrationEdges", edges},
                 {"targetIndex", p.targetIndex},
                 {"outcome", g.nodes[p.outcome].id},
                 {"exposure", g.nodes[p.exposure].id},
                 {"family", p.family}});
  if (p.method == "aipw" || p.method == "tmle") {
    out.set("propensityTerms", terms_json(g, p.propensity));
    out.set("outcomeModel", p.outcomeModel);
    out.set("learner", p.learner);
    if(p.bounded) out.set("outcomeBounds",array({p.lower,p.upper}));
  }
  if (p.method == "gformula") {
    J models = J::array();
    for (const auto &c : p.chain)
      models.push(object(
          {{"node", g.nodes[c.node].id}, {"terms", terms_json(g, c.terms)}}));
    out.set("models", models);
  }
  J notes = J::array();
  for (const auto &text : p.notes)
    notes.push(text);
  out.set("notes", notes);
  out.set("analysis", p.specification);
  out.set("folds", p.folds);
  J tested=J::array(); for(int i:p.tested) tested.push(i);
  if(!p.tested.empty()) out.set("testedIndices",tested);
  return out;
}
// Point-treatment AIPW / TMLE. Nuisance fits exclude held-out independent
// units. TMLE solves both arm-specific efficient-score equations and uses the
// substitution estimate. The same EIF is summed within independent clusters.
Analysis analyze_aipw(const J &e, const Data &d, const Plan &p) {
  const int n=d.n;
  const bool clustered=!d.cluster.empty();
  Vec z(n), y=d.values[p.outcome];
  std::vector<int> rows[2];
  for(int i=0;i<n;++i) {
    z[i]=d.values[p.exposure][i]==e.get("active").as_number() ? 1 : 0;
    rows[static_cast<int>(z[i])].push_back(i);
    if(p.bounded) {
      if(y[i]<p.lower-1e-12 || y[i]>p.upper+1e-12)
        throw Error("outcome_bounds","An observed outcome summary is outside its prespecified bounds.");
      y[i]=(y[i]-p.lower)/(p.upper-p.lower);
    }
  }
  if(rows[0].size()<2 || rows[1].size()<2) throw Error("empty_arm","An analysis arm has fewer than two observations.");
  const bool learned=p.learner=="superlearner";
  const auto pm=learned ? Matrix{} : model_matrix(d,p.propensity,true), om=learned ? Matrix{} : model_matrix(d,p.terms,true);
  const auto variables=[&](const std::vector<Term> &terms) {
    std::vector<int> out;
    for(const auto &t:terms) for(int v:{t.left,t.right}) if(v>=0 && std::find(out.begin(),out.end(),v)==out.end()) out.push_back(v);
    std::sort(out.begin(),out.end());
    return out;
  };
  const auto pv=variables(p.propensity);
  Vec propensityWeights(LEARNERS.size()), outcomeWeights(LEARNERS.size());
  // The Super Learner's outcome model is one regression of Y on the exposure
  // and the covariates over both arms, evaluated with the exposure set to each
  // level: the data are stacked as observed, at the reference level and at the
  // active level, and only observed rows are fitted.
  Data stacked;
  auto ov=variables(p.terms);
  if(learned) {
    // the exposure leads, so its products with the covariates (the effect
    // modification) are the ones kept when the product learner is capped
    ov.erase(std::remove(ov.begin(),ov.end(),p.exposure),ov.end());
    ov.insert(ov.begin(),p.exposure);
    stacked.n=3*n; stacked.values.resize(d.values.size());
    for(int v:ov) {
      auto &column=stacked.values[v]; column.resize(3*n);
      for(int s=0;s<3;++s) for(int i=0;i<n;++i) column[s*n+i]=v==p.exposure ? (s==0 ? z[i] : s-1) : d.values[v][i];
    }
  }
  Vec stackedY(learned ? 3*n : 0);
  for(int i=0;i<n && learned;++i) stackedY[i]=y[i];
  // the Super Learner's own folds keep clusters together, as the cross-fitting does
  std::vector<int> stackedClusters;
  if(learned && clustered) { stackedClusters.resize(3*n); for(int s=0;s<3;++s) for(int i=0;i<n;++i) stackedClusters[s*n+i]=d.cluster[i]; }
  const std::vector<int> *groups=clustered ? &d.cluster : nullptr, *stackedGroups=clustered ? &stackedClusters : nullptr;
  const int units=clustered ? d.units : n;
  if(units<=2 || (p.folds>1 && units<p.folds)) throw Error("few_units","Too few independent units for this analysis or fold count.");
  std::vector<int> unitArm(units,-1), unitFold(units), fold(n);
  for(int i=0;i<n;++i) {
    const int u=clustered ? d.cluster[i] : i, a=static_cast<int>(z[i]);
    if(unitArm[u]==-1) unitArm[u]=a; else if(unitArm[u]!=a) unitArm[u]=2;
  }
  int ranks[3]={0,0,0};
  for(int u=0;u<units;++u) unitFold[u]=ranks[unitArm[u]]++%p.folds;
  for(int i=0;i<n;++i) fold[i]=unitFold[clustered ? d.cluster[i] : i];
  Vec pi(n), m0(n), m1(n);
  for(int k=0;k<p.folds;++k) {
    std::vector<int> train, arms[2];
    for(int i=0;i<n;++i) if(p.folds==1 || fold[i]!=k) {
      train.push_back(i); arms[static_cast<int>(z[i])].push_back(i);
    }
    Vec prob,q0,q1;
    if(learned) {
      const std::string family=p.outcomeModel=="bernoulli" ? "bernoulli" : "gaussian";
      const auto a=super_learner(d,pv,z,train,"bernoulli",groups), b=super_learner(stacked,ov,stackedY,train,family,stackedGroups);
      prob=a.prediction;
      q0.assign(b.prediction.begin()+n,b.prediction.begin()+2*n); q1.assign(b.prediction.begin()+2*n,b.prediction.end());
      for(std::size_t l=0;l<LEARNERS.size();++l) {
        propensityWeights[l]+=a.weights[l]/p.folds; outcomeWeights[l]+=b.weights[l]/p.folds;
      }
    } else {
      prob=fitted_means(pm,z,train,"bernoulli");
      q0=fitted_means(om,y,arms[0],p.outcomeModel); q1=fitted_means(om,y,arms[1],p.outcomeModel);
    }
    for(int i=0;i<n;++i) if(p.folds==1 || fold[i]==k) { pi[i]=prob[i]; m0[i]=q0[i]; m1[i]=q1[i]; }
  }
  const double bound=p.specification.get("propensityBound").as_number();
  int truncated=0, clipped=0;
  double minProp=1,maxProp=0,maxWeight=0,maxScore=0;
  for(int i=0;i<n;++i) {
    minProp=std::min(minProp,pi[i]); maxProp=std::max(maxProp,pi[i]);
    const double raw=pi[i]; pi[i]=std::clamp(pi[i],bound,1-bound); truncated+=raw!=pi[i];
    if(!(pi[i]>0 && pi[i]<1)) throw Error("positivity","A fitted propensity equals 0 or 1. Inspect overlap or explicitly set a propensity bound.");
    maxWeight=std::max(maxWeight,std::max(1/pi[i],1/(1-pi[i])));
  }
  if(p.method=="tmle") {
    const auto expit=[](double x) { return x>=0 ? 1/(1+std::exp(-x)) : std::exp(x)/(1+std::exp(x)); };
    for(int arm=0;arm<2;++arm) {
      auto &q=arm ? m1 : m0;
      Vec offsets(n),h(n);
      for(int i=0;i<n;++i) {
        h[i]=1/(arm ? pi[i] : 1-pi[i]);
        if(p.bounded) {
          const double raw=q[i]; q[i]=std::clamp(q[i],1e-6,1-1e-6); clipped+=raw!=q[i];
          offsets[i]=std::log(q[i])-std::log1p(-q[i]);
        }
      }
      if(p.bounded) {
        const auto score=[&](double eps) {
          double value=0; for(int i:rows[arm]) value+=h[i]*(y[i]-expit(offsets[i]+eps*h[i])); return value;
        };
        double lo=-1,hi=1;
        for(int j=0;j<50 && score(lo)<0;++j) lo*=2;
        for(int j=0;j<50 && score(hi)>0;++j) hi*=2;
        if(score(lo)<0 || score(hi)>0) throw Error("nonconvergence","The TMLE targeting step did not bracket a solution.");
        for(int j=0;j<80;++j) { const double mid=(lo+hi)/2; if(score(mid)>0) lo=mid; else hi=mid; }
        for(int i=0;i<n;++i) q[i]=expit(offsets[i]+(lo+hi)/2*h[i]);
      } else {
        double score=0,information=0;
        for(int i:rows[arm]) { score+=h[i]*(y[i]-q[i]); information+=h[i]*h[i]; }
        const double epsilon=score/information;
        for(int i=0;i<n;++i) q[i]+=epsilon*h[i];
      }
      double score=0; for(int i:rows[arm]) score+=h[i]*(y[i]-q[i]);
      maxScore=std::max(maxScore,std::fabs(score)/n);
    }
  }
  if(p.bounded) for(int i=0;i<n;++i) {
    m0[i]=p.lower+(p.upper-p.lower)*m0[i]; m1[i]=p.lower+(p.upper-p.lower)*m1[i];
  }
  y=d.values[p.outcome];
  Vec phi0(n),phi1(n),influence(n);
  for(int i=0;i<n;++i) {
    phi1[i]=m1[i]+z[i]*(y[i]-m1[i])/pi[i];
    phi0[i]=m0[i]+(1-z[i])*(y[i]-m0[i])/(1-pi[i]);
  }
  const double psi0=moments(p.method=="tmle" ? m0 : phi0).mean,
               psi1=moments(p.method=="tmle" ? m1 : phi1).mean;
  const auto measure=e.get("measure").as_string();
  Analysis out;
  if(!positive(measure)) {
    out.estimate=out.testEstimate=psi1-psi0;
    for(int i=0;i<n;++i) influence[i]=phi1[i]-phi0[i]-out.estimate;
  } else {
    const bool rr=measure=="risk_ratio";
    if(!(psi0>0 && psi1>0) || (!rr && !(psi0<1 && psi1<1)))
      throw Error("boundary_response","An estimated arm mean is outside the probability range for this ratio.");
    const double d0=rr ? psi0 : psi0*(1-psi0), d1=rr ? psi1 : psi1*(1-psi1);
    for(int i=0;i<n;++i) influence[i]=(phi1[i]-psi1)/d1-(phi0[i]-psi0)/d0;
    out.estimate=rr ? psi1/psi0 : (psi1/(1-psi1))/(psi0/(1-psi0)); out.testEstimate=std::log(out.estimate);
  }
  if(clustered) {
    Vec sums(units);
    for(int i=0;i<n;++i) sums[d.cluster[i]]+=influence[i]/n;
    double variance=0; for(double v:sums) variance+=v*v;
    out.se=std::sqrt(variance*units/(units-1.0)); out.df=units-2;
  } else out.se=std::sqrt(moments(influence).variance/n);
  out.diagnostics=object({{"propensityLearners",learned ? array(propensityWeights) : J()},{"outcomeLearners",learned ? array(outcomeWeights) : J()},{"minPropensity",minProp},{"maxPropensity",maxProp},{"propensityTruncated",truncated},
    {"maxInverseWeight",maxWeight},{"initialPredictionsClipped",clipped},{"targetingScore",maxScore},
    {"independentUnits",units},{"subjects",n},{"referenceMean",psi0},{"activeMean",psi1}});
  return out;
}
Data complete_summary(const Graph &g, const J &e, const Data &d, const Plan &p) {
  const auto w=study_weights(g,e); const int visits=static_cast<int>(w.size());
  Data out; out.n=out.subjects=d.subjects; out.units=d.units;
  out.values.resize(g.nodes.size(),Vec(out.n)); out.cluster.resize(out.n);
  for(int i=0;i<out.n;++i) {
    const int first=i*visits; out.cluster[i]=d.cluster[first];
    for(int j=0;j<static_cast<int>(g.nodes.size());++j)
      if(g.nodes[j].scope=="subject" || g.nodes[j].scope=="cluster") out.values[j][i]=d.values[j][first];
    double value=0;
    for(int v=0;v<visits;++v) if(w[v]!=0) {
      if(!d.observed[p.outcome][first+v]) throw Error("missing_outcome","Trajectory TMLE requires every outcome contributing to the summary to be observed.");
      value+=w[v]*d.values[p.outcome][first+v];
    }
    out.values[p.outcome][i]=value;
  }
  return out;
}
Analysis analyze(const Graph &g, const J &e, const Data &d, const Plan &p) {
  Analysis out;
  const auto measure = e.get("measure").as_string();
  if(g.study.enabled && p.method=="tmle") out=analyze_aipw(e,complete_summary(g,e,d,p),p);
  else if(g.study.enabled) {
    const auto f=analyze_study(g,e,d,p.specification);
    out.estimate=out.testEstimate=f.coefficient; out.se=f.se; out.df=f.df;
  }
  else if (p.method == "aipw" || p.method == "tmle")
    out = analyze_aipw(e, d, p);
  else if(p.method=="joint") {
    const auto m=model_matrix(d,p.terms,true);
    ModelFit model;
    if(g.nodes[p.outcome].kind=="survival")
      model=fit_weibull(m,d.times[p.outcome],d.events[p.outcome],p.family=="exponential");
    else {
      Vec y=d.values[p.outcome];
      if(p.family=="lognormal") for(double &v:y) v=std::log(v);
      model=fit_model(m,y,model_family(p.family),true);
    }
    out.statistic=wald_statistic(model,p.tested);
    out.jointP=stats::pchisq(out.statistic,static_cast<double>(p.tested.size()),false);
    return out;
  }
  else if (p.method == "gformula" && p.chain.empty()) {
    // Zero by construction: an estimate without spread, never rejected.
    out.estimate = positive(measure) ? 1 : 0;
    return out;
  } else if (p.method == "gformula") {
    const auto r = gformula(
        g, d, p.chain, p.exposure, e.get("reference").as_number(),
        e.get("active").as_number(), measure,
        e.get("horizon").is_null() ? 0 : e.get("horizon").as_number());
    out.estimate = r.estimate;
    out.testEstimate = r.testEstimate;
    out.se = r.se;
  }
  else if (e.get("kind").as_string() == "coefficient") {
    const auto m = model_matrix(d, p.terms, p.method != "cox");
    const int target = p.targetIndex + (p.method == "cox" ? 0 : 1);
    Fit fit;
    if (p.method == "cox")
      fit = fit_cox(m, d.times[p.outcome], d.events[p.outcome], target);
    else if (p.method == "glm" && quasi(p.family)) {
      const auto model = fit_model(m, d.values[p.outcome], model_family(p.family), true);
      fit.coefficient = model.beta[target] / m.scales[target];
      fit.se = std::sqrt(model.covariance[target * m.p + target]) / m.scales[target];
    } else if (p.method == "glm")
      fit = fit_glm(m, d.values[p.outcome], target, model_family(p.family));
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
  return object({{"estimate", a.jointP >= 0 ? J() : J(a.estimate)},
                 {"testEstimate", a.jointP >= 0 ? J() : J(a.testEstimate)},
                 {"se", a.jointP >= 0 ? J() : J(a.se)},
                 {"df", std::isfinite(a.df) ? J(a.df) : J()},
                 {"statistic", a.statistic},{"diagnostics",a.diagnostics}});
}
J population(const Graph &g, const J &e, int n, std::uint64_t seed,
             const Token *token) {
  validate_estimand(g, e);
  if(g.study.enabled) return study_population(g,e,n,seed,token);
  const auto kind = e.get("kind").as_string(),
             measure = e.get("measure").as_string();
  const int a = g.index(e.get("exposure").as_string()),
            y = g.index(e.get("outcome").as_string());
  if(kind=="joint") {
    J coefficients=J::array(); bool null=true;
    for(const auto &t:g.nodes[y].terms)
      if(joint_term(t,a,measure)) {
        coefficients.push(t.coefficient); null &= t.coefficient==0;
      }
    return object({{"estimate",J()},{"mcSE",0},{"n",0},{"exact",true},
      {"jointNull",null},{"coefficients",coefficients},
      {"description","Generating coefficients of the tested predictor block; fitted-model hypotheses depend on the specified analysis basis."}});
  }
  if (kind == "coefficient") {
    double beta = 0;
    for (const auto &term : g.nodes[y].terms)
      if (term.left == a && term.right < 0 && term.transform.kind == "identity")
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
  J json, analysis;
  std::vector<int> sizes;
  int repetitions = 500, populationSize = 20000;
  std::uint64_t seed = 2026;
  double alpha = .05, targetPower = .8, lower = -3, upper = 3;
  std::string side = "two-sided", edge, nullMode = "no_effect";
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
  s.analysis=parse_analysis(raw.get("analysis"));
  if(raw.find("nullMode")) s.nullMode=raw.get("nullMode").as_string();
  if(s.nullMode!="no_effect" && s.nullMode!="average_effect") fail("Choose a no-effect or average-effect null.");
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
                   {"targetPower", s.targetPower}, {"analysis",s.analysis}, {"nullMode",s.nullMode}});
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
  if (chosen < 0) {
    // No arrow moves the target: the run keeps the drawn DGP and says so.
    auto out = answer(g.json, effect(g, verification), J(), 0);
    out.set("target", J());
    out.set("unreachable", target);
    return out;
  }
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
  const auto plan = compile_plan(g, e, s.analysis);
  if(e.get("kind").as_string()=="joint" && s.side!="two-sided") fail("Joint tests use an omnibus two-sided hypothesis.");
  if (!s.edge.empty() &&
      std::none_of(plan.edges.begin(), plan.edges.end(), [&](int i) {
        return edge_key(g.json.get("edges").at(i)) == s.edge;
      }))
    fail("The selected calibration arrow no longer belongs to this estimand. "
         "Choose it again.");
  double rows=0;
  for(int n:s.sizes) rows+=g.study.enabled ? study_rows(g,n) : n;
  const double draws = rows *
                       s.repetitions * g.nodes.size() * (s.checkNull ? 2 : 1);
  if (draws > 200000000)
    fail("This run exceeds the browser work limit. Reduce sample sizes, "
         "repetitions, or the number of variables.");
  // Bound concurrent workspaces, especially for wide regressions. The
  // calculator's pool is still shared.
  const double largest = g.study.enabled ? study_rows(g,s.sizes.back()) : s.sizes.back();
  // A Super Learner fit holds, for the data stacked three times over: the
  // main-term, product (at most 40 terms) and spline design matrices, every
  // learner's cross-validated and refitted predictions, 48 kernel features,
  // and the booster's bins and gradients, about 110 + 5 values per variable.
  double learnerBytes = 0;
  if(plan.learner=="superlearner") {
    std::set<int> variables={plan.exposure};
    for(const auto *terms:{&plan.terms,&plan.propensity}) for(const auto &t:*terms) { variables.insert(t.left); if(t.right>=0) variables.insert(t.right); }
    learnerBytes = 8.0 * 3 * largest * (110 + 5.0 * variables.size());
  }
  const double terms = plan.terms.size() + 1,
               bytes = 8 * (largest * (g.nodes.size() * 4 + terms + 10) +
                            8 * terms * terms) + (plan.method=="lmm" ? 2000.0*s.sizes.back() : 0) + learnerBytes;
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
    if(e.get("kind").as_string()=="joint") {
      J ng=g.json, nodes=J::array(), edges=J::array();
      const auto a=e.get("exposure").as_string(), y=e.get("outcome").as_string();
      const bool overall=e.get("measure").as_string()=="overall_effect";
      for(J edge:g.json.get("edges").elements()) {
        if(overall && edge.get("source").as_string()==a && edge.get("target").as_string()==y) edge.set("coefficient",0);
        edges.push(edge);
      }
      for(J node:g.json.get("nodes").elements()) {
        if(node.get("id").as_string()==y) for(const auto &key:{"interactions","functions"}) {
          if(!node.find(key)) continue;
          J terms=J::array();
          for(J term:node.get(key).elements()) {
            const auto t=parse_term(g,term);
            if(joint_term(t,plan.exposure,e.get("measure").as_string())) term.set("coefficient",0);
            terms.push(term);
          }
          node.set(key,terms);
        }
        nodes.push(node);
      }
      ng.set("nodes",nodes); ng.set("edges",edges);
      nullScenario=object({{"graph",ng},{"truth",population(parse_graph(ng),e,s.populationSize,truthSeed,token)},
        {"exactNull",true},{"construction","The tested generating coefficient block is zeroed; other outcome terms and paths remain. Under misspecified analysis, this need not be a null for the working-model coefficients."}});
    } else if (s.nullMode=="no_effect" && e.get("kind").as_string() == "marginal" &&
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
              : "One arrow is Monte Carlo calibrated to the null contrast while other terms, including heterogeneity, remain. The verified effect and its Monte Carlo uncertainty must be considered when interpreting rejection.");
    }
  }
  J warnings = J::array(), results = J::array(), first, notes = J::array();
  for (const auto &text : plan.notes) {
    warnings.push(text);
    notes.push(text);
  }
  for (const auto *scenario : {&alternative, &nullScenario})
    if (!scenario->get("unreachable").is_null())
      warnings.push(
          "No arrow lies on a path from " + e.get("exposure").as_string() +
          " to " + e.get("outcome").as_string() + ", so the " +
          (scenario == &alternative ? "planning effect " : "null value ") +
          scenario->get("unreachable").to_text() +
          " cannot be set; the run keeps the drawn model, whose effect is " +
          J(scenario->get("truth").get("estimate").as_number()).to_text() +
          ".");
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
    if(g.study.enabled) {
      row.set("measurementRows",study_rows(g,n)); row.set("subjects",study_rows(g,n)/static_cast<int>(g.study.times.size()));
      row.set("sizeUnit",g.study.clusterSizes.empty()?"subjects":"clusters");
    }
    for (int scenario = 0; scenario < scenarios; ++scenario) {
      const auto runGraph =
          parse_graph((scenario ? nullScenario : alternative).get("graph"));
      sim::Accumulator accumulator;
      std::map<std::string, int> reasons;
      int coverageCount=0, truncatedCount=0, clippedCount=0, learnedCount=0;
      double largestWeight=0;
      Vec propensityLearners(LEARNERS.size()), outcomeLearners(LEARNERS.size());
      const bool scalar=e.get("kind").as_string()!="joint";
      const double truth=(scenario ? nullScenario : alternative).get("truth").get("estimate").as_number();
      // Batches keep progress responsive and aggregation in replicate order,
      // independent of thread scheduling.
      constexpr int batchSize = 64;
      for (int begin = 0; begin < s.repetitions; begin += batchSize) {
        const int count = std::min(batchSize, s.repetitions - begin);
        std::vector<sim::Outcome> outcomes(count);
        std::vector<std::string> failures(count);
        std::vector<int> covered(count);
        std::vector<J> diagnostics(count);
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
                out.reject = (fit.jointP>=0 ? fit.jointP : p_value(fit.statistic, fit.df, s.side)) < s.alpha;
                out.estimate = fit.estimate;
                out.se = fit.se;
                out.df = fit.df; diagnostics[i]=fit.diagnostics;
                if(scalar) {
                  const double center=positive(e.get("measure").as_string()) ? std::log(truth) : truth;
                  const double critical=std::isfinite(fit.df) ? stats::qt(1-s.alpha/2,fit.df) : stats::qnorm(1-s.alpha/2);
                  covered[i]=std::isfinite(center) && std::fabs(fit.testEstimate-center)<=critical*fit.se;
                }
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
          coverageCount+=covered[i];
          truncatedCount+=static_cast<int>(diagnostics[i].get("propensityTruncated").as_number());
          clippedCount+=static_cast<int>(diagnostics[i].get("initialPredictionsClipped").as_number());
          largestWeight=std::max(largestWeight,diagnostics[i].get("maxInverseWeight").as_number());
          if(diagnostics[i].is_object() && diagnostics[i].get("outcomeLearners").is_array()) {
            ++learnedCount;
            for(std::size_t l=0;l<LEARNERS.size();++l) {
              propensityLearners[l]+=diagnostics[i].get("propensityLearners").at(l).as_number();
              outcomeLearners[l]+=diagnostics[i].get("outcomeLearners").at(l).as_number();
            }
          }
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
               {"meanEstimate", scalar ? J(summary.mean_estimate) : J()},
               {"estimateSD", scalar ? J(summary.sd_estimate) : J()},
               {"meanSE", scalar ? J(summary.mean_se) : J()},
               {"failureReasons", failures},
               {"targetedDiagnostics",(plan.method=="tmle" || plan.method=="aipw") ? [&] {
                 J d=object({{"maxInverseWeight",largestWeight},{"propensityTruncated",truncatedCount},{"initialPredictionsClipped",clippedCount}});
                 if(learnedCount) {
                   J names=J::array(); for(const auto &name:LEARNERS) names.push(name);
                   for(auto &w:propensityLearners) w/=learnedCount;
                   for(auto &w:outcomeLearners) w/=learnedCount;
                   d.set("learners",object({{"names",names},{"propensity",array(propensityLearners)},{"outcome",array(outcomeLearners)}}));
                 }
                 return d; }() : J()},
                {"bias",scalar && summary.reps>summary.failures ? J(summary.mean_estimate-truth) : J()},
                {"coverage",scalar && summary.reps>summary.failures ? J(static_cast<double>(coverageCount)/(summary.reps-summary.failures)) : J()},
                {"coverageValid",summary.reps-summary.failures}}));
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
                            {"testScale", plan.testScale},
                            {"notes", notes}, {"specification",s.analysis}})},
       {"alternative", alternative},
       {"nullScenario", nullScenario},
       {"results", results},
       {"firstTestedN", first},
       {"sizeUnit",g.study.enabled && !g.study.clusterSizes.empty()?"clusters":"subjects"},
       {"assumptions","Power is conditional on the specified data-generating scenario. Monte Carlo intervals quantify simulation error only; they do not quantify uncertainty about the generating model."},
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
    return plan_json(graph, compile_plan(graph, e, request.get("analysis")));
  if (op == "analyze") {
    const auto fit = analyze(graph, e, read_data(graph, request.get("data")),
                             compile_plan(graph, e, request.get("analysis")));
    J out = analysis_json(fit);
    const auto side = request.get("side").as_string();
    out.set("pValue",
            fit.jointP>=0 ? fit.jointP : p_value(fit.statistic, fit.df, side.empty() ? "two-sided" : side));
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

