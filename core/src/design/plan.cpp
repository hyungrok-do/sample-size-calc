#include "plan.hpp"
#include "program.hpp"
#include "stats/stats.hpp"
#include "util/numtext.hpp"
#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <string>

namespace ssc::design {
namespace {
using J = json::Value;
J object(std::initializer_list<std::pair<std::string, J>> values) {
  J out=J::object(); for (const auto& value:values) out.set(value.first,value.second); return out;
}
[[noreturn]] void fail(const std::string& message) {
  throw std::invalid_argument("Design plan: " + message);
}
void keys(const J& value, std::initializer_list<const char*> allowed) {
  if (!value.is_object()) fail("expected an object.");
  for (const auto& kv : value.items()) {
    bool found = false;
    for (const char* key : allowed) found |= kv.first == key;
    if (!found) fail("unsupported field '" + kv.first + "'.");
  }
}
std::string label(const J& value, const std::string& key) {
  const auto& item = value.get(key);
  if (!item.is_string() || item.as_string().empty() || item.as_string().size() > 4096)
    fail("'" + key + "' must be a nonempty string (at most 4096 bytes).");
  return item.as_string();
}
// Lean exports rationals as strings to avoid an implicit decimal rounding step.
// The numerical backend still evaluates them in binary64, not interval arithmetic.
double scalar(const J& value, const std::string& key) {
  const auto& item = value.get(key);
  double result = std::numeric_limits<double>::quiet_NaN();
  if (item.is_number()) result = item.as_number();
  else if (item.is_string()) {
    const auto& s = item.as_string();
    const auto slash = s.find('/');
    const auto integer = [](const std::string& x, bool signedValue) {
      std::size_t start = signedValue && !x.empty() && x[0] == '-' ? 1 : 0;
      return x.size() > start && x.size() <= 17 &&
        x.find_first_not_of("0123456789", start) == std::string::npos;
    };
    if (slash != std::string::npos && integer(s.substr(0, slash), true) &&
        integer(s.substr(slash + 1), false)) {
      double numerator, denominator;
      if (text::parse_double(s.substr(0, slash), numerator) &&
          text::parse_double(s.substr(slash + 1), denominator) &&
          std::fabs(numerator) <= 9007199254740991.0 &&
          denominator > 0 && denominator <= 9007199254740991.0)
        result = numerator / denominator;
    }
  }
  if (!std::isfinite(result)) fail("'" + key + "' must be finite numeric data or an integer/positive-integer rational.");
  return result;
}
int count(const J& value, const char* key, int lo, int hi) {
  const auto& raw = value.get(key);
  const double x = raw.as_number(-1);
  if (!raw.is_number() || !std::isfinite(x) || x != std::floor(x) || x < lo || x > hi)
    fail(std::string(key) + " must be an integer in [" + std::to_string(lo) + ", " + std::to_string(hi) + "].");
  return static_cast<int>(x);
}
J strings(std::initializer_list<const char*> values) {
  J a = J::array(); for (const char* s : values) a.push(s); return a;
}
J verification(const char* basis) {
  return object({{"basis", basis}, {"numerics", "binary64_not_certified"},
    {"inputProofAuthenticated", false},
    {"formalScope", "See lean/DesignCompiler: typed lowering and finite-search theorems. JSON alone does not carry a checked proof."},
    {"notVerified", strings({"C++/WASM implementation refinement", "floating-point error bounds",
                             "data-generating assumptions", "estimator asymptotic regularity"})}});
}
} // namespace

J execute(const J& p, Token* token) {
  if (p.get("schema").as_string() == "ssc.expression-plan") return execute_program(p, token);
  if (token) token->check();
  if (p.get("schema").as_string() != "ssc.design-plan" ||
      !p.get("version").is_number() || p.get("version").as_number() != 1)
    fail("expected ssc.design-plan version 1.");
  const auto method = label(p,"method");
  if (method != "normal_known" && method != "influence_wald" && method != "ate_eif")
    fail("unsupported method '" + method + "'.");
  keys(p, {"schema", "version", "method", "name", "parameters", "search", "solve"});
  const auto name = label(p,"name");
  const auto& x = p.get("parameters");
  if (method == "normal_known")
    keys(x, {"alternative", "nullValue", "sd", "alpha", "targetPower"});
  else if (method == "influence_wald")
    keys(x, {"alternative", "nullValue", "variance", "alpha", "targetPower", "independentUnit", "estimand", "estimator", "varianceSource"});
  else keys(x, {"nullValue", "strata", "alpha", "targetPower", "estimator"});
  const double alpha = scalar(x,"alpha"), target = scalar(x,"targetPower");
  if (!(alpha >= 1e-12 && alpha <= 1 - 1e-12 && target > alpha && target < 1))
    fail("require 1e-12 <= alpha <= 1-1e-12 and alpha < targetPower < 1.");
  const double nullValue = scalar(x,"nullValue");
  double alternative = 0, variance = 0;
  std::string unit = "participants_per_arm", estimator = "known_variance_z", estimand = "mean_difference";
  J assumptions;
  if (method == "normal_known") {
    alternative = scalar(x,"alternative");
    const double sd = scalar(x,"sd");
    if (!(sd > 0)) fail("sd must be strictly positive.");
    variance = 2 * sd * sd;
    assumptions = strings({"Independent normal outcomes", "Common known SD in both arms", "Equal fixed allocation", "Two-sided z test"});
  } else if (method == "influence_wald") {
    alternative = scalar(x,"alternative"); variance = scalar(x,"variance");
    unit = label(x,"independentUnit");
    if (unit != "participant" && unit != "cluster") fail("independentUnit must be participant or cluster.");
    estimator = label(x,"estimator"); estimand = label(x,"estimand");
    label(x,"varianceSource");
    assumptions = strings({"Supplied variance is per independent unit for the planned estimator", "Asymptotic linearity and a negligible remainder", "Consistent standard error on the supplied effect scale", "Normal approximation has no finite-sample guarantee"});
  } else {
    unit = "participant"; estimand = "ate"; estimator = label(x,"estimator");
    if (estimator != "aipw" && estimator != "tmle") fail("ATE EIF supports aipw or tmle.");
    const auto& strata = x.get("strata");
    if (!strata.is_array() || strata.size() == 0 || strata.size() > 10000)
      fail("ATE EIF needs 1–10000 strata.");
    double mass = 0;
    for (const auto& s : strata.elements()) {
      keys(s, {"weight", "propensity", "mean0", "mean1", "variance0", "variance1"});
      const double w = scalar(s,"weight"), g = scalar(s,"propensity");
      const double v0 = scalar(s,"variance0"), v1 = scalar(s,"variance1");
      if (!(w > 0 && g > 0 && g < 1 && v0 >= 0 && v1 >= 0))
        fail("stratum weights must be positive, propensities in (0,1), and conditional variances nonnegative.");
      mass += w;
      alternative += w * (scalar(s,"mean1") - scalar(s,"mean0"));
    }
    if (std::fabs(mass - 1) > 1e-12) fail("stratum weights must sum to one.");
    for (const auto& s : strata.elements()) {
      const double w = scalar(s,"weight"), g = scalar(s,"propensity");
      const double d = scalar(s,"mean1") - scalar(s,"mean0") - alternative;
      variance += w * (scalar(s,"variance1") / g + scalar(s,"variance0") / (1 - g) + d*d);
    }
    assumptions = strings({"Finite baseline strata with positivity", "Supplied conditional outcome means and variances", "Independent participants", "Planned AIPW/TMLE attains this influence function under its nuisance and remainder conditions", "Causal interpretation additionally requires consistency and exchangeability", "Normal approximation has no finite-sample guarantee"});
  }
  const double effect = alternative - nullValue;
  if (!(variance > 0) || !std::isfinite(variance) || !std::isfinite(effect))
    fail("the derived effect and positive variance must be finite; rescale inputs.");
  const auto& search = p.get("search");
  keys(search, {"minN", "maxN"});
  const int minimum = count(search,"minN",1,10000000), maximum = count(search,"maxN",1,10000000);
  if (minimum > maximum) fail("minN must not exceed maxN.");
  const double z = stats::qnorm(alpha / 2, false);
  const auto power = [&](int n) {
    const double lambda = std::fabs(effect) / std::sqrt(variance) * std::sqrt(static_cast<double>(n));
    return stats::pnorm(-z - lambda) + stats::pnorm(z - lambda, false);
  };
  const auto& solve = p.get("solve");
  std::string solveFor = "sample_size";
  if (!solve.is_null()) {
    keys(solve, {"for", "n"}); solveFor = label(solve,"for");
  }
  if (solveFor != "sample_size" && solveFor != "power" && solveFor != "effect")
    fail("solve.for must be sample_size, power, or effect.");
  if (solveFor == "sample_size" && solve.find("n")) fail("n is an output when solving for sample_size.");
  if (solveFor != "sample_size") {
    if (method == "ate_eif" && solveFor == "effect")
      fail("change the ATE planning strata explicitly; an EIF variance can change with the effect.");
    const int n = count(solve,"n",1,10000000);
    J out=object({{"name",name},{"method",method},{"solveFor",solveFor},{"n",n},{"sizeUnit",unit},
      {"alpha",alpha},{"targetPower",target},{"varianceCoefficient",variance},{"assumptions",assumptions},
      {"verification",verification(method=="normal_known" ? "exact_distributional_formula" : "asymptotic_approximation")}});
    if (method=="normal_known") out.set("totalParticipants",2*n);
    if (solveFor=="power") {
      const double value=power(n);
      if(!std::isfinite(value)) fail("non-finite power evaluation.");
      out.set("status","evaluated"); out.set("power",value);
      out.set("alternative",alternative); out.set("nullValue",nullValue);
    } else {
      // Bisect the noncentrality, whose upper endpoint has effectively unit power
      // over the accepted alpha range. Report a bracket, not a certified exact root.
      const auto at=[&](double lambda) {return stats::pnorm(-z-lambda)+stats::pnorm(z-lambda,false);};
      double lo=0,hi=64;
      for(int i=0;i<100;++i) {
        const double mid=lo+(hi-lo)/2;
        if(mid==lo || mid==hi) break;
        if(at(mid)>=target) hi=mid; else lo=mid;
      }
      const double scale=std::sqrt(variance)/std::sqrt(static_cast<double>(n));
      const double upper=hi*scale,lower=lo*scale;
      if(!std::isfinite(upper) || !std::isfinite(nullValue+upper) || !std::isfinite(nullValue-upper))
        fail("detectable effect overflowed; rescale inputs.");
      out.set("status","approximated"); out.set("detectableMagnitude",upper);
      out.set("effectLower",lower); out.set("effectUpper",upper);
      out.set("powerAtLower",at(lo)); out.set("powerAtUpper",at(hi));
      out.set("alternativeAbove",nullValue+upper); out.set("alternativeBelow",nullValue-upper);
      out.set("effectConvention","Absolute alternative-minus-null difference, holding the supplied variance fixed.");
    }
    return out;
  }
  constexpr double comparisonGuard = 1e-12; // ambiguity guard, NOT a proved error bound
  int selected = 0, unresolved = 0, evaluated = 0;
  double lastPower = 0, previousPower = 0;
  if (token) token->set_total(maximum - minimum + 1);
  for (int n = minimum; n <= maximum; ++n) {
    if (token && (n - minimum) % 1024 == 0) token->check();
    previousPower = lastPower; lastPower = power(n); ++evaluated;
    if (!std::isfinite(lastPower)) fail("non-finite power evaluation.");
    if (token) token->add_done();
    if (std::fabs(lastPower - target) <= comparisonGuard) { unresolved = n; break; }
    if (lastPower > target) { selected = n; break; }
    if (effect == 0) break; // constant power alpha: no solution in the domain
  }
  J out = object({{"name", name}, {"method", method},
    {"solveFor","sample_size"},
    {"status", selected ? "found" : unresolved ? "numerically_unresolved" : "not_found_in_domain"},
    {"selectedN", selected ? J(selected) : J()}, {"sizeUnit", unit},
    {"power", lastPower}, {"targetPower", target}, {"alpha", alpha},
    {"alternative", alternative}, {"nullValue", nullValue}, {"varianceCoefficient", variance},
    {"estimator", estimator}, {"estimand", estimand}, {"assumptions", assumptions},
    {"verification", verification(method == "normal_known" ? "exact_distributional_formula" : "asymptotic_approximation")},
    {"search", object({{"minN", minimum}, {"maxN", maximum}, {"evaluated", evaluated},
      {"algorithm", "ascending_exhaustive"}, {"comparisonGuard", comparisonGuard},
      {"unresolvedN", unresolved ? J(unresolved) : J()},
      {"minimality", selected ? "all_smaller_admissible_counts_fail_in_binary64" : "not_established"},
      {"predecessorN", selected > minimum ? J(selected - 1) : J()},
      {"predecessorPower", selected > minimum ? J(previousPower) : J()}})}});
  if (method == "normal_known") out.set("totalParticipants", selected ? J(2 * selected) : J());
  if (method == "influence_wald") out.set("varianceSource", x.get("varianceSource"));
  return out;
}
} // namespace ssc::design
