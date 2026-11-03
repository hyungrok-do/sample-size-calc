#include "design/plan.hpp"
#include "testing.hpp"
#define SSC_STATIC
#include "ssc/ssc.h"
#include <cmath>

namespace {
using J = ssc::json::Value;
J request(const std::string &source) {
  J r = J::object();
  r.set("op", "compile");
  r.set("source", source);
  return r;
}
J api(const J &request) {
  char *out = ssc_design_json(request.dump().c_str(), nullptr);
  if (!out)
    throw std::runtime_error("Missing compiler response");
  std::string text(out);
  ssc_free(out);
  return ssc::json::parse(text);
}
J compile(const std::string &source) {
  auto a = api(request(source));
  if (!a.get("ok").as_bool())
    throw std::runtime_error(a.get("error").as_string());
  return a.get("result");
}
J solve(const std::string &source) {
  auto a = api(compile(source).get("plan"));
  if (!a.get("ok").as_bool())
    throw std::runtime_error(a.get("error").as_string());
  return a.get("result");
}
void invalid(const std::string &source) {
  auto a = api(request(source));
  if (a.get("ok").as_bool())
    a = api(a.get("result").get("plan"));
  CHECK(!a.get("ok").as_bool());
}
} // namespace

TEST_CASE("program: typed functions, scopes, arrays and headerless queries") {
  auto r = solve(R"(
-- NormalInput in a comment must not select the old compiler.
def f (x : Real) : Real := let y := x + 1; y * y
def g (x : Real) : Real := f(x) + f(x + 1)
def a : Array Real := map([1, 2, 3], fun x => g(x))
evaluate a[1] + mean(a) + variance([1, 2, 3])
)");
  CHECK_CLOSE(r.get("value").as_number(), 25 + 79.0 / 3 + 1, 0, 1e-12);
  CHECK_EQ(solve("evaluate -2^2 + 2^3^2").get("value").as_number(), 508);
  CHECK_EQ(solve("import DesignCompiler.Study\nopen DesignCompiler\nevaluate "
                 "Real.sqrt 81")
               .get("value")
               .as_number(),
           9);
  CHECK_EQ(solve("def μ : ℝ := 2\nsolve for n : ℕ ∈ [1, 5] where n≥μ+.1")
               .get("solution")
               .as_number(),
           3);
  CHECK_EQ(
      solve("evaluate if false then 1 / 0 else 3").get("value").as_number(), 3);
}
TEST_CASE(
    "program: algebraic roots, numerical residuals and honest failure states") {
  auto r = solve("def f (x : Real) := x^2 + 2*x\nsolve for x : Real in [0,10] "
                 "where f(x) = 8");
  CHECK_EQ(r.get("solution").as_number(), 2);
  CHECK_EQ(r.get("method").as_string(), "algebraic_quadratic");
  CHECK_EQ(solve("solve for x : Real in [-10,10] where x*x = 4")
               .get("solution")
               .as_number(),
           -2);
  r = solve("solve for x : Real in [0,2] where exp(x) = 2");
  CHECK_CLOSE(r.get("solution").as_number(), std::log(2.0), 0, 1e-8);
  CHECK(std::abs(r.get("residual").as_number()) <= 1e-8);
  r = solve("solve for x : Real in [-1,1] where (if x < 0 then -1 else 1) = 0");
  CHECK_EQ(r.get("status").as_string(), "numerically_unresolved");
  CHECK_EQ(solve("solve for x : Real in [0,1] where exp(x) = 4")
               .get("status")
               .as_string(),
           "root_not_bracketed");
  CHECK_EQ(solve("solve for x : Grid in [0,1] step .1 where x > .35")
               .get("solution")
               .as_number(),
           .4);
  CHECK_EQ(solve("solve for x : Nat in [1,3] where x >= 2")
               .get("solution")
               .as_number(),
           2);
  CHECK_EQ(
      solve("solve for x : Nat in [1,3] where x > 4").get("status").as_string(),
      "not_found_in_domain");
  r = solve("evaluate integrate(fun x => x*x, 0, 1, 100)");
  CHECK_CLOSE(r.get("value").as_number(), 1.0 / 3, 0, 1e-14);
}
TEST_CASE("program: Gaussian lowering preserves shared and independent draws") {
  const std::string power = R"(
def power (n : Real) : Real :=
  let z := normal(.5 * sqrt(n/2), 1)
  probability(abs(z) > normalQuantile(.975))
solve for n : Nat in [2,1000] where power(n) >= .8
)";
  auto r = solve(power);
  CHECK_EQ(r.get("solution").as_number(), 63);
  CHECK_CLOSE(r.get("value").as_number(), .8013023941055797, 0, 2e-13);
  CHECK_EQ(r.get("evaluation").as_string(), "analytic_gaussian_probability");
  r = solve("evaluate let x := normal(0,1); probability(x + x > 2)");
  CHECK_CLOSE(r.get("value").as_number(), .5 * std::erfc(1 / std::sqrt(2.0)), 0,
              1e-14);
  r = solve("evaluate probability(normal(0,1) + normal(0,1) > 2)");
  CHECK_CLOSE(r.get("value").as_number(), .5 * std::erfc(1.0), 0, 1e-14);
  CHECK_EQ(solve("evaluate let x := normal(0,1); probability(x-x > 0)")
               .get("value")
               .as_number(),
           0);
  CHECK_EQ(
      solve("evaluate probability(normal(2,0) >= 2)").get("value").as_number(),
      1);
}
TEST_CASE("program: Monte Carlo preserves call-by-value and reports "
          "simultaneous bounds") {
  const std::string config = "simulation { replicates := 1000, seed := 12, "
                             "confidence := .95, method := simulation }\n";
  auto r = solve("def f (x : Real) : Bool := x-x = 0\n" + config +
                 "evaluate probability(f(normal(0,1)))");
  CHECK_EQ(r.get("value").as_number(), 1);
  r = solve("def f (n : Real) := let x := normalSample(n,0,1); "
            "probability(mean(x)-mean(x)=0)\n" +
            config + "evaluate f(3)");
  CHECK_EQ(r.get("value").as_number(), 1);
  auto source =
      config +
      "solve for n : Nat in [1,4] where probability(bernoulli(n/5) > 0) >= .5";
  auto first = solve(source), again = solve(source);
  CHECK_EQ(first.dump(), again.dump());
  CHECK_EQ(first.get("solution").as_number(), 3);
  double radius = std::sqrt(std::log(2.0 * 4 / .05) / (2 * 1000));
  CHECK_CLOSE(first.get("upper").as_number() - first.get("value").as_number(),
              radius, 0, 1e-14);
  r = solve(
      config +
      "solve for n : Nat in [1,3] where probability(bernoulli(.5) > 0) >= .5");
  CHECK_EQ(r.get("status").as_string(), "unresolved");
  CHECK(r.get("solution").is_null());
  r = solve(config +
            "solve for n : Nat in [1,3] where probability(bernoulli(if n = 1 "
            "then .5 else 1) > 0) >= .5");
  CHECK_EQ(r.get("status").as_string(), "bounded");
  CHECK_EQ(r.get("solutionLower").as_number(), 1);
  CHECK_EQ(r.get("solutionUpper").as_number(), 2);
  r = solve(config +
            "solve for n : Nat in [1,2] where probability(false) >= .8");
  CHECK_EQ(r.get("status").as_string(), "not_found_in_domain");
}
TEST_CASE(
    "program: invalid syntax, types, domains and tampered plans fail closed") {
  for (const char *source :
       {"evaluate sqrt(true)", "evaluate missing(2)",
        "def f (x : Real) := f(x)\nevaluate f(2)",
        "solve for n : Real in [0,1] where n >= .5",
        "solve for n : Nat in [1.1,5] where n > 2",
        "evaluate probability(normal(0,-1) > 0)", "evaluate normal(0,1)",
        "evaluate probability(probability(normal(0,1)>0) > .2)",
        "evaluate log(0)", "evaluate 1 / 0", "evaluate [1,2][2]",
        "evaluate variance([1])", "evaluate integrate(fun x => x,0,1,3)",
        "import IO\nevaluate 2",
        "simulation { method := analytic }\nevaluate "
        "probability(exp(normal(0,1)) > 2)",
        "simulation { method := simulation }\nsolve for x : Real in [0,1] "
        "where probability(normal(x,1)>0) = .8",
        "evaluate 2\nevaluate 3"})
    invalid(source);
  auto c = compile("evaluate 2 + 3"), plan = c.get("plan"), ir = plan.get("ir");
  ir.set("expression", J());
  plan.set("ir", ir);
  CHECK(!api(plan).get("ok").as_bool());
  CHECK(!c.get("leanKernelChecked").as_bool());
  CHECK(c.get("cpp").as_string().find("int main()") != std::string::npos);
}

TEST_CASE("program: exact derivations preserve decimal source and all roots") {
  auto c = compile("def a : Real := 0.1 + 0.2\nsolve for x : Real in [0,10] "
                   "where a*x + a = 1.4");
  const auto &d = c.get("derivation");
  CHECK_EQ(d.get("kind").as_string(), "exact_algebraic");
  CHECK(!d.get("leanKernelChecked").as_bool());
  CHECK(d.get("leanSource").as_string().find("(11 / 3)") != std::string::npos);
  CHECK(d.get("leanSource").as_string().find("theorem domain_solution") !=
        std::string::npos);
  c = compile("solve for x : Real in [0,10] where x*x = 2");
  CHECK(
      c.get("derivation").get("leanSource").as_string().find("Real.sqrt (8)") !=
      std::string::npos);
  c = compile("solve for x : Real in [0,10] where exp(x) = 2");
  CHECK_EQ(c.get("derivation").get("kind").as_string(), "numerical");
  CHECK_EQ(c.get("derivation").get("leanSource").as_string(), "");
}
