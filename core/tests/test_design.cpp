#include "design/plan.hpp"
#include "dag/model.hpp"
#include "stats/stats.hpp"
#define SSC_STATIC
#include "ssc/ssc.h"
#include "testing.hpp"
#include <fstream>
#include <sstream>

namespace {
using J = ssc::json::Value;
using ssc::dag::object;
J normal() {
  return ssc::json::parse(R"({"schema":"ssc.design-plan","version":1,"method":"normal_known","name":"normal",
    "parameters":{"alternative":"1/2","nullValue":0,"sd":1,"alpha":"1/20","targetPower":"4/5"},
    "search":{"minN":1,"maxN":10000}})");
}
J run(const J& plan) {
  char* answer = ssc_design_json(plan.dump().c_str(),nullptr);
  if (!answer) throw std::runtime_error("missing design response");
  const std::string copy = answer; ssc_free(answer);
  return ssc::json::parse(copy);
}
void bad(const J& p) { CHECK(!run(p).get("ok").as_bool()); }
J stratum(double w, double g, double m0, double m1, double v0, double v1) {
  return object({{"weight",w},{"propensity",g},{"mean0",m0},{"mean1",m1},{"variance0",v0},{"variance1",v1}});
}
J ate() {
  J p=normal(), strata=J::array(); strata.push(stratum(1,.5,0,.5,1,1));
  p.set("method","ate_eif");
  p.set("parameters",object({{"strata",strata},{"nullValue",0},{"alpha",.05},{"targetPower",.8},{"estimator","tmle"}}));
  return p;
}
}

TEST_CASE("design: exact normal example and predecessor") {
  const auto answer=run(normal()); CHECK(answer.get("ok").as_bool());
  const auto& r=answer.get("result");
  CHECK_EQ(r.get("selectedN").as_number(),63);
  CHECK_EQ(r.get("totalParticipants").as_number(),126);
  CHECK_CLOSE(r.get("power").as_number(),.8013023941055797,0,2e-14);
  CHECK_CLOSE(r.get("search").get("predecessorPower").as_number(),.7950080284018122,0,2e-14);
  CHECK(!r.get("verification").get("inputProofAuthenticated").as_bool());
  CHECK_EQ(r.get("verification").get("basis").as_string(),"exact_distributional_formula");
}
TEST_CASE("design: sign, scale and nonzero null invariance") {
  const auto reference=run(normal()).get("result");
  for (double scale : {.00001,1.0,100000.0}) for (double sign : {-1.0,1.0}) {
    J p=normal(), x=p.get("parameters");
    x.set("alternative",scale*(3+sign*.5)); x.set("nullValue",3*scale); x.set("sd",scale); p.set("parameters",x);
    auto r=run(p).get("result");
    CHECK_EQ(r.get("selectedN").as_number(),63);
    CHECK_CLOSE(r.get("power").as_number(),reference.get("power").as_number(),0,2e-13);
  }
}
TEST_CASE("design: domain limits and ambiguous power comparisons") {
  J p=normal(); p.set("search",object({{"minN",64},{"maxN",64}}));
  auto r=run(p).get("result"); CHECK_EQ(r.get("selectedN").as_number(),64);
  CHECK(r.get("search").get("predecessorN").is_null());
  p.set("search",object({{"minN",1},{"maxN",62}}));
  r=run(p).get("result"); CHECK_EQ(r.get("status").as_string(),"not_found_in_domain"); CHECK(r.get("selectedN").is_null());
  p=normal(); auto x=p.get("parameters"); x.set("targetPower",run(p).get("result").get("power")); p.set("parameters",x);
  r=run(p).get("result"); CHECK_EQ(r.get("status").as_string(),"numerically_unresolved");
  CHECK_EQ(r.get("search").get("unresolvedN").as_number(),63); CHECK(r.get("selectedN").is_null());
  x.set("alternative",0); x.set("targetPower",.8); p.set("parameters",x);
  CHECK_EQ(run(p).get("result").get("status").as_string(),"not_found_in_domain");
}
TEST_CASE("design: strict validation and no claimed proof through JSON") {
  for(const char* key : {"sd","alpha","targetPower","alternative","nullValue"}) {
    for (const J& invalid : {J(),J(true),J("nan"),J("1/0"),J("1/-2"),J("1.5/2"),J("1/2/3"),J("9007199254740992/1")}) {
      J p=normal(),x=p.get("parameters"); x.set(key,invalid);p.set("parameters",x);bad(p);
    }
  }
  for(double alpha : {0.0,-.1,1.0,1e-15}) { J p=normal(),x=p.get("parameters");x.set("alpha",alpha);p.set("parameters",x);bad(p); }
  J p=normal(); p.set("formalVerified",true);bad(p);
  p=normal();p.set("version",2);bad(p);
  p=normal();p.set("method","normal_unknown_variance");bad(p);
  p=normal();auto x=p.get("parameters");x.set("allocation","unequal");p.set("parameters",x);bad(p);
  p=normal();p.set("search",object({{"minN",2},{"maxN",1}}));bad(p);
  p=normal();p.set("search",object({{"minN",1.5},{"maxN",100}}));bad(p);
  p=normal();p.set("search",object({{"minN",1},{"maxN",10000001}}));bad(p);
  p=normal();x=p.get("parameters");x.set("sd",0);p.set("parameters",x);bad(p);
  char* empty=ssc_design_json(nullptr,nullptr);CHECK(empty!=nullptr); if(empty) {CHECK(!ssc::json::parse(empty).get("ok").as_bool());ssc_free(empty);}
}
TEST_CASE("design: influence variance counts the declared independent unit") {
  J p=normal();p.set("method","influence_wald");
  J x=object({{"alternative",.5},{"nullValue",0},{"variance",4},{"alpha",.05},{"targetPower",.8},
    {"independentUnit","participant"},{"estimand","ATE"},{"estimator","AIPW"},{"varianceSource","Known planning model"}});
  p.set("parameters",x);auto r=run(p).get("result");
  CHECK_EQ(r.get("selectedN").as_number(),126);CHECK(r.get("totalParticipants").is_null());
  CHECK_EQ(r.get("verification").get("basis").as_string(),"asymptotic_approximation");
  x.set("independentUnit","cluster");p.set("parameters",x);r=run(p).get("result");
  CHECK_EQ(r.get("selectedN").as_number(),126);CHECK_EQ(r.get("sizeUnit").as_string(),"cluster");
  x.set("independentUnit","measurement_row");p.set("parameters",x);bad(p);
}
TEST_CASE("design: ATE EIF variance includes overlap and effect heterogeneity") {
  J p=ate();auto r=run(p).get("result");CHECK_EQ(r.get("selectedN").as_number(),126);CHECK_CLOSE(r.get("varianceCoefficient").as_number(),4,0,1e-14);
  J s=J::array();s.push(stratum(.25,.25,0,1,2,1));s.push(stratum(.75,.75,1,1,1,3));
  auto x=p.get("parameters");x.set("strata",s);p.set("parameters",x);r=run(p).get("result");
  CHECK_CLOSE(r.get("alternative").as_number(),.25,0,1e-14);
  CHECK_CLOSE(r.get("varianceCoefficient").as_number(),.25*(4+2/.75+.75*.75)+.75*(4+4+.25*.25),0,1e-13);
  for(double g:{0.0,1.0}) {s=J::array();s.push(stratum(1,g,0,.5,1,1));x.set("strata",s);p.set("parameters",x);bad(p);}
  s=J::array();s.push(stratum(.9,.5,0,.5,1,1));x.set("strata",s);p.set("parameters",x);bad(p);
}

TEST_CASE("design: solve power and detectable effect with an explicit count") {
  J p=normal();p.set("solve",object({{"for","power"},{"n",63}}));
  auto r=run(p).get("result");CHECK_CLOSE(r.get("power").as_number(),.8013023941055797,0,2e-14);
  p.set("solve",object({{"for","effect"},{"n",63}}));r=run(p).get("result");
  CHECK(r.get("effectLower").as_number()<=r.get("effectUpper").as_number());
  CHECK(r.get("effectUpper").as_number()<.5);
  CHECK(r.get("powerAtLower").as_number()<.8);
  CHECK(r.get("powerAtUpper").as_number()>=.8);
  auto x=p.get("parameters");x.set("alternative",r.get("detectableMagnitude"));p.set("parameters",x);
  p.set("solve",object({{"for","power"},{"n",63}}));CHECK_CLOSE(run(p).get("result").get("power").as_number(),.8,0,1e-14);
  for(const auto& solve : {object({{"for","power"},{"n",0}}),object({{"for","sample_size"},{"n",63}}),object({{"for","other"}})}) {p=normal();p.set("solve",solve);bad(p);}
}

TEST_CASE("design: Lean declaration compiler accepts actual examples and rejects code outside its subset") {
  const std::string source=R"(import DesignCompiler.Input
open DesignCompiler
/- nested /- comment -/ -/
def study : NormalInput := {
  name := "normal"
  alternative := (1 / 2)
  sd := 1 -- known variance
  alpha := 0.05
  targetPower := 0.80
})";
  J request=object({{"op","compile"},{"source",source},{"solveFor","sample_size"}});
  const auto answer=run(request);CHECK(answer.get("ok").as_bool());
  const auto& c=answer.get("result");CHECK(!c.get("leanKernelChecked").as_bool());
  CHECK_EQ(run(c.get("plan")).get("result").get("selectedN").as_number(),63);
  CHECK(c.get("cpp").as_string().find("int main()")!=std::string::npos);
  for(const std::string& suffix : {"\n#eval 1+1", "\ndef extra := 4", "\naxiom falseClaim : False"}) {request.set("source",source+suffix);bad(request);}
  for(const std::string& fields : {"sd := 0", "sd := 1, sd := 2", "variance := 2", "alternative := 1 + 2", "sd := (1 / 0)", "sd := 1 alpha := 0.05"}) {
    request.set("source","import DesignCompiler.Input\nopen DesignCompiler\ndef study : NormalInput := { "+fields+" }");bad(request);
  }
  request.set("source",source);request.set("solveFor","power");request.set("n",63);CHECK(run(request).get("ok").as_bool());
}
