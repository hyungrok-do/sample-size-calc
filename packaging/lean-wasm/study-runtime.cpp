#include "study-runtime.hpp"
#include "util/json.hpp"
#include <emscripten/emscripten.h>
extern "C" EMSCRIPTEN_KEEPALIVE const char* ssc_evaluate_study(const char* input) {
  static std::string output;
  try {
    auto x=ssc::json::parse(input);
    ssc_lean::Plan p{x.get("variance").as_number(),x.get("effect").as_number(),
      x.get("alpha").as_number(),x.get("target").as_number(),
      int64_t(x.get("minimum").as_number()),int64_t(x.get("maximum").as_number()),
      int64_t(x.get("fixed").as_number()),x.get("mode").as_string(),{}};
    for(const auto &g:x.get("groups").elements())
      p.groups.push_back({g.get("label").as_string(),int64_t(g.get("units").as_number()),int64_t(g.get("participants").as_number())});
    output=ssc_lean::evaluate(p);
  } catch(const std::exception &e) { output="{\"error\":"+ssc_lean::quote(e.what())+"}"; }
  return output.c_str();
}
