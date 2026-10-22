// extern "C" interface (see core/include/ssc/ssc.h).
#include <cstdlib>
#include <cstring>
#include <exception>
#include <new>
#include <string>

#include "engine/procedure.hpp"
#include "engine/runner.hpp"
#include "ssc/ssc.h"
#include "stats/dispatch.hpp"
#include "util/memory.hpp"
#include "util/parallel.hpp"

struct ssc_token {
  ssc::Token token;
};

namespace {
char* dup_string(const std::string& s) {
  char* p = static_cast<char*>(std::malloc(s.size() + 1));
  if (!p) return nullptr;
  std::memcpy(p, s.c_str(), s.size() + 1);
  return p;
}

std::string error_json(const std::string& msg) {
  ssc::json::Value o = ssc::json::Value::object();
  o.set("ok", false);
  o.set("error", msg);
  return o.dump();
}
}  // namespace

extern "C" {

SSC_API const char* ssc_version(void) { return "0.1.0"; }

SSC_API int ssc_hardware_threads(void) { return ssc::hardware_threads(); }

SSC_API char* ssc_catalog_json(void) {
  try {
    return dup_string(ssc::catalog_json().dump());
  } catch (const std::exception& e) {
    return dup_string(error_json(e.what()));
  }
}

SSC_API char* ssc_procedure_json(const char* procedure_id) {
  if (!procedure_id) return nullptr;
  const ssc::Procedure* p = ssc::find_procedure(procedure_id);
  if (!p) return nullptr;
  try {
    return dup_string(ssc::procedure_json(*p).dump());
  } catch (const std::exception& e) {
    return dup_string(error_json(e.what()));
  }
}

SSC_API char* ssc_search_json(const char* query, int limit) {
  try {
    return dup_string(ssc::search_json(query ? query : "", limit > 0 ? static_cast<std::size_t>(limit) : 0).dump());
  } catch (const std::exception& e) {
    return dup_string(error_json(e.what()));
  } catch (...) {
    return dup_string(error_json("unknown error"));
  }
}

SSC_API char* ssc_run_json(const char* request_json, ssc_token* token) {
  if (!request_json) return dup_string(error_json("empty request"));
  try {
    const ssc::json::Value req = ssc::json::parse(request_json);
    ssc::json::Value res = ssc::run_request(req, token ? &token->token : nullptr);
    // Stop while the answer was put together: writing out a large one takes seconds too
    if (token && token->token.cancelled() && res.get("ok").as_bool(false))
      return dup_string(R"({"ok":false,"cancelled":true,"error":"Calculation cancelled."})");
    const std::string out = res.dump();
    char* result = dup_string(out);
    // a full simulation, or a large result (many scenarios): give the memory it used back to the system
    const ssc::Procedure* proc = ssc::find_procedure(req.get("procedure").as_string());
    const bool simulated = proc && proc->simulation && !req.get("options").get("approximate").as_bool(false);
    if (simulated || out.size() > (256u << 10)) ssc::memory::release_free_memory();
    return result;
  } catch (const std::bad_alloc&) {
    ssc::memory::release_free_memory();
    return dup_string(error_json("There is not enough memory for this calculation: use fewer combinations of input "
                                 "values or shorter lists."));
  } catch (const std::exception& e) {
    return dup_string(error_json(e.what()));
  } catch (...) {
    return dup_string(error_json("unknown error"));
  }
}

SSC_API void ssc_free(char* p) { std::free(p); }

SSC_API ssc_token* ssc_token_new(void) { return new (std::nothrow) ssc_token(); }

SSC_API void ssc_token_cancel(ssc_token* token) {
  if (token) token->token.cancel();
}

SSC_API double ssc_token_progress(const ssc_token* token) { return token ? token->token.progress() : 0.0; }

SSC_API void ssc_token_free(ssc_token* token) { delete token; }

SSC_API int ssc_eval(const char* name, const double* args, int nargs, double* out) {
  if (!name || !out || (nargs > 0 && !args)) return 2;
  try {
    return ssc::stats::eval(name, args, nargs, *out) ? 0 : 1;
  } catch (...) {
    return 3;
  }
}

}  // extern "C"
