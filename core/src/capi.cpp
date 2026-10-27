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
#include "dag/model.hpp"
#ifdef __EMSCRIPTEN__
#include <emscripten.h>
// Called only by the engine's owning worker, between parallel replicate batches.
EM_JS(void, ssc_dag_progress, (const char* text), {
  if (Module['onDagProgress']) Module['onDagProgress'](JSON.parse(UTF8ToString(text)));
});
#endif

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

// error_json(msg) copied for the caller, or NULL when there is not even memory for that: no exception may leave
// the C interface (it would end the process)
char* error_answer(const char* msg) noexcept {
  try {
    return dup_string(error_json(msg));
  } catch (...) {
    return nullptr;
  }
}

// The answer to a calculation that runs out of memory, as error_json writes it, ready so that only its copy needs
// memory
constexpr char NO_MEMORY[] = R"({"ok":false,"error":"There is not enough memory for this calculation: use fewer )"
                             R"(combinations of input values or shorter lists."})";

char* no_memory_answer() noexcept {
  char* p = static_cast<char*>(std::malloc(sizeof NO_MEMORY));
  if (p) std::memcpy(p, NO_MEMORY, sizeof NO_MEMORY);
  return p;
}
}  // namespace

extern "C" {

SSC_API const char* ssc_version(void) { return "0.1.0"; }

SSC_API int ssc_hardware_threads(void) { return ssc::hardware_threads(); }

// The catalog, a schema and search results: NULL when there is not enough memory for them (see ssc.h)
SSC_API char* ssc_catalog_json(void) {
  try {
    return dup_string(ssc::catalog_json().dump());
  } catch (const std::bad_alloc&) {
    return nullptr;
  } catch (const std::exception& e) {
    return error_answer(e.what());
  } catch (...) {
    return error_answer("unknown error");
  }
}

SSC_API char* ssc_procedure_json(const char* procedure_id) {
  if (!procedure_id) return nullptr;
  try {
    const ssc::Procedure* p = ssc::find_procedure(procedure_id);
    return p ? dup_string(ssc::procedure_json(*p).dump()) : nullptr;
  } catch (const std::bad_alloc&) {
    return nullptr;
  } catch (const std::exception& e) {
    return error_answer(e.what());
  } catch (...) {
    return error_answer("unknown error");
  }
}

SSC_API char* ssc_search_json(const char* query, int limit) {
  try {
    return dup_string(ssc::search_json(query ? query : "", limit > 0 ? static_cast<std::size_t>(limit) : 0).dump());
  } catch (const std::bad_alloc&) {
    return nullptr;
  } catch (const std::exception& e) {
    return error_answer(e.what());
  } catch (...) {
    return error_answer("unknown error");
  }
}

SSC_API char* ssc_run_json(const char* request_json, ssc_token* token) {
  if (!request_json) return error_answer("empty request");
  try {
    const ssc::json::Value req = ssc::json::parse(request_json);
    ssc::json::Value res = ssc::run_request(req, token ? &token->token : nullptr);
    // Stop while the answer was put together: writing out a large one takes seconds too
    if (token && token->token.cancelled() && res.get("ok").as_bool(false))
      return dup_string(R"({"ok":false,"cancelled":true,"error":"Calculation cancelled."})");
    const std::string out = res.dump();
    char* result = dup_string(out);
    if (!result) throw std::bad_alloc();  // no room for the copy: the message below, rather than no answer at all
    // a full simulation, or a large result (many scenarios): give the memory it used back to the system
    const ssc::Procedure* proc = ssc::find_procedure(req.get("procedure").as_string());
    const bool simulated = proc && proc->simulation && !req.get("options").get("approximate").as_bool(false);
    if (simulated || out.size() > (256u << 10)) ssc::memory::release_free_memory();
    return result;
  } catch (const std::bad_alloc&) {
    ssc::memory::release_free_memory();
    return no_memory_answer();
  } catch (const std::exception& e) {
    return error_answer(e.what());
  } catch (...) {
    return error_answer("unknown error");
  }
}

SSC_API void ssc_free(char* p) { std::free(p); }

SSC_API char* ssc_dag_json(const char* request_json, ssc_token* token) {
  if (!request_json) return error_answer("empty request");
  try {
    ssc::dag::Progress progress;
#ifdef __EMSCRIPTEN__
    progress = [](const ssc::json::Value& value) { ssc_dag_progress(value.dump().c_str()); };
#endif
    const auto result = ssc::dag::run_request(ssc::json::parse(request_json), token ? &token->token : nullptr, progress);
    if (token) token->token.check();
    char* answer = dup_string(ssc::dag::object({{"ok",true},{"result",result}}).dump());
    if (!answer) throw std::bad_alloc();
    ssc::memory::release_free_memory();
    return answer;
  } catch (const std::bad_alloc&) {
    ssc::memory::release_free_memory();
    return no_memory_answer();
  } catch (const std::exception& e) {
    return error_answer(e.what());
  } catch (...) {
    return error_answer("unknown DAG engine error");
  }
}

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
