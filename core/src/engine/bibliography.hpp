#pragma once

#include "engine/procedure.hpp"
#include "util/json.hpp"

namespace ssc::bibliography {
// Resolve the procedure's references from the retrieved, provenance-checked data.
void attach(Procedure& procedure);
json::Value to_json(const Procedure& procedure);
}  // namespace ssc::bibliography
