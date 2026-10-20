#pragma once

#include "engine/procedure.hpp"
#include "util/json.hpp"

namespace ssc {

/// {"categories": [...], "procedures": [{id, name, category, ...}, ...]}
json::Value catalog_json();

/// Full schema of one procedure for the UI.
json::Value procedure_json(const Procedure& p);

/// Procedures ranked for a free-text query (see search.hpp):
/// {"query", "corrected", "partial", "hits": [{id, name, category, score, terms}, ...]}.
json::Value search_json(const std::string& query, std::size_t limit);

/// Execute a calculation request (see core/include/ssc/ssc.h).
json::Value run_request(const json::Value& request, Token* token);

}  // namespace ssc
