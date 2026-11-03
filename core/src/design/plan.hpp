#pragma once
#include "util/json.hpp"
#include "util/parallel.hpp"

namespace ssc::design {
// Versioned plans emitted by lean/ or constructed by another frontend.
// Parsing a plan is not authentication of a Lean proof.
json::Value execute(const json::Value& plan, Token* token = nullptr);
json::Value compile_source(const json::Value& request);
}
